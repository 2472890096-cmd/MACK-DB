#!/usr/bin/env python3
"""Run bounded synthetic scientific tests; never downloads or installs software."""
from pathlib import Path
import argparse,collections,csv,datetime,hashlib,importlib.util,json,math,os,shutil,subprocess,sys,time,traceback
ROOT=Path(__file__).resolve().parents[1]
def load(relative,name):
 spec=importlib.util.spec_from_file_location(name,ROOT/relative);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
def read(p,delimiter='\t'):
 with Path(p).open() as f:return list(csv.DictReader(f,delimiter=delimiter))
def table(p,rows,fields):
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
class Blocked(Exception):pass
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--rscript',default=os.environ.get('MACKDB_RSCRIPT','Rscript'));ap.add_argument('--only',nargs='*');opts=ap.parse_args()
 rscript=shutil.which(opts.rscript);stamp=datetime.datetime.now().strftime('%Y%m%dT%H%M%S');run=ROOT/'tests/results'/stamp;run.mkdir(parents=True,exist_ok=False)
 env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
 results=[];log=None
 def cmd(argv,timeout=60):
  log.write('COMMAND '+json.dumps([str(x) for x in argv])+'\n');log.flush()
  q=subprocess.run(list(map(str,argv)),env=env,cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=timeout)
  log.write(q.stdout);log.flush()
  if q.returncode:raise AssertionError(f'command exit {q.returncode}; see log')
  return q.stdout
 def require_R(packages):
  if not rscript:raise Blocked('Rscript not found: '+opts.rscript)
  expr='p<-c('+','.join(json.dumps(x) for x in packages)+'); for(x in p)cat(x,if(requireNamespace(x,quietly=TRUE))as.character(packageVersion(x))else "MISSING","\\n");'
  s=cmd([rscript,'-e',expr])
  if 'MISSING' in s:raise Blocked('Required installed R package missing; see version probe in log')
 def test_profile(out):
  require_R(['tidyverse','optparse']);data=out/'input';data.mkdir();(out/'samples.txt').write_text('SYNTHETIC1\n')
  values=[0,1,2,3,4,5,6,10]
  lines=[f'Q\tTaxon_{i}\t{12 if i%2==0 else 13}\t0\t{x}' for i,x in enumerate(values)]+['Q\tExcluded_X3_11\t11\t0\t100']
  (data/'SYNTHETIC1_all.txt').write_text('\n'.join(lines)+'\n')
  cmd([rscript,ROOT/'pipeline/metakssd_profile/profile.R','-l',out/'samples.txt','-i',data,'-o',out/'profile'])
  rows=[x.split('\t') for x in (out/'profile/SYNTHETIC1_profile.txt').read_text().splitlines()]
  depths=[0,.02,.08,.69,1.57,2.49,3,7];expected={f'Taxon_{i}':round(x/sum(depths)*100,4) for i,x in enumerate(depths)}
  assert dict((k,float(v)) for k,v in rows)==expected
  assert all(len(v.partition('.')[2])<=4 for k,v in rows)
  assert abs(sum(float(v) for k,v in rows)-100)<.001
  (out/'EXPECTED_PROFILE.json').write_text(json.dumps(expected,indent=2)+'\n')
  log.write('PASS actual profile CLI: default sharekmer=12, X3=11 excluded / 12 included; X5<=5 lookup and >5 subtract3; percentage round4\n')
 def test_epi(out):
  require_R(['vegan']);cmd([rscript,ROOT/'tests/unit/epi.R',ROOT,out])
 def test_zipi(out):
  require_R(['microeco','igraph']);cmd([rscript,ROOT/'tests/unit/zipi.R',ROOT,out])
  m=load('pipeline/zipi/extract_signed_zipi_keys_by_environment.py','test_zipi_extract')
  header,rr,role=m.extract_key_rows(out/'ZIPI_ROLE_BOUNDARIES.tsv')
  assert {x['name'] for x in rr}=={'connector','module','network'}
  _,rr,_=m.extract_key_rows(out/'ZIPI_ACTUAL_NODE_TABLE.tsv');assert all(x['name'] not in {'i','j','k','l'} for x in rr)
  log.write('PASS actual Python candidate exporter excludes peripheral/NA and preserves three hub/connector classes\n')
 def test_dki(out):
  require_R(['dplyr','data.table']);cmd([rscript,ROOT/'tests/unit/dki.R',ROOT,out])
  (out/'projects.txt').write_text('PRJTOYDKI\n')
  cmd([rscript,ROOT/'pipeline/dki/calculate_structural_keystoneness_rk4_fast.R','--project-list',out/'projects.txt','--data-base',out/'data_base','--pred-base',out/'pred_base','--threads','1','--force','TRUE'])
  status=read(out/'pred_base/reports/ks_rk4_fast_status.tsv')[0];assert status['status']=='DONE' and status['n_valid']=='3'
  stream=out/'stream';stream.mkdir();(stream/'LOWMEM_REPRODUCTION.json').write_text(json.dumps(dict(status='PASS',scope='SYNTHETIC_ONLY: test independently checks the legacy full-qtrn counterexample; not production evidence'))+'\n')
  (stream/'taxon_name_map.tsv').write_text('species_id\tname\n1\tSynthetic_A\n2\tSynthetic_B\n3\tSynthetic_C\n')
  cmd([rscript,ROOT/'pipeline/dki/r_scoring_revision/score_R_streaming.R','PRJTOYDKI',out/'data_base/PRJTOYDKI/data',out/'pred_base/PRJTOYDKI/results',stream,'6'])
  scores=read(stream/'R_taxon_summary.tsv');assert len(scores)==1 and abs(float(scores[0]['median_Ks_pred'])-.16)<=1e-15
  st=read(stream/'R_SCORING_STATUS.tsv');st={r['field']:r['value'] for r in st};assert st['status']=='PASS' and st['valid']=='3'
  log.write('PASS actual whole R scorer CLI plus original streaming adapter and real gate on labelled synthetic inputs; no training/prediction\n')
 def test_dki_selection(out):
  m=load('pipeline/dki/prepare_dki_keys_for_mysql.py','test_dki_selection')
  mu,sd,threshold=m.calculate_threshold([.8]+[0.]*7,'sample')
  assert abs(mu-.1)<1e-15 and abs(sd-math.sqrt(.08))<1e-15 and abs(threshold-(.1+2*math.sqrt(.08)))<1e-15
  scores=out/'scores';dest=out/'export';dest.mkdir()
  for project,values in [('PRJTOYHIGH',[.8]+[0.]*7),('PRJTOYZERO',[0.]*8)]:
   f=scores/project/m.RESULT_RELATIVE_PATH;f.parent.mkdir(parents=True)
   with f.open('w') as h:
    w=csv.writer(h);w.writerow(['species','Ks_pred']);w.writerows((f'Synthetic_genus_{i}',v) for i,v in enumerate(values))
  m.process_level('species',scores,dest,{p:['synthetic_env'] for p in ['PRJTOYHIGH','PRJTOYZERO']},'sample')
  # Read actual scientific TSV output, independent of internal counters.
  candidates=read(dest/'key_species_dki_new.tsv')
  assert candidates,'Candidate export table missing'
  sets={p:{r['name'] for r in candidates if r['bioproject']==p} for p in ['PRJTOYHIGH','PRJTOYZERO']}
  assert sets['PRJTOYHIGH']=={'Synthetic_genus_0'},sets
  assert sets['PRJTOYZERO']=={'Synthetic_genus_0','Synthetic_genus_1','Synthetic_genus_2'},sets
  log.write('PASS real DKI selection: full-background sample SD, TopK3, >= threshold and lexical ties; zero-score taxa are valid, not empty inputs\n')
 def test_figure3(out):
  m=load('analysis/figures/figure3/figure3_release.py','test_figure3');rec={(l,method):set() for l in ['species','genus'] for method in m.methods}
  scope={'P0','P1','P2','P3'};envs={'P1':'E1','P2':'E1','P3':'E2','PX':'E3'}
  assignments=[('P1','A','Z'),('P1','B','E'),('P1','C','D'),('P1','D','ZE'),('P1','E','ZD'),('P1','F','ED'),('P1','G','ZED'),('P1','H','Z'),('P2','H','ED'),('P3','G','ZED'),('PX','I','ZED')]
  for p,n,ms in assignments:
   if p in scope:
    for x in ms:rec['species',dict(Z='ZiPi',E='EPI',D='DKI')[x]].add((p,n,envs[p]))
  for p,ms in [('P1','Z'),('P2','ED')]:
   for x in ms:rec['genus',dict(Z='ZiPi',E='EPI',D='DKI')[x]].add((p,'G',envs[p]))
  vals,ctx=m.calc(rec)
  assert [vals['A',cat,'species'] for cat in m.cats.values()]==[1,1,1,1,1,1,3]
  assert [vals[p,c,l] for p,c,l in [('C','Broad >=2 methods','species'),('C','Broad All three methods','species'),('D','Strict >=2 methods','species'),('D','Strict All three methods','species')]]==[5,2,5,1]
  assert sum(v['mask'].bit_count()>=2 for v in ctx['species','Strict'].values())==6
  assert sum(v['mask']==7 for v in ctx['species','Strict'].values())==2
  assert [vals[p,c,'genus'] for p,c in [('C','Broad >=2 methods'),('C','Broad All three methods'),('D','Strict >=2 methods'),('D','Strict All three methods')]]==[1,1,1,0]
  assert len(scope)==4 and all(p!='PX' for records in rec.values() for p,_,_ in records)
  # A second direct set formulation checks species/ge2/All3 without production membership masks.
  for level in ['species','genus']:
   for mode,panel in [('Broad','C'),('Strict','D')]:
    a,b,c=[{(e if mode=='Broad' else p,n) for p,n,e in rec[level,method]} for method in m.methods]
    assert vals[panel,mode+' >=2 methods',level]==len({n for _,n in (a&b)|(a&c)|(b&c)})
    assert vals[panel,mode+' All three methods',level]==len({n for _,n in a&b&c})
  table(out/'FIGURE3_TOY_22_VALUES.tsv',[dict(panel=p,category=c,level=l,value=v) for (p,c,l),v in vals.items()],['panel','category','level','value'])
  reduced={k:{r for r in v if r[0]!='P2'} for k,v in rec.items()};v2,_=m.calc(reduced);assert v2['C','Broad All three methods','species']==1
  log.write('PASS actual calc(): all seven partitions, Broad vs Strict, >=2 includes All3, cross-environment/project unique names, separate ranks, scope filtering and zero-candidate project retained in caller scope\n')
 def test_epi_worker(out):
  require_R(['vegan']);a=out/'abundance.matrix';mp=out/'taxon.matrix'
  saved={k:env.get(k) for k in ['TMPDIR','MACKDB_RSCRIPT']}
  env.update(TMPDIR=str(out),MACKDB_RSCRIPT=rscript)
  try:
   cmd(['bash',ROOT/'pipeline/epi_v3/genus/scripts/timed_rscript.sh','-e','stopifnot(1+1==2)'])
   assert (out/'resource_usage.txt').is_file()
  finally:
   for k,v in saved.items():
    if v is None:env.pop(k,None)
    else:env[k]=v

  # All taxa have prevalence1, so real worker must succeed with three empty keys and a complete summary.
  a.write_text('OTU_id\tS1\tS2\tS3\tS4\n0\t1\t2\t1\t2\n1\t2\t1\t2\t1\n')
  mp.write_text('0\tSynthetic_A\n1\tSynthetic_B\n')
  for level in ['species','genus']:
   job=out/level;argv=[sys.executable,ROOT/'tools/epi_project.py','--level',level,'--project','PRJTOYZERO','--abundance',a,'--mapping',mp,'--run-root',job,'--rscript',rscript]
   cmd(argv);st=read(job/'status/project_status.tsv')[0];assert st['status']=='SUCCESS' and st['exit_code']=='0'
   assert all((job/f'keys/PRJTOYZERO_EPI-{k}.key').stat().st_size==0 for k in ['D1','D2','Q'])
   assert (job/'final_2of3/PRJTOYZERO_EPI-2of3.key').stat().st_size==0
   before={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (job/'keys').iterdir()}
   cmd(argv+['--resume']);assert before=={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (job/'keys').iterdir()}
   cmd(['bash',job/'scripts/run_one_project.sh','PRJTOYZERO',a,mp,'4'])
   assert before=={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (job/'keys').iterdir()}
  log.write('PASS both compatibility workers: real EPI successful zero candidates, status+summary+hash evidence, 2of3 empty valid, resume skips intact results\n')
 def test_epi_votes(out):
  # An explicitly synthetic committed-output fixture tests the real postprocess() function; no scientific-success claim for this fixture.
  os.environ['MACKDB_EPI_RUN_ROOT']=str(out)
  w=load('pipeline/epi_v3/species/scripts/workflow.py','test_epi_votes')
  for d in ['audit','keys','summaries','status/projects','tmp/PRJTOYVOTES/synthetic','final_2of3']:(out/d).mkdir(parents=True,exist_ok=True)
  row=dict(project='PRJTOYVOTES',sample_count=8,taxa_count=5,status='SUCCESS',attempt='synthetic',D1_count=3,D2_count=3,Q_count=2,exit_code=0)
  w.write_tsv(w.status_path(row['project']),[row],w.STATUS_FIELDS)
  for k,values in [('D1',['A','B','C']),('D2',['B','C','D']),('Q',['C','E'])]:(out/f'keys/PRJTOYVOTES_EPI-{k}.key').write_text(''.join(x+'\n' for x in values))
  w.write_tsv(out/'summaries/PRJTOYVOTES_EPI_summary.tsv',[dict(D1_key_count=3,D2_key_count=3,Q_key_count=2)],['D1_key_count','D2_key_count','Q_key_count'])
  w.write_tsv(out/'tmp/PRJTOYVOTES/synthetic/committed_output_hashes.tsv',[dict(path=str(p),sha256=w.sha(p)) for p in w.outputs(row['project'])],['path','sha256'])
  w.write_tsv(out/'audit/protected_fingerprints.tsv',[],['path','sha256']);w.postprocess([row])
  assert w.key_set(out/'final_2of3/PRJTOYVOTES_EPI-2of3.key')=={'B','C'}
  del os.environ['MACKDB_EPI_RUN_ROOT'];log.write('PASS real EPI postprocess on declared synthetic candidate fixture: >=2 includes triple C; final B,C\n')
 def test_fastspar_post(out):
  require_R(['tidyverse','optparse']);names=['A','B','C','D'];corr=[[1,.3,-.3,.2999],[.3,1,.4,.7],[-.3,.4,1,float('nan')],[.2999,.7,float('nan'),1]];pv=[[0,.0001,0,0],[.0001,0,.0001001,.00005],[0,.0001001,0,float('nan')],[0,.00005,float('nan'),0]]
  for name,mat in [('corr',corr),('pv',pv)]:
   with (out/(name+'.csv')).open('w') as f:
    w=csv.writer(f);w.writerow(['']+names);w.writerows([n]+['NA' if math.isnan(x) else x for x in row] for n,row in zip(names,mat))
  (out/'taxa.tsv').write_text(''.join(f'{i}\t{n}\n' for i,n in enumerate(names)))
  cmd([rscript,ROOT/'pipeline/associations/fastspar/sparcc_corr_pvalue.R','-c',out/'corr.csv','-p',out/'pv.csv','-s',out/'taxa.tsv','--corr_num','0.3','--pvalue_num','0.0001','-o',out/'selected'])
  rr=read(out/'selected.csv',',');actual={(r['from'],r['to']):float(r['corr']) for r in rr};assert actual=={('A','B'):.3,('A','C'):-.3,('B','D'):.7},actual
  for name,mat in [('corr',corr),('pv',pv)]:
   with (out/(name+'.tsv')).open('w') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['OTU_id']+list(map(str,range(4))));w.writerows([i]+['NA' if math.isnan(x) else str(x) for x in row] for i,row in enumerate(mat))
  cmd([sys.executable,ROOT/'tools/network_matrix_io.py','--corr',out/'corr.tsv','--pvalue',out/'pv.tsv','--mapping',out/'taxa.tsv','--output',out/'converted'])
  cmd([rscript,ROOT/'pipeline/associations/fastspar/sparcc_corr_pvalue.R','-c',out/'converted/correlation.csv','-p',out/'converted/pvalues.csv','-s',out/'taxa.tsv','--corr_num','0.3','--pvalue_num','0.0001','-o',out/'converted_selected'])
  assert read(out/'converted_selected.csv',',')==rr
  io=load('tools/network_matrix_io.py','matrix_io_test')
  (out/'bad_taxa.tsv').write_text('1\tB\n0\tA\n2\tC\n3\tD\n')
  try:io.convert(out/'corr.tsv',out/'pv.tsv',out/'bad_taxa.tsv',out/'bad_output')
  except ValueError:pass
  else:raise AssertionError('Misaligned axes accepted')
  log.write('PASS labelled TSV to CSV IO adapter preserves values/order, rejects mismatched IDs, and feeds identical scientific postprocessing\n')
  log.write('PASS actual FastSpar postprocessing: three signed upper-triangle edges, inclusive correlation/P boundaries, matching NA masks\n')
 tests=[('PROFILE',test_profile),('EPI_CORE',test_epi),('EPI_WORKERS',test_epi_worker),('EPI_2OF3',test_epi_votes),('ZIPI_ROLES',test_zipi),('DKI_R_SCORING',test_dki),('DKI_SELECTION',test_dki_selection),('FIGURE3_SETS',test_figure3),('FASTSPAR_POSTPROCESS',test_fastspar_post)]
 for name,fn in tests:
  if opts.only and name not in opts.only:continue
  out=run/name;out.mkdir();start=time.monotonic();state='PASS';detail=''
  with (out/'execution.log').open('w') as log:
   try:fn(out)
   except Blocked as e:state='BLOCKED';detail=str(e);log.write(detail+'\n')
   except Exception as e:state='FAIL';detail=f'{type(e).__name__}: {e}';log.write(traceback.format_exc())
  results.append(dict(test=name,status=state,elapsed_seconds=f'{time.monotonic()-start:.3f}',log=str((out/'execution.log').relative_to(ROOT)),detail=detail,test_type='ACTUAL_SMALL_SCIENTIFIC_OR_ENTRY_TEST'))
  print(name,state,detail,flush=True)
 table(run/'TEST_EXECUTION_RESULTS.tsv',results,list(results[0]))
 # Latest complete run stays at the root; partial reruns retain their own separate results.
 if not opts.only:table(ROOT/'TEST_EXECUTION_RESULTS.tsv',results,list(results[0]))
 (run/'RUN_CONTEXT.json').write_text(json.dumps(dict(python=sys.version,rscript=rscript,run=str(run.relative_to(ROOT)),network_used=False,installations=False,full_training=False,full_bootstrap=False),indent=2)+'\n')
 print('RESULTS',run.relative_to(ROOT),flush=True)
 return int(any(x['status']!='PASS' for x in results))
if __name__=='__main__':raise SystemExit(main())
