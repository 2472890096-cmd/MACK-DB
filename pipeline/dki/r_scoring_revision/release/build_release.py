#!/usr/bin/env python3

# v0.1: repository-local IO configuration; scientific functions below are unchanged.
import sys as _path_sys
from pathlib import Path as _RepoPath
_repo_root = next(p for p in _RepoPath(__file__).resolve().parents if (p / "lib/mackdb_paths.py").is_file())
_path_sys.path.insert(0, str(_repo_root / "lib"))
from mackdb_paths import resolve_path
import csv,json,hashlib,math,shutil,datetime,importlib.util,sys,collections,re
from pathlib import Path
sys.dont_write_bytecode=True
P=Path(str(resolve_path('project_data')));O=Path(__file__).resolve().parents[1];A=P/'mackdb_dki_lowmem_impact_audit_20261012';version='MACKDB_DKI_R_20261012_v1'
def read(p,sep='\t'):return list(csv.DictReader(Path(p).open(),delimiter=sep))
def write(p,rr,fields=None):
 rr=list(rr);p=O/p;p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w') as h:
  w=csv.DictWriter(h,fieldnames=fields or list(rr[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
# Verify every indexed audit output, including all saved R case RDS, without rerunning scoring.
checks=[]
for r in read(A/'ARTIFACT_INDEX.tsv'):
 p=Path(r['path'])
 if p.is_symlink():
  checks.append(dict(path=str(p),status='LINK_REFERENCE_ONLY',SHA256=''));continue
 if r['SHA256']:
  actual=sha(p);assert actual==r['SHA256'],str(p)+' hash mismatch'
  assert p.stat().st_size==int(r['bytes']);checks.append(dict(path=str(p),status='SHA256_PASS',SHA256=actual))
write('evidence/PRIOR_ARTIFACT_SHA256_VALIDATION.tsv',checks)
assert json.loads((A/'evidence/FINAL_VALIDATION.json').read_text())['status']=='PASS'
# Also validate current original small sources against previous protection record.
protection=[]
for r in read(A/'evidence/PROTECTED_INPUTS_AFTER.tsv'):
 p=Path(r['path']);actual=sha(p) if re.fullmatch(r'[0-9a-f]{64}',r['after_SHA256']) else ''
 assert p.stat().st_size==int(r['bytes'])
 if actual:assert actual==r['after_SHA256'],str(p)+' source changed'
 protection.append(dict(path=str(p),bytes=p.stat().st_size,SHA256=actual or 'NOT_REHASHED_LARGE_SOURCE',mtime_ns=p.stat().st_mtime_ns))
write('evidence/ORIGINAL_SOURCE_PROTECTION.tsv',protection)
spec=importlib.util.spec_from_file_location('selector',_repo_root/'pipeline/dki/prepare_dki_keys_for_mysql.py');sel=importlib.util.module_from_spec(spec);sys.modules['selector']=sel;spec.loader.exec_module(sel)
configs=json.loads((A/'evidence/PROJECT_CONFIGS.json').read_text());prior={ (r['project'],r['name']):r for r in read(A/'FOUR_PROJECT_SCORE_COMPARISON.tsv')};cc={r['project']:r for r in read(A/'FOUR_PROJECT_CANDIDATE_COMPARISON.tsv')}
full=[];cand=[];diff=[];originals=[];prov=[];gates=[]
for c in configs:
 p=c['project'];src=Path(c['output']);orig=Path(c['result']);dst=O/'release/new_DKI_rk4_fast'/p/'results';dst.mkdir(parents=True,exist_ok=True);old=O/'archive/lowmem'/p/'results';old.mkdir(parents=True,exist_ok=True)
 gate=json.loads((src/'LOWMEM_REPRODUCTION.json').read_text());assert gate['status']=='PASS'
 status={r['field']:r['value'] for r in read(src/'R_SCORING_STATUS.tsv')};assert status['status']=='PASS' and int(status['valid'])==c['case_count']
 mp={int(r['species_id']):r['name'] for r in read(src/'taxon_name_map.tsv')};summary=read(src/'R_taxon_summary.tsv');assert len(summary)==c['species_count'] and len(mp)==len(set(mp.values()))==len(summary)
 assert set(mp)=={int(r['species']) for r in summary}==set(range(1,len(mp)+1))
 assert sum(int(r['n_valid']) for r in summary)==c['case_count']
 for r in summary:
  assert all(math.isfinite(float(v)) for v in r.values())
 # Check each R summary rank uses min ties on all taxa, independently.
 for vcol,rcol in [('median_Ks_pred','rank_median_Ks_pred'),('q90_Ks_pred','rank_q90_Ks_pred'),('max_Ks_pred','rank_max_Ks_pred')]:
  vs=sorted((float(r[vcol]) for r in summary),reverse=True);ranks={v:vs.index(v)+1 for v in set(vs)}
  assert all(int(r[rcol])==ranks[float(r[vcol])] for r in summary)
 scores={mp[int(r['species'])]:float(r['median_Ks_pred']) for r in summary};order=sorted(scores,key=lambda n:(-scores[n],n));vals=[scores[n] for n in order];mu,sd,thr=sel.calculate_threshold(vals,'sample');K=min(len(vals),max(3,math.ceil(.01*len(vals))));names={n for n in order[:K] if scores[n]>=thr};assert names==set(cc[p]['R_candidates'].split(';'))
 assert len(names)=={'PRJEB80533':17,'PRJNA1345363':6,'PRJNA1357645':19,'PRJEB72274':19}[p]
 byid={int(r['species']):r for r in summary}
 for sid,n in sorted(mp.items()):
  oldr=prior[p,n];r=byid[sid];assert float(oldr['R_Ks'])==scores[n];rank=order.index(n)+1;assert int(oldr['R_final_rank'])==rank
  row=dict(version=version,project=p,environment=c['environment'],species_id=sid,name=n,median_Ks=scores[n],rank=rank,R_summary_rank=int(r['rank_median_Ks_pred']),Top_K=K,background_mean=mu,sample_SD=sd,threshold=thr,candidate=n in names,n_valid=int(r['n_valid']),within_taxon_mean=float(r['mean_Ks_pred']),q90_Ks=float(r['q90_Ks_pred']),max_Ks=float(r['max_Ks_pred']),median_BC=float(r['median_BC_pred']),mean_BC=float(r['mean_BC_pred']))
  full.append(row)
  if n in names:cand.append(dict(name=n,bioproject=p,environment=c['environment'],Ks_pred=scores[n],rank=rank,Top_K=K,mean=mu,sample_SD=sd,threshold=thr,version=version))
  if oldr['candidate_change']!='NEITHER':diff.append(dict(project=p,name=n,environment=c['environment'],change=oldr['candidate_change'],old_Ks=oldr['Lowmem_saved_Ks'],new_Ks=scores[n],old_rank=oldr['Lowmem_final_rank'],new_rank=rank,old_threshold=oldr['Lowmem_threshold'],new_threshold=thr))
 # Production candidate consumer requires exactly species,Ks_pred; summary retains original R schema.
 with (dst/'Keystoneness_structural_median_final.csv').open('w') as h:
  w=csv.writer(h,lineterminator='\n');w.writerow(['species','Ks_pred']);w.writerows((mp[sid],format(float(byid[sid]['median_Ks_pred']),'.17g')) for sid in sorted(mp))
 shutil.copy2(src/'R_taxon_summary.tsv',dst/'R_taxon_summary.full_precision.tsv')
 for fn in ['Keystoneness_structural_summary.csv','Keystoneness_structural_median.csv','Keystoneness_structural_median_final.csv']:
  q=orig/fn;assert q.exists();shutil.copy2(q,old/fn);originals.append(dict(project=p,original=str(q),archive=str(old/fn),SHA256=sha(q),bytes=q.stat().st_size))
 for fn in ['R_taxon_summary.tsv','taxon_name_map.tsv','R_SCORING_STATUS.tsv','LOWMEM_REPRODUCTION.json','INPUT_CORRESPONDENCE.json','EXACT_R_SOURCE_EXPRESSIONS.R.txt']:
  q=O/'evidence'/p/fn;q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src/fn,q)
 # Saved case scores are versioned by checksum and read-only reference; no duplicate large predictions.
 for q,kind in [(src/'R_case_scores.rds','VERIFIED_R_CASE_SCORES'),(orig/'qtrn.csv','EXISTING_PREDICTION'),(orig/'qtst.csv','EXISTING_PREDICTION'),(orig/'best_model.pth','EXISTING_CHECKPOINT'),(orig/'model_metrics.json','EXISTING_MODEL_METADATA')]:
  prov.append(dict(project=p,kind=kind,path=str(q),exists=q.exists(),bytes=q.stat().st_size if q.exists() else '',SHA256=sha(q) if q.exists() and q.stat().st_size<100000000 else 'LARGE_REFERENCE_NO_NEW_HASH',historical_model_path=c['original_model_path']))
 gates.append(dict(project=p,taxa=len(summary),candidate_count=len(names),case_count=c['case_count'],lowmem_reproduction='PASS',R_status='PASS',mapping='BIJECTIVE',all_finite=True,all_R_ranks_checked=True,all_selector_ranks_checked=True,prior_candidates_exact=True))
write('FOUR_PROJECT_R_SCORE_FINAL.tsv',full);write('FOUR_PROJECT_R_CANDIDATES_FINAL.tsv',cand);write('FOUR_PROJECT_CANDIDATE_DIFF.tsv',diff);write('evidence/ORIGINAL_LOWMEM_FILE_SHA256.tsv',originals);write('evidence/PREDICTION_MODEL_PROVENANCE.tsv',prov);write('evidence/RELEASE_ACCEPTANCE.tsv',gates)
shutil.copy2(A/'FOUR_PROJECT_SCORE_COMPARISON.tsv',O/'FOUR_PROJECT_SCORE_DIFF.tsv')
for label in ['GAINED','LOST','RETAINED']:write('evidence/CANDIDATES_'+label+'.tsv',[r for r in diff if r['change']==label])
assert len(full)==5910 and len(cand)==61 and collections.Counter(r['change'] for r in diff)==dict(GAINED=31,LOST=31,RETAINED=30)
write('release/FOUR_PROJECT_MYSQL_ROWS.tsv',[{k:r[k] for k in ['name','bioproject','environment']} for r in cand])
(O/'evidence/RELEASE_VERSION.json').write_text(json.dumps(dict(version=version,generated=datetime.datetime.now().astimezone().isoformat(),status='VALIDATED_CANDIDATE_RELEASE_NOT_DEPLOYED',taxa=len(full),candidates=len(cand),scoring_reused=True,training=False,prediction=False,scoring_reexecuted=False),indent=2)+'\n')
print('PASS: 5910 taxa,61 candidates,31 gained/31 lost/30 retained; all indexed artifacts validated')
