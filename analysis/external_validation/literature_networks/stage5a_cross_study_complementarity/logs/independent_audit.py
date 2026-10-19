"""Independent audit. Does not import the builder or null generator."""

# v0.1: repository-local IO configuration; scientific functions below are unchanged.
import sys as _path_sys
from pathlib import Path as _RepoPath
_repo_root = next(p for p in _RepoPath(__file__).resolve().parents if (p / "lib/mackdb_paths.py").is_file())
_path_sys.path.insert(0, str(_repo_root / "lib"))
from mackdb_paths import resolve_path
from pathlib import Path
from itertools import combinations
import csv,json,hashlib,datetime,math,sys
from collections import defaultdict
import numpy as np
from scipy.stats import false_discovery_control
R=Path(__file__).resolve().parents[1];EXT=R.parent;B4=EXT/'stage4b_recurrence_recovery';A4=EXT/'stage4a_recurrence';P=['Chen','Wallen','Zhou','Guo','Lacruz','Li'];AUDIT=[];DB=Path(str(resolve_path('environment_network')))
def load(p):
 with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def dump(p,rows):
 with Path(p).open('w',newline='') as f:
  w=csv.DictWriter(f,rows[0].keys(),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def test(scope,check,ok,details=''):
 AUDIT.append(dict(scope=scope,check=check,status='PASS' if ok else 'FAIL',details=str(details)))
 if not ok:
  dump(R/'qc/stage5A_independent_recalculation.tsv',AUDIT);raise AssertionError((scope,check,details))
def eq(a,b):
 if a=='NA' or b=='NA':return a==b
 return math.isclose(float(a),float(b),rel_tol=1e-10,abs_tol=1e-12)
def vals(x):return x.split(';') if x else []
def pair(r):return (r['taxon1'],r['taxon2'])
def truth(x):return x=='True'
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 sets={}
 for n in P+['Lacruz_sensitivity']:
  fn=n+('_primary' if n=='Lacruz' else '')+'_paper_pairs.tsv';rs=load(A4/'paper_level'/fn);sets[n]={(x['canonical_species_a'],x['canonical_species_b']) for x in rs}
  test(n,'frozen_set_copy',sets[n]=={pair(r) for r in load(R/'01_study_sets'/f'{n}_paper_pairs.tsv')})
 nodes={p:set(sum(([a,b] for a,b in es),[])) for p,es in sets.items()}
 allnodes=set.union(*nodes.values());decode={r['database_label']:r['gtdb_r226_species'] for r in load(EXT/'stage3a/qc/03a_database_label_decoding.tsv')}
 # Fresh original-file scan: no Stage 5A support cache used for the independent truth.
 records={};projects=set();seen=set();bad=0
 with DB.open() as f:
  for rownum,r in enumerate(csv.DictReader(f,delimiter='\t'),2):
   k=tuple(sorted((decode[r['species_a']],decode[r['species_b']])));vv=[vals(r[f]) for f in ['proj_ID','proj_corr','proj_pv','proj_sam']]
   if len({len(x) for x in vv})!=1 or len(vv[0])!=int(r['proj_count']) or len(set(vv[0]))!=len(vv[0]):bad+=1
   projects.update(vv[0]);seen.add(k)
   if set(k)<=allnodes:records[k]=(vv,rownum)
 test('ALL','independent_original_database_full_scan',len(seen)==653071 and len(projects)==859 and bad==0,dict(rows=len(seen),projects=len(projects),retained_pairs=len(records)))
 sources={p:set() for p in P}
 for r in load(B4/'source_sets/04b_network_source_config.tsv'):
  if r['network']!='Lacruz_MHO':sources[r['paper']].update(set(vals(r['reported_source_projects']))&projects)
 lookup={}
 names=sorted(allnodes);gi={t:i for i,t in enumerate(names)}
 for A,B in combinations(P,2):
  mat=np.zeros((len(names),len(names)),bool);exclude=sources[A]|sources[B]
  for (a,b),(v,_) in records.items():
   if len(set(v[0])-exclude)>=2:mat[gi[a],gi[b]]=mat[gi[b],gi[a]]=True
  lookup[(A,B)]=mat;z=np.load(R/'00_input_audit'/f'{A}_{B}_source_excluded_lookup.npz')
  test(A+'__'+B,'source_excluded_lookup_against_original_DB',np.array_equal(z['taxa'],names) and np.array_equal(z['recovered'],mat) and set(z['excluded_projects'])==exclude)
 summaries=load(R/'04_source_excluded_recovery/stage5A_directional_recovery_summary.tsv')+load(R/'04_source_excluded_recovery/stage5A_strict_directional_recovery_summary.tsv')
 testsummaries=summaries+load(R/'07_statistics/stage5A_shared_pair_reference.tsv')
 strict={p:{r['gtdb_r226_species'] for r in load(R/'02_universe_audit'/f'{p}_strict_tested_taxa.tsv')} for p in P if (R/'02_universe_audit'/f'{p}_strict_tested_taxa.tsv').exists()}
 by_direction={};allrows={};support_spots=[]
 for s in testsummaries:
  A,B=s['study_A'],s['study_B'];typ=s['analysis'];scope=A+'_to_'+B+'_'+typ
  filename=f'{A}_and_{B}_shared_pair_level.tsv' if typ=='shared_reference' else f'{A}_to_{B}_strict_pair_level.tsv' if typ=='strict' else f'{A}_to_{B}_pair_level.tsv';rs=load(R/'04_source_excluded_recovery'/filename);vis=strict[B] if typ=='strict' else nodes[B]
  expected=sets[A]&sets[B] if typ=='shared_reference' else {e for e in sets[A]-sets[B] if set(e)<=vis}
  test(scope,'independent_directional_definition',len(rs)==len(expected) and {pair(r) for r in rs}==expected)
  before=after=removed=0;known=[];ok=True;metadata_ok=True
  for r in rs:
   k=pair(r);rec=records.get(k);exclude=sources[A]|sources[B];metadata_ok&=set(vals(r['AB_source_projects_to_exclude']))==exclude and set(vals(r['A_source_projects_to_exclude']))==sources[A] and set(vals(r['B_source_projects_to_exclude']))==sources[B] and truth(r['in_A']) and truth(r['in_B'])==(k in sets[B]) and truth(r['B_endpoint1_visible'])==(k[0] in vis) and truth(r['B_endpoint2_visible'])==(k[1] in vis)
   if rec:
    v,line=rec;keep=[i for i,x in enumerate(v[0]) if x not in exclude];delete=[i for i,x in enumerate(v[0]) if x in exclude];nr=len(keep);before+=1;after+=nr>=2;removed+=len(delete);known.append(nr)
    ok&=truth(r['baseline_mackdb_present']) and int(r['baseline_support_count'])==len(v[0]) and int(r['remaining_support_count'])==nr and truth(r['recovered_after_AB_exclusion'])==(nr>=2) and int(r['database_row_number'])==line
    for j,f in enumerate(['support_projects','correlations','pvalues','sample_counts']):
     ok&=vals(r['baseline_'+f])==v[j] and vals(r['remaining_'+f])==[v[j][i] for i in keep] and vals(r['removed_'+f])==[v[j][i] for i in delete]
    ok&=int(r['n_A_supports_removed'])==sum(x in sources[A] for x in v[0]) and int(r['n_B_supports_removed'])==sum(x in sources[B] for x in v[0]) and int(r['n_AB_unique_supports_removed'])==len(delete)
   else:ok&=not truth(r['baseline_mackdb_present']) and r['baseline_support_count']=='NA' and r['remaining_support_count']=='NA' and not truth(r['recovered_after_AB_exclusion']) and not r['remaining_support_projects']
  test(scope,'all_pair_original_project_vectors_and_synchronous_removal',bool(ok),f'{len(rs)} pairs directly checked against fresh original DB scan')
  test(scope,'all_pair_source_union_and_visibility_flags',bool(metadata_ok))
  test(scope,'independent_summary_counts_rates_losses',int(s['directional_pairs'])==len(rs) and int(s['baseline_recovered'])==before and int(s['AB_source_excluded_recovered'])==after and int(s['pairs_lost_after_AB_exclusion'])==before-after and eq(s['baseline_recovery_rate'],before/len(rs) if rs else 'NA') and eq(s['AB_source_excluded_recovery_rate'],after/len(rs) if rs else 'NA') and int(s['total_project_support_entries_removed'])==removed,dict(E=len(rs),baseline=before,recovered=after,lost=before-after,removed_supports=removed))
  if known:
   q=np.quantile(known,[.25,.5,.75]);test(scope,'independent_support_distribution',all(eq(s[f],v) for f,v in zip(['remaining_support_q25','median_remaining_support','remaining_support_q75'],q)) and int(s['remaining_support_count_known_pairs'])==len(known))
  if typ=='primary':by_direction[(A,B)]=s;allrows[(A,B)]=rs
 for A,B in [('Chen','Wallen'),('Wallen','Chen'),('Zhou','Guo'),('Guo','Zhou'),('Lacruz','Li'),('Li','Lacruz')]:
  rs=allrows[(A,B)];sel=np.random.default_rng(20261009+P.index(A)*6+P.index(B)).choice(len(rs),min(10,len(rs)),replace=False)
  for i in sel:
   r=rs[i];rec=records.get(pair(r));support_spots.append(dict(study_A=A,study_B=B,taxon1=r['taxon1'],taxon2=r['taxon2'],sample_method='seeded random without replacement; all pairs when n<10',original_DB_row=rec[1] if rec else 'NA_AUDITED_ABSENT',original_projects=';'.join(rec[0][0]) if rec else '',AB_excluded=';'.join(sorted(sources[A]|sources[B])),remaining_projects=r['remaining_support_projects'],recovered=r['recovered_after_AB_exclusion'],status='PASS'))
  test(A+'_to_'+B,'required_random_original_DB_spotcheck',len(sel)>0,f'{len(sel)} pairs')
 dump(R/'qc/stage5A_original_DB_random_spotchecks.tsv',support_spots)
 # Sensitivity and primary heatmap values independently reconstructed.
 for s in load(R/'07_statistics/stage5A_lacruz_MHO_sensitivity.tsv'):
  A,B=s['study_A'],s['study_B'];ea=sets['Lacruz_sensitivity'] if A=='Lacruz' else sets[A];eb=sets['Lacruz_sensitivity'] if B=='Lacruz' else sets[B];vis={t for e in eb for t in e};d={e for e in ea-eb if set(e)<=vis};rs=load(R/'04_source_excluded_recovery'/f'{A}_to_{B}_MHO_sensitivity_pair_level.tsv');rc=sum(e in records and len(set(records[e][0][0])-(sources[A]|sources[B]))>=2 for e in d)
  test(A+'_to_'+B,'MHO_sensitivity_independent_definition_and_recovery',d=={pair(r) for r in rs} and len(d)==int(s['sensitivity_directional_pair_count']) and rc==int(s['sensitivity_recovered']) and eq(s['sensitivity_recovery'],rc/len(d)))
 for r in load(R/'07_statistics/stage5A_primary_directional_matrix.tsv'):
  A=r['study_A']
  for B in P:
   if A==B:wanted='NA'
   else:
    s=by_direction[(A,B)];wanted=f"{s['AB_source_excluded_recovered']}/{s['directional_pairs']} ({100*float(s['AB_source_excluded_recovery_rate']):.1f}%)" if int(s['directional_pairs']) else '0/0 (NA)'
   test(A+'__'+B,'heatmap_matrix_consistency',r[B]==wanted)
 if '--data-only' in sys.argv:
  dump(R/'qc/stage5A_independent_data_recalculation.tsv',AUDIT)
  print('INDEPENDENT_DATA_AUDIT_PASS',len(AUDIT),flush=True)
  return
 # Stored graphs and actual sampled pairs independently evaluated using original-DB membership.
 ns=load(R/'07_statistics/stage5A_null_summary.tsv')
 for s in ns:
  A,B=s['study_A'],s['study_B'];scope=A+'_to_'+B;model=s['null_model'];E=int(s['observed_pairs']);obs=int(s['observed_recovered']);ordered=tuple(sorted([A,B],key=P.index));mat=lookup[ordered];directory='05_null_a' if model=='NullA' else '06_null_b';its=load(R/directory/f'{scope}_{model[0].lower()+model[1:]}_iterations.tsv')
  if s['status']!='VALID':
   test(scope+'_'+model,'invalid_null_statistics_NA',all(s[k]=='NA' for k in ['enrichment','z_score','empirical_p','BH_adjusted_p']))
   if s['status']=='NA_SWAP_CAP':
    z=np.load(R/directory/f'{scope}_nullB_graphs.npz');local=list(z['taxa']);original={pair(r) for r in allrows[(A,B)]};expected_degree=defaultdict(int)
    for x,y in original:expected_degree[x]+=1;expected_degree[y]+=1
    ok=len(its)==len(z['graphs']) and len(its)>0
    for g,it,st in zip(z['graphs'],its,z['stats']):
     edges={tuple(sorted((local[x],local[y]))) for x,y in g};degree=defaultdict(int)
     for x,y in edges:degree[x]+=1;degree[y]+=1
     ok &= len(edges)==E and dict(degree)==dict(expected_degree) and not edges&sets[B] and all(x!=y for x,y in edges) and list(st)==[int(it['successful_swaps']),int(it['attempts']),int(it['recovered_pairs'])]
    last=its[-1];ok &= int(last['successful_swaps'])<50*E and int(last['attempts'])==int(last['maximum_attempts'])==max(10000,100*50*E) and int(s['null_valid_replicates'])==sum(int(it['successful_swaps'])==50*E for it in its)
    test(scope,'independent_NullB_attempt_cap_and_partial_graph_audit',bool(ok),last)
   if s['status']=='NA_NO_LEGAL_SWAP':
    d={pair(r) for r in allrows[(A,B)]};legal=False
    for (a,b),(c,d0) in combinations(d,2):
     if len({a,b,c,d0})<4:continue
     for x,y in [((a,c),(b,d0)),((a,d0),(b,c))]:
      x=tuple(sorted(x));y=tuple(sorted(y))
      if x not in d|sets[B] and y not in d|sets[B]:legal=True
    test(scope,'independent_no_legal_swap_proof',not legal)
   continue
  test(scope+'_'+model,'1000_replicates',len(its)==int(s['null_valid_replicates'])==1000)
  if model=='NullA':
   z=np.load(R/directory/f'{scope}_nullA_sampled_pairs.npz');cand=z['candidate_edges'];sampled=z['sampled_indices'];wanted={(gi[a],gi[b]) for a,b in combinations(sorted(nodes[B]),2) if (a,b) not in sets[B]}
   test(scope,'independent_nullA_candidate_universe',set(map(tuple,cand))==wanted and np.array_equal(z['taxa'],names))
   good=sampled.shape==(1000,E) and all(len(set(x))==E for x in sampled);ev=cand[sampled];counts=mat[ev[:,:,0],ev[:,:,1]].sum(axis=1)
   for i in [0,499,999]:good&=np.array_equal(sampled[i],np.random.default_rng(int(its[i]['seed'])).choice(len(cand),E,replace=False))
   test(scope,'independent_NullA_samples_size_uniqueness_seed_replay',bool(good))
  else:
   z=np.load(R/directory/f'{scope}_nullB_graphs.npz');graphs=z['graphs'];local=list(z['taxa']);original={pair(r) for r in allrows[(A,B)]};ed={tuple(sorted((local[x],local[y]))) for x,y in z['original_edges']};test(scope,'independent_NullB_original_graph',ed==original)
   degree=defaultdict(int)
   for a,b in original:degree[a]+=1;degree[b]+=1
   counts=[];good=graphs.shape==(1000,E,2)
   for i,g in enumerate(graphs):
    re={tuple(sorted((local[x],local[y]))) for x,y in g};de=defaultdict(int)
    for a,b in re:de[a]+=1;de[b]+=1
    good&=len(re)==E and dict(de)==dict(degree) and not re&sets[B] and all(a!=b for a,b in re)
    counts.append(sum(mat[gi[a],gi[b]] for a,b in re));it=its[i];good&=int(it['requested_swaps'])==int(it['successful_swaps'])==50*E and int(it['successful_swaps'])<=int(it['attempts'])<=int(it['maximum_attempts']) and eq(it['swap_success_rate'],int(it['successful_swaps'])/int(it['attempts'])) and list(z['stats'][i])==[int(it['successful_swaps']),int(it['attempts']),int(it['recovered_pairs'])]
   counts=np.array(counts);test(scope,'independent_NullB_all_graph_degrees_nodes_edges_forbidden_and_50E',bool(good))
  rates=counts/E;mu=float(rates.mean());sd=float(rates.std(ddof=1));p0=(1+int((counts>=obs).sum()))/1001;p=max(.001,p0);test(scope+'_'+model,'independent_recovery_P_enrichment_z_quantiles',all(int(r['recovered_pairs'])==int(c) and eq(r['recovery_rate'],c/E) for r,c in zip(its,counts)) and eq(s['null_mean_recovery'],mu) and eq(s['null_sd'],sd) and eq(s['empirical_p_exact'],p0) and eq(s['raw_empirical_p'],p0) and eq(s['empirical_p'],p) and eq(s['enrichment'],(obs/E)/mu if mu else 'NA') and eq(s['z_score'],((obs/E)-mu)/sd if sd>1e-15 else 'NA') and eq(s['null_q025'],float(np.quantile(rates,.025))) and eq(s['null_q975'],float(np.quantile(rates,.975))))
 for model in ['NullA','NullB']:
  valid=[s for s in ns if s['null_model']==model and s['status']=='VALID'];p=np.array([float(s['raw_empirical_p']) for s in valid]);q=false_discovery_control(p,method='bh');test(model,'independent_BH_scipy_all_valid_primary',all(eq(s['BH_adjusted_p'],v) for s,v in zip(valid,q)),len(valid))
 integrity=[]
 for r in load(R/'00_input_audit/stage5A_input_manifest.tsv'):
  digest=sha(r['path']);ok=digest==r['sha256_before'];integrity.append(dict(**r,sha256_after=digest,status='PASS' if ok else 'FAIL'));test('FROZEN_INPUT','before_after_sha256',ok,r['path'])
 dump(R/'qc/stage5A_frozen_input_integrity.tsv',integrity);dump(R/'qc/stage5A_independent_recalculation.tsv',AUDIT)
 print('INDEPENDENT_AUDIT_PASS',len(AUDIT),'checks',len(support_spots),'random spotchecks',flush=True)
if __name__=='__main__':main()
