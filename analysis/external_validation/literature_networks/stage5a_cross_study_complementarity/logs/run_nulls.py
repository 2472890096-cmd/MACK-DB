"""Specified Stage 5A nulls: 1000 independent restarts; 50E successful swaps."""
from pathlib import Path
import csv,json,datetime,time,math,os
import numpy as np
from numba import njit
ROOT=Path(__file__).resolve().parents[1]
PAPERS=['Chen','Wallen','Zhou','Guo','Lacruz','Li'];REPS=1000;SEED=20261009;NA='NA'
def read(p):
 with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows,fields=None):
 rows=list(rows);p=Path(p);tmp=p.with_suffix(p.suffix+'.tmp')
 with tmp.open('w',newline='') as f:
  w=csv.DictWriter(f,fields or list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
 os.replace(tmp,p)
def qc(scope,name,ok,detail=''):
 p=ROOT/'qc/stage5A_qc_summary.tsv';rs=read(p);rs.append(dict(scope=scope,check=name,status='PASS' if ok else 'FAIL',details=detail,timestamp=datetime.datetime.now(datetime.timezone.utc).isoformat()));write(p,rs)
 if not ok:raise RuntimeError(scope+' '+name+' '+detail)
@njit(cache=False)
def feasible(edges,adj,forbidden):
 E=len(edges)
 for i in range(E):
  a,b=edges[i]
  for j in range(i+1,E):
   c,d=edges[j]
   if a==c or a==d or b==c or b==d:continue
   if not (adj[a,c] or adj[b,d] or forbidden[a,c] or forbidden[b,d]):return True
   if not (adj[a,d] or adj[b,c] or forbidden[a,d] or forbidden[b,c]):return True
 return False
@njit(cache=False)
def rewire(edges0,adj0,forbidden,lookup,seeds,requested,cap):
 reps=len(seeds);E=len(edges0)
 graphs=np.empty((reps,E,2),np.int32);stats=np.zeros((reps,3),np.int64);done=0
 for r in range(reps):
  np.random.seed(seeds[r]);edges=edges0.copy();adj=adj0.copy();success=0;attempts=0
  while success<requested and attempts<cap:
   attempts+=1;i=np.random.randint(E);j=np.random.randint(E)
   if i==j:continue
   a,b=edges[i];c,d=edges[j]
   if a==c or a==d or b==c or b==d:continue
   if np.random.randint(2):c,d=d,c
   if adj[a,c] or adj[b,d] or forbidden[a,c] or forbidden[b,d]:continue
   adj[a,b]=False;adj[b,a]=False;adj[c,d]=False;adj[d,c]=False
   adj[a,c]=True;adj[c,a]=True;adj[b,d]=True;adj[d,b]=True
   edges[i,0]=a;edges[i,1]=c;edges[j,0]=b;edges[j,1]=d;success+=1
  count=0
  for i in range(E):count+=int(lookup[edges[i,0],edges[i,1]])
  graphs[r]=edges;stats[r,0]=success;stats[r,1]=attempts;stats[r,2]=count;done+=1
  if success!=requested:break
 return graphs[:done],stats[:done]
def summarise(row,model,values,status,reason,candidates=NA):
 E=int(row['directional_pairs']);R=int(row['AB_source_excluded_recovered']);obs=R/E if E else NA
 result=dict(study_A=row['study_A'],study_B=row['study_B'],universe_type=row['universe_type'],observed_pairs=E,observed_recovered=R,observed_recovery_rate=obs,null_model=model,null_mean_recovery=NA,null_sd=NA,null_q025=NA,null_q975=NA,enrichment=NA,z_score=NA,empirical_p=NA,empirical_p_exact=NA,raw_empirical_p=NA,BH_adjusted_p=NA,null_valid_replicates=len(values),status=status,reason=reason,candidate_pairs=candidates)
 if status=='VALID':
  rates=np.asarray(values,dtype=float)/E;mean=float(rates.mean());sd=float(rates.std(ddof=1));nge=int((np.asarray(values)>=R).sum());exact=(1+nge)/(REPS+1);p=max(.001,exact)
  result.update(null_mean_recovery=mean,null_sd=sd,null_q025=float(np.quantile(rates,.025)),null_q975=float(np.quantile(rates,.975)),enrichment=obs/mean if mean else NA,z_score=(obs-mean)/sd if sd>1e-15 else NA,empirical_p=p,raw_empirical_p=exact,empirical_p_exact=exact)
  if sd<=1e-15:result['reason']='Valid degree-constrained null with zero variance; z-score NA' if model=='NullB' else 'Valid zero-variance sampled null; z-score NA'
  if mean==0:result['reason']+='; zero null mean, enrichment NA'
 return result
def main():
 sets={p:{(r['taxon1'],r['taxon2']) for r in read(ROOT/'01_study_sets'/f'{p}_paper_pairs.tsv')} for p in PAPERS};summaries=[];qa=[];primary=[r for r in read(ROOT/'04_source_excluded_recovery/stage5A_directional_recovery_summary.tsv') if r['analysis']=='primary']
 settings=dict(replicates=REPS,base_seed=SEED,nullA='sample E unique pairs without replacement from all_pairs(B published network nodes) minus E_B',nullB='each replicate resets original D graph; exactly 50E successful double-edge swaps; forbidden E_B; reject self-loops, duplicates and shared-endpoint proposals',nullB_max_attempts='max(10000,100*50E)',empirical_p='(1+n_null>=observed)/1001; exact column retained; minimum displayed p=0.001; BH uses exact raw_empirical_p',BH='separate NullA and NullB families, all valid primary comparisons',sd='sample standard deviation (ddof=1); zero variance z=NA',enrichment='observed/null mean; zero mean=NA',optional_fisher='not performed; shared-edge analysis remains descriptive')
 (ROOT/'00_input_audit/stage5A_null_protocol.json').write_text(json.dumps(settings,indent=2)+'\n')
 for di,row in enumerate(primary):
  A,B=row['study_A'],row['study_B'];scope=A+'_to_'+B;E=int(row['directional_pairs']);start=time.time();aa,bb=sorted([A,B],key=PAPERS.index)
  data=np.load(ROOT/'00_input_audit'/f'{aa}_{bb}_source_excluded_lookup.npz');names=list(data['taxa']);ix={t:i for i,t in enumerate(names)};lookup=data['recovered'];bnodes=sorted({x for e in sets[B] for x in e});bi=np.array([ix[x] for x in bnodes]);a,b=np.triu_indices(len(bnodes),1);cand=np.column_stack((bi[a],bi[b])).astype(np.int32);bidx={tuple(sorted((ix[x],ix[y]))) for x,y in sets[B]};cand=np.array([x for x in cand if tuple(x) not in bidx],np.int32)
  # taxa globally sorted, so node order and pair orientation are canonical.
  qc(scope,'null_candidate_universe',len(cand)==len(bnodes)*(len(bnodes)-1)//2-len(sets[B]) and E<=len(cand),str(len(cand)))
  seedsA=np.random.SeedSequence([SEED,di,1]).generate_state(REPS,dtype=np.uint32);seedsB=np.random.SeedSequence([SEED,di,2]).generate_state(REPS,dtype=np.uint32)
  afields=['iteration','seed','random_pairs','recovered_pairs','recovery_rate','status'];bfields=afields+['requested_swaps','successful_swaps','attempts','maximum_attempts','swap_success_rate','degree_QC','edge_count_QC','node_set_QC','forbidden_B_edge_QC','simple_graph_QC','graph_index']
  if not E:
   write(ROOT/'05_null_a'/f'{scope}_nullA_iterations.tsv',[],afields);write(ROOT/'06_null_b'/f'{scope}_nullB_iterations.tsv',[],bfields)
   summaries.extend([summarise(row,m,[],'NA_EMPTY_DIRECTIONAL_SET','No directional pairs; recovery rate undefined',len(cand) if m=='NullA' else NA) for m in ['NullA','NullB']]);qa.append(dict(study_A=A,study_B=B,E=0,node_count=0,requested_swaps=0,maximum_attempts=0,valid_replicates=0,status='NA_EMPTY_DIRECTIONAL_SET',reason='Empty observed set'));continue
  values=[];ar=[];samples=np.empty((REPS,E),np.int32)
  for rep,seed in enumerate(seedsA):
   selected=np.random.default_rng(int(seed)).choice(len(cand),E,replace=False);samples[rep]=selected;edges=cand[selected];count=int(lookup[edges[:,0],edges[:,1]].sum());values.append(count)
   ar.append(dict(iteration=rep+1,seed=int(seed),random_pairs=E,recovered_pairs=count,recovery_rate=count/E,status='VALID'))
  np.savez_compressed(ROOT/'05_null_a'/f'{scope}_nullA_sampled_pairs.npz',candidate_edges=cand,sampled_indices=samples,taxa=np.asarray(names),seeds=seedsA)
  write(ROOT/'05_null_a'/f'{scope}_nullA_iterations.tsv',ar,afields);summaries.append(summarise(row,'NullA',values,'VALID','',len(cand)));qc(scope,'NullA_1000_unique_fixed_size_samples',len(ar)==REPS and all(len(set(s))==E for s in samples),'1000 replicates; same A+B-excluded lookup')
  dr=read(ROOT/'03_directional_pairs'/f'{scope}_node_visible_pairs.tsv');local=sorted({r[k] for r in dr for k in ['taxon1','taxon2']});li={x:i for i,x in enumerate(local)};gi=np.array([ix[x] for x in local]);edges=np.array([[li[r['taxon1']],li[r['taxon2']]] for r in dr],np.int32);n=len(local);adj=np.zeros((n,n),np.bool_);adj[edges[:,0],edges[:,1]]=True;adj|=adj.T;forbidden=np.zeros_like(adj)
  for x,y in sets[B]:
   if x in li and y in li:forbidden[li[x],li[y]]=forbidden[li[y],li[x]]=True
  ml=lookup[np.ix_(gi,gi)];degree=np.bincount(edges.ravel(),minlength=n);requested=50*E;cap=max(10000,100*requested)
  qc(scope,'observed_lookup_consistency',int(ml[edges[:,0],edges[:,1]].sum())==int(row['AB_source_excluded_recovered']))
  if not feasible(edges,adj,forbidden):
   reason='Exhaustive two-edge/two-orientation search proves no legal first double-edge swap under forbidden B edges'
   write(ROOT/'06_null_b'/f'{scope}_nullB_iterations.tsv',[],bfields);summaries.append(summarise(row,'NullB',[],'NA_NO_LEGAL_SWAP',reason));qa.append(dict(study_A=A,study_B=B,E=E,node_count=n,requested_swaps=requested,maximum_attempts=cap,valid_replicates=0,status='NA_NO_LEGAL_SWAP',reason=reason));print(scope,'NullA complete; NullB infeasible',flush=True);continue
  graphs,stats=rewire(edges,adj,forbidden,ml,seedsB,requested,cap);br=[];valid=[];allok=True
  for rep,(graph,stat) in enumerate(zip(graphs,stats)):
   success,attempts,count=map(int,stat);deg_ok=np.array_equal(np.bincount(graph.ravel(),minlength=n),degree);edge_ok=len(graph)==E;node_ok=set(graph.ravel())==set(range(n));forbid_ok=not forbidden[graph[:,0],graph[:,1]].any();simple_ok=not (graph[:,0]==graph[:,1]).any() and len(set(map(tuple,np.sort(graph,axis=1))))==E;cnt_ok=count==int(ml[graph[:,0],graph[:,1]].sum());ok=bool(deg_ok and edge_ok and node_ok and forbid_ok and simple_ok and cnt_ok);allok&=ok
   isvalid=ok and success==requested
   if isvalid:valid.append(count)
   br.append(dict(iteration=rep+1,seed=int(seedsB[rep]),random_pairs=E,recovered_pairs=count,recovery_rate=count/E,status='VALID' if isvalid else 'FAILED_SWAP_TARGET',requested_swaps=requested,successful_swaps=success,attempts=attempts,maximum_attempts=cap,swap_success_rate=success/attempts if attempts else NA,degree_QC='PASS' if deg_ok else 'FAIL',edge_count_QC='PASS' if edge_ok else 'FAIL',node_set_QC='PASS' if node_ok else 'FAIL',forbidden_B_edge_QC='PASS' if forbid_ok else 'FAIL',simple_graph_QC='PASS' if simple_ok else 'FAIL',graph_index=rep))
  qc(scope,'NullB_all_graph_invariants_and_lookup',allok,f'{len(graphs)} graphs; degree, exact node set, E edges, simple, no E_B')
  np.savez_compressed(ROOT/'06_null_b'/f'{scope}_nullB_graphs.npz',graphs=graphs,original_edges=edges,taxa=np.asarray(local),forbidden_B=forbidden,seeds=seedsB[:len(graphs)],stats=stats,degree=degree)
  write(ROOT/'06_null_b'/f'{scope}_nullB_iterations.tsv',br,bfields);status='VALID' if len(valid)==REPS else 'NA_SWAP_CAP';reason='' if status=='VALID' else 'Could not generate all 1000 valid replicates at fixed 50E successes within the attempt cap; no replacement null';summaries.append(summarise(row,'NullB',valid,status,reason));qa.append(dict(study_A=A,study_B=B,E=E,node_count=n,requested_swaps=requested,maximum_attempts=cap,valid_replicates=len(valid),status=status,reason=reason));print(scope,'E',E,'NullB',status,'elapsed',round(time.time()-start,2),'sec',flush=True)
  write(ROOT/'07_statistics/stage5A_null_summary_partial.tsv',summaries)
 for model in ['NullA','NullB']:
  valid=[r for r in summaries if r['null_model']==model and r['status']=='VALID'];ordered=sorted(valid,key=lambda r:r['raw_empirical_p']);m=len(ordered);running=1.
  for j in range(m-1,-1,-1):running=min(running,ordered[j]['raw_empirical_p']*m/(j+1));ordered[j]['BH_adjusted_p']=running
  qc(model,'BH_full_valid_directional_family',m>0,f'{m} valid primary comparisons; excluded NA only')
 write(ROOT/'07_statistics/stage5A_null_summary.tsv',summaries);write(ROOT/'06_null_b/stage5A_nullB_feasibility.tsv',qa)
 (ROOT/'logs/state.json').write_text(json.dumps(dict(status='RUNNING',step='nulls_complete_independent_audit_pending',timestamp=datetime.datetime.now(datetime.timezone.utc).isoformat()))+'\n')
 print('NULLS_COMPLETE',flush=True)
if __name__=='__main__':main()
