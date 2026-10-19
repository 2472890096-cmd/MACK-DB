"""Independent-seed degree-preserving double-edge swaps, restarted from observed graph."""
from pathlib import Path
import csv,datetime,hashlib,json,os,time
import numpy as np
from numba import njit
ROOT=Path(__file__).resolve().parents[1]
REPS=1000;BASE_SEED=20260930

def read(p):
 with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rs):
 p=Path(p);assert ROOT in p.resolve().parents
 with p.with_suffix(p.suffix+'.tmp').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
 os.replace(p.with_suffix(p.suffix+'.tmp'),p)
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
@njit(cache=False)
def rewire(original,n,seed,requested,max_attempts,fs):
 np.random.seed(seed)
 edges=original.copy();adj=np.zeros((n,n),np.bool_)
 for a,b in edges:adj[a,b]=True;adj[b,a]=True
 m=edges.shape[0];success=0;attempts=0
 while success<requested and attempts<max_attempts and m>=2:
  attempts+=1;i=np.random.randint(m);j=np.random.randint(m-1)
  if j>=i:j+=1
  a,b=edges[i];c,d=edges[j]
  if np.random.randint(2)==1:c,d=d,c
  if a==c or a==d or b==c or b==d:continue
  if adj[a,c] or adj[b,d]:continue
  adj[a,b]=False;adj[b,a]=False;adj[c,d]=False;adj[d,c]=False
  adj[a,c]=True;adj[c,a]=True;adj[b,d]=True;adj[d,b]=True
  edges[i,0]=a;edges[i,1]=c;edges[j,0]=b;edges[j,1]=d
  success+=1
 shared=0
 for a,b in edges:
  if fs[a,b]:shared+=1
 return edges,success,attempts,shared

def main():
 summaries=read(ROOT/'07_statistics/stage5B2_primary_summary.tsv');qc=read(ROOT/'qc/stage5B2_qc_summary.tsv')
 qc=[r for r in qc if not (r['check'].startswith('null_') or r['check']=='1000_replicates_50E_successful_swaps')]
 seeds=np.random.SeedSequence(BASE_SEED).generate_state(len(summaries)*REPS,dtype=np.uint32);summary=[]
 # Compile using a separate throwaway graph, not an analysis replicate.
 test=np.array([[0,1],[2,3]],dtype=np.int32);rewire(test,4,np.uint32(1),1,100,np.zeros((4,4),np.bool_))
 for idx,s in enumerate(summaries):
  network=s['network'];start=time.monotonic()
  nodes=[r['gtdb_r226_species'] for r in read(ROOT/'03_common_universe'/f'{network}_common_taxa.tsv')];index={t:i for i,t in enumerate(nodes)};n=len(nodes)
  pub=read(s['published_eligible_edge_file']);fsrows=read(s['fastspar_eligible_edge_file']);E=len(pub)
  original=np.array([[index[r['taxon1']],index[r['taxon2']]] for r in pub],dtype=np.int32);degree=np.bincount(original.ravel(),minlength=n)
  fs=np.zeros((n,n),np.bool_)
  for r in fsrows:a,b=index[r['taxon1']],index[r['taxon2']];fs[a,b]=fs[b,a]=True
  observed=int(sum(fs[a,b] for a,b in original));assert observed==int(s['shared_edges'])
  requested=50*E;max_attempts=max(10000,100*requested);iterations=[];graphs=np.zeros((REPS,E,2),dtype=np.int32)
  print('START_NULL',network,'N',n,'E',E,'requested',requested,'max_attempts',max_attempts,flush=True)
  for rep in range(REPS):
   seed=seeds[idx*REPS+rep]
   edges,success,attempts,shared=rewire(original,n,seed,requested,max_attempts,fs)
   # Verify using Python/NumPy code independent of the rewiring implementation.
   canonical=np.sort(edges,axis=1);codes=canonical[:,0].astype(np.int64)*n+canonical[:,1]
   degok=np.array_equal(np.bincount(edges.ravel(),minlength=n),degree)
   edgeok=len(codes)==E and len(np.unique(codes))==E
   loopfree=bool(np.all(canonical[:,0]<canonical[:,1]));sharedcheck=int(fs[edges[:,0],edges[:,1]].sum())
   valid=success==requested and degok and edgeok and loopfree and shared==sharedcheck
   graphs[rep]=canonical
   iterations.append(dict(network=network,replicate=rep+1,seed=int(seed),node_count=n,edge_count=E,requested_swaps=requested,successful_swaps=success,attempts=attempts,max_attempts=max_attempts,success_rate=success/attempts if attempts else 0,degree_sequence_preserved=bool(degok),edge_count_preserved=bool(edgeok),self_loop_free=loopfree,null_shared_count=shared,randomized_graph_sha256=hashlib.sha256(canonical[np.argsort(codes)].tobytes()).hexdigest(),status='PASS' if valid else 'FAIL'))
   if (rep+1)%100==0:
    write(ROOT/'06_null_model'/f'{network}_degree_null_iterations.tsv',iterations)
    (ROOT/'logs/state.json').write_text(json.dumps(dict(status='RUNNING',step='degree_preserving_null',network=network,replicates_completed=rep+1,timestamp=now()))+'\n')
    print('NULL_PROGRESS',network,rep+1,'seconds',round(time.monotonic()-start,2),flush=True)
  np.savez_compressed(ROOT/'06_null_model'/f'{network}_randomized_graphs.npz',nodes=np.asarray(nodes),original_edges=original,randomized_edges=graphs,degree_sequence=degree,seeds=seeds[idx*REPS:(idx+1)*REPS])
  valid=[r for r in iterations if r['status']=='PASS'];values=np.array([r['null_shared_count'] for r in valid]);ok=len(valid)==REPS
  if ok:
   mean=float(values.mean());sd=float(values.std(ddof=1));q=np.quantile(values,[.025,.975]);p=(1+int((values>=observed).sum()))/(REPS+1)
   enrichment=observed/mean if mean>0 else 'NA';z=(observed-mean)/sd if sd>0 else 'NA'
  else:mean=sd=p=enrichment=z='NA';q=['NA','NA']
  stats=dict(network=network,observed_shared_count=observed,published_eligible_edges=E,node_count=n,requested_replicates=REPS,valid_replicates=len(valid),failed_replicates=REPS-len(valid),requested_swaps_per_replicate=requested,min_successful_swaps=min(r['successful_swaps'] for r in iterations),max_attempts_per_replicate=max_attempts,null_mean_shared=mean,null_sd_shared=sd,null_sd_ddof=1,null_2_5_percentile=q[0],null_97_5_percentile=q[1],null_enrichment=enrichment,null_z=z,null_count_ge_observed=int((values>=observed).sum()) if ok else 'NA',empirical_p=p,empirical_p_display=f'{p:.3f}' if ok else 'NA',observed_above_null_97_5=bool(observed>q[1]) if ok else 'NA',base_seed=BASE_SEED,status='PASS' if ok else 'FAIL',interpretation='unsigned simple-graph degree-preserving null on eligible published network; isolated common nodes retained; each replicate restarts from original graph; no connectedness/sign-degree constraint')
  summary.append(stats);write(ROOT/'06_null_model/stage5B2_degree_null_summary.tsv',summary)
  s.update(null_mean_shared=mean,null_enrichment=enrichment,null_z=z,empirical_p=p)
  for test,result in [('null_degree_preservation',all(r['degree_sequence_preserved'] for r in iterations)),('null_edge_count_preservation',all(r['edge_count_preserved'] for r in iterations)),('null_no_self_loops_or_duplicate_edges',all(r['self_loop_free'] and r['edge_count_preserved'] for r in iterations)),('1000_replicates_50E_successful_swaps',ok)]:
   qc.append(dict(network=network,check=test,status='PASS' if result else 'FAIL',details=f'{len(valid)}/1000 valid; {requested} successful swaps required per replicate; cap={max_attempts}',timestamp=now()))
  write(ROOT/'qc/stage5B2_qc_summary.tsv',qc)
  print('NULL_SUMMARY',json.dumps(stats),flush=True)
 write(ROOT/'07_statistics/stage5B2_primary_summary.tsv',summaries)
 (ROOT/'logs/state.json').write_text(json.dumps(dict(status='RUNNING' if all(r['status']=='PASS' for r in summary) else 'FAILED',step='null_complete_delivery_pending',timestamp=now()))+'\n')
 assert all(r['status']=='PASS' for r in summary),'Null replicate QC failure'
 print('ALL_FIVE_NULLS_COMPLETE',flush=True)
if __name__=='__main__':main()
