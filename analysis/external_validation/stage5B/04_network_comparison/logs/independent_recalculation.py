"""Independent readback implementation: does not import the analysis or rewiring code."""
from pathlib import Path
import csv,datetime,hashlib,json,math,os
from collections import Counter
import numpy as np
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[1]
def read(path):
 with Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(path,rs):
 with Path(path).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
def pairs(rows):
 result={}
 for x in rows:
  a,b=x['taxon1'],x['taxon2'];assert a<b and (a,b) not in result;result[(a,b)]=x
 return result
def ratio(a,b):return a/b if b else 'NA'
def equal(a,b):
 if a=='NA' or b=='NA':return a==b
 return math.isclose(float(a),float(b),rel_tol=1e-12,abs_tol=1e-12)
def main():
 primary=read(ROOT/'07_statistics/stage5B2_primary_summary.tsv')
 sensitivity=[r for r in read(ROOT/'07_statistics/stage5B2_wallen_threshold_sensitivity.tsv') if 'sensitivity' in r['analysis_label']]
 rows=[];qcrecords=[];sg=read(ROOT/'05_sign_concordance/stage5B2_sign_concordance_summary.tsv')
 def verify(n,metric,expected,actual):
  ok=equal(expected,actual);rows.append(dict(network=n,metric=metric,reported_value=expected,recalculated_value=actual,status='PASS' if ok else 'FAIL'));assert ok,(n,metric,expected,actual)
 for s in primary+sensitivity:
  n=s['network'];tag=n+':'+s['analysis_label'];pub=pairs(read(s['published_eligible_edge_file']));fs=pairs(read(s['fastspar_eligible_edge_file']));shared=pub.keys() & fs.keys()
  allpub=pairs(read(ROOT/'01_published_networks'/f'{n}_gtdb_r226_edges.tsv'))
  common={r['gtdb_r226_species'] for r in read(ROOT/'03_common_universe'/f'{n}_common_taxa.tsv')}
  assert set(pub)=={p for p in allpub if set(p)<=common}
  verify(tag,'published_final_gtdb_edges',s['published_final_gtdb_edges'],len(allpub))
  verify(tag,'published_eligible_edges',s['published_eligible_edges'],len(pub));verify(tag,'fastspar_edges',s['fastspar_edges'],len(fs));verify(tag,'shared_edges',s['shared_edges'],len(shared))
  verify(tag,'published_edge_recovery',s['published_edge_recovery'],ratio(len(shared),len(pub)))
  verify(tag,'fastspar_to_published_overlap',s['fastspar_to_published_overlap'],ratio(len(shared),len(fs)))
  verify(tag,'jaccard',s['jaccard'],ratio(len(shared),len(set(pub)|set(fs))))
  same=sum(pub[p]['published_sign'] in ['positive','negative'] and pub[p]['published_sign']==fs[p]['sign'] for p in shared)
  verify(tag,'sign_concordance',s['sign_concordance'],ratio(same,len(shared)))
  signrow=next(r for r in sg if r['network']==n and r['analysis']==('sensitivity' if 'sensitivity' in s['analysis_label'] else 'primary'))
  for a,b,col in [('positive','positive','positive_positive'),('negative','negative','negative_negative'),('positive','negative','positive_negative'),('negative','positive','negative_positive')]:verify(tag,col,signrow[col],sum(pub[p]['published_sign']==a and fs[p]['sign']==b for p in shared))
  assert set(pairs(read(s['shared_pair_file'])))==shared
  restored=read(ROOT/'00_inputs'/f'{n}_restored_taxa.tsv');names=[r['gtdb_r226_species'] for r in restored]
  def matrix(kind):
   with (ROOT.parent/'03_fastspar_result'/n/f'{n}_{kind}.tsv').open() as f:rr=list(csv.reader(f,delimiter='\t'))
   assert rr[0][1:]==[r['count_species_id'] for r in restored]==[r[0] for r in rr[1:]]
   return [[float(x) for x in row[1:]] for row in rr[1:]]
  corr=matrix('correlation');pv=matrix('pvalues');extracted=set();total=0
  paper='paper_matched' in s['analysis_label']
  for i in range(len(names)):
   for j in range(i+1,len(names)):
    rval,p=corr[i][j],pv[i][j]
    passes=abs(rval)>0.2 and p<0.05 if paper else abs(rval)>=0.3 and p<=0.0001
    if passes:
     total+=1
     if names[i] in common and names[j] in common:
      edge=tuple(sorted([names[i],names[j]]));extracted.add(edge)
      assert float(fs[edge]['r'])==rval and float(fs[edge]['pvalue'])==p
      assert fs[edge]['sign']==('positive' if rval>0 else 'negative')
  assert extracted==set(fs);verify(tag,'fastspar_edges_all_input_taxa',s['fastspar_edges_all_input_taxa'],total)
  rows.append(dict(network=tag,metric='independent_matrix_threshold_reconstruction',reported_value=len(fs),recalculated_value=len(extracted),status='PASS'))
  if n.startswith('Wallen') and paper:
   valid=[p for p in shared if pub[p]['published_r']!='NA'];result=spearmanr([float(pub[p]['published_r']) for p in sorted(valid)],[float(fs[p]['r']) for p in sorted(valid)])
   old=read(ROOT/'07_statistics'/f'{n}_shared_r_spearman.tsv')[0]
   verify(tag,'Spearman_pair_count',old['pairs_with_unique_reliable_published_r'],len(valid));verify(tag,'Spearman_rho',old['spearman_rho'],result.statistic);verify(tag,'Spearman_p',old['spearman_p'],result.pvalue)
 for s in primary:
  n=s['network'];it=read(ROOT/'06_null_model'/f'{n}_degree_null_iterations.tsv');null=read(ROOT/'06_null_model/stage5B2_degree_null_summary.tsv');ns=next(r for r in null if r['network']==n)
  arr=np.load(ROOT/'06_null_model'/f'{n}_randomized_graphs.npz',allow_pickle=False);nodes=list(arr['nodes']);index={v:i for i,v in enumerate(nodes)};fs=pairs(read(s['fastspar_eligible_edge_file']));fscode={min(index[a],index[b])*len(nodes)+max(index[a],index[b]) for a,b in fs};E=len(arr['original_edges']);counts=[]
  assert len(it)==1000==arr['randomized_edges'].shape[0]
  assert np.array_equal(np.bincount(arr['original_edges'].ravel(),minlength=len(nodes)),arr['degree_sequence'])
  for k,edges in enumerate(arr['randomized_edges']):
   codes=edges[:,0].astype(np.int64)*len(nodes)+edges[:,1];assert len(edges)==len(set(codes))==E and np.all(edges[:,0]<edges[:,1]);assert np.array_equal(np.bincount(edges.ravel(),minlength=len(nodes)),arr['degree_sequence'])
   assert int(it[k]['successful_swaps'])==50*E and int(it[k]['requested_swaps'])==50*E and it[k]['status']=='PASS'
   assert int(it[k]['seed'])==arr['seeds'][k]
   assert hashlib.sha256(edges[np.argsort(codes)].tobytes()).hexdigest()==it[k]['randomized_graph_sha256']
   count=len(set(codes)&fscode);assert count==int(it[k]['null_shared_count']);counts.append(count)
  v=np.array(counts);obs=int(s['shared_edges']);mu=float(np.mean(v));sd=float(np.std(v,ddof=1));p=(1+sum(x>=obs for x in counts))/1001
  verify(n,'null_mean_shared',ns['null_mean_shared'],mu);verify(n,'null_SD',ns['null_sd_shared'],sd);verify(n,'null_enrichment',s['null_enrichment'],obs/mu if mu else 'NA');verify(n,'null_z',s['null_z'],(obs-mu)/sd if sd else 'NA');verify(n,'empirical_p',s['empirical_p'],p)
  rows.append(dict(network=n,metric='independent_null_1000_graph_degree_edge_swap_validation',reported_value=1000,recalculated_value=1000,status='PASS'))
 write(ROOT/'qc/stage5B2_independent_recalculation.tsv',rows)
 qc=read(ROOT/'qc/stage5B2_qc_summary.tsv')
 for n in [r['network'] for r in primary]:qc.append(dict(network=n,check='independent_recalculation_and_null_archive_readback',status='PASS',details='Separate implementation reloaded final TSVs and original matrices; checked counts, thresholds, ratios, signs, Spearman and 1000 randomized graph archives',timestamp=datetime.datetime.now(datetime.timezone.utc).isoformat()))
 write(ROOT/'qc/stage5B2_qc_summary.tsv',qc)
 print('INDEPENDENT_RECALCULATION_PASS',len(rows))
if __name__=='__main__':main()
