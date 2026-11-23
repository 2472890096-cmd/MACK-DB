from common import *
from collections import Counter
from scipy.stats import spearmanr
PIDS=['PRJNA751448','PRJNA834801']
def project(pid,ref):
 a=ABUNDANCE[pid+'_abundance.matrix'];cnt=PROJECT/'GDTB226_species_abundance_result/count'/f'{pid}_count.matrix';sp=PROJECT/'GDTB226_species_abundance_result/species'/f'{pid}_species.matrix';heads,ids,ab=matrix(a);ch,ci,counts=matrix(cnt)
 with reg(sp,'production_species_row_dictionary').open() as f:dictionary=list(csv.reader(f,delimiter='\t'))
 check(pid,'abundance_count_dictionary_alignment',ids==ci==[r[0] for r in dictionary] and heads==ch and len(set(ids))==len(ids) and np.isfinite(ab).all() and np.isfinite(counts).all())
 vocab={r['production_encoded_label']:r['gtdb_r226_species'] for r in rr(N5/'00_inputs/gtdb_r226_species_vocabulary.tsv')};taxa=[];taxmap=[]
 for i,(rid,label) in enumerate(dictionary):
  prefix,payload=label.split('_',1);assert prefix.isdigit();name=vocab[payload];taxa.append(name);taxmap.append(dict(project=pid,numeric_row_id=rid,row_position=i+1,production_species_label=label,gtdb_r226_species=name,rule='reuse exact frozen R226 serialization inverse; preserve GTDB suffix'))
 check(pid,'unique_exact_R226_taxon_identity',len(set(taxa))==len(taxa))
 old=rr(N5/'00_inputs'/f'{ref}_restored_taxa.tsv');check(pid,'matches_frozen_Stage5B_restored_taxa',taxa==[r['gtdb_r226_species'] for r in old] and ids==[r['original_species_id'] for r in old]);write(ROOT/'data'/f'{pid}_taxon_row_mapping.tsv',taxmap)
 mats={};audit=[]
 for typ in ['correlation','pvalues','covariance']:
  p=FASTSPAR[pid+'_'+typ+'.tsv'];cols,rows,x=matrix(p);asym=float(np.abs(x-x.T).max());tol=1.000001e-4 if typ=='covariance' else 1e-10;ok=cols==rows==ids and x.shape==(len(ids),len(ids)) and np.isfinite(x).all() and asym<=tol;check(pid,typ+'_identities_dimensions_symmetry_finite',ok,f'max asymmetry={asym}; tolerance={tol}; covariance allows one last printed decimal; original unchanged')
  if typ=='correlation':check(pid,'correlation_bounds',np.abs(x).max()<=1.00001)
  if typ=='pvalues':check(pid,'pvalue_bounds',x.min()>=0 and x.max()<=1)
  mats[typ]=x;audit.append(dict(project=pid,matrix_type=typ,path=str(p),rows=x.shape[0],columns=x.shape[1],abundance_run_count=ab.shape[1],taxon_identity_status='EXACT_FROZEN_R226',matrix_QC='PASS_WITH_COVARIANCE_ROUNDING_NOTE' if asym>1e-10 else 'PASS',maximum_asymmetry=asym,symmetry_tolerance=tol,sha256=sha(p),input_linkage='project filename + exact count/abundance/species row identities + frozen history; original outputs do not embed seed/run manifest'))
 verified=[r for r in rr(RES/'evidence/verified_local_run_columns.tsv') if r['source_bioproject']==pid];verified.sort(key=lambda r:int(r['column_position_1based']));runfile=reg(verified[0]['run_list_file'],'original_ordered_run_list');runs=runfile.read_text().split();check(pid,'complete_project_ordered_run_identity',len(runs)==len(set(runs))==len(heads) and runs==[r['run_accession'] for r in verified] and heads==[r['column_index'] for r in verified] and all(r['profile_values_match_every_retained_taxon']=='Yes' for r in verified))
 write(ROOT/'data'/f'{pid}_run_columns.tsv',[dict(project=pid,column_id=h,run_accession=r) for h,r in zip(heads,runs)])
 np.savez_compressed(ROOT/'data'/f'{pid}_project_matrices.npz',taxa=np.array(taxa),correlation=mats['correlation'],pvalues=mats['pvalues'])
 return dict(pid=pid,taxa=taxa,ids=ids,abundance=ab,count=counts,runs=runs,matrices=mats,audit=audit)

def published():
 comps={};raws={};audits={};universes={};full=rr(EXT/'stage2b/04_all_external_edges_gtdb_r226_stage2b.tsv');audit=rr(EXT/'stage2b/04_edge_mapping_audit_stage2b.tsv')
 for name in ['Guo_AggControl','Guo_AggMECFS']:
  nid=name.replace('Guo_','Guo2023_');rs=[r for r in full if r['network_id']==nid];ar=[r for r in audit if r['network_id']==nid];byid={r['raw_record_id']:r for r in ar};raw=rr(EXT/'parsed'/('guo2023_'+name.split('_',1)[1]+'_final_edges.tsv'));raws[name]=[dict(raw_edge_id=str(i+1),taxon1=r['taxon1_raw'],taxon2=r['taxon2_raw'],published_sign=r['sign'],value=r['association_value']) for i,r in enumerate(raw)]
  comps[name]=[dict(network=name,taxon1=r['species_a'],taxon2=r['species_b'],published_sign=r['sign'] if r['sign'] in ['positive','negative'] else NA,sign_status=r['sign_status'],raw_edge_ids=r['raw_record_ids'],coefficients=json.dumps([byid[i]['association_value'] for i in json.loads(r['raw_record_ids'])]),frozen_file=str(EXT/'stage2b/04_all_external_edges_gtdb_r226_stage2b.tsv')) for r in rs]
  audits[name]=[dict(network=name,raw_edge_id=r['raw_record_id'],raw_taxon1=r['taxon1_raw'],raw_taxon2=r['taxon2_raw'],gtdb_taxon1=r['endpoint1_gtdb_species'],gtdb_taxon2=r['endpoint2_gtdb_species'],mapping_status=r['primary_fate'],both_mapped=r['both_endpoints_accepted']=='True',source_file=str(EXT/'stage2b/04_edge_mapping_audit_stage2b.tsv')) for r in ar]
 for name in ['Guo_Control','Guo_MECFS','Guo_Common','Wallen_PD','Wallen_NHC']:
  path=N5/'01_published_networks'/f'{name}_gtdb_r226_edges.tsv';rs=rr(path);comps[name]=[dict(network=name,taxon1=r['taxon1'],taxon2=r['taxon2'],published_sign=r['published_sign'],sign_status=r['sign_status'],raw_edge_ids=r['raw_edge_ids'],coefficients=r['published_weights_or_r'],frozen_file=str(path)) for r in rs];ar=rr(N5/'01_published_networks'/f'{name}_mapping_audit.tsv');raws[name]=[dict(raw_edge_id=r['raw_edge_id'],taxon1=r['species1'],taxon2=r['species2'],published_sign=r['published_sign'],value=r['published_weight_or_r']) for r in ar];audits[name]=[dict(network=name,raw_edge_id=r['raw_edge_id'],raw_taxon1=r['species1'],raw_taxon2=r['species2'],gtdb_taxon1=r['gtdb_taxon1'],gtdb_taxon2=r['gtdb_taxon2'],mapping_status=r['mapping_status'],both_mapped=r['gtdb_taxon1'] not in ['',NA] and r['gtdb_taxon2'] not in ['',NA],source_file=str(N5/'01_published_networks'/f'{name}_mapping_audit.tsv')) for r in ar]
  universes[name]=rr(N5/'03_common_universe'/f'{name}_published_analysis_taxa.tsv')
 universes['Guo_AggControl']=universes['Guo_Common'];universes['Guo_AggMECFS']=universes['Guo_Common'];rr(N5/'00_inputs/published_universe_evidence.tsv');reg(N5/'00_inputs/Guo_analysis_universe_methods_excerpt.txt','analysis_universe_evidence')
 for n in comps:
  check(n,'frozen_simple_canonical_component',len({edgekey(r) for r in comps[n]})==len(comps[n]) and all(r['taxon1']<r['taxon2'] for r in comps[n]));write(ROOT/'data'/f'{n}_frozen_mapped_edges.tsv',comps[n]);write(ROOT/'data'/f'{n}_raw_edge_mapping_audit.tsv',audits[n])
 return comps,raws,audits,universes

def compare(cid,study,components,pp,comps,raws,audits,universes,threshold,role,sampleclass,fitclass):
 u=defaultdict(list)
 for n in components:
  for r in comps[n]:u[edgekey(r)].append(r)
 union=[]
 for (a,b),rs in sorted(u.items()):
  signs={r['published_sign'] for r in rs};known=signs&{'positive','negative'};conflict=len(known)>1 or any('conflict' in r['sign_status'] for r in rs);state='conflict' if conflict else 'consistent' if len(known)==1 and signs<=known else 'unknown';pubsign=next(iter(known)) if state=='consistent' else NA
  union.append(dict(comparison_id=cid,taxon1=a,taxon2=b,published_sign=pubsign,sign_status=state,cross_component_sign_conflict=len(known)>1,component_count=len(rs),source_networks=pack(r['network'] for r in rs),component_signs=json.dumps({r['network']:r['published_sign'] for r in rs}),component_coefficients=json.dumps({r['network']:json.loads(r['coefficients']) for r in rs}),component_raw_edge_ids=json.dumps({r['network']:json.loads(r['raw_edge_ids']) for r in rs}),coefficient_rule='no averaging; retain every component coefficient; ambiguous sign excluded from sign-concordance denominator'))
 pu={r['gtdb_r226_species'] for n in components for r in universes[n] if r['gtdb_r226_species'] not in ['',NA]};rawtax={r['raw_taxon'] for n in components for r in universes[n]};pt=set(pp['taxa']);common=pu&pt;unionkeys=set(u);eligible={e for e in unionkeys if set(e)<=common};check(cid,'all_mapped_endpoints_in_published_analysis_universe',{t for e in u for t in e}<=pu)
 r,p=pp['matrices']['correlation'],pp['matrices']['pvalues'];ix,jx=np.triu_indices(len(pp['taxa']),1);mask=(np.abs(r[ix,jx])>=.3)&(p[ix,jx]<=1e-4) if threshold=='MACKDB_STANDARD' else (np.abs(r[ix,jx])>.2)&(p[ix,jx]<.05);fs={};allproject=0
 for i,j in zip(ix[mask],jx[mask]):
  allproject+=1;a,b=sorted([pp['taxa'][i],pp['taxa'][j]])
  if a in common and b in common:fs[(a,b)]=dict(taxon1=a,taxon2=b,project_r=float(r[i,j]),project_p=float(p[i,j]),project_sign=sign(r[i,j]))
 shared=eligible&fs.keys();byedge={edgekey(s):s for s in union};signable=[e for e in shared if byedge[e]['published_sign'] in ['positive','negative']];agree=sum(byedge[e]['published_sign']==fs[e]['project_sign'] for e in signable);rawnets=[{tuple(sorted((r['taxon1'],r['taxon2']))) for r in raws[n]} for n in components];rawunion=set.union(*rawnets);unknown_raw=set();loops_raw=set();reasons=[]
 for n in components:
  for ar in audits[n]:
   a,b=ar['gtdb_taxon1'],ar['gtdb_taxon2'];reason='ELIGIBLE' if ar['both_mapped'] and a!=b and {a,b}<=common else 'UNMAPPED_TAXON' if not ar['both_mapped'] else 'SELF_LOOP_AFTER_MAPPING' if a==b else 'MAPPED_ENDPOINT_ABSENT_FROM_PROJECT_INPUT'
   if reason=='UNMAPPED_TAXON':unknown_raw.add(tuple(sorted((ar['raw_taxon1'],ar['raw_taxon2']))))
   if reason=='SELF_LOOP_AFTER_MAPPING':loops_raw.add(tuple(sorted((ar['raw_taxon1'],ar['raw_taxon2']))))
   reasons.append(dict(comparison_id=cid,**ar,eligibility_reason=reason))
 for x in union:x['eligible']=edgekey(x) in eligible;x['in_project_network']=edgekey(x) in fs;x['shared']=edgekey(x) in shared;x['exclusion_reason']='' if edgekey(x) in eligible else 'MAPPED_ENDPOINT_ABSENT_FROM_PROJECT_INPUT'
 write(ROOT/'data'/f'{cid}_published_union.tsv',union);write(ROOT/'data'/f'{cid}_published_eligible.tsv',[x for x in union if x['eligible']],list(union[0]));write(ROOT/'data'/f'{cid}_project_eligible.tsv',fs.values(),['taxon1','taxon2','project_r','project_p','project_sign']);write(ROOT/'data'/f'{cid}_shared_edges.tsv',[dict(**byedge[e],**{k:v for k,v in fs[e].items() if k not in ['taxon1','taxon2']},sign_eligible=e in signable,sign_agreement=byedge[e]['published_sign']==fs[e]['project_sign'] if e in signable else NA) for e in sorted(shared)]);write(ROOT/'data'/f'{cid}_raw_edge_eligibility.tsv',reasons);write(ROOT/'data'/f'{cid}_common_taxa.tsv',[dict(taxon=t) for t in sorted(common)])
 rows=[]
 for n in components:rows.append(dict(comparison_id=cid,component=n,raw_unique_pairs=len({tuple(sorted((r['taxon1'],r['taxon2']))) for r in raws[n]}),mapped_unique_pairs=len(comps[n])))
 if len(components)==2:
  rows.append(dict(comparison_id=cid,component='INTERSECTION',raw_unique_pairs=len(rawnets[0]&rawnets[1]),mapped_unique_pairs=len({edgekey(x) for x in comps[components[0]]}&{edgekey(x) for x in comps[components[1]]})))
 rows.append(dict(comparison_id=cid,component='UNION',raw_unique_pairs=len(rawunion),mapped_unique_pairs=len(union)));write(ROOT/'data'/f'{cid}_component_union_counts.tsv',rows)
 if cid in ['Guo_aggregate_union','Wallen_union_published_matched']:
  frozen=rr(EXT/'stage4a_recurrence/paper_level'/f'{study}_paper_pairs.tsv');check(cid,'union_exactly_matches_Stage4_publication_set',unionkeys=={(x['canonical_species_a'],x['canonical_species_b']) for x in frozen})
 check(cid,'shared_union_and_denominators',shared<=eligible and shared<=fs.keys() and len(shared)<=len(signable)+sum(byedge[e]['published_sign']==NA for e in shared))
 return dict(comparison_id=cid,study=study,analysis_role=role,project=pp['pid'],project_runs=len(pp['runs']),sample_classification=sampleclass,network_design=fitclass,published_components=pack(components),threshold=threshold,published_union_edges=len(rawunion),published_mapped_edges=len(union),eligible_published_edges=len(eligible),published_analysis_raw_taxa=len(rawtax),published_analysis_mapped_taxa=len(pu),project_input_taxa=len(pt),comparable_taxa=len(common),project_edges_all_input_taxa=allproject,project_edges_within_common_taxa=len(fs),shared_edges=len(shared),published_edge_recovery=ratio(len(shared),len(eligible)),jaccard=ratio(len(shared),len(eligible|fs.keys())),jaccard_scope='common_mapped_published_analysis_taxa_and_project_input',sign_eligible_shared_edges=len(signable),sign_agree_edges=agree,sign_concordance=ratio(agree,len(signable)),union_sign_conflict_pairs=sum(x['sign_status']=='conflict' for x in union),cross_component_sign_conflict_pairs=sum(x['cross_component_sign_conflict'] for x in union),unmapped_raw_unique_pairs=len(unknown_raw),self_loop_raw_unique_pairs=len(loops_raw),mapped_untestable_pairs=len(union)-len(eligible),external_independent_validation=False)

def samples_and_rerun(g,w):
 exact=rr(B5/'02_prepare_fastspar_input/qc/stage5B1A_resolved_exact_samples.tsv');sets={n:{r['run_accession'] for r in exact if r['network']==n} for n in ['Guo_Control','Guo_MECFS','Guo_Common','Wallen_PD','Wallen_NHC']};check('Guo','197_equals_disjoint_91_plus_106',len(sets['Guo_Control'])==91 and len(sets['Guo_MECFS'])==106 and not sets['Guo_Control']&sets['Guo_MECFS'] and sets['Guo_Control']|sets['Guo_MECFS']==sets['Guo_Common']==set(g['runs']));ws=sets['Wallen_PD']|sets['Wallen_NHC'];check('Wallen','724_of_725_exact_extra_run',len(sets['Wallen_PD'])==490 and len(sets['Wallen_NHC'])==234 and not sets['Wallen_PD']&sets['Wallen_NHC'] and not ws-set(w['runs']) and set(w['runs'])-ws=={'SRR19064765'})
 sa=[]
 for study,obj,ps,typ,fit in [('Guo',g,sets['Guo_Common'],'EXACT','complete parent sample set; aggregate union is not a single fit; subsampling draws unresolved'),('Wallen',w,ws,'NEAR_COMPLETE_PROJECT','PD/NHC union is not a pooled 724-sample fit; complete project has one additional run')]:
  sa.append(dict(study=study,project=obj['pid'],published_runs=len(ps),local_runs=len(obj['runs']),shared_runs=len(ps&set(obj['runs'])),missing_runs=len(ps-set(obj['runs'])),extra_runs=len(set(obj['runs'])-ps),extra_run_ids=pack(set(obj['runs'])-ps),run_coverage=ratio(len(ps&set(obj['runs'])),len(ps)),sample_classification=typ,missing_percent=100*len(ps-set(obj['runs']))/len(ps),extra_percent=100*len(set(obj['runs'])-ps)/len(ps),analysis_design=fit))
 write(ROOT/'data/stage6B_Guo_Wallen_sample_audit.tsv',sa)
 # Compare historical full-project output with completed 197-run rerun, without restarting FastSpar.
 ch,ids,cnt=matrix(B5/'02_prepare_fastspar_input/count_matrix/Guo_Common_count.matrix');ah,ai,ab=matrix(B5/'02_prepare_fastspar_input/abundance_subset/Guo_Common_abundance.matrix');equalcount=ids==g['ids'] and np.array_equal(cnt,g['count']);ordered=[r['run_accession'] for r in exact if r['network']=='Guo_Common'];equalab=ai==g['ids'] and ordered==g['runs'] and np.array_equal(ab,g['abundance']);check('Guo_Common','same_197_order_and_count_abundance_values',equalcount and equalab)
 mats={}
 for typ in ['correlation','pvalues']:
  p=B5/'03_fastspar_result/Guo_Common'/f'Guo_Common_{typ}.tsv';cols,rows,arr=matrix(p);check('Guo_Common',typ+'_rerun_taxa_symmetry_finite',cols==rows==g['ids'] and np.isfinite(arr).all() and np.allclose(arr,arr.T));mats[typ]=arr
 i,j=np.triu_indices(len(g['taxa']),1);r0=g['matrices']['correlation'][i,j];r1=mats['correlation'][i,j];p0=g['matrices']['pvalues'][i,j];p1=mats['pvalues'][i,j];s0=set(np.where((np.abs(r0)>=.3)&(p0<=1e-4))[0]);s1=set(np.where((np.abs(r1)>=.3)&(p1<=1e-4))[0]);both=s0&s1;changes=[]
 for k in sorted(s0|s1):
  a,b=sorted((g['taxa'][i[k]],g['taxa'][j[k]]));changes.append(dict(taxon1=a,taxon2=b,original_project_r=r0[k],rerun_r=r1[k],original_project_p=p0[k],rerun_p=p1[k],in_original_project=k in s0,in_completed_rerun=k in s1,sign_agreement=sign(r0[k])==sign(r1[k]) if k in both else NA))
 write(ROOT/'data/Guo_Common_original_vs_rerun_edge_audit.tsv',changes)
 out=dict(comparison_id='Guo_Common_original_project_vs_completed_rerun',samples=197,taxa=len(g['taxa']),abundance_values_identical=equalab,count_values_identical=equalcount,all_pairs=len(i),correlation_identical=np.array_equal(r0,r1),pvalues_identical=np.array_equal(p0,p1),correlation_Pearson=float(np.corrcoef(r0,r1)[0,1]),correlation_Spearman=float(spearmanr(r0,r1).statistic),correlation_mean_absolute_difference=float(np.abs(r0-r1).mean()),correlation_max_absolute_difference=float(np.abs(r0-r1).max()),pvalue_mean_absolute_difference=float(np.abs(p0-p1).mean()),original_project_standard_edges=len(s0),completed_rerun_standard_edges=len(s1),shared_standard_edges=len(both),edge_jaccard=ratio(len(both),len(s0|s1)),original_only_edges=len(s0-s1),rerun_only_edges=len(s1-s0),shared_sign_concordance=ratio(sum(sign(r0[k])==sign(r1[k]) for k in both),len(both)),interpretation='Inputs numerically identical; output differences are run-to-run variation compatible with stochastic FastSpar/bootstrap. Original run seeds/version metadata unavailable, so do not attribute exclusively to RNG.')
 out['interpretation']='Inputs numerically identical; existing correlation and pvalue outputs numerically identical; no observed variation. This file comparison does not establish independent stochastic replication.' if out['correlation_identical'] and out['pvalues_identical'] else out['interpretation']
 write(ROOT/'stage6B_guo_fullproject_rerun_audit.tsv',[out]);reg(PROJECT/'script/run_fastspar.sh','original_production_command_rule')

def main():
 QC.clear();ensure_manifests();reg(N5/'STAGE5B2_REPORT.md','frozen_processing_rules');g=project(PIDS[0],'Guo_Common');w=project(PIDS[1],'Wallen_PD');write(ROOT/'data/stage6B_project_matrix_audit.tsv',g['audit']+w['audit']);samples_and_rerun(g,w);comps,raws,audits,universes=published();summ=[]
 specs=[('Guo_aggregate_union','Guo',['Guo_AggControl','Guo_AggMECFS'],g,'MACKDB_STANDARD','PRIMARY','EXACT_PARENT_PROJECT','aggregate_final_network_union_not_single_fit'),('Guo_singlefit_union','Guo',['Guo_Control','Guo_MECFS'],g,'MACKDB_STANDARD','SENSITIVITY_A','EXACT_PARENT_PROJECT','union_of_91_and_106_separate_fits_not_pooled_fit'),('Guo_Common','Guo',['Guo_Common'],g,'MACKDB_STANDARD','SENSITIVITY_B','EXACT','single_fit_same_197_samples_cross_framework'),('Wallen_union_published_matched','Wallen',['Wallen_PD','Wallen_NHC'],w,'PUBLISHED_MATCHED','PRIMARY','NEAR_COMPLETE_PROJECT','union_of_490_and_234_separate_fits_vs_725_run_project'),('Wallen_union_mackdb_standard','Wallen',['Wallen_PD','Wallen_NHC'],w,'MACKDB_STANDARD','PRIMARY','NEAR_COMPLETE_PROJECT','union_of_490_and_234_separate_fits_vs_725_run_project')]
 for cid,study,cs,pp,thr,role,sc,fit in specs:
  s=compare(cid,study,cs,pp,comps,raws,audits,universes,thr,role,sc,fit);summ.append(s);print(json.dumps(s),flush=True)
 write(ROOT/'stage6B_guo_union_comparison.tsv',summ[:1]);write(ROOT/'stage6B_guo_singlefit_and_common_sensitivity.tsv',summ[1:3]);write(ROOT/'stage6B_wallen_union_comparison.tsv',summ[3:]);write(ROOT/'data/stage6B_all_comparison_summary.tsv',summ)
 prior=rr(N5/'07_statistics/stage5B2_primary_summary.tsv');write(ROOT/'data/stage6B_preserved_exact_sample_subgroup_results.tsv',prior)
 finish_inputs();print('PROJECT_COMPARISONS_COMPLETE',flush=True)
if __name__=='__main__':main()
