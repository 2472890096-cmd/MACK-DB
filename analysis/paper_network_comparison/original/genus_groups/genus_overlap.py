"""Coarsen BOTH complete side-sets of existing species edges; never infer genus correlations."""
from common import *

def main():
 datasets=json.loads((OUT/'data/comparison_datasets.json').read_text());vpath=N5/'00_inputs/gtdb_r226_species_vocabulary.tsv';vrows=rr(vpath,'frozen_accepted_R226_serialization');vocab={r['gtdb_r226_species']:r['production_encoded_label'] for r in vrows};check('GENUS','frozen_species_vocabulary_unique',len(vocab)==len(vrows))
 labels=defaultdict(set)
 for pid in ['PRJNA751448','PRJNA834801']:
  for r in rr(S6/'data'/f'{pid}_taxon_row_mapping.tsv'):labels[r['gtdb_r226_species']].add(r['production_species_label'])
 used=set();sets={}
 for d in datasets:
  for field in ['published_all_path','project_all_path','published_common_path','project_common_path']:
   p=d[field]
   if p not in sets:_,sets[p]=edge_rows(p);used|={t for e in sets[p] for t in e}
 mapping={};maprows=[]
 for t in sorted(used):
  payload=vocab.get(t,'');group=payload.split('_',1)[0] if payload else '';state='EXACT_FROZEN_R226_VOCABULARY' if payload and group else 'MAPPING_FAILED';actual=sorted(labels[t]);legacy={s.split(' ')[0].split('_')[1] for s in actual}
  if actual:check(t,'actual_production_rule_equals_vocabulary_payload_group',legacy=={group} and all(s.split('_',1)[1]==payload for s in actual))
  mapping[t]=group if state!='MAPPING_FAILED' else None
  maprows.append(dict(gtdb_r226_species=t,frozen_production_encoded_payload=payload or NA,actual_production_labels=';'.join(actual) or NA,species_name_genus_token=t.split(' ')[0],mackdb_genus_group=group or NA,genus_suffix_collapsed=t.split(' ')[0]!=group if group else NA,mapping_status=state,mapping_source=str(vpath),rule='exact frozen species-to-encoded-payload lookup; first payload token; equivalent to standard numeric-prefixed production rule; no invented prefix ID'))
 write(OUT/'species_to_mackdb_genus_group.tsv',maprows)
 collapsed=defaultdict(set)
 for r in maprows:
  if r['mapping_status']!='MAPPING_FAILED':collapsed[r['mackdb_genus_group']].add(r['species_name_genus_token'])
 write(OUT/'genus_group_suffix_collapse_audit.tsv',[dict(mackdb_genus_group=g,distinct_species_name_genus_tokens=len(ts),merged_genus_tokens=';'.join(sorted(ts)),strict_complete_GTDB_genus_grouping=False) for g,ts in sorted(collapsed.items()) if len(ts)>1],['mackdb_genus_group','distinct_species_name_genus_tokens','merged_genus_tokens','strict_complete_GTDB_genus_grouping'])
 rulepath=reg(PROJECT/'script/SpeciesAbundanceMatrix_TO_GenusAbundanceMatrix.py','accepted_genus_rule_source');write(OUT/'genus_rule_evidence.tsv',[dict(source_file=str(rulepath),line_number=j,code=line) for j,line in enumerate(rulepath.read_text().splitlines(),1) if "split('_')[1]" in line])
 def convert(cid,scope,side,edges):
  groups=defaultdict(set);trace=[];fail=within=inter=0
  for a,b in sorted(edges):
   ga,gb=mapping[a],mapping[b];g=tuple(sorted((ga,gb))) if ga and gb else (NA,NA);status='MAPPING_FAILED' if not ga or not gb else 'INTRA_GENUS_GROUP_REMOVED' if ga==gb else 'RETAINED_INTERGENUS_GROUP_PAIR'
   if status=='MAPPING_FAILED':fail+=1
   elif status=='INTRA_GENUS_GROUP_REMOVED':within+=1
   else:inter+=1;groups[g].add((a,b))
   trace.append(dict(comparison_id=cid,scope=scope,side=side,species_pair_id=edgeid((a,b)),species_a=a,species_b=b,genus_group_for_a=ga or NA,genus_group_for_b=gb or NA,canonical_genus_a=g[0],canonical_genus_b=g[1],status=status))
  tracefile=OUT/'genus'/f'{cid}__{scope}__{side}_species_to_genus_pairs.tsv';pairfile=OUT/'genus'/f'{cid}__{scope}__{side}_unique_genus_pairs.tsv';write(tracefile,trace)
  write(pairfile,[dict(genus_pair_id=edgeid(e),genus_a=e[0],genus_b=e[1],supporting_species_pairs=len(ss),supporting_species_pair_ids=';'.join(edgeid(x) for x in sorted(ss))) for e,ss in sorted(groups.items())],['genus_pair_id','genus_a','genus_b','supporting_species_pairs','supporting_species_pair_ids'])
  check(cid+'_'+scope+'_'+side,'conversion_count_conservation',len(edges)==fail+within+inter and inter-len(groups)>=0 and len(edges)-fail-within-(inter-len(groups))==len(groups))
  return groups,dict(input_species_pairs=len(edges),mapping_failed_species_pairs=fail,intra_genus_species_pairs_removed=within,intergenus_species_pairs_before_dedup=inter,species_pair_collapses=inter-len(groups),unique_genus_pairs=len(groups),trace_file=str(tracefile.relative_to(OUT)),genus_pair_file=str(pairfile.relative_to(OUT)))
 results=[];diagnostics=[]
 for d in datasets:
  for scope,pf,df in [('COMMON_SPECIES','published_common_path','project_common_path'),('ALL_MAPPED_SPECIES','published_all_path','project_all_path')]:
   P,D=sets[d[pf]],sets[d[df]];pg,pdiag=convert(d['comparison_id'],scope,'published',P);dg,ddiag=convert(d['comparison_id'],scope,'mackdb',D);GP,GD=set(pg),set(dg);shared=GP&GD;union=GP|GD;SS=P&D;only_from_species_shared=set()
   for a,b in SS:
    ga,gb=mapping[a],mapping[b]
    if ga and gb and ga!=gb:only_from_species_shared.add(tuple(sorted((ga,gb))))
   check(d['comparison_id']+scope,'whole_side_sets_converted_before_intersection',only_from_species_shared<=shared)
   sharedfile=OUT/'genus'/f'{d["comparison_id"]}__{scope}__shared_genus_pairs.tsv'
   write(sharedfile,[dict(genus_pair_id=edgeid(g),genus_a=g[0],genus_b=g[1],published_supporting_species_pairs=len(pg[g]),mackdb_supporting_species_pairs=len(dg[g]),has_exact_shared_species_support=bool(pg[g]&dg[g]),coarse_match_without_exact_species_support=not bool(pg[g]&dg[g])) for g in sorted(shared)],['genus_pair_id','genus_a','genus_b','published_supporting_species_pairs','mackdb_supporting_species_pairs','has_exact_shared_species_support','coarse_match_without_exact_species_support'])
   row=dict(comparison_id=d['comparison_id'],family=d['family'],network=d['network'],analysis_role=d['analysis_role'],project_threshold=d['threshold'],scope=scope,species_published_pairs=len(P),species_mackdb_pairs=len(D),species_shared_pairs=len(SS),species_union_pairs=len(P|D),species_published_pair_overlap=ratio(len(SS),len(P)),species_Jaccard=ratio(len(SS),len(P|D)),genus_published_pairs=len(GP),genus_mackdb_pairs=len(GD),genus_shared_pairs=len(shared),genus_union_pairs=len(union),genus_published_pair_overlap=ratio(len(shared),len(GP)),genus_Jaccard=ratio(len(shared),len(union)),shared_numerator_change=len(shared)-len(SS),published_denominator_change=len(GP)-len(P),overlap_change_percentage_points=(len(shared)/len(GP)-len(SS)/len(P))*100 if GP and P else NA,Jaccard_change_percentage_points=(len(shared)/len(union)-len(SS)/len(P|D))*100 if union and (P|D) else NA,genus_shared_from_species_shared_only=len(only_from_species_shared),additional_coarse_matches_missed_by_shared_only_conversion=len(shared-only_from_species_shared),published_mapping_failed_pairs=pdiag['mapping_failed_species_pairs'],mackdb_mapping_failed_pairs=ddiag['mapping_failed_species_pairs'],published_intra_genus_removed=pdiag['intra_genus_species_pairs_removed'],mackdb_intra_genus_removed=ddiag['intra_genus_species_pairs_removed'],published_species_pair_collapses=pdiag['species_pair_collapses'],mackdb_species_pair_collapses=ddiag['species_pair_collapses'],published_species_edge_source=d[pf],mackdb_species_edge_source=d[df],published_genus_file=pdiag['genus_pair_file'],mackdb_genus_file=ddiag['genus_pair_file'],shared_genus_file=str(sharedfile.relative_to(OUT)),interpretation='legacy MACK-DB genus grouping of existing species edges; not genus-level FastSpar; no coefficient averaging and no sign-concordance condition')
   results.append(row)
   for side,diag in [('published',pdiag),('mackdb',ddiag)]:diagnostics.append(dict(comparison_id=d['comparison_id'],scope=scope,side=side,**diag))
 write(OUT/'species_genus_overlap_comparison.tsv',results);write(OUT/'genus_conversion_diagnostics.tsv',diagnostics);finish();print('GENUS_OVERLAP_COMPLETE',len(results),'condition-scope comparisons',len(maprows),'species mappings',flush=True)
 print('MAIN_COMMON_GENUS',json.dumps([{k:r[k] for k in ['comparison_id','genus_published_pairs','genus_mackdb_pairs','genus_shared_pairs','genus_union_pairs','genus_published_pair_overlap','genus_Jaccard']} for r in results if r['family']=='MAIN_COMPLETE_PROJECT' and r['scope']=='COMMON_SPECIES']),flush=True)
if __name__=='__main__':main()
