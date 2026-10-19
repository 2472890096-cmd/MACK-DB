"""Stage 5A: reuse frozen identities/support logic; no taxonomy inference or FastSpar."""

# v0.1: repository-local IO configuration; scientific functions below are unchanged.
import sys as _path_sys
from pathlib import Path as _RepoPath
_repo_root = next(p for p in _RepoPath(__file__).resolve().parents if (p / "lib/mackdb_paths.py").is_file())
_path_sys.path.insert(0, str(_repo_root / "lib"))
from mackdb_paths import resolve_path
from pathlib import Path
from collections import Counter,defaultdict
import csv,datetime,hashlib,json,math,os,re,traceback
import numpy as np
ROOT=Path(__file__).resolve().parents[1];EXT=ROOT.parent;BASE=EXT.parent
A4=EXT/'stage4a_recurrence';B4=EXT/'stage4b_recurrence_recovery';B3=EXT/'stage3b_source_exclusion'
B5=BASE/'9_22_stage5B_same_data_cross_method/04_network_comparison';CACHE=BASE/'9_22_stage5B_same_data_cross_method/01_sample_rescue/cache'
DB=Path(str(resolve_path('environment_network')))
PAPERS=['Chen','Wallen','Zhou','Guo','Lacruz','Li'];NODE='PUBLISHED_NETWORK_NODE_UNIVERSE';STRICT='STRICT_TESTED_UNIVERSE';NA='NA';QC=[];INPUT={}
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def reg(p,role='frozen_input'):
 p=Path(p).resolve()
 if ROOT not in p.parents and str(p) not in INPUT:INPUT[str(p)]=dict(path=str(p),role=role,bytes=p.stat().st_size,sha256_before=sha(p))
 return p
def read(p):
 with reg(p).open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rs,fields=None):
 p=Path(p);assert ROOT in p.resolve().parents;rs=list(rs);fields=fields or list(rs[0]);tmp=p.with_name(p.name+'.tmp')
 with tmp.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
 os.replace(tmp,p)
def check(name,ok,details='',scope='ALL',level='FAIL'):
 QC.append(dict(scope=scope,check=name,status='PASS' if ok else level,details=str(details),timestamp=now()));write(ROOT/'qc/stage5A_qc_summary.tsv',QC)
 if not ok and level=='FAIL':raise RuntimeError(name+': '+str(details))
def info(name,details,scope='ALL',status='INFO'):
 QC.append(dict(scope=scope,check=name,status=status,details=str(details),timestamp=now()));write(ROOT/'qc/stage5A_qc_summary.tsv',QC)
def tokens(v):return v.split(';') if v else []
def key(r):return (r['canonical_species_a'],r['canonical_species_b'])
def packed(s):return ';'.join(sorted(set(s)))
def ratio(a,b):return a/b if b else NA

def study_sets():
 manifest={r['absolute_path']:r for r in read(A4/'qc/04a_formal_output_manifest.tsv')}
 overlaps={r['paper']:r for r in read(A4/'summary/04a_primary_paper_pair_overlap_counts.tsv')};sensover={r['paper']:r for r in read(A4/'summary/04a_sensitivity_paper_pair_overlap_counts.tsv')}
 vocab={r['gtdb_r226_species'] for r in read(B5/'00_inputs/gtdb_r226_species_vocabulary.tsv')};sets={};audit=[]
 for n in PAPERS+['Lacruz_sensitivity']:
  fname=n+('_primary' if n=='Lacruz' else '')+'_paper_pairs.tsv';p=A4/'paper_level'/fname;rows=read(p);es=[key(r) for r in rows];nodes={x for e in es for x in e}
  check('frozen_pair_file_hash',sha(p)==manifest[str(p)]['sha256'],p,scope=n)
  expected=int(sensover['Lacruz']['Lacruz']) if n=='Lacruz_sensitivity' else int(overlaps[n][n])
  loops=sum(x==y for x,y in es);dups=len(es)-len(set(es));canon=all(x<y for x,y in es)
  check('pair_count_consistency',len(es)==expected,len(es),scope=n)
  check('canonical_GTDB_R226_pair_integrity',not loops and not dups and canon and nodes<=vocab and all(r['pair_key']==r['canonical_species_a']+'\t'+r['canonical_species_b'] for r in rows),scope=n)
  sets[n]=set(es);audit.append(dict(study=n,pair_file=str(p),pair_count=len(es),unique_taxa=len(nodes),duplicate_pairs=dups,self_loops=loops,canonical_pair_status='PASS',sha256=sha(p),status='PASS'))
  write(ROOT/'01_study_sets'/f'{n}_paper_pairs.tsv',[dict(taxon1=x,taxon2=y) for x,y in sorted(es)],['taxon1','taxon2'])
  write(ROOT/'02_universe_audit'/f'{n}_published_network_nodes.tsv',[dict(study=n,gtdb_r226_species=t,universe_type=NODE) for t in sorted(nodes)])
  if n=='Lacruz':check('Lacruz_primary_excludes_MHO',all(r['MHO_included_in_analysis']=='0' and 'Lacruz_MHO' not in r['subnetwork_sources'].split('|') for r in rows),scope=n)
 for a in PAPERS:
  for b in PAPERS:check('frozen_intersection_consistency',len(sets[a]&sets[b])==int(overlaps[a][b]),scope=a+'__'+b)
 write(ROOT/'00_input_audit/stage5A_frozen_pairset_audit.tsv',audit)
 return sets,vocab

def strict_universes(sets,vocab):
 strict={};audit=[];raw_evidence=[]
 # Guo/Wallen full analysis identities and frozen mapping decisions already audited in Stage 5B-2.
 for study,networks in [('Guo',['Guo_Common']),('Wallen',['Wallen_PD','Wallen_NHC'])]:
  records={};sources=[]
  for n in networks:
   p=B5/'03_common_universe'/f'{n}_published_analysis_taxa.tsv';sources.append(str(p));rows=read(p)
   for r in rows:
    raw=r['raw_taxon'];t=r['gtdb_r226_species']
    if raw in records:check('strict_frozen_decision_consistency',records[raw]['gtdb_r226_species']==t,raw,study)
    records[raw]=dict(study=study,raw_taxon=raw,gtdb_r226_species=t,mapping_status=r['mapping_status'],decision_source=str(p))
  ev=read(B5/'00_inputs/published_universe_evidence.tsv');relevant=[r for r in ev if r['network'] in networks]
  check('strict_universe_direct_evidence',len(relevant)==len(networks) and all(r['universe_type']=='complete_published_analysis_taxa' for r in relevant),scope=study)
  strict[study]={r['gtdb_r226_species'] for r in records.values() if r['gtdb_r226_species'] not in ['',NA]}
  check('strict_includes_frozen_network_nodes',{t for e in sets[study] for t in e}<=strict[study]<=vocab,scope=study)
  write(ROOT/'02_universe_audit'/f'{study}_strict_tested_taxa.tsv',[dict(study=study,gtdb_r226_species=t,universe_type=STRICT,scope='mapped portion of complete author analysis taxa; no new mapping') for t in sorted(strict[study])])
  write(ROOT/'02_universe_audit'/f'{study}_strict_raw_identity_audit.tsv',records.values())
  audit.append(dict(study=study,strict_available='YES',strict_raw_tested_taxa=len(records),strict_mapped_unique_GTDB_taxa=len(strict[study]),strict_unresolved_raw_taxa=sum(r['gtdb_r226_species'] in ['',NA] for r in records.values()),published_network_nodes=len({t for e in sets[study] for t in e}),evidence=' | '.join(r['evidence'] for r in relevant),evidence_files=';'.join(sources),limitation='Paper-level union of known inference taxa; does not assert every taxon pair was statistically evaluated in every component; incomplete frozen mapping coverage'))
 # Li supplied an exhaustive upper triangle including diagonal, which proves identity coverage.
 lp=EXT/'stage2c/li2024/raw/sparcc_all_species.txt';reg(lp,'full_published_pairwise_table')
 taxa=set();observed=set();rows_count=0;loops=0
 with lp.open() as f:
  for r in csv.DictReader(f,delimiter='\t'):
   a,b=r['Var1'],r['Var2'];taxa.update([a,b]);observed.add(tuple(sorted([a,b])));loops+=a==b;rows_count+=1
 check('Li_exhaustive_upper_triangle_including_diagonal',rows_count==len(observed)==len(taxa)*(len(taxa)+1)//2 and loops==len(taxa)==1164,dict(rows=rows_count,unique_taxa=len(taxa),diagonal=loops),scope='Li')
 del observed
 specific=read(EXT/'stage2c/li2024/taxonomy/li2024_taxonomy_mapping.tsv');spec={r['original_taxon']:r for r in specific}
 pool=defaultdict(list)
 gpath=EXT/'stage2b/03_external_taxa_to_gtdb_r226_stage2b.tsv';wpath=BASE/'figure56/PMID_36376318_taxonomy_rescue_overlap_v44/02_taxon_mapping_audit.tsv'
 def norm(s):return ' '.join(s.replace('_',' ').split()).casefold()
 for r in read(gpath):
  if r['mapping_status'].startswith('mapped_'):
   pool[norm(r['raw_taxon'])].append((r['gtdb_r226_species'],str(gpath),r['raw_taxon_id']))
 for r in read(wpath):
  if r['gtdb_status'].startswith('accepted'):pool[norm(r['external_taxon'])].append((r['gtdb_species'],str(wpath),r['external_taxon']))
 mapped=[]
 for t in sorted(taxa):
  if t in spec:
   q=spec[t];name=q['gtdb_r226_species'] if q['mapping_status'].startswith('mapped_') else NA;status=q['mapping_status'];source=str(EXT/'stage2c/li2024/taxonomy/li2024_taxonomy_mapping.tsv')
  else:
   hits=pool.get(norm(t),[]);names={x[0] for x in hits};name=next(iter(names)) if len(names)==1 else NA;status='reused_frozen_unique_exact_name_decision' if len(names)==1 else 'unresolved_no_unique_frozen_decision';source=';'.join(sorted({x[1]+'#'+x[2] for x in hits}))
  mapped.append(dict(study='Li',raw_taxon=t,gtdb_r226_species=name,mapping_status=status,decision_source=source or NA))
 strict['Li']={r['gtdb_r226_species'] for r in mapped if r['gtdb_r226_species']!=NA}
 check('Li_strict_frozen_labels_and_network_node_coverage',{t for e in sets['Li'] for t in e}<=strict['Li']<=vocab,scope='Li')
 write(ROOT/'02_universe_audit/Li_strict_tested_taxa.tsv',[dict(study='Li',gtdb_r226_species=t,universe_type=STRICT,scope='mapped part of complete 1164-species pairwise table') for t in sorted(strict['Li'])])
 write(ROOT/'02_universe_audit/Li_strict_raw_identity_audit.tsv',mapped)
 audit.append(dict(study='Li',strict_available='YES',strict_raw_tested_taxa=len(taxa),strict_mapped_unique_GTDB_taxa=len(strict['Li']),strict_unresolved_raw_taxa=sum(r['gtdb_r226_species']==NA for r in mapped),published_network_nodes=len({t for e in sets['Li'] for t in e}),evidence=f'Exhaustive N(N+1)/2 upper triangle: {rows_count} distinct pairs including {loops} diagonals, N=1164; no species IDs missing',evidence_files=str(lp),limitation='Only existing exact taxon decisions reused; unmapped taxa excluded. Identity set proves inference input membership, not significance testing success for every pair.'))
 unavailable=[('Chen',[CACHE/'Chen_author_network.R',CACHE/'Chen_Source_Data.xlsx'],'Author code references sp abundance object without a complete taxa identity table; available Source Data sheets are sample-level Figure 5C/2A/2B data, not the full inference input.'),('Zhou',[EXT/'raw/zhou2022/NIHMS1880430-supplement-Supplementary_Table_S6.xlsx',CACHE.parent.parent/'00_metadata_cache/Zhou_author_network.Rmd',CACHE/'Zhou_taxa_list.rda'],'Supplementary edge tables are not exhaustive pairwise matrices; author code applies cohort-specific filtering and separate analysis objects. Available taxa.list is not certified as the exact post-filter network input.'),('Lacruz',[CACHE/'Lacruz_author_v1.0.zip',EXT/'raw/lacruz2026/whole_cohort_data.csv',EXT/'raw/lacruz2026/Lacruz2026_Supplementary_Data3.xlsx'],'Author code requires physeq_adj.rds and phenotype-specific taxon prevalence filtering; repository archive lacks the required RDS. Metadata and final edges do not certify complete inference taxa.')]
 for study,paths,reason in unavailable:
  for p in paths:reg(p,'strict_universe_evidence_audit')
  audit.append(dict(study=study,strict_available='NO',strict_raw_tested_taxa=NA,strict_mapped_unique_GTDB_taxa=NA,strict_unresolved_raw_taxa=NA,published_network_nodes=len({t for e in sets[study] for t in e}),evidence=reason,evidence_files=';'.join(map(str,paths)),limitation='Strict analysis NA; primary remains A-positive / B-node-visible / B-edge-absent'))
 audit.sort(key=lambda r:PAPERS.index(r['study']));write(ROOT/'02_universe_audit/stage5A_universe_audit.tsv',audit)
 for r in audit:info('strict_universe_evidence',r,scope=r['study'],status='WARN' if r['strict_available']=='YES' and int(r['strict_unresolved_raw_taxa']) else 'INFO')
 return strict

def database(sets):
 reg(B4/'scripts/04b_recurrence_source_excluded_recovery.py','frozen_support_logic')
 expected={r['input_file']:r['sha256'] for r in read(B4/'qc/04b_input_manifest.tsv')};reg(DB,'original_frozen_f2_database')
 check('official_database_hash',sha(DB)==expected[str(DB)])
 decoderrows=read(EXT/'stage3a/qc/03a_database_label_decoding.tsv');decoder={r['database_label']:r['gtdb_r226_species'] for r in decoderrows}
 check('exact_frozen_decoder',len(decoder)==len(decoderrows) and all(r['decoding_status']=='exact_unique_R226_serialization' for r in decoderrows))
 universe={t for es in sets.values() for e in es for t in e};seen=set();records={};support=Counter();alias={r['reported_accession']:r['canonical_project_accession'] for r in read(B3/'provenance/03b_accession_alias_map.tsv')};bad=0
 with DB.open(newline='') as f:
  reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames
  check('database_support_fields',fields==['species_a','species_b','aver_corr','aver_pv','proj_count','proj_ID','proj_corr','proj_pv','proj_sam','super_env','sub_env'])
  for line,r in enumerate(reader,2):
   vs=[tokens(r[k]) for k in ['proj_ID','proj_corr','proj_pv','proj_sam']];ids=vs[0];count=int(r['proj_count'])
   if not (count>=2 and all(len(v)==count for v in vs) and len(set(ids))==count and all(alias.get(x,x)==x for x in ids)):bad+=1
   if not all(math.isfinite(float(x)) for v in vs[1:] for x in v):bad+=1
   support.update(ids);a,b=sorted([decoder[r['species_a']],decoder[r['species_b']]])
   if a==b or (a,b) in seen:raise RuntimeError('Duplicate/self-loop database row '+str(line))
   seen.add((a,b))
   if a in universe and b in universe:records[(a,b)]=dict(taxon1=a,taxon2=b,database_row_number=line,**r)
 check('database_global_scan_complete',len(seen)==653071 and len(support)==859 and not bad,dict(rows=len(seen),projects=len(support),malformed=bad,study_node_pairs_retained=len(records)))
 prioractual={r['project']:int(r['supporting_database_association_count']) for r in read(B4/'qc/04b_actual_project_support_universe.tsv')}
 check('actual_support_universe_matches_frozen_Stage4B',dict(support)==prioractual)
 configs=read(B4/'source_sets/04b_network_source_config.tsv');source={p:set() for p in PAPERS};reported={p:set() for p in PAPERS}
 for r in configs:
  assert r['status']=='PASS'
  if r['network']=='Lacruz_MHO':continue
  source[r['paper']].update(tokens(r['overlapping_source_projects']));reported[r['paper']].update(tokens(r['reported_source_projects']))
 for p in PAPERS:check('source_project_intersection',source[p]==reported[p]&set(support),source[p],p)
 check('MHO_adds_no_overlapping_source',all(not r['overlapping_source_projects'] for r in configs if r['network']=='Lacruz_MHO'))
 write(ROOT/'00_input_audit/stage5A_study_source_sets.tsv',[dict(study=p,reported_source_projects=packed(reported[p]),actual_overlapping_sources=packed(source[p]),source_evidence=str(B4/'source_sets/04b_network_source_config.tsv'),rule='paper-level union across primary subnetworks; exclude both A and B even when B did not report edge') for p in PAPERS])
 write(ROOT/'00_input_audit/stage5A_mackdb_project_support_cache.tsv',records.values())
 # Reconcile all available Stage4B baseline vectors, without using its differently conditioned post-exclusion result.
 prior=read(B4/'pair_level/04b_sensitivity_pair_recovery.tsv');mismatches=0
 for r in prior:
  q=records.get(key(r));present=q is not None
  if present!=(r['baseline_mackdb_present']=='Yes'):mismatches+=1;continue
  if present and any(q[k]!=r['baseline_'+k] for k in ['proj_ID','proj_corr','proj_pv','proj_sam']):mismatches+=1
 check('frozen_Stage4B_baseline_support_vectors_identical',mismatches==0,len(prior))
 read(B4/'source_sets/04b_primary_pair_source_sets.tsv');read(B4/'source_sets/04b_sensitivity_pair_source_sets.tsv')
 (ROOT/'00_input_audit/stage5A_support_lookup_policy.txt').write_text('The frozen f>=2 database was fully scanned (653071 unique pairs; 859 projects), and baseline vectors reconciled against Stage 4B. An absent pair is not present in the audited f>=2 database, so cannot satisfy f>=2 after source removal. Its exact raw/sub-f2 support count is unknown, not assumed zero: baseline_support_count and remaining_support_count are NA. Recovery=0 follows audited f>=2 absence, not imputation of support count. Median/IQR are calculated only for known counts, with coverage counts reported.\n')
 return records,source

def pairrows(A,B,es,EB,visible,typ,records,sources,analysis):
 result=[];sa,sb=sources[A],sources[B];remove=sa|sb
 for a,b in sorted(es):
  r=records.get((a,b));vs=[tokens(r[k]) for k in ['proj_ID','proj_corr','proj_pv','proj_sam']] if r else [[],[],[],[]];ids=vs[0]
  ia=[i for i,p in enumerate(ids) if p in sa];ib=[i for i,p in enumerate(ids) if p in sb];removed=sorted(set(ia)|set(ib));keep=[i for i,p in enumerate(ids) if p not in remove];retained=[[v[i] for i in keep] for v in vs];deleted=[[v[i] for i in removed] for v in vs]
  assert list(zip(*retained))==[t for t in zip(*vs) if t[0] not in remove] and not set(retained[0])&remove
  recovered=bool(r) and len(keep)>=2
  out=dict(study_A=A,study_B=B,analysis=analysis,taxon1=a,taxon2=b,in_A=True,in_B=(a,b) in EB,B_universe_type=typ,B_endpoint1_visible=a in visible,B_endpoint2_visible=b in visible,baseline_mackdb_present=bool(r),baseline_support_count=len(ids) if r else NA,support_lookup_status='PRESENT_AUDITED_F2_ROW' if r else 'ABSENT_FROM_AUDITED_F2_DATABASE_SUBF2_COUNT_UNKNOWN',A_source_projects_to_exclude=packed(sa),B_source_projects_to_exclude=packed(sb),AB_source_projects_to_exclude=packed(remove),A_source_projects_removed=';'.join(ids[i] for i in ia),B_source_projects_removed=';'.join(ids[i] for i in ib),n_A_supports_removed=len(ia) if r else NA,n_B_supports_removed=len(ib) if r else NA,n_AB_unique_supports_removed=len(removed) if r else NA,remaining_support_count=len(keep) if r else NA,recovered_after_AB_exclusion=recovered,remaining_support_projects=';'.join(retained[0]),baseline_support_projects=';'.join(ids),baseline_correlations=';'.join(vs[1]),baseline_pvalues=';'.join(vs[2]),baseline_sample_counts=';'.join(vs[3]),remaining_correlations=';'.join(retained[1]),remaining_pvalues=';'.join(retained[2]),remaining_sample_counts=';'.join(retained[3]),removed_support_projects=';'.join(deleted[0]),removed_correlations=';'.join(deleted[1]),removed_pvalues=';'.join(deleted[2]),removed_sample_counts=';'.join(deleted[3]),removed_indices_zero_based=';'.join(map(str,removed)),database_row_number=r['database_row_number'] if r else NA,source_alignment_QC='PASS')
  result.append(out)
 return result
PAIR_FIELDS=list(pairrows('A','B',[('a','b')],set(),{'a','b'},NODE,{}, {'A':set(),'B':set()},'primary')[0])
def summary(A,B,rows,typ,analysis):
 known=[int(r['remaining_support_count']) for r in rows if r['remaining_support_count']!=NA];before=sum(r['baseline_mackdb_present'] for r in rows);after=sum(r['recovered_after_AB_exclusion'] for r in rows);q=np.quantile(known,[.25,.5,.75]) if known else [NA]*3
 return dict(study_A=A,study_B=B,analysis=analysis,universe_type=typ,directional_pairs=len(rows),baseline_recovered=before,baseline_recovery_rate=ratio(before,len(rows)),AB_source_excluded_recovered=after,AB_source_excluded_recovery_rate=ratio(after,len(rows)),pairs_lost_after_AB_exclusion=before-after,pairs_with_any_source_support_removed=sum(r['n_AB_unique_supports_removed']!=NA and r['n_AB_unique_supports_removed']>0 for r in rows),total_project_support_entries_removed=sum(r['n_AB_unique_supports_removed'] for r in rows if r['n_AB_unique_supports_removed']!=NA),median_remaining_support=q[1],IQR_remaining_support=NA if not known else f'{q[0]}–{q[2]}',remaining_support_q25=q[0],remaining_support_q75=q[2],remaining_support_count_known_pairs=len(known),remaining_support_count_unknown_pairs=len(rows)-len(known),support_distribution_scope='known recorded counts only; no support-count zero imputation',classification='A-positive / B-node-visible / B-edge-absent' if typ==NODE else 'A-positive / B-testable / B-edge-absent' if typ==STRICT else 'shared published reference')

def main():
 (ROOT/'logs/state.json').write_text(json.dumps(dict(status='RUNNING',step='input_and_directional_audit',timestamp=now()))+'\n')
 sets,vocab=study_sets();strict=strict_universes(sets,vocab);records,sources=database(sets)
 counts=[];summaries=[];sens=[];references=[];nodes={p:{t for e in es for t in e} for p,es in sets.items()}
 for A in PAPERS:
  for B in PAPERS:
   if A==B:continue
   shared=sets[A]&sets[B];only=sets[A]-sets[B];d={e for e in only if set(e)<=nodes[B]};ds={e for e in only if set(e)<=strict[B]} if B in strict else None
   counts.append(dict(study_A=A,study_B=B,A_total_edges=len(sets[A]),B_total_edges=len(sets[B]),A_B_shared_edges=len(shared),A_only_edges=len(only),A_only_B_node_visible_edges=len(d),A_only_B_strict_testable_edges=len(ds) if ds is not None else NA))
   check('directional_set_definition',len(sets[A])==len(shared)+len(only) and not d&sets[B] and all(set(e)<=nodes[B] for e in d),scope=A+'__'+B)
   for typ,es,vis,label in [(NODE,d,nodes[B],'node_visible')]+([(STRICT,ds,strict[B],'strict_testable')] if ds is not None else []):
    pairs=[dict(study_A=A,study_B=B,taxon1=x,taxon2=y,universe_type=typ) for x,y in sorted(es)]
    write(ROOT/'03_directional_pairs'/f'{A}_to_{B}_{label}_pairs.tsv',pairs,['study_A','study_B','taxon1','taxon2','universe_type'])
    pr=pairrows(A,B,es,sets[B],vis,typ,records,sources,'primary' if typ==NODE else 'strict')
    fn=f'{A}_to_{B}_pair_level.tsv' if typ==NODE else f'{A}_to_{B}_strict_pair_level.tsv';write(ROOT/'04_source_excluded_recovery'/fn,pr,PAIR_FIELDS)
    summaries.append(summary(A,B,pr,typ,'primary' if typ==NODE else 'strict'))
    check('aligned_AB_source_removal_and_f2_recovery',all(r['source_alignment_QC']=='PASS' and (not r['recovered_after_AB_exclusion'] or r['remaining_support_count']>=2) for r in pr),f'{len(pr)} rows; exclusion={packed(sources[A]|sources[B])}',A+'__'+B+'_'+label)
   if A=='Lacruz' or B=='Lacruz':
    ea=sets['Lacruz_sensitivity'] if A=='Lacruz' else sets[A];eb=sets['Lacruz_sensitivity'] if B=='Lacruz' else sets[B];vis={x for e in eb for x in e};ed={e for e in ea-eb if set(e)<=vis};pr=pairrows(A,B,ed,eb,vis,NODE,records,sources,'Lacruz_MHO_sensitivity');write(ROOT/'04_source_excluded_recovery'/f'{A}_to_{B}_MHO_sensitivity_pair_level.tsv',pr,PAIR_FIELDS)
    sr=summary(A,B,pr,NODE,'Lacruz_MHO_sensitivity');old=next(r for r in summaries if r['study_A']==A and r['study_B']==B and r['analysis']=='primary')
    sens.append(dict(study_A=A,study_B=B,primary_directional_pair_count=old['directional_pairs'],sensitivity_directional_pair_count=sr['directional_pairs'],primary_recovered=old['AB_source_excluded_recovered'],sensitivity_recovered=sr['AB_source_excluded_recovered'],primary_recovery=old['AB_source_excluded_recovery_rate'],sensitivity_recovery=sr['AB_source_excluded_recovery_rate'],recovery_rate_difference=sr['AB_source_excluded_recovery_rate']-old['AB_source_excluded_recovery_rate'] if old['directional_pairs'] and sr['directional_pairs'] else NA,source_set_changed=False))
   if PAPERS.index(A)<PAPERS.index(B):
    pr=pairrows(A,B,shared,sets[B],nodes[B],'SHARED_PUBLISHED_PAIRS',records,sources,'shared_reference');write(ROOT/'04_source_excluded_recovery'/f'{A}_and_{B}_shared_pair_level.tsv',pr,PAIR_FIELDS);references.append(summary(A,B,pr,'SHARED_PUBLISHED_PAIRS','shared_reference'))
 write(ROOT/'03_directional_pairs/stage5A_directional_pair_counts.tsv',counts);write(ROOT/'04_source_excluded_recovery/stage5A_directional_recovery_summary.tsv',[x for x in summaries if x['analysis']=='primary']);write(ROOT/'04_source_excluded_recovery/stage5A_strict_directional_recovery_summary.tsv',[x for x in summaries if x['analysis']=='strict']);write(ROOT/'07_statistics/stage5A_lacruz_MHO_sensitivity.tsv',sens);write(ROOT/'07_statistics/stage5A_shared_pair_reference.tsv',references)
 primary=[r for r in summaries if r['analysis']=='primary'];matrix=[]
 for A in PAPERS:
  row={'study_A':A}
  for B in PAPERS:
   q=next((r for r in primary if r['study_A']==A and r['study_B']==B),None)
   row[B]=NA if not q else f"{q['AB_source_excluded_recovered']}/{q['directional_pairs']} ({100*q['AB_source_excluded_recovery_rate']:.1f}%)" if q['directional_pairs'] else '0/0 (NA)'
  matrix.append(row)
 write(ROOT/'07_statistics/stage5A_primary_directional_matrix.tsv',matrix)
 # Compact source-excluded membership matrices for both null models; absent f2 entries are non-recovered.
 universe=sorted({t for es in sets.values() for e in es for t in e});index={t:i for i,t in enumerate(universe)};n=len(universe)
 for ia,A in enumerate(PAPERS):
  for B in PAPERS[ia+1:]:
   recovered=np.zeros((n,n),dtype=np.bool_);remove=sources[A]|sources[B]
   for (a,b),r in records.items():
    if sum(p not in remove for p in tokens(r['proj_ID']))>=2:recovered[index[a],index[b]]=recovered[index[b],index[a]]=True
   np.savez_compressed(ROOT/'00_input_audit'/f'{A}_{B}_source_excluded_lookup.npz',taxa=np.asarray(universe),recovered=recovered,excluded_projects=np.asarray(sorted(remove)))
 write(ROOT/'00_input_audit/stage5A_input_manifest.tsv',INPUT.values())
 (ROOT/'logs/state.json').write_text(json.dumps(dict(status='RUNNING',step='directional_recovery_complete_null_pending',primary_comparisons=30,strict_comparisons=sum(r['analysis']=='strict' for r in summaries),timestamp=now()))+'\n')
 for r in primary:print(r['study_A'],r['study_B'],r['directional_pairs'],r['AB_source_excluded_recovered'],r['AB_source_excluded_recovery_rate'],flush=True)
 print('DIRECTIONAL_RECOVERY_COMPLETE',flush=True)
if __name__=='__main__':
 try:main()
 except Exception as ex:
  (ROOT/'logs/failure.json').write_text(json.dumps(dict(error=repr(ex),traceback=traceback.format_exc(),timestamp=now()),indent=2)+'\n')
  if INPUT:write(ROOT/'00_input_audit/stage5A_input_manifest.tsv',INPUT.values())
  raise
