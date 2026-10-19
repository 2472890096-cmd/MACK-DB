"""Build provenance tables only. Never fit networks or modify frozen inputs."""
from rescue_io import *
import collections, copy, re, xml.etree.ElementTree as ET
import openpyxl
from bs4 import BeautifulSoup

DOI={'Chen':'10.1038/s41467-020-17840-y','Zhou':'10.1016/j.cell.2022.08.021','Guo':'10.1016/j.chom.2023.01.004','Lacruz':'10.1038/s41467-026-72588-1','Li':'10.1002/advs.202308313'}
REPO={'Chen':'https://github.com/GRONINGEN-MICROBIOME-CENTRE/Groningen-Microbiome/tree/master/Projects/Microbial%20co-abundance%20network','Zhou':'https://github.com/BaranziniLab/iMSMS_study','Guo':'No study-specific sample manifest repository identified; author-deposited BioSample records used','Lacruz':'https://zenodo.org/records/19484218','Li':'https://github.com/junhuili/IBS_co-occurrence_network'}
METHOD={'Chen':'SparCC AND SPIEC-EASI v1.0.6 glasso; four separate cohort fits','Zhou':'SparCC (R SpiecEasi); corpcor adjustment for age/sex/BMI','Guo':'FlashWeave-S; covariate-adjusted conditional independence','Lacruz':'SPIEC-EASI v1.1.3; Meinshausen-Buhlmann neighborhood selection','Li':'SparCC and bootstrap (R SpiecEasi), cross-cohort discovery analysis'}
MF='study network publication_sample_id participant_id group timepoint biosample run_accession bioproject study_accession evidence_level evidence_source evidence_url_or_local_file found_in_local_abundance local_abundance_file notes status exact_network_assignment input_role column_index library_strategy public_identity_evidence_level public_sample_alias collection_timestamp wet_dry source_cohort publication_row sample_accession_namespace'.split()
SF='study network doi published_method published_sample_count published_participant_count source_projects_reported source_projects_discovered exact_samples_recovered exact_runs_recovered samples_found_in_local_abundance samples_missing_from_local_abundance recovery_fraction sample_set_status classification main_evidence_source repository_url notes audit_scope sample_count_unit exact_biosamples_recovered local_runs_recovered participants_with_any_local_run local_complete_fraction next_stage_ready can_rerun_fastspar candidate_public_runs candidate_local_runs'.split()
META={r['run_accession']:r for r in rows(E/'public_run_sample_crosswalk.tsv')}
LOCAL={r['run_accession']:r for r in rows(E/'verified_local_run_columns.tsv')}
assert len(LOCAL)==len(rows(E/'verified_local_run_columns.tsv'))
MASTER=[]; SUMMARY=[]; QC=[]; COVERAGE=[]; EXTRA=[]

def check(scope,name,ok,details,negative='FAIL'):
    QC.append(dict(scope=scope,check=name,status='PASS' if ok else negative,details=str(details)))
def info(scope,name,details,status='INFO'):
    QC.append(dict(scope=scope,check=name,status=status,details=str(details)))
def distinct(rs,key):return {r[key] for r in rs if r[key]}
def base(study,network,run='',**kw):
    m=META.get(run,{}); loc=LOCAL.get(run,{})
    r=dict.fromkeys(MF,'')
    r.update(study=study,network=network,run_accession=run,biosample=m.get('biosample_accession',''),
             bioproject=m.get('source_bioproject',''),study_accession=m.get('study_accession',''),
             library_strategy=m.get('library_strategy',''),public_sample_alias=m.get('sample_name',''),
             evidence_level='U',public_identity_evidence_level='E4' if m else 'U',
             status='unresolved',exact_network_assignment='NO',input_role='unresolved_project_candidate',
             found_in_local_abundance='YES' if loc else 'NO',local_abundance_file=loc.get('abundance_file',''),
             column_index=loc.get('column_index',''),sample_accession_namespace='BioSample' if m.get('biosample_accession') else '')
    r.update(kw);return r
def addsummary(study,network,n,projects,classification,source,notes,scope='frozen_external_input',**overrides):
    rs=[r for r in MASTER if r['network']==network];exact=[r for r in rs if r['exact_network_assignment']=='YES']
    by=collections.defaultdict(list)
    for r in exact:by[r['participant_id'] or r['publication_sample_id']].append(r)
    complete=sum(all(r['found_in_local_abundance']=='YES' for r in a) for a in by.values())
    anylocal=sum(any(r['found_in_local_abundance']=='YES' for r in a) for a in by.values())
    for pid,a in by.items():
        COVERAGE.append(dict(study=study,network=network,participant_id=pid,group=a[0]['group'],
            required_runs=len(distinct(a,'run_accession')),local_runs=sum(r['found_in_local_abundance']=='YES' for r in a),
            complete_local_input='YES' if all(r['found_in_local_abundance']=='YES' for r in a) else 'NO',
            publication_sample_ids=';'.join(sorted(distinct(a,'publication_sample_id'))),biosamples=';'.join(sorted(distinct(a,'biosample')))))
    sr=dict.fromkeys(SF,'')
    sr.update(study=study,network=network,doi=DOI[study],published_method=METHOD[study],published_sample_count=n,
        published_participant_count=n,source_projects_reported=projects,source_projects_discovered=projects,
        exact_samples_recovered=len(by),exact_runs_recovered=len(distinct(exact,'run_accession')),
        samples_found_in_local_abundance=complete,samples_missing_from_local_abundance=(n-complete if isinstance(n,int) and exact else 'unresolved'),
        recovery_fraction=(f'{len(by)/n:.8f}' if isinstance(n,int) else 'unresolved'),
        sample_set_status=('complete_public_identity_and_local_input' if classification=='PASS_EXACT' else 'complete_public_identity_incomplete_local_input' if exact else 'unresolved_network_membership'),
        classification=classification,main_evidence_source=source,repository_url=REPO[study],notes=notes,audit_scope=scope,
        sample_count_unit='biological_participant; complete means ALL documented input runs locally present',
        exact_biosamples_recovered=len(distinct(exact,'biosample')),local_runs_recovered=sum(r['found_in_local_abundance']=='YES' for r in exact),
        participants_with_any_local_run=anylocal,local_complete_fraction=(f'{complete/n:.8f}' if isinstance(n,int) else 'unresolved'),
        next_stage_ready='YES' if classification in ['PASS_EXACT','MULTI_COHORT_RECONSTRUCTABLE'] else 'NO',
        can_rerun_fastspar='YES_NEXT_STAGE_ONLY' if classification in ['PASS_EXACT','MULTI_COHORT_RECONSTRUCTABLE'] else 'NO',
        candidate_public_runs=len(distinct([r for r in rs if r['exact_network_assignment']!='YES'],'run_accession')),
        candidate_local_runs=sum(r['found_in_local_abundance']=='YES' for r in rs if r['exact_network_assignment']!='YES'))
    sr.update(overrides);SUMMARY.append(sr)
    return sr
def wbrows(p,sheet):
    w=openpyxl.load_workbook(p,read_only=True,data_only=True);v=iter(w[sheet].values);head=next(v)
    data=[dict(zip(head,row)) for row in v];w.close();return data

def zhou():
    p=C/'Zhou_publisher_Table_S1.xlsx'; seq=wbrows(p,'Sequenced samples');ph=wbrows(p,'Sample phenotype')
    url='https://ars.els-cdn.com/content/image/1-s2.0-S0092867422011151-mmc2.xlsx'
    subjects={r['iMSMS_ID']:r for r in ph};assert len(subjects)==1152
    groups={'MS':{k for k,r in subjects.items() if r['disease']=='MS' and r['treatment_status']=='Untreated'},
            'HHC':{k for k,r in subjects.items() if r['disease']=='Control' and r['treatment_control']=='Untreated_HHC'}}
    assert not groups['MS']&groups['HHC']
    for g,ids in groups.items():
        network='Zhou_'+g; chosen=[(i,r) for i,r in enumerate(seq,2) if r['iMSMS_ID'] in ids and r['target_gene']=='WGS' and r['paired_used']=='Yes']
        assert len(ids)==209 and {r['iMSMS_ID'] for _,r in chosen}==ids
        for i,r in chosen:
            run=r['run_accession'];m=META[run]
            assert r['sample_accession']==m['biosample_accession'] and r['experiment_accession']==m['experiment_accession']
            assert r['study_accession']==m['source_bioproject']=='PRJEB32762'
            MASTER.append(base('Zhou',network,run,publication_sample_id=r['sample_id'],participant_id=r['iMSMS_ID'],group=g,
                timepoint=str(r['collection_timestamp']),collection_timestamp=str(r['collection_timestamp']),wet_dry=str(r['wet_dry']),
                source_cohort=str(r['Cohort']),publication_row=f'Sequenced samples:{i}',evidence_level='E1',evidence_source='Publisher Table S1, Sequenced samples AND Sample phenotype; article WGS participant count summation',
                evidence_url_or_local_file=f'{p};{url};{PROV}/raw/Zhou_article.html',status='partial',exact_network_assignment='YES',input_role='published_untreated_network_participant_contributing_run',
                notes='Exact participant join, WGS AND paired_used=Yes; counts from repeated collections/preparations/technical sequencing are summed per participant, not treated as independent subjects. '+('Publication WGS conflicts with public library_strategy=OTHER; retained by explicit publisher designation. ' if m['library_strategy']=='OTHER' else '')+'Archived executed species.treatment.rda not available; this is the publication-defined input, not a recovered fitted object.'))
        rs=[r for r in MASTER if r['network']==network]
        addsummary('Zhou',network,209,'PRJEB32762 (ERP115476)','PARTIAL',str(p),
            'Complete 209 published untreated-group identities; run and BioSample mapping uniquely verified. Local participant count requires every contributing run, not merely any run. Final fitted group R object not deposited in inspected tree.')
        check(network,'publisher_run_sample_experiment_public_crosswalk',True,f'{len(rs)} unique linked runs; publisher WGS designation retained')
        info(network,'publication_vs_public_library_strategy',dict(collections.Counter(r['library_strategy'] for r in rs)))
    chosenruns={r['run_accession'] for r in MASTER if r['study']=='Zhou'}
    lookup=collections.defaultdict(list)
    for r in seq:lookup[r['run_accession']].append(r)
    for run,l in LOCAL.items():
        if l['source_bioproject']!='PRJEB32762' or run in chosenruns:continue
        matched=lookup.get(run,[])
        reason='Publisher record outside target untreated WGS paired-used participant selection' if matched else 'Not located in publisher Table S1; excluded rather than inferred'
        EXTRA.append(dict(study='Zhou',run_accession=run,project=l['source_bioproject'],reason=reason,evidence=str(p),status='EXCLUDED',details=json.dumps([{k:x[k] for k in ['iMSMS_ID','target_gene','paired_used']} for x in matched])))

def guo():
    p=C/'Guo_direct_BioSample_attributes.xml'; nodes=ET.parse(p).findall('.//BioSample')
    b={x.get('accession'):{a.get('attribute_name'):a.text or '' for a in x.findall('./Attributes/Attribute')} for x in nodes}
    assert len(b)==len(nodes)==197
    mr=[r for r in META.values() if r['source_bioproject']=='PRJNA751448'];assert len(mr)==197
    assert {r['biosample_accession'] for r in mr}==set(b)
    gc=[];attributes=[]
    for m in mr:
        a=b[m['biosample_accession']]; old=json.loads(m['available_sample_attributes'])
        assert a['SampleType'] in ['Control','Case'] and a['mecfs_status']=={'Control':'N','Case':'Y'}[a['SampleType']]
        assert a['Original_ID']==a['participant_id'] and a['STOOL_ID']
        for k in ['SampleType','mecfs_status','participant_id','STOOL_ID']:assert a[k]==old[k]
        group={'Control':'Control','Case':'MECFS'}[a['SampleType']]
        r=base('Guo','Guo_'+group,m['run_accession'],publication_sample_id=a['STOOL_ID'],participant_id=a['participant_id'],group=group,
            timepoint='first_study_visit (article-defined)',collection_timestamp=a.get('collection_date',''),
            evidence_level='E5',public_identity_evidence_level='E4',status='exact',exact_network_assignment='YES',input_role='published_single_fit_group',
            evidence_source='Author-submitted BioSample SampleType + mecfs_status + participant_id/Original_ID + STOOL_ID; paper first-visit group definition; Table S7',
            evidence_url_or_local_file=f'{p};https://www.ncbi.nlm.nih.gov/biosample/{m["biosample_accession"]};{PROV}/raw/Guo_fulltext.xml;{OLD}/raw/guo2023/1-s2.0-S193131282300029X-mmc8.xlsx',
            notes='Unique evidence chain: paper deposits all analyzed first-visit stool metagenomes in PRJNA751448 and fits all Control/Case samples; author BioSample explicitly labels SampleType, corroborated by mecfs_status; participant_id=Original_ID and unique STOOL_ID -> BioSample -> unique SRA/ENA run -> verified original-profile local column. Not inferred from disease text or matching counts.')
        gc.append(r);attributes.append(dict(run_accession=m['run_accession'],biosample=m['biosample_accession'],SampleType=a['SampleType'],mecfs_status=a['mecfs_status'],participant_id=a['participant_id'],Original_ID=a['Original_ID'],STOOL_ID=a['STOOL_ID'],attributes=json.dumps(a,sort_keys=True),evidence=str(p)))
    for key in ['participant_id','publication_sample_id','biosample','run_accession']:assert len(distinct(gc,key))==197
    assert all(r['found_in_local_abundance']=='YES' for r in gc)
    MASTER.extend(gc)
    for r in gc:MASTER.append(dict(r,network='Guo_Common',input_role='published_all197_single_fit'))
    for group,n in [('Control',91),('MECFS',106),('Common',197)]:
        assert len([r for r in MASTER if r['network']=='Guo_'+group])==n
        addsummary('Guo','Guo_'+group,n,'PRJNA751448','PASS_EXACT',str(p),
            'Published Table S7 intermediate network, not the frozen Agg network. Full unique first-visit membership and all local run columns recovered. Sample eligibility only; conditional FlashWeave weights/covariates are not equivalent to marginal SparCC correlations.',scope='published_intermediate_component')
    for group in ['Control','MECFS']:
        network='Guo_Agg'+group
        for r in gc:
            MASTER.append(dict(r,network=network,status='partial',input_role='Common197_and_group_specific_parent' if r['group']==group else 'Common197_parent_only',
                notes=r['notes']+' Aggregate edge table combines Common197 with validated group-specific edges. Individual validation subsampling draws were not recovered. This row identifies the parent input membership, not an exact subsampling draw.'))
        addsummary('Guo',network,197,'PRJNA751448','PARTIAL',str(p),
            'All 197 parent identities and local runs complete; no single fitted sample set corresponds to Agg. Published validation subsampling draw identities/seeds not recovered. PARTIAL refers to composite fitting provenance, not missing first-visit samples. Use separately published Table S7 Control/MECFS/Common networks for a single-fit comparison.',
            published_sample_count='Common197 + '+group+('91' if group=='Control' else '106')+' + validation subsampling',
            sample_set_status='complete_parent_memberships;composite_not_single_fit;subsampling_memberships_unresolved',sample_count_unit='197 unique parent participants; overlapping component sets, not additive')
    write(E/'Guo_explicit_group_attribute_audit.tsv',attributes)
    info('Guo','aggregate_vs_single_fit','AggControl/AggMECFS are not single-group fits. Parent membership complete; subsampling provenance unresolved. Conservative PARTIAL retained separately from the three PASS_EXACT intermediate objects.','WARN')
    check('Guo','explicit_author_labels_crossvalidated',True,'197 unique participants/STOOL_ID/BioSamples/runs; SampleType matches mecfs_status and separate SRA snapshot for every sample')

def lacruz():
    fp=PROV/'03b_lacruz_sample_project_evidence.tsv';mp=OLD/'parsed/lacruz2026_sample_source_manifest.tsv'
    register([fp,mp,OLD/'raw/lacruz2026/Lacruz2026_Supplementary_Data3.xlsx',OLD/'raw/lacruz2026/whole_cohort_data.csv'])
    frozen=rows(fp);original=rows(mp);assert len(frozen)==len(original)==931
    original={r['subject_id']:r for r in original};assert len(original)==931
    with (C/'Lacruz_data__whole_cohort_data.csv').open() as f:cols=next(csv.reader(f))[1:]
    assert len(cols)==len(set(cols));assert {r['abundance_sample_id'] for r in frozen}<=set(cols)
    assert sha(C/'Lacruz_data__whole_cohort_data.csv')==sha(OLD/'raw/lacruz2026/whole_cohort_data.csv')
    for r in frozen:
        o=original[r['subject_id']];assert o['run_accession']==r['run_accessions'] and o['phenotype']==r['phenotype']
        visit=re.search(r'_V(\d+)$',r['abundance_sample_id'])
        for run in r['run_accessions'].split(';'):
            m=META[run];assert m['source_bioproject']==r['canonical_project_accession']
            MASTER.append(base('Lacruz',r['network'],run,publication_sample_id=r['abundance_sample_id'],participant_id=r['subject_id'],group=r['phenotype'],
                timepoint='V'+visit.group(1) if visit else 'source_sample_selected_in_published_manifest',source_cohort=r['study_name'],publication_row=r['metadata_excel_row'],
                evidence_level='E5',status='partial',exact_network_assignment='YES',input_role='published_phenotype_network_input',
                evidence_source='Supplementary Data3 phenotype/run manifest + versioned author abundance columns + ENA sample alias/visit crosswalk',
                evidence_url_or_local_file=f'{fp};{mp};{OLD}/raw/lacruz2026/Lacruz2026_Supplementary_Data3.xlsx;{C}/Lacruz_data__whole_cohort_data.csv;{r["ena_evidence"]}',
                notes='Unique evidence chain: publication row specifies subject, phenotype and run(s); versioned author whole_cohort_data.csv contains exact mapped sample column; '+r['identity_rule']+'. Multiple sequencing runs retained as one selected biological subject; AI4Food visit preserved from author selection, never chosen anew.'))
    allprojects=sorted({r['canonical_project_accession'] for r in frozen})
    for g,n in [('MHNO',103),('MHO',19),('MUNO',377),('MUO',432)]:
        network='Lacruz_'+g;rs=[r for r in MASTER if r['network']==network]
        assert len(distinct(rs,'participant_id'))==n
        addsummary('Lacruz',network,n,';'.join(sorted(distinct(rs,'bioproject'))),'MULTI_COHORT_NOT_RECONSTRUCTABLE',str(fp),
            'Complete publication subject/run selection recovered; only PRJEB7774 has local abundance among the five actual source projects. Related PRJEB41311 is a different dataset and supplies no verified target sample.')
    target=[r for r in MASTER if r['study']=='Lacruz'];truns=distinct(target,'run_accession');tbs=distinct(target,'biosample');aliases=distinct(target,'public_sample_alias')|distinct(target,'participant_id')|distinct(target,'publication_sample_id')
    cand=[r for r in META.values() if r['source_bioproject']=='PRJEB41311']; candlocal=[r for r in cand if r['run_accession'] in LOCAL]
    interruns={r['run_accession'] for r in cand}&truns;interbs={r['biosample_accession'] for r in cand}&tbs;interaliases={r['sample_name'] for r in cand}&aliases
    check('Lacruz','related_MetaCardis_project_is_not_exact_primary_input',not interruns and not interbs and not interaliases,
        f'PRJEB41311: {len(cand)} public runs / {len(candlocal)} local columns; exact primary intersections runs={len(interruns)}, BioSamples={len(interbs)}, aliases={len(interaliases)}')
    for run,l in LOCAL.items():
        if l['source_bioproject'] not in ['PRJEB7774','PRJEB41311'] or run in truns:continue
        EXTRA.append(dict(study='Lacruz',run_accession=run,project=l['source_bioproject'],reason='Absent from complete published 931-subject run manifest; not eligible',evidence=str(fp),status='EXCLUDED',details='Related project/cohort naming is insufficient for sample identity'))
    cohort=[]
    for project in allprojects:
        rs=[r for r in target if r['bioproject']==project]
        cohort.append(dict(project=project,cohorts=';'.join(sorted(distinct(rs,'source_cohort'))),subjects=len(distinct(rs,'participant_id')),runs=len(distinct(rs,'run_accession')),local_target_runs=sum(r['found_in_local_abundance']=='YES' for r in rs),evidence=str(fp)))
    write(E/'Lacruz_all_source_cohorts.tsv',cohort)

def chen():
    # Public accession relationships are rescued; final network selection is not inferred from dataset size.
    cohort_defs=[('LLD',1135,'EGAD00001001991','EGAS00001001704'),('IBD',496,'EGAD00001004194','EGAS00001002702'),('300OB',298,'EGAD00001005083','EGAS00001003508')]
    ca=[]
    for name,n,d,study in cohort_defs:
        sp=C/f'Chen_{d}_samples.json';rp=C/f'Chen_{d}_runs.json'
        samples={r['accession_id']:r for r in json.loads(sp.read_text())};runs=json.loads(rp.read_text())
        for run in runs:
            a=samples[run['sample_accession_id']]
            MASTER.append(base('Chen','Chen_'+name,run['accession_id'],participant_id=a.get('subject_id',''),biosample=a['accession_id'],
                bioproject='',study_accession=study,sample_accession_namespace='EGA_sample (EGAN, not an NCBI BioSample)',source_cohort=name,
                public_identity_evidence_level='E4',evidence_source='Public EGA dataset run -> EGAN sample -> subject metadata; final published network inclusion unresolved',
                evidence_url_or_local_file=f'{sp};{rp};https://metadata.ega-archive.org/datasets/{d}/samples',
                input_role='public_cohort_candidate_not_network_selected',notes=f'Public dataset {d}; {len(samples)} public subjects versus {n} published network inputs. Matching counts do not prove final selection. No publication-to-EGA final fitting manifest or local run crosswalk recovered.'))
        ca.append(dict(cohort=name,initial_or_source_count={'LLD':'1500 recruited; EGA description 1200','IBD':'EGA description 461; public visible 355','300OB':'302 recruited; public visible 297'}[name],published_network_n=n,archive_dataset=d,archive_study=study,public_samples=len(samples),public_runs=len(runs),final_selection='unresolved',local_exact_network_runs='unresolved',notes='Open metadata retrieved; sequencing data controlled access; EGA aliases not treated as NCBI runs'))
        addsummary('Chen','Chen_'+name,n,d,'FAIL',str(sp),'Public EGA subject/sample/run crosswalk recovered, but final network membership cannot be certified from count or cohort identity alone.',scope='published_cohort_component',source_projects_discovered=f'{d};{study}',samples_found_in_local_abundance='unresolved',local_complete_fraction='unresolved')
    for m in META.values():
        if m['source_bioproject'] not in ['PRJNA319574','PRJNA942468']:continue
        resequence=m['source_bioproject']=='PRJNA942468'
        MASTER.append(base('Chen','Chen_500FG',m['run_accession'],participant_id=m.get('participant_id_public',''),source_cohort='500FG',
            evidence_source='Public project metadata; author network script supplies cohort N but no sample IDs',
            evidence_url_or_local_file=';'.join(x for x in [m['ncbi_evidence_file'],m['ena_evidence_file'],str(C/'Chen_author_network.R'),str(STAGE/'00_metadata_cache/Chen_fulltext.xml'),str(C/'Chen_2024_ENA_study.xml') if resequence else ''] if x),
            input_role='later_resequencing_not_original_2020_runs' if resequence else 'original_project_candidate_final450_unknown',
            notes='PRJNA942468 explicitly describes resequencing using QIAamp Fast DNA Stool Mini Kit; the 2020 paper specifies AllPrep DNA/RNA Mini Kit and deposits 500FG at PRJNA319574. Thus not certified as original network run input or a replacement. First-public date2023 is supporting context, not proof of sequencing date.' if resequence else 'Original 500FG project; final 450 samples not uniquely enumerated by recovered author materials. Not all BioProject runs are network inputs.'))
    ca.append(dict(cohort='500FG',initial_or_source_count='534 recruited; 875 original-project public runs; 517 later resequencing runs',published_network_n=450,archive_dataset='PRJNA319574',archive_study='PRJNA319574;PRJNA942468 (later resequencing)',public_samples='see public metadata',public_runs=1392,final_selection='unresolved',local_exact_network_runs='unresolved',notes='Local original-project runs466; later resequencing517; neither count certifies final450 membership'))
    addsummary('Chen','Chen_500FG',450,'PRJNA319574','FAIL',str(C/'Chen_author_network.R'),'Final450 selection unavailable; 466 local original-project candidate runs. PRJNA942468 has517 local later-resequenced runs, not original2020 run input.',scope='published_cohort_component',source_projects_discovered='PRJNA319574;PRJNA942468 (later resequencing)',samples_found_in_local_abundance='unresolved',local_complete_fraction='unresolved')
    for r in list(MASTER):
        if r['study']=='Chen':MASTER.append(dict(r,network='Chen',input_role=r['network']+':'+r['input_role']))
    addsummary('Chen','Chen','LLD1135;500FG450;300OB298;IBD496','PRJNA319574;EGAD00001001991;EGAD00001004194;EGAC00001001143','MULTI_COHORT_NOT_RECONSTRUCTABLE',str(C/'Chen_author_network.R'),
        'Frozen Chen pair union spans four independent cohort fits. Public EGA candidates and original/resequenced500FG candidates are preserved as U, not exact network samples.',source_projects_discovered='PRJNA319574;PRJNA942468 (later resequencing);EGAD00001001991;EGAD00001004194;EGAD00001005083;EGAS00001001704;EGAS00001002702;EGAS00001003508',published_participant_count=2379,samples_found_in_local_abundance='unresolved',local_complete_fraction='unresolved')
    write(E/'Chen_cohort_selection_audit.tsv',ca)
    info('Chen','private_archive_public_metadata_rescue','LLD1135 / IBD355 / 300OB297 public sample/run identities retrieved; none automatically promoted to final-network exact membership.','WARN')
    for run,l in LOCAL.items():
        if l['source_bioproject'] in ['PRJNA319574','PRJNA942468']:
            EXTRA.append(dict(study='Chen',run_accession=run,project=l['source_bioproject'],reason='Final450 inclusion unresolved' if l['source_bioproject']=='PRJNA319574' else 'Later resequencing; not original2020 run',evidence=str(C/'Chen_author_network.R'),status='UNRESOLVED' if l['source_bioproject']=='PRJNA319574' else 'EXCLUDED',details='Not certified as a published network input'))

def li():
    p=PROV/'03b_li_table_s1_cohort_audit.tsv';register([p,PROV/'raw/ADVS-11-2308313-s002.xlsx',PROV/'raw/Li_supplementary.zip'])
    co=rows(p);discovery={r['canonical_project_accession'] for r in co if r['role']=='Discovery'}
    assert len(discovery)==13
    audit=rows(E/'local_project_and_alias_audit.tsv');availability={r['project']:r for r in audit}
    details=[]
    for r in co:
        a=availability.get(r['canonical_project_accession'],{})
        details.append(dict(r,primary_discovery_project='YES' if r['role']=='Discovery' else 'NO',local_project_columns=a.get('local_run_columns','0'),local_project_files=a.get('abundance_files',''),exact_network_sample_selection='unresolved; cohort membership is not final run selection'))
    write(E/'Li_discovery_validation_treatment_audit.tsv',details)
    for m in META.values():
        if m['source_bioproject'] not in discovery:continue
        MASTER.append(base('Li','Li_IBS',m['run_accession'],participant_id=m.get('participant_id_public',''),timepoint=m.get('timepoint_public',''),
            evidence_source='Published Table S1 discovery project plus public SRA/ENA accession crosswalk; final IBS network sample selection unresolved',
            evidence_url_or_local_file=';'.join(x for x in [str(p),m['ncbi_evidence_file'],m['ena_evidence_file']] if x),
            notes='U denotes network membership, not accession identity. Candidate project may contain case/control, treatment, repeated or different assays. Do not use all project runs; no unique paper fitting manifest/QC selection recovered. Public pooled-sample experiment aliases resolved using exact SRA RUN/Pool/Member where needed.'))
    # GSA public project page identifies BioSamples; archive pages are challenge responses, not a run manifest.
    gp=C/'Li_GSA_project_PRJCA002555.html'; soup=BeautifulSoup(gp.read_text(),'html.parser');seen=set()
    for a in soup.find_all('a',href=True):
        match=re.search(r'SAMC\d+',a['href'])
        if not match or match.group() in seen:continue
        bs=match.group();seen.add(bs);tr=a.find_parent('tr');txt=tr.get_text(' ',strip=True) if tr else a.get_text(' ',strip=True)
        MASTER.append(base('Li','Li_IBS',biosample=bs,bioproject='PRJCA002555',public_identity_evidence_level='E4',sample_accession_namespace='GSA_BioSample',
            evidence_source='Canonical public NGDC BioProject BioSample list; no final publication selection or run list',evidence_url_or_local_file=str(gp),
            input_role='public_GSA_BioSample_candidate',notes='Public project record: '+txt+'; exact archive links CRA002557 and CRA002559. Archive run pages blocked by security challenge. Network membership U.'))
    assert len(seen)==154
    addsummary('Li','Li_IBS','unresolved',';'.join(sorted(discovery)),'MULTI_COHORT_NOT_RECONSTRUCTABLE',str(p),
        'Discovery:14 assay cohorts/13 projects. PRJEB42304 contributes reported139 WGS and141 amplicon subjects; overlap between assays not assumed independent. Exact network N and final run/sample selection unavailable. Local139 PRJEB42304 WGS columns are candidates. Newly found PRJEB11419 local467 columns belong to validation, not primary discovery.',
        source_projects_discovered='Discovery:'+','.join(sorted(discovery))+';Validation:PRJEB11419,PRJEB37924,PRJEB19857,PRJNA812699;Treatment-only:PRJNA279031,PRJNA391915;NGDC:CRA002557,CRA002559',samples_found_in_local_abundance='unresolved',local_complete_fraction='unresolved')
    for run,l in LOCAL.items():
        if l['source_bioproject'] not in ['PRJEB42304','PRJEB11419']:continue
        EXTRA.append(dict(study='Li',run_accession=run,project=l['source_bioproject'],reason='Discovery-project run; final network membership unresolved' if l['source_bioproject']=='PRJEB42304' else 'Validation cohort; excluded from discovery network',evidence=str(p),status='UNRESOLVED' if l['source_bioproject']=='PRJEB42304' else 'EXCLUDED',details='Local availability is not evidence of published inclusion'))

def qc_and_save():
    for sr in SUMMARY:
        network=sr['network'];rs=[r for r in MASTER if r['network']==network];exact=[r for r in rs if r['exact_network_assignment']=='YES']
        rr=[r['run_accession'] for r in rs if r['run_accession']]
        check(network,'duplicate_run_within_network',len(rr)==len(set(rr)),f'rows with run={len(rr)}; unique runs={len(set(rr))}')
        for field,name in [('run_accession','run_maps_to_multiple_BioSamples'),('biosample','BioSample_ambiguous_participant'),('participant_id','participant_conflicting_groups'),('publication_sample_id','sample_conflicting_groups')]:
            mapping=collections.defaultdict(set)
            for r in exact:
                if r[field]:mapping[r[field]].add(r['biosample'] if field=='run_accession' else r['participant_id'] if field=='biosample' else r['group'])
            bad={k:sorted(v) for k,v in mapping.items() if len(v)>1}
            check(network,name,not bad,json.dumps(bad) if bad else f'{len(mapping)} unambiguous identities; multiple technical runs allowed')
        if exact:
            check(network,'reported_participant_count_vs_recovered',len(distinct(exact,'participant_id'))==int(sr['published_participant_count']),f'published={sr["published_participant_count"]}; recovered={len(distinct(exact,"participant_id"))}')
        else:info(network,'reported_participant_count_vs_recovered',f'Published={sr["published_participant_count"]}; final membership unresolved; public candidates are not exact','WARN')
        check(network,'exact_samples_vs_local_complete_input',sr['next_stage_ready']!='YES' or sr['samples_missing_from_local_abundance']==0,
            f'Complete participants={sr["samples_found_in_local_abundance"]}; any local run={sr["participants_with_any_local_run"]}; exact local runs={sr["local_runs_recovered"]}')
        missing=[r['run_accession'] for r in exact if r['found_in_local_abundance']=='NO']
        info(network,'exact_sample_runs_not_found_locally',f'{len(missing)} missing; see individual sample list',status='WARN' if missing else 'INFO')
        check(network,'provenance_evidence_completeness',all(r['evidence_source'] and r['evidence_url_or_local_file'] and r['evidence_level']!='U' and (r['evidence_level']!='E5' or 'chain:' in r['notes']) for r in exact),f'{len(exact)} exact assignment rows; E5 chains explicit')
        check(network,'unresolved_candidates_excluded_from_exact_sets',all(r['evidence_level']!='U' for r in exact),f'U rows={len(rs)-len(exact)} retained only in master/candidate lists')
        # Files named exact_samples contain either verified membership rows or one explicit unresolved sentinel, never U candidate runs.
        if exact:out=exact
        else:out=[base(sr['study'],network,evidence_source=sr['main_evidence_source'],evidence_url_or_local_file=sr['main_evidence_source'],input_role='unresolved_no_certified_exact_samples',notes='No certified exact network sample rows. Public identities are in master/candidate list; do not use them for inference.')]
        write(S/(network+'_exact_samples.tsv'),out,MF)
        candidates=[r for r in rs if r['exact_network_assignment']!='YES']
        if candidates:write(S/(network+'_unresolved_candidates.tsv'),candidates,MF)
    for study,networks in [('Zhou',['Zhou_MS','Zhou_HHC']),('Guo',['Guo_Control','Guo_MECFS']),('Lacruz',['Lacruz_MHNO','Lacruz_MHO','Lacruz_MUNO','Lacruz_MUO'])]:
        for field in ['participant_id','biosample','run_accession']:
            sets=[distinct([r for r in MASTER if r['network']==n],field) for n in networks]
            intersections=sum(len(sets[i]&sets[j]) for i in range(len(sets)) for j in range(i+1,len(sets)))
            check(study,'subgroup_mutual_exclusivity_'+field,intersections==0,f'{intersections} conflicting identities across disjoint groups')
    cols=rows(E/'verified_local_run_columns.tsv')
    check('local','duplicate_abundance_columns',len({(r['abundance_file'],r['column_index']) for r in cols})==len(cols),f'{len(cols)} indexed columns across {len({r["abundance_file"] for r in cols})} matrices')
    check('local','exact_run_to_local_column_identity',all(r['profile_values_match_every_retained_taxon']=='Yes' for r in cols),'Frozen verified profile proof reused for seven projects; two additional local projects fully compared to original profiles')
    for study in ['Zhou','Guo','Lacruz','Chen','Li']:
        counts=collections.Counter(r['status'] for r in EXTRA if r['study']==study)
        info(study,'extra_local_run_disposition',json.dumps(counts)+'; every extra row explicitly excluded or unresolved; see local_extra_run_disposition.tsv',status='WARN' if counts['UNRESOLVED'] else 'INFO')
    info('metadata','pooled_SRA_experiment_resolution',str(sum(r['ncbi_ena_accession_consistency']=='NCBI_run_member_resolves_ENA_pool' for r in META.values()))+' runs resolved using exact RUN/Pool/Member; not by taking first BioSample in an experiment package')
    for r in rows(Q/'public_metadata_qc.tsv'):QC.append(dict(scope='public_metadata',**r))
    MASTER.sort(key=lambda r:(r['study'],r['network'],r['participant_id'],r['publication_sample_id'],r['run_accession']))
    SUMMARY.sort(key=lambda r:(r['study'],r['network']))
    write(R/'stage5B_05_sample_master.tsv',MASTER,MF)
    write(R/'stage5B_05_sample_rescue_summary.tsv',SUMMARY,SF)
    write(E/'participant_complete_input_coverage.tsv',COVERAGE)
    write(E/'local_extra_run_disposition.tsv',EXTRA)
    write(Q/'stage5B_05_sample_rescue_qc.tsv',QC)
    event('build_rescue.py: join publication-defined memberships, public identities and verified local columns; no inference',E,R/'stage5B_05_sample_master.tsv',int(any(r['status']=='FAIL' for r in QC)))
    print('MASTER_ROWS',len(MASTER),'SUMMARY_ROWS',len(SUMMARY),'QC',dict(collections.Counter(r['status'] for r in QC)))
    for r in SUMMARY:print(r['network'],r['classification'],'exact participants',r['exact_samples_recovered'],'runs',r['exact_runs_recovered'],'complete local participants',r['samples_found_in_local_abundance'],'local runs',r['local_runs_recovered'])
    if any(r['status']=='FAIL' for r in QC):raise RuntimeError('QC failures must be resolved')

if __name__=='__main__':
    zhou();guo();lacruz();chen();li();qc_and_save()
