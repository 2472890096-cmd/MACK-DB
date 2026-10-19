from rescue_io import *
import xml.etree.ElementTree as ET
import collections
def main():
 p=STAGE/'00_audit/stage5B_public_run_metadata.tsv';register([p]);data={r['run_accession']:r for r in rows(p)};issues=[]
 for p in sorted(C.glob('*_ncbi_sra_*.xml')):
  project=p.name.split('_ncbi')[0]
  try:
   root=ET.parse(p)
   for package in root.findall('.//EXPERIMENT_PACKAGE'):
    sample=package.find('SAMPLE');ex=package.find('EXPERIMENT');study=package.find('STUDY')
    if sample is None or ex is None:continue
    attrs={a.findtext('TAG'):a.findtext('VALUE','') for a in sample.findall('.//SAMPLE_ATTRIBUTE')}
    bio=sample.find(".//EXTERNAL_ID[@namespace='BioSample']");bio=bio.text if bio is not None else ''
    for run in package.findall('.//RUN_SET/RUN'):
     samples={s.get('accession'):s for s in package.findall('SAMPLE')}
     members={s.get('accession') for s in run.findall('./Pool/Member')}
     if len(members)==1 and next(iter(members)) in samples:sample=samples[next(iter(members))]
     elif len(samples)==1:sample=next(iter(samples.values()))
     else:
      issues.append(dict(run=run.get('accession'),type='unresolved_multiple_sample_package',evidence=str(p)));continue
     attrs={a.findtext('TAG'):a.findtext('VALUE','') for a in sample.findall('.//SAMPLE_ATTRIBUTE')}
     bio=sample.find(".//EXTERNAL_ID[@namespace='BioSample']");bio=bio.text if bio is not None else ''
     accession=run.get('accession');r=dict(run_accession=accession,source_bioproject=project,biosample_accession=bio,
      sra_sample_accession=sample.get('accession',''),sample_name=sample.get('alias',''),experiment_accession=ex.get('accession',''),
      study_accession=study.get('accession','') if study is not None else '',library_strategy=ex.findtext('.//LIBRARY_STRATEGY',''),library_source=ex.findtext('.//LIBRARY_SOURCE',''),
      host=attrs.get('host',''),disease=attrs.get('host_disease',attrs.get('disease','')),phenotype=attrs.get('phenotype',''),condition=attrs.get('condition',''),
      group=attrs.get('SampleType',attrs.get('Case_status','')),treatment=attrs.get('treatment',''),collection_metadata=attrs.get('collection_date',''),
      available_sample_attributes=json.dumps(attrs,ensure_ascii=False,sort_keys=True),ncbi_evidence_file=str(p),ena_evidence_file='',ncbi_ena_accession_consistency='NA')
     if accession in data and data[accession]['biosample_accession'] not in ['',bio]:issues.append(dict(run=accession,type='conflicting_BioSample',evidence=str(p)))
     data[accession]=r
  except Exception as ex:issues.append(dict(run='',type='parse_error:'+str(ex),evidence=str(p)));event('parse NCBI XML',p,repr(ex),1)
 for row in rows(E/'project_metadata_availability.tsv'):
  if not row['ena_file']:continue
  p=Path(row['ena_file']);register([p])
  for r in rows(p):
   run=r['run_accession'];sample_values=set(r.get('sample_accession','').split(';')+r.get('secondary_sample_accession','').split(';'))
   if run not in data:
    data[run]=dict(run_accession=run,source_bioproject=row['project'],biosample_accession=r.get('sample_accession',''),sra_sample_accession=r.get('secondary_sample_accession',''),sample_name=r.get('sample_alias',''),
     experiment_accession=r.get('experiment_accession',''),study_accession=r.get('study_accession',''),library_strategy=r.get('library_strategy',''),library_source=r.get('library_source',''),
     host='',disease='',phenotype='',condition='',group='',treatment='',collection_metadata='',available_sample_attributes='{}',ncbi_evidence_file='',ena_evidence_file=str(p),ncbi_ena_accession_consistency='ENA_only')
   else:
    n=data[run];same=bool(sample_values & {n['biosample_accession'],n['sra_sample_accession']}) and n['experiment_accession']==r.get('experiment_accession',n['experiment_accession'])
    n.update(ena_evidence_file=str(p),ncbi_ena_accession_consistency=('NCBI_run_member_resolves_ENA_pool' if ';' in r.get('sample_accession','') else 'Yes') if same else 'No')
    if not same:issues.append(dict(run=run,type='NCBI_ENA_disagree',evidence=str(p)))
   data[run]['ena_biosample_candidates']=r.get('sample_accession','')
 for r in data.values():
  r.setdefault('ena_biosample_candidates','')
  a=json.loads(r['available_sample_attributes']);r['participant_id_public']=a.get('participant_id',a.get('host_subject_id',a.get('subject_id','')))
  r['timepoint_public']=a.get('timepoint',a.get('time_point',''))
 write(E/'public_run_sample_crosswalk.tsv',[data[k] for k in sorted(data)])
 write(Q/'public_metadata_qc.tsv',[dict(check='NCBI_ENA_run_BioSample_identity',status='PASS' if not issues else 'FAIL',details=json.dumps(issues)),dict(check='unique_public_Run_records',status='PASS',details=str(len(data)))])
 event('compile_rescue_metadata.py',C,E/'public_run_sample_crosswalk.tsv',1 if issues else 0)
 print('Normalized public runs',len(data),'issues',len(issues))
if __name__=='__main__':main()
