
# v0.1: repository-local IO configuration; scientific functions below are unchanged.
import sys as _path_sys
from pathlib import Path as _RepoPath
_repo_root = next(p for p in _RepoPath(__file__).resolve().parents if (p / "lib/mackdb_paths.py").is_file())
_path_sys.path.insert(0, str(_repo_root / "lib"))
from mackdb_paths import resolve_path
from rescue_io import *
import collections,xml.etree.ElementTree as ET
AB=Path(str(resolve_path('species_input')))
def main():
 inventory={r['project']:Path(r['path']) for r in rows(E/'local_abundance_root_inventory.tsv')}
 prior=rows(STAGE/'00_audit/stage5B_mackdb_sample_index.tsv');oldprojects={r['source_bioproject'] for r in prior}
 local=list(prior);audit=[]
 for r in rows(E/'project_metadata_availability.tsv'):
  project=r['project'];aliases={project}
  p=C/(project+'_study_aliases.xml')
  if p.exists():
   t=ET.parse(p);aliases|={e.text for e in t.findall('.//IDENTIFIERS/SECONDARY_ID') if e.text}
  elif r['ena_file']:
   for e in rows(r['ena_file']):aliases|={e.get('study_accession',''),e.get('secondary_study_accession','')}
  aliases.discard('');found=[a for a in sorted(aliases) if a in inventory]
  count=0;status='UNAVAILABLE';note='No canonical or exact public alias filename in 1577-file abundance root'
  if found:
   for a in found:
    p=inventory[a];register([p]);header=p.open().readline().rstrip('\r\n').split('\t')[1:];count+=len(header)
    assert len(set(header))==len(header)
    if a in oldprojects:
     index=[x for x in prior if x['source_bioproject']==a]
     assert header==[x['column_index'] for x in index]
     status='PASS';note='Reused frozen Stage5B-0 per-column original-profile identity proof; header count/order independently checked'
    else:
     runlist=Path(str(resolve_path('project_data', 'all_sample')))/(a+'_uniq.txt');species=AB/'species'/(a+'_species.matrix');register([runlist,species])
     runs=runlist.read_text().splitlines();assert len(runs)==len(header)==len(set(runs)) and header==[str(i) for i in range(len(runs))]
     with p.open() as f:next(f);values=[l.rstrip('\r\n').split('\t') for l in f]
     taxa=[l.split('\t') for l in species.read_text().splitlines()];assert [r[0] for r in taxa]==[r[0] for r in values]
     profiles=[Path(str(resolve_path('project_data', 'profile_result')))/(run+'_profile.txt') for run in runs];register(profiles)
     for j,(run,pf) in enumerate(zip(runs,profiles)):
      d=dict(l.split('\t')[:2] for l in pf.read_text().splitlines() if l)
      assert all(float(v[j+1])==float(d.get(t[1],'0')) for v,t in zip(values,taxa)),(a,run)
      local.append(dict(source_bioproject=a,column_index=header[j],column_position_1based=j+1,run_accession=run,abundance_file=str(p),run_list_file=str(runlist),source_profile_file=str(pf),profile_values_match_every_retained_taxon='Yes',identity_rule='ordered run list AND all retained original profile taxon values equal column'))
     status='PASS';note='New local candidate project: every retained profile value in every run column checked'
  audit.append(dict(project=project,study_role=r['studies'],exact_public_aliases=';'.join(sorted(aliases)),matching_abundance_filename_projects=';'.join(found),
    abundance_files=';'.join(str(inventory[a]) for a in found),local_run_columns=count,status=status,notes=note))
 write(E/'local_project_and_alias_audit.tsv',audit);write(E/'verified_local_run_columns.tsv',local)
 event('audit_local_rescue.py: exact canonical/secondary alias inventory and column validation',AB,E/'local_project_and_alias_audit.tsv')
 print('Projects checked',len(audit),'available',sum(r['status']=='PASS' for r in audit),'total indexed run columns including frozen Wallen',len(local))
if __name__=='__main__':main()
