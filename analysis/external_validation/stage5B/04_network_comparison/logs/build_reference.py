
# v0.1: repository-local IO configuration; scientific functions below are unchanged.
import sys as _path_sys
from pathlib import Path as _RepoPath
_repo_root = next(p for p in _RepoPath(__file__).resolve().parents if (p / "lib/mackdb_paths.py").is_file())
_path_sys.path.insert(0, str(_repo_root / "lib"))
from mackdb_paths import resolve_path
from pathlib import Path
import csv,gzip,hashlib,json,time
ROOT=Path(__file__).resolve().parents[1]
REF=Path(str(resolve_path('gtdb_r226')))
prior=Path(str(resolve_path('literature', 'stage4a_recurrence/qc/04a_input_manifest.tsv')))
with prior.open() as f: expected={r['input_file']:r['sha256'] for r in csv.DictReader(f,delimiter='\t')}
names=set();manifest=[]
for p in [REF/'bac120_metadata_r226.tsv.gz',REF/'ar53_metadata_r226.tsv.gz']:
 h=hashlib.file_digest(p.open('rb'),'sha256').hexdigest();assert h==expected[str(p)]
 with gzip.open(p,'rt') as f:
  header=f.readline().rstrip('\n').split('\t');idx=header.index('gtdb_taxonomy')
  for line in f:
   species=line.split('\t')[idx].split(';')[-1]
   if species.startswith('s__') and len(species)>3:names.add(species[3:])
 manifest.append(dict(path=str(p),sha256=h,bytes=p.stat().st_size));print('REFERENCE_READ',p.name,len(names),flush=True)
encoded={}
for n in sorted(names):encoded.setdefault(n.replace(' ','_'),[]).append(n)
assert all(len(x)==1 for x in encoded.values()),'GTDB serialized label not globally unique'
with (ROOT/'00_inputs/gtdb_r226_species_vocabulary.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['gtdb_r226_species','production_encoded_label']);w.writerows((n,n.replace(' ','_')) for n in sorted(names))
(ROOT/'00_inputs/reference_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('VOCABULARY_COMPLETE',len(names),flush=True)
import numba,numpy,scipy,matplotlib,openpyxl
print('VERSIONS',json.dumps({m.__name__:m.__version__ for m in [numba,numpy,scipy,matplotlib,openpyxl]}),flush=True)
