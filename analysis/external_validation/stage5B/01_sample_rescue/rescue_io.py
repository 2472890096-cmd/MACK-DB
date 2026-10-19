"""Scoped read-only provenance rescue utilities. No network-inference code."""

# v0.1: repository-local IO configuration; scientific functions below are unchanged.
import sys as _path_sys
from pathlib import Path as _RepoPath
_repo_root = next(p for p in _RepoPath(__file__).resolve().parents if (p / "lib/mackdb_paths.py").is_file())
_path_sys.path.insert(0, str(_repo_root / "lib"))
from mackdb_paths import resolve_path
from pathlib import Path
import csv,json,hashlib,datetime,urllib.request,urllib.parse,sys,threading
R=Path(str(resolve_path('stage5b', '01_sample_rescue')))
STAGE=R.parent;OLD=STAGE.parent/'9_16_external_validation';PROV=OLD/'stage3b_source_exclusion/provenance'
C=R/'cache';E=R/'evidence';L=R/'logs';Q=R/'qc';S=R/'sample_lists'
csv.field_size_limit(100000000)
lock=threading.Lock()
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,data,fields=None):
 p=Path(p);assert R in p.parents;data=list(data)
 with p.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
def event(cmd,inp='',out='',code=0):
 with lock:
  p=L/'execution_log.tsv';new=not p.exists()
  with p.open('a',newline='') as f:
   w=csv.writer(f,delimiter='\t',lineterminator='\n')
   if new:w.writerow(['command','input','output','exit_code','timestamp'])
   w.writerow([cmd,inp,out,code,now()])
def evidence(study,kind,title,url,file,support,status='available',error='',**extra):
 p=Path(file) if file else None
 r=dict(study=study,resource_type=kind,title=title,url=url,access_date=now(),downloaded_file=str(file),sha256=sha(p) if p and p.is_file() else '',what_it_supports=support,status=status,error=error,**extra)
 with lock:
  with (E/'resources.jsonl').open('a') as f:f.write(json.dumps(r,ensure_ascii=False)+'\n')
 return r
def register(paths):
 p=Q/'input_snapshot.json';old=json.loads(p.read_text()) if p.exists() else {}
 for f in map(Path,paths):
  if R in f.parents:continue
  if str(f) not in old:old[str(f)]=dict(exists=f.is_file(),sha256=sha(f) if f.is_file() else '',size=f.stat().st_size if f.is_file() else 0)
 p.write_text(json.dumps(old,ensure_ascii=False,indent=2)+'\n')
def get(study,url,name,kind='metadata',support='',timeout=40,max_mb=200):
 p=C/name
 if p.is_file() and p.stat().st_size:
  evidence(study,kind,name,url,p,support,'cache_reuse');event('reuse cached public resource',url,p);return p
 try:
  req=urllib.request.Request(url,headers={'User-Agent':'MACKDB-public-sample-provenance-audit/0.5','Accept':'*/*'})
  with urllib.request.urlopen(req,timeout=timeout) as response:
   b=response.read(max_mb*1024*1024+1);resolved=response.url
   if len(b)>max_mb*1024*1024:raise ValueError('Resource exceeds bounded metadata download size')
  if name.endswith(('.xlsx','.zip')) and not b.startswith(b'PK'):raise ValueError('Not ZIP/XLSX: '+repr(b[:80]))
  if name.endswith('.json'):json.loads(b)
  if name.endswith('.tsv') and b'\t' not in b[:5000]:raise ValueError('Not a TSV: '+repr(b[:80]))
  if name.endswith('.pdf') and not b.startswith(b'%PDF'):raise ValueError('Not a PDF')
  p.write_bytes(b);evidence(study,kind,name,url,p,support,'downloaded',resolved_url=resolved,size_bytes=len(b));event('HTTPS GET',url,p)
  print('OK',study,name,len(b),flush=True);return p
 except Exception as ex:
  msg=type(ex).__name__+': '+str(ex);evidence(study,kind,name,url,'',support,'failed',msg);event('HTTPS GET',url,p,1)
  print('FAILED',study,name,msg,flush=True);return None
if __name__=='__main__':
 if sys.argv[1]=='get':get(*sys.argv[2:])
 elif sys.argv[1]=='snapshot':
  paths=[STAGE/'STAGE5B_0_REPORT.md']+list((STAGE/'00_audit').glob('*'))+list((STAGE/'00_metadata_cache').glob('*'))
  register([p for p in paths if p.is_file()]);event('snapshot all frozen Stage5B-0 files',STAGE,Q/'input_snapshot.json')
  print('Frozen',len(paths),'Stage5B-0 paths')
