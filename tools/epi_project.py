#!/usr/bin/env python3
"""Portable single-project setup for the unchanged EPI v3 production worker.
Only input validation, runtime paths, and entry orchestration live here.
"""
from pathlib import Path
import argparse,csv,hashlib,importlib.util,json,math,os,re,shutil,sys
REPO=Path(__file__).resolve().parents[1]
CORE_SHA='98d832663fea152e0d63c481106ff24fb5faca4e34ce332d0187814cbead5cc2'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def table(p,rows,fields):
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def inspect_input(a,m):
 with a.open() as f:
  rr=csv.reader(f,delimiter='\t');head=next(rr);n=len(head)-1
  if n<2 or len(set(head[1:]))!=n:raise ValueError('Need at least two unique sample columns')
  ids=[]
  for row in rr:
   if len(row)!=n+1:raise ValueError('Abundance matrix is not rectangular')
   v=[float(x) for x in row[1:]]
   if not all(math.isfinite(x) and x>=0 for x in v):raise ValueError('Nonfinite or negative abundance')
   ids.append(row[0])
 if not ids or ids!=list(map(str,range(len(ids)))):raise ValueError('Taxon IDs must be 0-based ordered rows')
 with m.open() as f: pairs=list(csv.reader(f,delimiter='\t'))
 if len(pairs)!=len(ids) or any(len(x)!=2 for x in pairs) or [x[0] for x in pairs]!=ids:raise ValueError('Taxon mapping must correspond exactly to abundance rows; no automatic transpose')
 names=[x[1] for x in pairs]
 if len(set(names))!=len(names) or any(not x or x!=x.strip() for x in names):raise ValueError('Taxon names empty/duplicate/whitespace padded')
 return n,len(ids)
def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--level',choices=['species','genus'],required=True);p.add_argument('--project',required=True)
 p.add_argument('--abundance',type=Path,required=True);p.add_argument('--mapping',type=Path,required=True)
 p.add_argument('--run-root',type=Path,required=True);p.add_argument('--rscript',default=os.environ.get('MACKDB_RSCRIPT','Rscript'));p.add_argument('--resume',action='store_true')
 a=p.parse_args()
 if not re.fullmatch('PRJ[A-Z0-9]+',a.project):p.error('project must match PRJ[A-Z0-9]+')
 rscript=shutil.which(a.rscript)
 if not rscript:p.error('Rscript executable unavailable; no installation is attempted')
 root=a.run_root.resolve(); abundance=a.abundance.resolve();mapping=a.mapping.resolve();n,t=inspect_input(abundance,mapping)
 if root.exists() and not a.resume:p.error('run-root exists; explicitly use --resume after inspecting status')
 for d in ['scripts','audit','manifests','status/projects','logs/projects','logs/master','keys','summaries','tmp','final_2of3']: (root/d).mkdir(parents=True,exist_ok=True)
 origin=REPO/'pipeline/epi_v3'/a.level/'scripts'
 workername='workflow.py' if a.level=='species' else 'species_production_workflow.py'
 core=REPO/'pipeline/epi_v3/core/get_keystone_EPI_single_safe_qecdf_v3_fast.R'
 assert digest(core)==CORE_SHA,'Frozen EPI core hash mismatch'
 if not a.resume:
  for name in [workername,'run_one_project.sh']:shutil.copy2(origin/name,root/'scripts'/name)
  destcore=root/'scripts'/core.name;shutil.copy2(core,destcore)
  row=dict(project=a.project,sample_count=n,taxa_count=t,abundance_path=str(abundance),species_path=str(mapping))
  for kind,file in [('abundance',abundance),('species',mapping)]:row.update({kind+'_sha256':digest(file),kind+'_bytes':file.stat().st_size,kind+'_mtime_ns':file.stat().st_mtime_ns})
  fields=list(row)
  table(root/'manifests/all_valid_projects.tsv',[row],fields)
  for label,include in [('le5000',n<=5000),('gt5000',n>5000)]:table(root/f'manifests/projects_{label}.tsv',[row] if include else [],fields)
  cfg=dict(core_script=str(destcore),core_sha256=CORE_SHA,rscript=rscript,valid_projects=1,production_root=str(root),scientific_version='EPI_qecdf_v3',level=a.level)
  (root/'audit/config.json').write_text(json.dumps(cfg,indent=2)+'\n')
  protected=[dict(path=str(x),sha256=digest(x)) for x in [abundance,mapping,core]]
  table(root/'audit/protected_fingerprints.tsv',protected,['path','sha256'])
  files=list((root/'scripts').iterdir())+list((root/'manifests').iterdir())+[root/'audit/config.json']
  table(root/'audit/workflow_fingerprints.tsv',[dict(relative_path=str(x.relative_to(root)),sha256=digest(x)) for x in files],['relative_path','sha256'])
 os.environ['MACKDB_EPI_RUN_ROOT']=str(root)
 spec=importlib.util.spec_from_file_location('mackdb_epi_worker',root/'scripts'/workername);w=importlib.util.module_from_spec(spec);sys.modules[spec.name]=w;spec.loader.exec_module(w)
 rows=w.verify_setup();row=rows[0]
 if [row['project'],row['abundance_path'],row['species_path']]!=[a.project,str(abundance),str(mapping)]:raise ValueError('Resume arguments disagree with frozen manifest')
 rc=w.worker(a.project,str(abundance),str(mapping),str(n))
 if rc==0:w.postprocess(rows)
 return rc
if __name__=='__main__':raise SystemExit(main())
