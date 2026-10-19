#!/usr/bin/env python3
"""Read-only dependency probe; no installation, downloads, or scientific jobs."""
from pathlib import Path
import argparse,csv,importlib.metadata,importlib.util,json,os,shutil,subprocess,sys
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--rscript',default=os.environ.get('MACKDB_RSCRIPT','Rscript'));p.add_argument('--dki-python');p.add_argument('--output',type=Path,default=ROOT/'audit/DEPENDENCY_PROBE.json');a=p.parse_args()
 result={'python':sys.version,'python_executable':sys.executable,'python_packages':{},'executables':{},'R':{},'configured_external_inputs':{}}
 for name in ['numpy','pandas','scipy','statsmodels','openpyxl','matplotlib','torch','torchdiffeq']:
  try:result['python_packages'][name]=importlib.metadata.version(name)
  except importlib.metadata.PackageNotFoundError:result['python_packages'][name]='MISSING_IN_THIS_PYTHON'
 for name in ['metakssd','fastspar','fastspar_bootstrap','fastspar_pvalues','parallel','Rscript','bash']:
  var={'metakssd':'METAKSSD_BIN','fastspar':'FASTSPAR_BIN','fastspar_bootstrap':'FASTSPAR_BOOTSTRAP_BIN','fastspar_pvalues':'FASTSPAR_PVALUES_BIN','parallel':'PARALLEL_BIN'}.get(name,'')
  result['executables'][name]=shutil.which(os.environ.get(var,name))
 rs=shutil.which(a.rscript)
 if rs:
  packages=['vegan','microeco','igraph','tidyverse','optparse','dplyr','data.table','htmlwidgets','networkD3','gridExtra','ggrepel','RColorBrewer','future.apply']
  expr='cat("R_VERSION\\t",R.version.string,"\\n",sep=""); for(p in c('+','.join(json.dumps(x) for x in packages)+')) { ok<-requireNamespace(p,quietly=TRUE); cat(p,if(ok)as.character(packageVersion(p))else "MISSING",if(ok)packageDescription(p)$License else "NA",sep="\\t");cat("\\n") }'
  q=subprocess.run([rs,'-e',expr],capture_output=True,text=True,timeout=45)
  result['R']={'executable':rs,'returncode':q.returncode,'versions_and_license':q.stdout,'stderr':q.stderr}
 else:result['R']={'status':'BLOCKED','reason':'Rscript executable not found'}
 if a.dki_python:
  q=subprocess.run([a.dki_python,'-c','import sys,importlib.metadata as m,json; print(json.dumps(dict(python=sys.version,torch=m.version("torch"),torchdiffeq=m.version("torchdiffeq"),numpy=m.version("numpy"))))'],capture_output=True,text=True,timeout=20)
  result['dki_python_probe']={'executable':a.dki_python,'returncode':q.returncode,'stdout':q.stdout,'stderr':q.stderr}
 core=os.environ.get('FASTCCLASSO_CORE',str(ROOT/'external/fastCCLasso/RCode/fastCCLasso.R'))
 result['configured_external_inputs']['fastCCLasso_core']={'path':core,'exists':Path(core).is_file(),'status':'AVAILABLE_NOT_EXECUTED' if Path(core).is_file() else 'BLOCKED_EXTERNAL_CORE_NOT_BUNDLED_OR_CONFIGURED'}
 result['scientific_algorithms_executed']=False;result['installed_or_downloaded']=False
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n');print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
