from pathlib import Path
from collections import defaultdict
import csv,json,hashlib,datetime,os
import numpy as np
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent.parent;EXT=BASE/'9_16_external_validation';B5=BASE/'9_22_stage5B_same_data_cross_method';N5=B5/'04_network_comparison';RES=B5/'01_sample_rescue';PROJECT=Path('/home/wangmengyang/project');NA='NA'
def read(p):
 with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rs,fields=None):
 p=Path(p);assert ROOT in p.resolve().parents;rs=list(rs);fields=fields or list(rs[0]);temp=p.with_suffix(p.suffix+'.tmp')
 with temp.open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
 os.replace(temp,p)
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
HIST=defaultdict(list)
HFILES=[(B5/'00_audit/stage5B_artifact_manifest.tsv','path','sha256'),(RES/'qc/delivery_artifact_manifest.tsv','path','sha256'),(B5/'02_prepare_fastspar_input/qc/stage5B1A_artifact_manifest.tsv','path','sha256'),(N5/'00_inputs/input_manifest.tsv','path','sha256_before'),(N5/'qc/stage5B2_delivery_manifest.tsv','path','sha256'),(EXT/'stage4a_recurrence/qc/04a_input_manifest.tsv','input_file','sha256'),(EXT/'stage4a_recurrence/qc/04a_formal_output_manifest.tsv','absolute_path','sha256'),(EXT/'stage5a_cross_study_complementarity/00_input_audit/stage5A_input_manifest.tsv','path','sha256_before')]
for p,k,h in HFILES:
 for r in read(p):HIST[r[k]].append((str(p),r[h]))
sp=RES/'qc/input_snapshot.json'
for path,r in json.loads(sp.read_text()).items():
 if r.get('exists') and r.get('sha256'):HIST[path].append((str(sp),r['sha256']))
MAN={r['input_path']:r for r in read(ROOT/'stage6B_input_manifest.tsv')} if (ROOT/'stage6B_input_manifest.tsv').exists() else {}
QC=read(ROOT/'stage6B_taxon_mapping_and_edge_qc.tsv') if (ROOT/'stage6B_taxon_mapping_and_edge_qc.tsv').exists() else []
def check(scope,name,ok,detail='',warn=False):
 QC.append(dict(scope=scope,check=name,status='PASS' if ok else 'WARN' if warn else 'FAIL',details=str(detail)))
 if not ok and not warn:
  save();raise AssertionError((scope,name,detail))
def reg(p,role='frozen_input'):
 p=Path(p).resolve();s=str(p)
 if s not in MAN:
  h=sha(p);expected=HIST.get(s,[]);check(role,'historical_SHA256_match',all(h==x for _,x in expected),s)
  MAN[s]=dict(input_path=s,role=role,size_bytes=p.stat().st_size,sha256_before=h,sha256_after=NA,historical_sha256=expected[0][1] if expected else NA,historical_evidence=';'.join(sorted({m for m,h in expected})) or NA,status='VERIFIED_BEFORE')
 return p
def rr(p,role='frozen_input'):return read(reg(p,role))
def save():
 if MAN:write(ROOT/'stage6B_input_manifest.tsv',MAN.values())
 if QC:write(ROOT/'stage6B_taxon_mapping_and_edge_qc.tsv',QC)
def finish_inputs():
 for p,r in MAN.items():r['sha256_after']=sha(p);r['status']='PASS' if r['sha256_after']==r['sha256_before'] else 'FAIL';check('FROZEN','input_SHA256_unchanged',r['status']=='PASS',p)
 save()
def matrix(p):
 with reg(p,'numeric_matrix').open() as f:
  rd=csv.reader(f,delimiter='\t');head=next(rd);rows=list(rd)
 return head[1:],[r[0] for r in rows],np.asarray([r[1:] for r in rows],dtype=float)
def ratio(a,b):return a/b if b else NA
def pack(xs):return ';'.join(sorted(set(xs)))
def sign(x):return 'positive' if x>0 else 'negative' if x<0 else 'zero'
def edgekey(r):return r['taxon1'],r['taxon2']
def ensure_manifests():
 for p,k,h in HFILES:reg(p,'historical_manifest')
 reg(sp,'historical_input_snapshot')
ABUNDANCE={p.name:p for p in (PROJECT/'GDTB226_species_abundance_result/abundance').iterdir() if p.is_file()}
FASTSPAR={p.name:p for p in (PROJECT/'r226_fastspar_species_result').iterdir() if p.is_file()}
