from pathlib import Path
from collections import defaultdict,Counter
import csv,json,hashlib,os,math
import numpy as np
OUT=Path(__file__).resolve().parents[1];FOLLOW=OUT.parent;BASE=FOLLOW.parent;EXT=BASE/'9_16_external_validation';B5=BASE/'9_22_stage5B_same_data_cross_method';N5=B5/'04_network_comparison';S6=FOLLOW/'02_stage6B_project_network_comparison';PROJECT=Path('/home/wangmengyang/project');NA='NA'
CURRENT=FOLLOW/'04_figure7_revision_20261023/final_clean/figures/revised_common_species';CURRENT_SCRIPT=CURRENT/'plot_figure7_revised_common_species.py'
def read(p):
 with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows,fields=None):
 p=Path(p);assert OUT in p.resolve().parents;rows=list(rows);fields=fields or list(rows[0]);tmp=p.with_suffix(p.suffix+'.tmp')
 with tmp.open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
 os.replace(tmp,p)
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
MAN={r['path']:r for r in read(OUT/'input_manifest.tsv')};HIST=defaultdict(list)
for mp,pk,hk in [(N5/'00_inputs/input_manifest.tsv','path','sha256_before'),(N5/'qc/stage5B2_delivery_manifest.tsv','path','sha256'),(S6/'stage6B_input_manifest.tsv','input_path','sha256_before')]:
 for r in read(mp):HIST[r[pk]].append(r[hk])
Q=read(OUT/'qc/builder_qc.tsv') if (OUT/'qc/builder_qc.tsv').exists() else []
def check(scope,name,ok,detail=''):
 Q.append(dict(scope=scope,check=name,status='PASS' if bool(ok) else 'FAIL',details=str(detail)))
 if not ok:save();raise AssertionError((scope,name,detail))
def reg(p,role='scientific_input'):
 p=Path(p).resolve();s=str(p)
 if s not in MAN:
  h=sha(p);expected=HIST.get(s,[]);check(role,'historical_hash_matches',all(h==v for v in expected),s);MAN[s]=dict(path=s,role=role,size_bytes=p.stat().st_size,sha256_before=h,sha256_after=NA,historical_sha256=expected[0] if expected else NA,status='BEFORE')
 elif role not in MAN[s]['role'].split(';'):MAN[s]['role']+=';'+role
 return p
def rr(p,role='scientific_input'):return read(reg(p,role))
def save():
 write(OUT/'input_manifest.tsv',MAN.values())
 if Q:write(OUT/'qc/builder_qc.tsv',Q)
def finish():
 for p,r in MAN.items():r['sha256_after']=sha(p);r['status']='PASS' if r['sha256_before']==r['sha256_after'] else 'FAIL';check('INPUT','unchanged_SHA256',r['status']=='PASS',p)
 save()
def matrix(p):
 with reg(p,'original_existing_matrix').open() as f:
  rd=csv.reader(f,delimiter='\t');head=next(rd);rows=list(rd)
 return head[1:],[r[0] for r in rows],np.array([r[1:] for r in rows],dtype=float)
def edge_rows(p):
 rs=rr(p,'existing_species_edge_table');es={(r['taxon1'],r['taxon2']) for r in rs};check(str(p),'canonical_unique_species_pairs',len(es)==len(rs) and all(a<b for a,b in es));return rs,es
def ratio(a,b):return a/b if b else NA
def dump_edges(p,es):write(p,[dict(taxon1=a,taxon2=b) for a,b in sorted(es)],['taxon1','taxon2'])
def key(r):return r['taxon1'],r['taxon2']
def edgeid(e):return hashlib.sha256(('\t'.join(e)).encode()).hexdigest()
def find_name(root,name):
 xs=[p for p in Path(root).rglob(name) if p.is_file()];assert len(xs)==1,(root,name,xs);return xs[0]
