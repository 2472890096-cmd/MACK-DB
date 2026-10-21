#!/usr/bin/env python3
import csv,json,collections,hashlib
from pathlib import Path
O = None
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def write(fn,rr,fields=None):
 rr=list(rr);p=O/fn;p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w') as h:
  w=csv.DictWriter(h,fieldnames=fields or list(rr[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
cats={1:'Zi-Pi only',2:'EPI only',4:'DKI only',3:'Zi-Pi + EPI only',5:'Zi-Pi + DKI only',6:'EPI + DKI only',7:'All three methods'};methods=['ZiPi','EPI','DKI'];bits=dict(zip(methods,[1,2,4]));projects={'PRJEB80533','PRJNA1345363','PRJNA1357645','PRJEB72274'}
scopes={1571:'COMMON_PROJECTS_SCIENTIFICALLY_USABLE.txt',1570:'COMMON_PROJECTS_SCIENTIFICALLY_USABLE_EXCLUDING_UNVERIFIED_DKI.txt',1568:'COMMON_PROJECTS_STRICT_SUCCESS.txt'}
def calc(records):
 vals={};contexts={}
 for level in ['species','genus']:
  for mode in ['Broad','Strict']:
   d={}
   for m in methods:
    for p,n,e in records[level,m]:
     key=(e,n) if mode=='Broad' else (p,n);v=d.setdefault(key,dict(environment=e,mask=0,sources={x:set() for x in methods}));v['mask']|=bits[m];v['sources'][m].add(p)
   contexts[level,mode]=d
   sets=[{(e,n) if mode=='Broad' else (p,n) for p,n,e in records[level,m]} for m in methods];z,e,k=sets
   partitions=[z-e-k,e-z-k,k-z-e,(z&e)-k,(z&k)-e,(e&k)-z,z&e&k];assert sum(map(len,partitions))==len(z|e|k)==len(d)
   # independent membership masks vs direct set differences for every context
   for mask,st in zip(cats,partitions):assert st=={key for key,v in d.items() if v['mask']==mask}
   if mode=='Broad':
    for mask,st in zip(cats,partitions):vals['A' if level=='species' else 'B',cats[mask],level]=len(st)
   ge2=(z&e)|(z&k)|(e&k);a3=z&e&k
   for label,st,check in [('>=2 methods',ge2,lambda v:v['mask'].bit_count()>=2),('All three methods',a3,lambda v:v['mask']==7)]:
    unique={n for _,n in st};assert unique=={n for (_,n),v in d.items() if check(v)};vals['C' if mode=='Broad' else 'D',mode+' '+label,level]=len(unique)
 return vals,contexts

def main():
 import argparse
 global O
 parser=argparse.ArgumentParser(description="Accepted unified-R Figure3 statistics; no database access")
 parser.add_argument('--base-analysis',type=Path,required=True,help='Frozen analysis directory: inputs/ six candidate tables and old22 archive')
 parser.add_argument('--scopes-dir',type=Path,required=True)
 parser.add_argument('--r-candidates',type=Path,required=True)
 parser.add_argument('--impact-audit',type=Path,required=True)
 parser.add_argument('--output',type=Path,required=True)
 args=parser.parse_args()
 F=args.base_analysis.resolve(); S=args.scopes_dir.resolve(); A=args.impact_audit.resolve(); O=args.output.resolve()
 O.mkdir(parents=True,exist_ok=True); (O/'evidence').mkdir(exist_ok=True)
 source={};envs=collections.defaultdict(set)
 for l in ['species','genus']:
  for m,suff in zip(methods,['','_epi','_dki']):
   rr=read(F/'inputs'/f'key_{l}{suff}.tsv');st={(r['bioproject'],r['name'],r['environment']) for r in rr};assert len(st)==len(rr);source[l,m]=st
   for p,n,e in st:envs[p].add(e)
 assert all(len(x)==1 for x in envs.values())
 newrows=read(args.r_candidates);newset={(r['bioproject'],r['name'],r['environment']) for r in newrows};assert len(newset)==61
 arch={(int(r['scope']),r['panel'],r['category'],r['taxonomic_level']):int(r['value']) for r in read(F/'FIGURE3_ALL_FOUR_SCOPES_22_VALUES.tsv')}
 allvals=[];trans=[];support=[];scopeaudit=[];qa={}
 for count,fn in scopes.items():
  text=(S/fn).read_text();items=text.split();scope=set(items);assert len(items)==len(scope)==count and 'PRJEB51353' not in scope and projects<=scope
  target=O/'evidence'/fn;target.write_text(text);scopeaudit.append(dict(scope=count,path=str(S/fn),SHA256=hashlib.sha256((S/fn).read_bytes()).hexdigest(),unique_projects=len(scope),duplicates=0,PRJEB51353_present=False))
  base={k:{r for r in v if r[0] in scope} for k,v in source.items()};new={k:set(v) for k,v in base.items()};new['species','DKI']={r for r in base['species','DKI'] if r[0] not in projects}|newset
  bv,bc=calc(base);nv,nc=calc(new)
  for key,v in bv.items():assert arch[(count,*key)]==v,(count,key)
  for key,v in nv.items():
   panel,cat,level=key;allvals.append(dict(scope=count,panel=panel,category=cat,taxonomic_level=level,old_mixed_value=bv[key],unified_R_value=v,difference=v-bv[key],unit='environment_taxon_records' if panel in ['A','B'] else 'unique_taxa'))
  for level in ['species','genus']:
   for mode in ['Broad','Strict']:
    b=bc[level,mode];a=nc[level,mode];rows=[]
    for key,v in sorted(a.items()):
     ctx,n=key;row=dict(scope=count,taxonomic_level=level,mode=mode,name=n,environment=v['environment'],bioproject=ctx if mode=='Strict' else '',combination=cats[v['mask']],method_count=v['mask'].bit_count(),**{m+'_projects':';'.join(sorted(v['sources'][m])) for m in methods});rows.append(row)
    write(f'details/{count}_{level}_{mode}_ALL_COMBINATIONS.tsv',rows)
    for key in sorted(set(a)|set(b)):
     av=a.get(key,dict(mask=0,sources={m:set() for m in methods}));bv_=b.get(key,dict(mask=0,sources={m:set() for m in methods}));ctx,n=key
     if av['mask']!=bv_['mask'] or av['sources']!=bv_['sources']:
      trans.append(dict(scope=count,taxonomic_level=level,mode=mode,name=n,environment=(a.get(key) or b[key])['environment'],bioproject=ctx if mode=='Strict' else '',old_combination=cats.get(bv_['mask'],'NONE'),new_combination=cats.get(av['mask'],'NONE'),**{m+'_'+label+'_projects':';'.join(sorted(v['sources'][m])) for m in methods for label,v in [('old',bv_),('new',av)]}))
   bo={(p,n,v['environment']) for (p,n),v in bc[level,'Strict'].items() if v['mask']==7};an={(p,n,v['environment']) for (p,n),v in nc[level,'Strict'].items() if v['mask']==7};assert bo==an
   for p,n,e in sorted(an):support.append(dict(scope=count,taxonomic_level=level,bioproject=p,name=n,environment=e,old_and_new='IDENTICAL'))
  qa[count]=dict(baseline_22_exact=True,independent_sets=True,all3_names_and_support_identical=True,seven_exclusive_partition=True)
 expected={(r['panel'],r['category'],r['taxonomic_level']):int(r['simulated_R_value']) for r in read(A/'FIGURE3_POTENTIAL_IMPACT.tsv') if r['scenario']=='ALL_FOUR_R'}
 assert {(r['panel'],r['category'],r['taxonomic_level']):r['unified_R_value'] for r in allvals if r['scope']==1571}==expected
 write('FIGURE3_UNIFIED_R_22_VALUES.tsv',[r for r in allvals if r['scope']==1571]);write('FIGURE3_UNIFIED_R_SENSITIVITY.tsv',[r for r in allvals if r['scope']!=1571]);write('FIGURE3_OLD_NEW_DIFFERENCES.tsv',allvals);write('details/ENVIRONMENT_AND_PROJECT_TRANSITIONS.tsv',trans);write('details/STRICT_ALL3_NAMES_AND_SUPPORT.tsv',support);write('details/HUMAN_GUT_STRICT_ALL3_NAMES_AND_SUPPORT.tsv',[r for r in support if r['environment']=='human_gut']);write('evidence/PROJECT_SCOPE_AUDIT.tsv',scopeaudit)
 (O/'evidence/FIGURE3_VALIDATION.json').write_text(json.dumps(qa,indent=2)+'\n')
 print('PASS: old22 reproduced for1571/1570/1568; independent fullcontext sets agree; main matches prior simulation; StrictAll3 supports unchanged')

if __name__ == '__main__':
 main()
