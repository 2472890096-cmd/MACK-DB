#!/usr/bin/env python3
"""Validate labelled tab-delimited correlation/P matrices and emit CSV for the unchanged R reader.
No reordering, transposition, inference, filtering, scaling, or numeric rewriting.
"""
from pathlib import Path
import argparse,csv,math

def convert(corr,pvalue,mapping,out):
 with mapping.open() as f:mp=list(csv.reader(f,delimiter='\t'))
 if any(len(r)!=2 for r in mp) or not mp:raise ValueError('Mapping must have two columns')
 ids=[r[0] for r in mp]
 if len(set(ids))!=len(ids):raise ValueError('Duplicate mapping IDs')
 matrices=[]
 for path in [corr,pvalue]:
  with path.open() as f:rows=list(csv.reader(f,delimiter='\t'))
  if len(rows)!=len(ids)+1 or rows[0][1:]!=ids:raise ValueError('Matrix must have a header with the exact ordered mapping IDs')
  if [r[0] for r in rows[1:]]!=ids or any(len(r)!=len(ids)+1 for r in rows):raise ValueError('Row/column IDs or dimensions mismatch; no automatic reorder')
  for row in rows[1:]:
   for value in row[1:]:
    if value not in ('NA','NaN','nan') and not math.isfinite(float(value)):raise ValueError('Unsupported nonfinite cell')
  matrices.append(rows)
 for i in range(1,len(ids)+1):
  for j in range(i+1,len(ids)+1):
   missing=[m[i][j] in ('NA','NaN','nan') for m in matrices]
   if missing[0]!=missing[1]:raise ValueError('Correlation/P upper-triangle NA masks differ')
 out.mkdir(parents=True,exist_ok=False)
 for name,rows in zip(['correlation.csv','pvalues.csv'],matrices):
  with (out/name).open('w') as f:
   w=csv.writer(f,lineterminator='\n');w.writerow(['']+ids);w.writerows(rows[1:])

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for flag in ['corr','pvalue','mapping','output']:p.add_argument('--'+flag,type=Path,required=True)
 a=p.parse_args();convert(a.corr,a.pvalue,a.mapping,a.output)
if __name__=='__main__':main()
