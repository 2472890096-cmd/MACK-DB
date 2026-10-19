from pathlib import Path
import csv,json
import numpy as np
from scipy.stats import false_discovery_control
R=Path(__file__).resolve().parents[1]
p=R/'07_statistics/stage5A_null_summary.tsv'
with p.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
for r in rows:
 if r['status']=='VALID':r['raw_empirical_p']=r['empirical_p_exact'];r['empirical_p']=max(.001,float(r['empirical_p_exact']))
for model in ['NullA','NullB']:
 valid=[r for r in rows if r['null_model']==model and r['status']=='VALID'];adjusted=false_discovery_control(np.array([float(r['raw_empirical_p']) for r in valid]),method='bh')
 for r,q in zip(valid,adjusted):r['BH_adjusted_p']=float(q)
with p.open('w') as f:
 w=csv.DictWriter(f,rows[0].keys(),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
p=R/'00_input_audit/stage5A_null_protocol.json';j=json.loads(p.read_text());j['empirical_p']='Exact (1+n_null>=observed)/1001 stored in raw_empirical_p and empirical_p_exact; empirical_p is presentation value with minimum 0.001';j['BH']='separate NullA and NullB families, all valid primary comparisons; use unrounded raw_empirical_p';p.write_text(json.dumps(j,indent=2)+'\n')
print('Exact empirical P values retained; separate BH corrections use unrounded values. Display minimum P = 0.001.')
