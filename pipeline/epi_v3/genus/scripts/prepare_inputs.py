#!/usr/bin/env python3
"""Validate real genus inputs and pin unmodified species production code."""
import csv
import datetime as dt
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import sys
import numpy as np

ROOT = Path(os.environ.get("MACKDB_EPI_RUN_ROOT", str(Path(__file__).resolve().parent.parent))).resolve()
REFERENCE = Path(os.environ.get('MACKDB_EPI_SPECIES_REFERENCE', str(Path(__file__).resolve().parents[2] / 'species'))).resolve()
INPUT = Path(os.environ.get('MACKDB_GENUS_INPUT', 'data/genus')).resolve()
FORMAL = Path(os.environ.get('MACKDB_PROJECT_LIST', str(Path(__file__).resolve().parents[4] / 'docs/scopes/mackdb_1575_bioprojects.txt'))).resolve()
ENVDIR = Path(os.environ.get('MACKDB_ENV_DIR', 'data/env')).resolve()
EXPECTED = '98d832663fea152e0d63c481106ff24fb5faca4e34ce332d0187814cbead5cc2'

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def tsv(p, rows, fields):
    with Path(p).open('w') as f:
        wr=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n',extrasaction='ignore')
        wr.writeheader();wr.writerows(rows)

def matrix(project, official, environments):
    a=INPUT/'abundance'/(project+'_abundance.matrix');g=INPUT/'genus'/(project+'_genus.matrix')
    r=dict(project=project,official=official,abundance_path=str(a),genus_path=str(g),species_path=str(g),
           sample_count=0,taxa_count=0,eligible_taxa=0,environment=';'.join(environments),
           input_status='VALID',validation_error='',mapping_reordered=False)
    try:
        if not a.is_file() or not g.is_file():
            r['input_status']='MISSING_INPUT';raise ValueError('Missing abundance or genus mapping')
        for kind,p in [('abundance',a),('species',g)]:
            before=p.stat();b=p.read_bytes();after=p.stat()
            if (before.st_mtime_ns,before.st_size)!=(after.st_mtime_ns,after.st_size):raise ValueError('Input changed while reading')
            r[kind+'_sha256']=hashlib.sha256(b).hexdigest();r[kind+'_bytes']=after.st_size;r[kind+'_mtime_ns']=after.st_mtime_ns
            if kind=='abundance': abundance=b.decode('utf-8')
            else: mapping=b.decode('utf-8')
        lines=abundance.splitlines()
        if not lines:raise ValueError('Empty abundance file')
        header=lines[0].split('\t');n=len(header)-1
        if n<=0 or len(lines)<=1:raise ValueError('Empty matrix: no samples/taxa')
        if len(set(header[1:]))!=n:raise ValueError('Duplicate sample column identifiers')
        ids=[];sums=np.zeros(n);eligible=0;max_value=0.0
        for line in lines[1:]:
            if not line.strip():continue
            parts=line.split('\t')
            if len(parts)!=n+1:raise ValueError('Nonrectangular abundance matrix')
            ids.append(parts[0]);v=np.asarray(parts[1:],dtype=np.float64)
            if not np.isfinite(v).all():raise ValueError('Nonfinite abundance values')
            if (v<0).any():raise ValueError('Negative abundance values')
            sums+=v;max_value=max(max_value,float(v.max()))
            prev=np.count_nonzero(v>0)/n;eligible+=int(.25<=prev<=.75)
        if not ids:raise ValueError('Empty abundance taxon dimension')
        if len(set(ids))!=len(ids):raise ValueError('Duplicate abundance taxon identifiers')
        if ids!=[str(i) for i in range(len(ids))]:raise ValueError('Abundance identifiers are not ordered 0..N-1')
        pairs=[x.split('\t') for x in mapping.splitlines() if x.strip()]
        if any(len(x)!=2 for x in pairs):raise ValueError('Mapping must have exactly ID and genus-name columns')
        if len(pairs)!=len(ids):raise ValueError('Genus mapping count differs from abundance rows')
        names={}
        for taxid,name in pairs:
            if taxid in names:raise ValueError('Duplicate mapping identifier')
            if not name or name!=name.strip():raise ValueError('Empty or whitespace-padded genus name')
            names[taxid]=name
        if set(names)!=set(ids):raise ValueError('Missing or extra mapping identifiers')
        if len(set(names.values()))!=len(names):raise ValueError('Duplicate genus names')
        if [x[0] for x in pairs]!=ids:
            adapted=ROOT/'adapted_mappings';adapted.mkdir(exist_ok=True)
            path=adapted/g.name;path.write_text(''.join(i+'\t'+names[i]+'\n' for i in ids))
            r.update(species_path=str(path),mapping_reordered=True,species_sha256=sha(path),
                     species_bytes=path.stat().st_size,species_mtime_ns=path.stat().st_mtime_ns)
        if official and (len(environments)!=1 or environments[0].lower()=='unknown'):
            raise ValueError('Formal project requires exactly one non-unknown environment')
        r.update(sample_count=n,taxa_count=len(ids),eligible_taxa=eligible,column_sum_min=float(sums.min()),
                 column_sum_median=float(np.median(sums)),column_sum_max=float(sums.max()),
                 zero_sum_samples=int((sums==0).sum()),maximum_abundance=max_value,
                 ids_zero_based=True,names_unique=True,numeric_finite_nonnegative=True)
    except (ValueError,OSError,UnicodeError) as exc:
        if r['input_status']!='MISSING_INPUT':r['input_status']='INVALID_INPUT'
        r['validation_error']=str(exc)
    return r

def main():
    # Prepare only an explicitly selected isolated runtime workspace.
    if 'MACKDB_EPI_RUN_ROOT' not in os.environ:
        raise RuntimeError('Set MACKDB_EPI_RUN_ROOT to a new isolated run directory')
    for directory in ['audit','scripts','manifests','status/projects','logs/projects','logs/master','keys','summaries','tmp','pilots','final_2of3']:
        (ROOT/directory).mkdir(parents=True,exist_ok=True)
    source_dir=Path(__file__).resolve().parent
    for name in ['genus_workflow.py','run_one_project.sh','timed_rscript.sh','capture_v3.R']:
        dst=ROOT/'scripts'/name
        if source_dir/name != dst:
            if dst.exists() and sha(dst)!=sha(source_dir/name): raise RuntimeError('Existing runtime script differs: '+name)
            if not dst.exists(): shutil.copy2(source_dir/name,dst)
    gate=ROOT/'audit/input_validation.json'
    if gate.exists():raise RuntimeError('Validation already exists; preserve current audit')
    formal=[x.strip() for x in FORMAL.read_text().splitlines() if x.strip()]
    if len(formal)!=1575 or len(set(formal))!=1575:raise RuntimeError('Formal list must contain 1575 unique projects')
    cfg=json.loads((REFERENCE/'audit/config.json').read_text());core=Path(cfg['core_script'])
    if not core.is_absolute(): core=(REFERENCE/core).resolve()
    cfg['core_script']=str(core)
    if cfg['core_sha256']!=EXPECTED or sha(core)!=EXPECTED:raise RuntimeError('Cannot establish verified V3 source')
    frozen=ROOT/'scripts/get_keystone_EPI_single_safe_qecdf_v3_fast.R'
    shutil.copy2(core,frozen)
    shutil.copy2(REFERENCE/'scripts/workflow.py',ROOT/'scripts/species_production_workflow.py')
    (ROOT/'audit/species_reference_config.json').write_text(json.dumps(cfg,indent=2)+'\n')
    protected=[dict(role='source_v3_core',path=str(core),sha256=sha(core)),
        dict(role='frozen_v3_core',path=str(frozen),sha256=sha(frozen)),
        dict(role='source_production_workflow',path=str(REFERENCE/'scripts/workflow.py'),sha256=sha(REFERENCE/'scripts/workflow.py')),
        dict(role='copied_production_workflow',path=str(ROOT/'scripts/species_production_workflow.py'),sha256=sha(ROOT/'scripts/species_production_workflow.py')),
        dict(role='formal_project_list',path=str(FORMAL),sha256=sha(FORMAL))]
    env={}
    for p in sorted(ENVDIR.glob('*.list')):
        protected.append(dict(role='environment_mapping',path=str(p),sha256=sha(p)))
        for name in set(x.strip() for x in p.read_text().splitlines() if x.strip()):env.setdefault(name,[]).append(p.stem)
    tsv(ROOT/'audit/protected_fingerprints.tsv',protected,['role','path','sha256'])
    aids={p.name[:-len('_abundance.matrix')] for p in (INPUT/'abundance').glob('*_abundance.matrix')}
    gids={p.name[:-len('_genus.matrix')] for p in (INPUT/'genus').glob('*_genus.matrix')}
    census=[]
    for i,p in enumerate(sorted(set(formal)|aids|gids),1):
        census.append(matrix(p,p in formal,env.get(p,[])))
        if i%250==0:print('Actual genus inputs validated',i,flush=True)
    fields=list(dict.fromkeys(k for r in census for k in r))
    tsv(ROOT/'audit/input_compatibility.tsv',census,fields)
    official=[r for r in census if r['official']]
    valid=sorted([r for r in official if r['input_status']=='VALID'],key=lambda r:(r['sample_count'],r['project']))
    invalid=[r for r in official if r['input_status']!='VALID'];extras=[r for r in census if not r['official']]
    tsv(ROOT/'manifests/formal_project_coverage.tsv',official,fields)
    tsv(ROOT/'manifests/excluded_nonformal_projects.tsv',extras,fields)
    tsv(ROOT/'manifests/missing_or_invalid_inputs.tsv',invalid,fields)
    tsv(ROOT/'manifests/all_valid_projects.tsv',valid,fields)
    for name,pred in [('le5000',lambda r:r['sample_count']<=5000),('gt5000',lambda r:r['sample_count']>5000)]:
        selected=[r for r in valid if pred(r)]
        tsv(ROOT/('manifests/projects_'+name+'.tsv'),selected,fields)
        (ROOT/('manifests/projects_'+name+'.list')).write_text(''.join(r['project']+'\n' for r in selected))
    # Required projects plus low-, middle-, high-genus projects with modest sample counts.
    candidates=sorted([r for r in valid if 30<=r['sample_count']<=150 and r['eligible_taxa']>0],key=lambda r:(r['taxa_count'],r['project']))
    pilot_ids=['PRJDB10010','PRJDB10829']
    for index in [0,len(candidates)//2,len(candidates)-1]:
        p=candidates[index]['project']
        if p not in pilot_ids:pilot_ids.append(p)
    pilots=[next(r for r in valid if r['project']==p) for p in pilot_ids]
    tsv(ROOT/'manifests/pilot_projects.tsv',pilots,fields)
    result=dict(checked_at=dt.datetime.now().astimezone().isoformat(),formal_projects=len(formal),abundance_files=len(aids),
        mapping_files=len(gids),formal_valid=len(valid),formal_invalid_or_missing=len(invalid),
        excluded_nonformal=[r['project'] for r in extras],mapping_reorders=sum(bool(r['mapping_reordered']) for r in valid),
        passed=not invalid,core_sha256=EXPECTED,pilot_projects=pilot_ids,
        le5000=sum(r['sample_count']<=5000 for r in valid),gt5000=sum(r['sample_count']>5000 for r in valid))
    gate.write_text(json.dumps(result,indent=2)+'\n')
    report=f'''# Genus EPI V3 输入兼容性报告

正式项目 {len(formal)}；原始 abundance {len(aids)}、genus mapping {len(gids)}。正式有效输入 {len(valid)}；缺失/无效 {len(invalid)}。排除非正式项目：{', '.join(r['project'] for r in extras)}。

核心来源为 species production audit/config.json 指向的脚本，SHA256 `{EXPECTED}`。本目录核心副本逐字节相同；完整 production workflow.py 同样原样保存为 species_production_workflow.py。核心没有硬编码 species 路径或要求 `_species.matrix` 后缀，`-s` 读取两列映射的第二列。因此 genus 输入通过 `-s <project>_genus.matrix` 使用同一读取函数；仅参数名及原始 summary 中 species_total/species_scored 字段保留旧称，实际代表 genus。

逐文件验证的实际方向为 genus 行 × sample 列：首行为 sample 列名，第一列 OTU_id；abundance 行号和 genus 映射号为从 0 开始连续编号，维度、唯一性和逐行对应关系均检查。没有自动转置。检查覆盖空文件、缺失/重复编号、重复 genus 名、非数值/非有限/负值；每个项目详细结果见 input_compatibility.tsv。必要映射重排仅写新目录，本次重排数量 {result['mapping_reorders']}。

丰度为输入提供的非负 abundance 数值；逐列总量的 min/median/max 和零和列数均记录在 TSV，用于辨识百分比丰度及舍入误差。保持原值输入，不预先改变尺度。与 species 正式流程完全一致，仅由固定核心在删除 focal genus 后对剩余群落按列归一化；D2 对 present/absent 组也按原实现处理。

prevalence 两端包含 0.25 和 0.75；D1 为删除 focal 后 present–absent 跨组 Bray–Curtis 距离均值；D2 复用原分组归一化及组均值 Bray–Curtis；Q 复用唯一非对角 sample pairs 的 ECDF、pQ=0.1、首跳回退及原 modularity。候选背景仅有限分数，使用 R sample SD 和严格 `score > mean + 2*sd`，其余边界行为不改。

原始 genus 名称从映射第二列原样输出，数字、连字符及 GTDB 后缀不删除。最终整理只接受映射中确实存在的名称；若遇编号前缀，仅按确切 `ID_名称` 映射解码，禁止通用正则剥离正常名称。

environment 来自 env/*.list 的文件名，正式项目要求唯一且非 unknown。最终仅正式名单进入结果；原始额外项目单独报告。试跑项目：{', '.join(pilot_ids)}。状态门禁：{'PASS' if not invalid else 'BLOCKED'}；出现正式无效输入不会启动全量。
'''
    (ROOT/'audit/INPUT_COMPATIBILITY_REPORT.md').write_text(report)
    print(json.dumps(result,indent=2),flush=True)
    print('PILOTS',[(r['project'],r['sample_count'],r['taxa_count'],r['eligible_taxa']) for r in pilots],flush=True)
    if invalid:raise SystemExit('Formal input validation failed; full run prohibited')

if __name__=='__main__':main()
