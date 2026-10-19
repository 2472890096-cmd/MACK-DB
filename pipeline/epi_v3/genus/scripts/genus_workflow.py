#!/usr/bin/env python3
"""Genus input binding, evidence checks and scheduling; all EPI math is pinned V3."""
import collections
import datetime as dt
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import time
import traceback
import species_production_workflow as w

ROOT = Path(os.environ.get("MACKDB_EPI_RUN_ROOT", str(Path(__file__).resolve().parent.parent))).resolve()
METRICS=['D1','D2','Q']

def load(p):return json.loads(Path(p).read_text())

def json_write(p,v):w.atomic_text(p,json.dumps(v,indent=2,ensure_ascii=False)+'\n')

def environment(work):
    env=os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',
               VECLIB_MAXIMUM_THREADS='1',BLIS_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',
               HOME=str(ROOT/'tmp/runtime_home'),XDG_CACHE_HOME=str(ROOT/'tmp/cache'),
               TMPDIR=str(work),TMP=str(work),TEMP=str(work),LC_ALL='C.UTF-8')
    return env

def mapping(row):
    with Path(row['species_path']).open() as f:
        return dict(line.rstrip('\r\n').split('\t') for line in f if line.strip())

def names_in_key(path,row):
    ids=mapping(row);valid=set(ids.values());encoded={i+'_'+n:n for i,n in ids.items()}
    result=set()
    for line in Path(path).read_text().splitlines():
        value=line.strip()
        if not value:continue
        if value in valid:result.add(value)
        elif value in encoded:result.add(encoded[value])
        else:raise RuntimeError('Candidate is not an exact input genus name: '+value)
    return result

def validate_files(row,paths):
    if not all(p.is_file() for p in paths):raise RuntimeError('Missing required outputs: '+row['project'])
    sr=w.read_tsv(paths[3])
    if len(sr)!=1:raise RuntimeError('Summary must contain one project')
    sr=sr[0]
    if sr['project']!=row['project'] or sr['sample_total']!=row['sample_count'] or sr['species_total']!=row['taxa_count']:
        raise RuntimeError('Summary dimensions/identity mismatch')
    if sr['mode']!='prevalence' or int(sr['species_scored'])!=int(row['eligible_taxa']):
        raise RuntimeError('Scored genus count differs from eligible genus count')
    counts=[len(names_in_key(p,row)) for p in paths[:3]]
    if counts!=[int(sr[k+'_key_count']) for k in METRICS]:raise RuntimeError('Key counts disagree with summary')
    return counts

def completed(row):
    s=w.status(row)
    if s['status']!='SUCCESS' or s.get('exit_code')!='0':return False
    try:
        validate_files(row,w.outputs(row['project']))
        committed=ROOT/'tmp'/row['project']/s['attempt']/'committed_output_hashes.tsv'
        records=w.read_tsv(committed)
        if {r['path'] for r in records}!={str(p) for p in w.outputs(row['project'])}:return False
        return all(w.sha(r['path'])==r['sha256'] for r in records)
    except (OSError,ValueError,KeyError,RuntimeError):return False

def resource_usage(work):
    p=work/'resource_usage.txt'
    if not p.exists():return {}
    text=p.read_text();result={}
    for field,pattern in [('max_rss_kb',r'Maximum resident set size \(kbytes\):\s*(\d+)'),
                          ('user_seconds',r'User time \(seconds\):\s*([0-9.]+)'),
                          ('system_seconds',r'System time \(seconds\):\s*([0-9.]+)')]:
        match=re.search(pattern,text)
        if match:result[field]=match.group(1)
    return result

def setup():
    check=load(ROOT/'audit/input_validation.json')
    if not check['passed']:raise RuntimeError('Input gate failed')
    cfg=dict(created_at=w.now(),production_root=str(ROOT),core_script=str(ROOT/'scripts/get_keystone_EPI_single_safe_qecdf_v3_fast.R'),
        core_sha256=check['core_sha256'],rscript=str(ROOT/'scripts/timed_rscript.sh'),
        real_rscript=shutil.which(os.environ.get('MACKDB_RSCRIPT','Rscript')) or os.environ.get('MACKDB_RSCRIPT','Rscript'),
        mode='prevalence',prevalence_min=.25,prevalence_max=.75,network_quantile=.1,
        valid_projects=check['formal_valid'],stage1_workers=8,stage2_workers=1,per_project_timeout=None,
        genus_input_adapter='-s points to validated genus mapping; no matrix transpose or scaling',
        memory_reserve_gib=8)
    if (ROOT/'audit/config.json').exists():raise RuntimeError('Setup already exists')
    json_write(ROOT/'audit/config.json',cfg)
    for p in (ROOT/'scripts').glob('*.sh'):p.chmod(0o755)
    for p in (ROOT/'scripts').glob('*.sh'):subprocess.run(['bash','-n',str(p)],check=True)
    for row in w.manifest():w.write_tsv(w.status_path(row['project']),[dict(row,status='PENDING')],w.STATUS_FIELDS)
    w.merge(w.manifest())
    fingerprint_files=list((ROOT/'scripts').iterdir())+list((ROOT/'manifests').iterdir())+[ROOT/'audit/config.json']
    w.write_tsv(ROOT/'audit/workflow_fingerprints.tsv',
        [dict(relative_path=str(p.relative_to(ROOT)),sha256=w.sha(p)) for p in fingerprint_files if p.is_file()],['relative_path','sha256'])
    w.verify_setup()
    w.atomic_text(ROOT/'status/current_stage.txt','INPUT_VALIDATED\n')
    print('GENUS_SETUP_PASS',flush=True)

def pilots():
    w.verify_setup()
    if (ROOT/'audit/pilot_validation.json').exists():raise RuntimeError('Pilot audit already exists; no duplicate runs')
    rows=w.read_tsv(ROOT/'manifests/pilot_projects.tsv');runs=[];promote=[]
    reference=load(ROOT/'audit/species_reference_config.json')['core_script']
    local=w.config()['core_script']
    for row in rows:
        project=row['project'];directories=[]
        for repeat,core in [(1,reference),(2,local)]:
            attempt='pilot_'+str(repeat);work=ROOT/'tmp'/project/attempt;work.mkdir(parents=True,exist_ok=False)
            log=ROOT/'pilots'/project;log.mkdir(parents=True,exist_ok=True)
            out=log/('repeat'+str(repeat)+'.stdout.log');err=log/('repeat'+str(repeat)+'.stderr.log')
            env=environment(work);env['EPI_CAPTURE_CORE']=core
            cmd=[str(ROOT/'scripts/timed_rscript.sh'),str(ROOT/'scripts/capture_v3.R'),
                '-a',row['abundance_path'],'-s',row['species_path'],'-o',project,'-d',str(work),
                '--mode','prevalence','--prevalence-min','0.25','--prevalence-max','0.75','--network-quantile','0.1']
            json_write(work/'command.json',dict(command=cmd,core_source=core,core_sha256=w.sha(core)))
            started=time.monotonic();start=w.now()
            with out.open('xb') as stdout,err.open('xb') as stderr:
                proc=subprocess.run(cmd,env=env,cwd=str(work),stdout=stdout,stderr=stderr)
            elapsed=time.monotonic()-started
            record=dict(project=project,repeat=repeat,sample_count=row['sample_count'],genus_count=row['taxa_count'],
                eligible_genera=row['eligible_taxa'],start_time=start,end_time=w.now(),runtime_seconds=elapsed,exit_code=proc.returncode,
                core_path=core,**resource_usage(work))
            runs.append(record)
            if proc.returncode!=0:raise RuntimeError('Pilot failed: '+project+' repeat '+str(repeat)+'; see '+str(err))
            paths=[work/p.name for p in w.outputs(project)];counts=validate_files(row,paths)
            if not (work/'audit_scores.rds').is_file():raise RuntimeError('Raw score audit missing')
            record.update(D1_count=counts[0],D2_count=counts[1],Q_count=counts[2]);directories.append(work)
            print(w.now(),'PILOT',project,'repeat',repeat,'PASS',round(elapsed,2),'seconds',flush=True)
        compared=[p.name for p in w.outputs(project)]+['audit_scores.rds','audit_scores.tsv']
        if any(w.sha(directories[0]/name)!=w.sha(directories[1]/name) for name in compared):
            raise RuntimeError('Pilot numerical/candidate reproducibility mismatch: '+project)
        promote.append((row,directories[0],runs[-2]))
    # Only after all pilots pass, commit their validated original-core outputs once.
    for row,work,run in promote:
        project=row['project'];paths=[work/p.name for p in w.outputs(project)];counts=validate_files(row,paths)
        for src,dst in zip(paths,w.outputs(project)):w.install_output(src,dst)
        w.write_tsv(work/'committed_output_hashes.tsv',[dict(path=str(p),sha256=w.sha(p)) for p in w.outputs(project)],['path','sha256'])
        for suffix in ['stdout','stderr']:
            shutil.copy2(ROOT/'pilots'/project/('repeat1.'+suffix+'.log'),ROOT/'logs/projects'/(project+'.'+suffix+'.log'))
        s=dict(row,status='SUCCESS',start_time=run['start_time'],end_time=run['end_time'],runtime_seconds=run['runtime_seconds'],
               exit_code=0,D1_count=counts[0],D2_count=counts[1],Q_count=counts[2],error_message='',attempt='pilot_1')
        w.write_tsv(w.status_path(project),[s],w.STATUS_FIELDS)
    fields=list(dict.fromkeys(k for r in runs for k in r));w.write_tsv(ROOT/'audit/pilot_runs.tsv',runs,fields)
    proof=dict(passed=True,checked_at=w.now(),projects=[r['project'] for r in rows],repeats_per_project=2,
        raw_scores_bitwise_equal=True,candidate_files_byte_equal=True,summary_byte_equal=True,
        original_core_vs_frozen_core_sha256_equal=w.sha(reference)==w.sha(local),promoted_verified_pilots=True)
    json_write(ROOT/'audit/pilot_validation.json',proof)
    text='# Genus V3 试跑报告\n\n5 个代表项目各运行两次：第一遍 source species 正式核心，第二遍 source 字节相同的冻结副本；仅在核心执行后导出审计对象。\n\n'
    text+='D1/D2/Q 原始 R 数值对象（RDS）、逐 genus 分数及 candidate 标记、三个 key 和 summary 均逐字节相同。映射名称、样本数、genus 数、eligible 数和候选数量检查通过。详细耗时、RSS、退出码见 pilot_runs.tsv。\n\n'
    for r in runs[::2]:text+=f"- {r['project']}: samples={r['sample_count']}, genera={r['genus_count']}, eligible={r['eligible_genera']}, D1/D2/Q={r['D1_count']}/{r['D2_count']}/{r['Q_count']}。\n"
    text+='\n第一遍已验证结果按正式提交流程保留，后续全量跳过这些项目。固定 V3 核心和 production worker 未作修改。\n'
    w.atomic_text(ROOT/'audit/PILOT_REPORT.md',text);w.merge(w.manifest());print('ALL_PILOTS_PASS',json.dumps(proof),flush=True)

def mem_available():
    values={line.split(':')[0]:int(line.split()[1])*1024 for line in Path('/proc/meminfo').read_text().splitlines() if ':' in line and line.split()[1].isdigit()}
    return values['MemAvailable']

def estimate_memory(row):
    n=int(row['sample_count']);g=int(row['taxa_count'])
    return max(600*1024**2,192*n*n+32*n*g)

def process_memory(pid):
    total=0;todo=[pid];seen=set()
    while todo:
        p=todo.pop()
        if p in seen:continue
        seen.add(p)
        try:
            status=Path('/proc/'+str(p)+'/status').read_text()
            match=re.search(r'^VmRSS:\s*(\d+)',status,re.M)
            if match:total+=int(match.group(1))*1024
            todo.extend(int(x) for x in Path('/proc/'+str(p)+'/task/'+str(p)+'/children').read_text().split())
        except OSError:pass
    return total

def phase(rows,workers,active):
    pending=collections.deque();skipped=0
    for row in rows:
        if completed(row):skipped+=1
        else:
            s=w.status(row)
            if s['status']=='RUNNING' and any(w.project_alive(s.get(k,''),row['project']) for k in ['worker_pid','r_pid']):
                raise RuntimeError('Live worker prevents overlapping resume: '+row['project'])
            if s['status']=='SUCCESS':raise RuntimeError('SUCCESS integrity failed; preserved files require review: '+row['project'])
            pending.append(row)
    print(w.now(),'PHASE_QUEUE',len(pending),'SKIP_VERIFIED_SUCCESS',skipped,'MAX_WORKERS',workers,flush=True)
    last_report=0;peak=0
    while pending or active:
        for pid,info in list(active.items()):
            proc,row=info['proc'],info['row'];rss=process_memory(pid);info['rss']=rss;info['peak']=max(info['peak'],rss)
            rc=proc.poll()
            if rc is None:continue
            s=w.status(row)
            if rc!=0 or not completed(row):
                try:os.killpg(pid,signal.SIGTERM)
                except ProcessLookupError:pass
                s.update(status='FAIL',end_time=w.now(),exit_code=rc,
                    error_message=s.get('error_message') or 'Worker exited without complete validated SUCCESS')
                w.write_tsv(w.status_path(row['project']),[s],w.STATUS_FIELDS)
            work=ROOT/'tmp'/row['project']/s.get('attempt','')
            usage=dict(project=row['project'],status=s['status'],runtime_seconds=s.get('runtime_seconds',''),
                process_tree_peak_rss_bytes=info['peak'],**resource_usage(work))
            json_write(ROOT/'status/projects'/(row['project']+'.resources.json'),usage)
            del active[pid]
            print(w.now(),'PROJECT_FINISHED',row['project'],s['status'],s.get('runtime_seconds',''),flush=True)
        while pending and len(active)<workers:
            row=pending[0];future=sum(max(0,estimate_memory(i['row'])-i['rss']) for i in active.values())
            if mem_available()<8*1024**3+future+estimate_memory(row):break
            pending.popleft()
            cmd=[str(ROOT/'scripts/run_one_project.sh'),row['project'],row['abundance_path'],row['species_path'],row['sample_count']]
            proc=subprocess.Popen(cmd,cwd=str(ROOT),env=environment(ROOT/'tmp'),start_new_session=True)
            active[proc.pid]=dict(proc=proc,row=row,rss=0,peak=0)
            print(w.now(),'PROJECT_START',row['project'],'worker_pid',proc.pid,flush=True)
        peak=max(peak,sum(i['rss'] for i in active.values()))
        if time.monotonic()-last_report>=60:
            json_write(ROOT/'status/progress.json',dict(at=w.now(),remaining_in_phase=len(pending),running=len(active),
                available_memory_bytes=mem_available(),observed_active_peak_rss_bytes=peak,memory_wait=bool(pending and len(active)<workers)))
            print(w.now(),'PROGRESS','pending',len(pending),'active',len(active),'available_GiB',round(mem_available()/1024**3,2),flush=True)
            last_report=time.monotonic()
        if pending or active:time.sleep(2)
    w.merge(w.manifest())

def finalize(rows):
    w.merge(rows,finalize=True)
    official=w.read_tsv(ROOT/'manifests/formal_project_coverage.tsv');valid={r['project']:r for r in rows}
    coverage=[];records=[];metric_records=[]
    for row in official:
        project=row['project'];s=w.status(row) if project in valid else dict(status=row['input_status'])
        item=dict(project=project,environment=row['environment'],input_status=row['input_status'],computation_status=s['status'],
            sample_count=row['sample_count'],genus_count=row['taxa_count'],eligible_genera=row['eligible_taxa'],
            D1_count='',D2_count='',Q_count='',final_2of3_count='',zero_all_metrics='',zero_final_candidates='',
            runtime_seconds=s.get('runtime_seconds',''),max_rss_kb='',error_message=s.get('error_message',row.get('validation_error','')))
        if s['status']=='SUCCESS':
            if not completed(row):raise RuntimeError('Postprocess integrity mismatch: '+project)
            sets=[names_in_key(p,row) for p in w.outputs(project)[:3]]
            votes=collections.Counter(n for group in sets for n in group);chosen=sorted(n for n,c in votes.items() if c>=2)
            for k,group in zip(METRICS,sets):item[k+'_count']=len(group)
            item.update(final_2of3_count=len(chosen),zero_all_metrics=not any(sets),zero_final_candidates=not chosen)
            usage=resource_usage(ROOT/'tmp'/project/s['attempt']);item['max_rss_kb']=usage.get('max_rss_kb','')
            env=row['environment']
            if not env or env.lower()=='unknown' or ';' in env:raise RuntimeError('Invalid final environment')
            records.extend(dict(name=n,bioproject=project,environment=env) for n in chosen)
            for n in sorted(votes):metric_records.append(dict(name=n,bioproject=project,environment=env,
                D1_support=n in sets[0],D2_support=n in sets[1],Q_support=n in sets[2],support_count=votes[n],final_candidate=n in chosen))
        coverage.append(item)
    failures=[r for r in coverage if r['computation_status']!='SUCCESS']
    records.sort(key=lambda r:(r['bioproject'],r['name']))
    if len({(r['name'],r['bioproject']) for r in records})!=len(records):raise RuntimeError('Duplicate project/genus final records')
    for label,items in [('project_coverage',coverage),('successful_projects',[r for r in coverage if r['computation_status']=='SUCCESS']),
                        ('failed_or_missing_projects',failures),('zero_candidate_projects',[r for r in coverage if r['computation_status']=='SUCCESS' and r['zero_final_candidates']])]:
        w.write_tsv(ROOT/'audit'/(label+'.tsv'),items,list(coverage[0]))
    final=ROOT/('key_genus_epi_partial_1575.tsv' if failures else 'key_genus_epi_final_1575.tsv')
    w.write_tsv(final,records,['name','bioproject','environment'])
    comparison=ROOT/'comparisons';comparison.mkdir(exist_ok=True)
    w.write_tsv(comparison/'key_genus_epi_v3_standardized.tsv',records,['name','bioproject','environment'])
    w.write_tsv(comparison/'genus_metric_support.tsv',metric_records,
        ['name','bioproject','environment','D1_support','D2_support','Q_support','support_count','final_candidate'])
    stats=dict(formal_projects=len(official),valid_inputs=len(rows),success_projects=len(coverage)-len(failures),failed_or_missing_projects=len(failures),
        zero_final_candidate_projects=sum(r['computation_status']=='SUCCESS' and r['zero_final_candidates'] for r in coverage),
        zero_all_metric_projects=sum(r['computation_status']=='SUCCESS' and r['zero_all_metrics'] for r in coverage),
        D1_candidates=sum(r['D1_count'] for r in coverage if r['computation_status']=='SUCCESS'),
        D2_candidates=sum(r['D2_count'] for r in coverage if r['computation_status']=='SUCCESS'),
        Q_candidates=sum(r['Q_count'] for r in coverage if r['computation_status']=='SUCCESS'),
        final_project_genus_records=len(records),projects_with_final_candidates=len({r['bioproject'] for r in records}),unique_genera=len({r['name'] for r in records}))
    json_write(ROOT/'audit/final_statistics.json',stats)
    comparison_rows=[dict(metric=k,old_value=old,new_value=stats[k],difference=stats[k]-old) for k,old in
        [('final_project_genus_records',3425),('projects_with_final_candidates',990),('unique_genera',1592)]]
    w.write_tsv(comparison/'legacy_count_comparison.tsv',comparison_rows,['metric','old_value','new_value','difference'])
    report='# Genus EPI V3 最终报告\n\n状态：'+('存在失败/缺失，仅提供 partial TSV。' if failures else '1575 个正式项目全部验证成功。')+'\n\n'
    report+='复用 species production 的固定 V3 核心和 worker；未修改任何 EPI 公式、prevalence 或阈值。完整参数与 SHA256 见 audit/config.json 和 protected_fingerprints.tsv；输入及试跑报告已保存。\n\n'
    report+='\n'.join('- '+k+': '+str(v) for k,v in stats.items())+'\n\n'
    report+='零候选只来自 exit=0、四份输出完整且经 summary/候选名称/提交哈希验证的 SUCCESS。失败和缺失输入另列，不视为零候选。最终名称按输入映射核验，数字、连字符和 GTDB 后缀保留；所有记录有唯一非 unknown 环境。\n\n'
    report+='旧版 3425 条、990 项目、1592 unique genera 仅作汇总对照，不作为目标。comparisons 中提供标准化结果和逐指标支持表；未读取或修改 AWS MySQL，也未用旧结果补缺。\n'
    w.atomic_text(ROOT/'audit/FINAL_REPORT.md',report)
    w.postprocess(rows)
    w.atomic_text(ROOT/'status/GENUS_DELIVERABLES_FINISHED.ok',w.now()+'\n')
    print(w.now(),'GENUS_FINAL_STATISTICS',json.dumps(stats),flush=True)

def master():
    lock=(ROOT/'status/master.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    active={}
    def stop(signum,frame):raise InterruptedError('Master signal '+str(signum))
    signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
    try:
        rows=w.verify_setup()
        if not load(ROOT/'audit/pilot_validation.json')['passed']:raise RuntimeError('Pilot gate not passed')
        for row in rows:
            s=w.status(row)
            if s['status']=='RUNNING' and any(w.project_alive(s.get(k,''),row['project']) for k in ['worker_pid','r_pid']):
                raise RuntimeError('Existing live project prevents resume: '+row['project'])
        w.atomic_text(ROOT/'status/master.pid',str(os.getpid())+'\n')
        json_write(ROOT/'status/master_identity.json',dict(pid=os.getpid(),started=w.now(),proc_stat=Path('/proc/self/stat').read_text()))
        archive=ROOT/'status/previous_runs'/dt.datetime.now().strftime('%Y%m%dT%H%M%S')
        for name in ['MASTER_FAILED.flag','ALL_FINISHED.ok','FINISHED_WITH_FAILURES.flag','GENUS_DELIVERABLES_FINISHED.ok',
                     'STAGE1_FINISHED','STAGE2_FINISHED','POSTPROCESS_FINISHED']:
            p=ROOT/'status'/name
            if p.exists():archive.mkdir(parents=True,exist_ok=True);p.rename(archive/name)
        w.set_stage('STAGE1_RUNNING');phase([r for r in rows if int(r['sample_count'])<=5000],8,active)
        (ROOT/'status/STAGE1_RUNNING').unlink(missing_ok=True);w.set_stage('STAGE1_FINISHED')
        w.set_stage('STAGE2_RUNNING');phase([r for r in rows if int(r['sample_count'])>5000],1,active)
        (ROOT/'status/STAGE2_RUNNING').unlink(missing_ok=True);w.set_stage('STAGE2_FINISHED')
        finalize(rows)
        return 0
    except BaseException:
        detail=traceback.format_exc()
        try:print(detail,file=sys.stderr,flush=True)
        except OSError:pass
        for pid,i in active.items():
            try:os.killpg(pid,signal.SIGTERM)
            except ProcessLookupError:pass
        try:
            w.atomic_text(ROOT/'status/MASTER_FAILED.flag',w.now()+'\n'+detail)
            w.atomic_text(ROOT/'status/current_stage.txt','MASTER_FAILED\n')
        except OSError:pass
        return 1
    finally:lock.close()

if __name__=='__main__':
    mode=sys.argv[1]
    if mode=='setup':setup()
    elif mode=='pilots':
        try:pilots()
        except BaseException:
            json_write(ROOT/'audit/pilot_validation.json',dict(passed=False,error=traceback.format_exc()))
            raise
    elif mode=='master':sys.exit(master())
    elif mode=='status':w.check_status()
    elif mode=='finalize':finalize(w.manifest())
    else:raise SystemExit('Unknown mode')
