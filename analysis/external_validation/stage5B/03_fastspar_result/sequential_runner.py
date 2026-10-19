"""Run the original MACK-DB shell script sequentially; never copy or edit it.

This supervisor survives a disconnected client, records child PIDs/exit codes,
and applies only the requested single-process scheduling via GNU parallel's
documented environment options. Statistical commands remain in the original script.
"""

# v0.1: repository-local IO configuration; scientific functions below are unchanged.
import sys as _path_sys
from pathlib import Path as _RepoPath
_repo_root = next(p for p in _RepoPath(__file__).resolve().parents if (p / "lib/mackdb_paths.py").is_file())
_path_sys.path.insert(0, str(_repo_root / "lib"))
from mackdb_paths import resolve_path
from pathlib import Path
import csv, datetime, fcntl, hashlib, json, os, shlex, signal, subprocess, sys, time
import numpy as np

ROOT=Path(__file__).resolve().parent
STAGE=ROOT.parent; PREP=STAGE/'02_prepare_fastspar_input'
INPUT=PREP/'count_matrix'; LOG=ROOT/'logs'; QC=ROOT/'qc'
SCRIPT=Path(str(resolve_path('fastspar_script')))
NETWORKS=[('Guo_Control',91,377),('Guo_MECFS',106,377),('Guo_Common',197,377),('Wallen_PD',490,598),('Wallen_NHC',234,598)]
STATUS_FIELDS='network input start_time end_time PID exit_code status correlation_file pvalue_file qc_status covariance_file output_prefix log_file runtime_seconds nproc main_threads command bootstrap_joblog failure_reason'.split()
STATE=[]; OUTPUT_QC=[]; PREFLIGHT=[]

def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def readtsv(p):
    with Path(p).open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def table(path,rs,fields=None):
    path=Path(path);assert ROOT in path.resolve().parents
    rs=list(rs);temp=path.with_name(path.name+'.tmp')
    with temp.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
    os.replace(temp,path)
def append(path,row):
    new=not path.exists()
    with path.open('a',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(row),delimiter='\t',lineterminator='\n')
        if new:w.writeheader()
        w.writerow(row)
def event(command,inp='',out='',code='',pid='',start='',end='',runtime=''):
    append(LOG/'execution_log.tsv',dict(command=command,PID=pid,start_time=start or now(),finish_time=end,input=str(inp),output=str(out),exit_code=code,runtime_seconds=runtime))
def save_status():table(QC/'stage5B1B_run_status.tsv',STATE,STATUS_FIELDS)
def processes():
    # Required process inspection; return only FastSpar binaries and actual original-pipeline bash processes.
    p=subprocess.run(['ps','-eo','pid=,ppid=,pgid=,stat=,pcpu=,pmem=,comm=,args='],capture_output=True,text=True,check=True)
    selected=[]
    for line in p.stdout.splitlines():
        parts=line.strip().split(None,7)
        if len(parts)<8:continue
        pid,ppid,pgid,state,cpu,mem,comm,args=parts
        tokens=args.split();exe=Path(tokens[0]).name if tokens else ''
        kind=''
        if exe in ['fastspar','fastspar_bootstrap','fastspar_pvalues']:kind=exe
        elif exe=='bash' and len(tokens)>1 and tokens[1]==str(SCRIPT):kind='original_pipeline'
        if kind and not state.startswith('Z'):selected.append(dict(pid=int(pid),ppid=int(ppid),pgid=int(pgid),state=state,cpu=cpu,mem=mem,kind=kind,args=args))
    return selected
def record_processes(network,items):
    append(LOG/'process_monitor.tsv',dict(timestamp=now(),network=network,main_fastspar_count=sum(p['kind']=='fastspar' for p in items),pipeline_count=sum(p['kind']=='original_pipeline' for p in items),processes=json.dumps(items,separators=(',',':'))))
def wait_idle(network):
    while True:
        p=processes();record_processes(network,p)
        if not p:return
        (LOG/'heartbeat.json').write_text(json.dumps(dict(timestamp=now(),network=network,state='WAITING_FOR_EXISTING_FASTSPAR',processes=p),indent=2)+'\n')
        time.sleep(10)
def count_info(p):
    with p.open() as f:
        header=f.readline().rstrip('\r\n').split('\t');data=[line.rstrip('\r\n').split('\t') for line in f]
    return header,data
def preflight(network,samples,taxa,baseline):
    p=INPUT/(network+'_count.matrix');cpus=int(subprocess.check_output(['nproc'],text=True).strip())
    event('nproc',network,cpus,0)
    header,data=count_info(p)
    ok=(p.stat().st_size>0 and header[0]=='OTU_id' and len(header)-1==samples and len(data)==taxa and all(len(r)==samples+1 for r in data)
        and len({r[0] for r in data})==taxa and len(set(header[1:]))==samples and sha(p)==baseline[str(p)]['sha256'])
    PREFLIGHT.append(dict(network=network,input=str(p),exists=p.is_file(),size_bytes=p.stat().st_size,samples=len(header)-1,expected_samples=samples,taxa=len(data),expected_taxa=taxa,first_header=header[0],nproc=cpus,status='PASS' if ok else 'FAIL',timestamp=now(),sha256=sha(p)))
    table(QC/'stage5B1B_preflight_qc.tsv',PREFLIGHT)
    assert ok,'Input count-matrix preflight mismatch'
    return cpus,[r[0] for r in data]
def bootstrap_summary(path):
    if not path.is_file():return dict(rows=0,failures=0,max_simultaneous=0,available=False)
    rs=readtsv(path);intervals=[]
    for r in rs:
        a=float(r['Starttime']);b=a+float(r['JobRuntime']);intervals.extend([(a,1),(b,-1)])
    live=maximum=0
    for _,delta in sorted(intervals,key=lambda x:(x[0],x[1])):live+=delta;maximum=max(maximum,live)
    return dict(rows=len(rs),failures=sum(r['Exitval']!='0' or r['Signal']!='0' for r in rs),max_simultaneous=maximum,available=True)
def check_outputs(rec,expected_taxa):
    job=bootstrap_summary(Path(rec['bootstrap_joblog']));errors=[];results=[]
    for kind,key in [('correlation','correlation_file'),('covariance','covariance_file'),('pvalues','pvalue_file')]:
        p=Path(rec[key]);r=dict(network=rec['network'],matrix_type=kind,file=str(p),exists=p.is_file(),bytes=p.stat().st_size if p.is_file() else 0,
            expected_taxa=len(expected_taxa),row_taxa='',column_taxa='',taxa_names_match='NO',duplicate_taxa='',NA_or_nonfinite='',max_symmetry_difference='',symmetric='NO',
            bootstrap_jobs=job['rows'],bootstrap_failed_jobs=job['failures'],bootstrap_max_simultaneous=job['max_simultaneous'],script_exit_code=rec['exit_code'],status='FAIL',error='')
        try:
            assert p.is_file() and p.stat().st_size>0,'Missing or empty output'
            with p.open() as f:
                head=f.readline().rstrip('\r\n').split('\t');rr=[l.rstrip('\r\n').split('\t') for l in f]
            names=[v[0] for v in rr];cols=head[1:];arr=np.array([v[1:] for v in rr],dtype=float)
            r.update(row_taxa=len(names),column_taxa=len(cols),duplicate_taxa=(len(names)-len(set(names)))+(len(cols)-len(set(cols))))
            assert names==expected_taxa and cols==expected_taxa,'Taxa identity/order mismatch with input count matrix'
            r['taxa_names_match']='YES'
            assert arr.shape==(len(expected_taxa),len(expected_taxa)),'Wrong matrix shape'
            nonfinite=int((~np.isfinite(arr)).sum());r['NA_or_nonfinite']=nonfinite;assert nonfinite==0,'NA/nonfinite values'
            delta=float(np.max(np.abs(arr-arr.T)));r['max_symmetry_difference']=delta;r['symmetric']='YES' if delta<=1e-12 else 'NO'
            # User requires numerical symmetry for correlation and pvalues, not covariance.
            # Covariance serialization can differ at printed precision; retain its measured
            # difference without rejecting a run for an unrequested additional criterion.
            if kind in ['correlation','pvalues']:assert delta<=1e-12,'Matrix not symmetric (absolute tolerance1e-12)'
            assert rec['exit_code']=='0','Original pipeline exited nonzero'
            assert job['available'] and job['rows']==10000 and job['failures']==0,'Bootstrap job log incomplete or contains failure'
            assert job['max_simultaneous']<=1,'Bootstrap concurrency exceeded1'
            r['status']='PASS'
        except Exception as ex:r['error']=str(ex);errors.append(kind+': '+str(ex))
        results.append(r)
    OUTPUT_QC.extend(results);table(QC/'stage5B1B_fastspar_output_qc.tsv',OUTPUT_QC)
    return not errors,'; '.join(errors)
def integrity(baseline):
    rs=[]
    for p,r in baseline.items():
        after=sha(p) if Path(p).is_file() else ''
        rs.append(dict(path=p,before_sha256=r['sha256'],after_sha256=after,status='PASS' if r['sha256']==after else 'FAIL'))
    table(QC/'stage5B1B_input_integrity.tsv',rs)
    return all(r['status']=='PASS' for r in rs)
def report(final=False):
    title='COMPLETE' if final and all(r['status']=='COMPLETE' for r in STATE) and len(STATE)==5 else 'FINISHED_WITH_FAILURES' if final else 'RUNNING'
    lines=[];by={r['network']:r for r in STATE}
    for n,s,t in NETWORKS:
        r=by.get(n,{})
        lines.append(f'| {n} | {s} | {t} | {r.get("status","尚未启动")} | {r.get("PID","")} | {r.get("runtime_seconds","")} | {r.get("qc_status","")} |')
    body=f'''# Stage 5B-1B — Sequential FastSpar exact-sample rerun

状态：**{title}**。更新：{now()}。

严格顺序：Guo_Control → Guo_MECFS → Guo_Common → Wallen_PD → Wallen_NHC。每组使用`nohup bash {SCRIPT} <input> <prefix>`启动；原脚本未复制、未修改。前组PID退出且完成输出检查后才启动后组。若前组失败，记录FAILED后继续，失败不伪装为完成。

| Network | Samples | Taxa | Status | PID | Runtime seconds | QC |
|---|---:|---:|---|---|---:|---|
{chr(10).join(lines)}

输入目录：`{INPUT}`。每个network使用`<network>_count.matrix`；输出前缀为`{ROOT}/<network>/<network>`。

## 原脚本参数记录

直接读取原脚本，未重设统计参数：主拟合`--iterations 50 --exclude_iterations 20 --threshold 0.2`；bootstrap `--number 10000`；每个bootstrap拟合`-i 5`；P值`--permutations 10000`。主拟合、bootstrap生成、P值步骤的线程参数按原脚本为`/proc/cpuinfo`的processor行数−3。本机启动前`nproc`记录在preflight及运行日志中。Bootstrap拟合没有显式设置exclude_iterations、threshold或threads，本报告不把它们解释为主拟合参数。

为落实“禁止同时运行多个FastSpar”，仅通过GNU parallel支持的`PARALLEL`环境变量设置`--jobs 1`，并添加`--joblog`记录每个bootstrap进程退出码及时间，`--will-cite`关闭交互提示。原脚本内统计命令保持原样。`TMPDIR`及joblog均在本阶段输出目录内；HOME未改动。GNU parallel自身已有的命令长度缓存只读沿用，未创建其他目录下的结果。

## 日志与QC

- [运行状态](qc/stage5B1B_run_status.tsv)：记录PID、start/end time、exit code、输入、输出、runtime及失败原因。
- [预检查](qc/stage5B1B_preflight_qc.tsv)：输入存在且非空、样本/物种数、首行OTU_id和输入SHA256。
- [输出QC](qc/stage5B1B_fastspar_output_qc.tsv)：三个结果文件、方阵维度、输入taxa身份/顺序、NA；correlation和pvalues按要求检查数值对称性（绝对容差1e-12，实际差值保留）。covariance对称差只记录，不作为任务未要求的额外失败条件。还检查10000条bootstrap job日志全部成功且最大并发≤1；不只依赖原脚本末尾rm命令的退出码。
- [原始输入完整性](qc/stage5B1B_input_integrity.tsv)：count矩阵、原脚本及登记的冻结源文件SHA256不变。
- `logs/<network>_fastspar.log`保存原脚本完整stdout/stderr；`logs/<network>_bootstrap_joblog.tsv`保留逐bootstrap运行记录。
- [执行记录](logs/execution_log.tsv)、[进程监控](logs/process_monitor.tsv)、[心跳](logs/heartbeat.json)。

## 各组文件与失败记录

'''
    for n,_,_ in NETWORKS:
        r=by.get(n)
        body+=f'### {n}\n\n输入：`{INPUT}/{n}_count.matrix`。输出目录：`{ROOT}/{n}`。\n\n'
        if r:
            body+=f'开始：{r["start_time"]}；结束：{r["end_time"] or "运行中"}；退出码：{r["exit_code"] or "待定"}。\n\n'
            body+=f'失败原因：{r["failure_reason"] or "无已记录失败"}。\n\n'
            for kind,key in [('correlation','correlation_file'),('covariance','covariance_file'),('pvalues','pvalue_file')]:
                if Path(r[key]).is_file():body+=f'- [{n}_{kind}.tsv]({n}/{n}_{kind}.tsv)\n'
            body+='\n'
    body+='本阶段未计算network overlap、Jaccard、sign concordance或null model；未过滤taxa、修改threshold或重新taxonomy mapping。\n'
    (ROOT/'STAGE5B1B_REPORT.md').write_text(body)

def main():
    for d in [LOG,QC,LOG/'tmp']+[ROOT/n for n,_,_ in NETWORKS]:d.mkdir(parents=True,exist_ok=True)
    lock=(LOG/'supervisor.lock').open('a+')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:raise SystemExit('An existing Stage5B-1B supervisor holds the lock; do not launch duplicates')
    (LOG/'supervisor_pid.json').write_text(json.dumps(dict(PID=os.getpid(),start_time=now(),command=shlex.join(sys.argv)),indent=2)+'\n')
    existing=QC/'stage5B1B_run_status.tsv'
    if existing.exists():STATE.extend(readtsv(existing))
    if (QC/'stage5B1B_fastspar_output_qc.tsv').exists():OUTPUT_QC.extend(readtsv(QC/'stage5B1B_fastspar_output_qc.tsv'))
    if (QC/'stage5B1B_preflight_qc.tsv').exists():PREFLIGHT.extend(readtsv(QC/'stage5B1B_preflight_qc.tsv'))
    freeze=QC/'stage5B1B_input_snapshot.json'
    if freeze.exists():baseline=json.loads(freeze.read_text())
    else:
        paths={SCRIPT,PREP/'qc/stage5B1A_artifact_manifest.tsv',PREP/'qc/stage5B1A_count_conversion_qc.tsv'}
        paths.update(Path(p) for p in json.loads((PREP/'qc/stage5B1A_input_snapshot.json').read_text()))
        paths.update(INPUT/(n+'_count.matrix') for n,_,_ in NETWORKS)
        baseline={str(p):dict(sha256=sha(p),size=p.stat().st_size) for p in sorted(paths)}
        prior={r['path']:r['sha256'] for r in readtsv(PREP/'qc/stage5B1A_artifact_manifest.tsv')}
        for p,r in baseline.items():
            if p in prior:assert r['sha256']==prior[p],'Input changed since Stage5B-1A: '+p
        freeze.write_text(json.dumps(baseline,indent=2)+'\n')
    assert integrity(baseline),'Original/frozen input hash mismatch'
    script=SCRIPT.read_text()
    for token in ['--iterations 50','--exclude_iterations 20','--threshold 0.2','--number 10000','-i 5','--permutations 10000']:assert token in script
    for network,expected_samples,expected_rows in NETWORKS:
        previous=next((r for r in STATE if r['network']==network),None)
        if previous and previous['status'] in ['COMPLETE','FAILED']:continue
        if previous:raise RuntimeError('Existing RUNNING record needs process/exit-code reconciliation before supervisor restart: '+network)
        rec=dict.fromkeys(STATUS_FIELDS,'');prefix=ROOT/network/network
        rec.update(network=network,input=str(INPUT/(network+'_count.matrix')),output_prefix=str(prefix),log_file=str(LOG/(network+'_fastspar.log')),
            correlation_file=str(prefix)+'_correlation.tsv',covariance_file=str(prefix)+'_covariance.tsv',pvalue_file=str(prefix)+'_pvalues.tsv',
            bootstrap_joblog=str(LOG/(network+'_bootstrap_joblog.tsv')),qc_status='PENDING')
        try:
            wait_idle(network)
            cpus,names=preflight(network,expected_samples,expected_rows,baseline)
            assert integrity(baseline),'Input integrity changed before launch'
            assert not any(Path(rec[k]).exists() for k in ['correlation_file','covariance_file','pvalue_file']),'Output already exists; refuse implicit overwrite'
            rec['nproc']=str(cpus);rec['main_threads']=str(sum(line.startswith('processor') for line in Path('/proc/cpuinfo').read_text().splitlines())-3)
            cmd=['nohup','bash',str(SCRIPT),rec['input'],str(prefix)]
            env=dict(os.environ,TMPDIR=str(LOG/'tmp'),PARALLEL='--jobs 1 --will-cite --joblog '+rec['bootstrap_joblog'])
            rec['command']='PARALLEL='+shlex.quote(env['PARALLEL'])+' TMPDIR='+shlex.quote(env['TMPDIR'])+' '+shlex.join(cmd)+' > '+shlex.quote(rec['log_file'])+' 2>&1 &'
            start=time.monotonic();rec['start_time']=now()
            with Path(rec['log_file']).open('ab') as out:
                proc=subprocess.Popen(cmd,cwd=ROOT,env=env,stdin=subprocess.DEVNULL,stdout=out,stderr=subprocess.STDOUT,start_new_session=True)
            rec.update(PID=str(proc.pid),status='RUNNING');STATE.append(rec);save_status();report()
            event(rec['command'],rec['input'],prefix,pid=proc.pid,start=rec['start_time'])
            print('START',network,'PID',proc.pid,rec['start_time'],flush=True)
            while proc.poll() is None:
                items=processes();record_processes(network,items)
                (LOG/'heartbeat.json').write_text(json.dumps(dict(timestamp=now(),network=network,state='RUNNING',PID=proc.pid,runtime_seconds=round(time.monotonic()-start,1),processes=items),indent=2)+'\n')
                time.sleep(10)
            code=proc.wait();rec.update(exit_code=str(code),end_time=now(),runtime_seconds=str(round(time.monotonic()-start,3)))
            # Do not start another network if any descendant inference process is still alive.
            wait_idle(network)
            ok,reason=check_outputs(rec,names)
            if not integrity(baseline):ok=False;reason+='; original/frozen input changed'
            rec.update(status='COMPLETE' if ok else 'FAILED',qc_status='PASS' if ok else 'FAIL',failure_reason=reason)
            event(rec['command'],rec['input'],prefix,code,proc.pid,rec['start_time'],rec['end_time'],rec['runtime_seconds'])
            print(rec['status'],network,rec['runtime_seconds'],reason,flush=True)
        except Exception as ex:
            rec.update(status='FAILED',qc_status='FAIL',end_time=now(),failure_reason=type(ex).__name__+': '+str(ex))
            if not rec['exit_code']:rec['exit_code']='not_launched' if not rec['PID'] else 'unknown'
            if not any(r is rec for r in STATE):STATE.append(rec)
            event('Network exception; mark FAILED and continue',rec['input'],rec['failure_reason'],rec['exit_code'],rec['PID'],rec['start_time'],rec['end_time'])
            print('FAILED',network,rec['failure_reason'],flush=True)
        save_status();report()
    integrity_ok=integrity(baseline);report(final=True)
    good=len(STATE)==5 and all(r['status']=='COMPLETE' for r in STATE) and integrity_ok
    artifacts=[Path(r[k]) for r in STATE for k in ['correlation_file','covariance_file','pvalue_file'] if Path(r[k]).is_file()]
    table(QC/'stage5B1B_output_manifest.tsv',[dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in artifacts],['path','bytes','sha256'])
    result=dict(stage='Stage5B-1B',status='COMPLETE' if good else 'FINISHED_WITH_FAILURES',timestamp=now(),network_status={r['network']:r['status'] for r in STATE},input_integrity='PASS' if integrity_ok else 'FAIL')
    (LOG/'completion.json').write_text(json.dumps(result,indent=2)+'\n');(LOG/'heartbeat.json').write_text(json.dumps(result,indent=2)+'\n')
    terminal=('Stage 5B-1B COMPLETE' if good else 'Stage 5B-1B FINISHED_WITH_FAILURES')+'\n'
    terminal+='\n'.join(r[k] for r in STATE for k in ['correlation_file','pvalue_file'])+'\n'+str(ROOT/'STAGE5B1B_REPORT.md')+'\n'
    (LOG/'completion_terminal.txt').write_text(terminal);print(terminal,flush=True)
    event(result['status'],INPUT,ROOT/'STAGE5B1B_REPORT.md',0 if good else 1)

if __name__=='__main__':
    try:main()
    except Exception as ex:
        if LOG.exists():(LOG/'supervisor_error.json').write_text(json.dumps(dict(timestamp=now(),error=repr(ex)),indent=2)+'\n')
        raise
