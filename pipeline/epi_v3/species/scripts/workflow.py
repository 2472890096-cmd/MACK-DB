#!/usr/bin/env python3
"""Production scheduling only. The externally fixed EPI R program is never edited."""
import argparse
import collections
import csv
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import shutil
import subprocess
import sys
import time
import traceback
import tempfile

ROOT = Path(os.environ.get("MACKDB_EPI_RUN_ROOT", str(Path(__file__).resolve().parent.parent))).resolve()
STATUS_FIELDS = ['project', 'sample_count', 'taxa_count', 'start_time', 'end_time',
                 'runtime_seconds', 'exit_code', 'status', 'D1_count', 'D2_count',
                 'Q_count', 'error_message', 'worker_pid', 'r_pid', 'attempt']

def now():
    return dt.datetime.now().astimezone().isoformat(timespec='seconds')

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def atomic_text(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.tmp.' + str(os.getpid()))
    with temp.open('w') as f:
        f.write(value)
        f.flush()
        os.fsync(f.fileno())
    os.replace(str(temp), str(path))

def write_tsv(path, rows, fields):
    import io
    out = io.StringIO()
    writer = csv.DictWriter(out, fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
    writer.writeheader()
    writer.writerows(rows)
    atomic_text(path, out.getvalue())

def read_tsv(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))

def config():
    cfg = json.loads((ROOT / 'audit/config.json').read_text())
    core = Path(cfg['core_script'])
    cfg['core_script'] = str(core if core.is_absolute() else (ROOT / core).resolve())
    cfg['rscript'] = shutil.which(cfg['rscript']) or str((ROOT / cfg['rscript']).resolve())
    return cfg

def manifest():
    return read_tsv(ROOT / 'manifests/all_valid_projects.tsv')

def status_path(project):
    return ROOT / 'status/projects' / (project + '.status.tsv')

def status(row):
    p = status_path(row['project'])
    return read_tsv(p)[0] if p.exists() else dict(row, status='PENDING')

def outputs(project):
    return [ROOT / 'keys' / (project + '_EPI-' + key + '.key') for key in ['D1', 'D2', 'Q']] + [
        ROOT / 'summaries' / (project + '_EPI_summary.tsv')]

def key_set(path):
    return {line.strip() for line in Path(path).read_text().splitlines() if line.strip()}

def counts(project):
    return [len(key_set(p)) for p in outputs(project)[:3]]

def verify_fingerprints():
    for row in read_tsv(ROOT / 'audit/protected_fingerprints.tsv'):
        if sha(row['path']) != row['sha256']:
            raise RuntimeError('Protected file fingerprint changed: ' + row['path'])

def verify_setup():
    cfg = config()
    verify_fingerprints()
    for row in read_tsv(ROOT / 'audit/workflow_fingerprints.tsv'):
        if sha(ROOT / row['relative_path']) != row['sha256']:
            raise RuntimeError('Production setup changed: ' + row['relative_path'])
    if not Path(cfg['rscript']).is_file():
        raise RuntimeError('Rscript missing')
    rows = manifest()
    if len(rows) != cfg['valid_projects'] or len({r['project'] for r in rows}) != len(rows):
        raise RuntimeError('Manifest project count/uniqueness mismatch')
    for row in rows:
        for kind in ['abundance', 'species']:
            p = Path(row[kind + '_path'])
            st = p.stat()
            # All inputs were freshly parsed and hashed during setup. On resume,
            # rehash changed metadata; each worker always rehashes its own inputs.
            if st.st_size != int(row[kind + '_bytes']) or st.st_mtime_ns != int(row[kind + '_mtime_ns']):
                if sha(p) != row[kind + '_sha256']:
                    raise RuntimeError('Input changed: ' + str(p))
    for group, select in [('le5000', lambda n: n <= 5000), ('gt5000', lambda n: n > 5000)]:
        expected = [r for r in rows if select(int(r['sample_count']))]
        actual = read_tsv(ROOT / ('manifests/projects_' + group + '.tsv'))
        if actual != expected:
            raise RuntimeError('Partition mismatch: ' + group)
    return rows

def install_output(src, dst):
    """Never overwrite an existing production output, including crash recovery."""
    if dst.exists():
        if sha(src) != sha(dst):
            raise RuntimeError('Existing output differs; retained both copies: ' + str(dst))
        return
    os.link(str(src), str(dst))  # exclusive destination creation, same filesystem
    src.unlink()

def worker(project, abundance, species, samples):
    if not re.fullmatch(r'PRJ[A-Z0-9]+', project):
        raise ValueError('Invalid project identifier')
    row = next(r for r in manifest() if r['project'] == project)
    if [abundance, species, str(samples)] != [row['abundance_path'], row['species_path'], row['sample_count']]:
        raise ValueError('Worker arguments disagree with frozen manifest')
    lock = (ROOT / 'status/projects' / (project + '.lock')).open('a')
    fcntl.flock(lock, fcntl.LOCK_EX)
    previous = status(row)
    if previous.get('status') == 'SUCCESS' and all(p.is_file() for p in outputs(project)):
        print(now(), 'SKIP_SUCCESS', project, flush=True)
        return 0
    attempt = dt.datetime.now().strftime('%Y%m%dT%H%M%S.%f') + '_' + str(os.getpid())
    work = ROOT / 'tmp' / project / attempt
    work.mkdir(parents=True, exist_ok=False)
    atomic_text(work / 'previous_status.json', json.dumps(previous, indent=2) + '\n')
    cfg = config()
    env = os.environ.copy()
    env.update(TMPDIR=str(work), TMP=str(work), TEMP=str(work))
    env['MACKDB_EPI_RESOURCE_LOG'] = str(work / 'resource_usage.txt')
    temporary_r_dir = None
    # R startup rejects whitespace in R_TempDir; output/commit paths stay unchanged.
    if any(c.isspace() for c in str(work)):
        temporary_r_dir = tempfile.TemporaryDirectory(prefix='mackdb_epi_r_', dir='/tmp')
        env.update(TMPDIR=temporary_r_dir.name, TMP=temporary_r_dir.name, TEMP=temporary_r_dir.name)
    child = None
    started = time.monotonic()
    s = dict(row, start_time=now(), end_time='', runtime_seconds='', exit_code='',
             status='RUNNING', D1_count='', D2_count='', Q_count='', error_message='',
             worker_pid=os.getpid(), r_pid='', attempt=attempt)
    write_tsv(status_path(project), [s], STATUS_FIELDS)
    rc = 1
    try:
        if sha(cfg['core_script']) != cfg['core_sha256']:
            raise RuntimeError('V3 core script hash changed')
        for kind in ['abundance', 'species']:
            if sha(row[kind + '_path']) != row[kind + '_sha256']:
                raise RuntimeError(kind + ' input hash changed')
        cmd = [cfg['rscript'], cfg['core_script'], '-a', abundance, '-s', species,
               '-o', project, '-d', str(work), '--mode', 'prevalence',
               '--prevalence-min', '0.25', '--prevalence-max', '0.75', '--network-quantile', '0.1']
        atomic_text(work / 'command.json', json.dumps(cmd, indent=2) + '\n')
        out_path = ROOT / 'logs/projects' / (project + '.stdout.log')
        err_path = ROOT / 'logs/projects' / (project + '.stderr.log')
        with out_path.open('ab', buffering=0) as out, err_path.open('ab', buffering=0) as err:
            header = ('\n=== attempt ' + attempt + ' ' + now() + ' ===\n').encode()
            out.write(header)
            err.write(header)
            child = subprocess.Popen(cmd, cwd=str(work), env=env, stdout=out, stderr=err)
            s['r_pid'] = child.pid
            write_tsv(status_path(project), [s], STATUS_FIELDS)
            rc = child.wait()
        s['exit_code'] = rc
        if rc != 0:
            raise RuntimeError('Rscript exit code ' + str(rc) + '; see project stderr log')
        files = [work / p.name for p in outputs(project)]
        if not all(p.is_file() for p in files):
            raise RuntimeError('Rscript exit 0 but one or more required outputs missing')
        summary = read_tsv(files[3])
        cs = [len(key_set(p)) for p in files[:3]]
        if len(summary) != 1 or summary[0]['project'] != project:
            raise RuntimeError('Invalid project summary')
        sr = summary[0]
        if int(sr['sample_total']) != int(row['sample_count']) or int(sr['species_total']) != int(row['taxa_count']):
            raise RuntimeError('Summary dimensions disagree with manifest')
        if cs != [int(sr[k + '_key_count']) for k in ['D1', 'D2', 'Q']]:
            raise RuntimeError('Key counts disagree with summary')
        for src, dst in zip(files, outputs(project)):
            install_output(src, dst)
        write_tsv(work / 'committed_output_hashes.tsv',
                  [{'path': str(p), 'sha256': sha(p)} for p in outputs(project)], ['path', 'sha256'])
        s.update(status='SUCCESS', exit_code=0, D1_count=cs[0], D2_count=cs[1], Q_count=cs[2])
    except BaseException as exc:
        if child is not None and child.poll() is None:
            child.terminate()
            child.wait()
        s.update(status='FAIL', error_message=str(exc).replace('\n', ' ')[:2000])
        if s['exit_code'] == '':
            s['exit_code'] = child.returncode if child is not None else 1
        atomic_text(work / 'failure.txt', traceback.format_exc())
    finally:
        if temporary_r_dir is not None:
            temporary_r_dir.cleanup()
        s.update(end_time=now(), runtime_seconds='{:.6f}'.format(time.monotonic() - started))
        write_tsv(status_path(project), [s], STATUS_FIELDS)
        lock.close()
    print(now(), s['status'], project, s['runtime_seconds'], s['error_message'], flush=True)
    return 0 if s['status'] == 'SUCCESS' else 1

def merge(rows, finalize=False):
    all_status = []
    for row in rows:
        s = status(row)
        if finalize and s['status'] != 'SUCCESS' and s['status'] != 'FAIL':
            s.update(status='FAIL', end_time=now(), exit_code=1,
                     error_message='Worker did not finish; check master and GNU parallel logs')
            write_tsv(status_path(row['project']), [s], STATUS_FIELDS)
        all_status.append(s)
    write_tsv(ROOT / 'status/project_status.tsv', all_status, STATUS_FIELDS)
    return all_status

def set_stage(stage):
    atomic_text(ROOT / 'status/current_stage.txt', stage + '\n')
    atomic_text(ROOT / 'status' / stage, now() + '\n')
    print(now(), stage, flush=True)

def postprocess(rows):
    states = merge(rows, finalize=True)
    final_rows = []
    output_hashes = []
    for s in states:
        if s['status'] != 'SUCCESS':
            continue
        project = s['project']
        paths = outputs(project)
        if not all(p.is_file() for p in paths):
            raise RuntimeError('SUCCESS output missing: ' + project)
        committed = ROOT / 'tmp' / project / s['attempt'] / 'committed_output_hashes.tsv'
        for record in read_tsv(committed):
            if sha(record['path']) != record['sha256']:
                raise RuntimeError('Committed production output changed: ' + record['path'])
        sets = [key_set(p) for p in paths[:3]]
        cs = [len(x) for x in sets]
        sr = read_tsv(paths[3])[0]
        if cs != [int(s[k + '_count']) for k in ['D1', 'D2', 'Q']] or cs != [int(sr[k + '_key_count']) for k in ['D1', 'D2', 'Q']]:
            raise RuntimeError('Final count integrity failure: ' + project)
        votes = collections.Counter(t for group in sets for t in group)
        final = sorted(t for t, n in votes.items() if n >= 2)
        content = ''.join(t + '\n' for t in final)
        path = ROOT / 'final_2of3' / (project + '_EPI-2of3.key')
        if path.exists():
            if path.read_text() != content:
                raise RuntimeError('Existing 2of3 differs: ' + project)
        else:
            temp = ROOT / 'tmp' / (project + '.2of3.' + str(os.getpid()))
            atomic_text(temp, content)
            install_output(temp, path)
        if key_set(path) != set(final):
            raise RuntimeError('2of3 integrity failure: ' + project)
        for p in paths + [path]:
            output_hashes.append({'project': project, 'path': str(p), 'sha256': sha(p)})
        final_rows.append(dict(project=project, D1_count=cs[0], D2_count=cs[1], Q_count=cs[2], final_2of3_count=len(final)))
    verify_fingerprints()
    write_tsv(ROOT / 'audit/final_output_sha256.tsv', output_hashes, ['project', 'path', 'sha256'])
    write_tsv(ROOT / 'status/final_2of3_summary.tsv', final_rows,
              ['project', 'D1_count', 'D2_count', 'Q_count', 'final_2of3_count'])
    failed = [s for s in states if s['status'] != 'SUCCESS']
    write_tsv(ROOT / 'status/failed_projects.tsv', failed, STATUS_FIELDS)
    totals = dict(input_valid_projects=len(rows), le5000_projects=sum(int(r['sample_count']) <= 5000 for r in rows),
                  gt5000_projects=sum(int(r['sample_count']) > 5000 for r in rows),
                  success_projects=len(final_rows), failed_projects=len(failed),
                  D1_total_candidates=sum(r['D1_count'] for r in final_rows),
                  D2_total_candidates=sum(r['D2_count'] for r in final_rows),
                  Q_total_candidates=sum(r['Q_count'] for r in final_rows),
                  two_of_three_total_candidates=sum(r['final_2of3_count'] for r in final_rows))
    if totals['success_projects'] + totals['failed_projects'] != len(rows):
        raise RuntimeError('Final project accounting failed')
    write_tsv(ROOT / 'status/FINAL_SUMMARY.tsv', [totals], list(totals))
    atomic_text(ROOT / 'audit/final_integrity.json', json.dumps(dict(checked_at=now(), passed=True,
                protected_fingerprints_unchanged=True, successful_outputs_checked=len(final_rows)), indent=2) + '\n')
    set_stage('POSTPROCESS_FINISHED')
    flag = 'FINISHED_WITH_FAILURES.flag' if failed else 'ALL_FINISHED.ok'
    atomic_text(ROOT / 'status' / flag, now() + '\n')

def alive(pid):
    try:
        os.kill(int(pid), 0)
        return Path('/proc/' + str(pid) + '/stat').read_text().split(') ', 1)[1].split()[0] != 'Z'
    except (OSError, ValueError):
        return False

def project_alive(pid, project):
    if not alive(pid):
        return False
    try:
        args = Path('/proc/' + str(pid) + '/cmdline').read_bytes().split(b'\0')
        return project.encode() in args
    except OSError:
        return False

def master_alive(pid):
    if not alive(pid):
        return False
    try:
        identity = json.loads((ROOT / 'status/master_identity.json').read_text())
        original = identity['proc_stat'].split(') ', 1)[1].split()[19]
        current = Path('/proc/' + str(pid) + '/stat').read_text().split(') ', 1)[1].split()[19]
        return int(pid) == identity['pid'] and original == current
    except (OSError, ValueError, KeyError, IndexError):
        return False

def check_status():
    pid_path = ROOT / 'status/master.pid'
    pid = pid_path.read_text().strip() if pid_path.exists() else ''
    rows = manifest()
    states = [status(r) for r in rows]
    count = collections.Counter(s['status'] for s in states)
    stage = ROOT / 'status/current_stage.txt'
    print('Current time:', now())
    print('Master PID:', pid or 'not started', 'alive:', bool(pid and master_alive(pid)))
    print('Current stage:', stage.read_text().strip() if stage.exists() else 'NOT_STARTED')
    print('Total projects:', len(rows))
    print('SUCCESS:', count['SUCCESS'], 'FAIL:', count['FAIL'])
    running = [s for s in states if s['status'] == 'RUNNING']
    print('Running projects:', sum(project_alive(s.get('worker_pid', ''), s['project']) for s in running), '(recorded RUNNING:', len(running), ')')
    print('Recent completed projects:')
    for s in sorted([s for s in states if s.get('end_time')], key=lambda s: s['end_time'])[-5:]:
        print(s['project'], s['status'], s['end_time'], s.get('runtime_seconds', ''))
    print('Stage 2 started:', (ROOT / 'status/STAGE2_RUNNING').exists() or (ROOT / 'status/STAGE2_FINISHED').exists())
    print('Final summary exists:', (ROOT / 'status/FINAL_SUMMARY.tsv').exists())

def master():
    lock = (ROOT / 'status/master.lock').open('a')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        print('A production master already holds the lock', file=sys.stderr)
        return 2
    # A master killed with SIGKILL can leave workers alive. Refuse overlapping resume.
    for p in (ROOT / 'status/projects').glob('*.status.tsv'):
        old = read_tsv(p)[0]
        if old['status'] == 'RUNNING' and (project_alive(old.get('worker_pid', ''), old['project']) or project_alive(old.get('r_pid', ''), old['project'])):
            raise RuntimeError('Existing live worker prevents overlapping master: ' + old['project'])
    atomic_text(ROOT / 'status/master.pid', str(os.getpid()) + '\n')
    atomic_text(ROOT / 'status/master_identity.json', json.dumps({'pid': os.getpid(), 'started': now(),
                'proc_stat': Path('/proc/self/stat').read_text()}, indent=2) + '\n')
    run_id = dt.datetime.now().strftime('%Y%m%dT%H%M%S') + '_' + str(os.getpid())
    archive = ROOT / 'status/previous_runs' / run_id
    for name in ['ALL_FINISHED.ok', 'FINISHED_WITH_FAILURES.flag', 'MASTER_FAILED.flag',
                 'STAGE1_RUNNING', 'STAGE1_FINISHED', 'STAGE2_RUNNING', 'STAGE2_FINISHED',
                 'POSTPROCESS_FINISHED', 'FINAL_SUMMARY.tsv']:
        p = ROOT / 'status' / name
        if p.exists():
            archive.mkdir(parents=True, exist_ok=True)
            p.rename(archive / name)
    child = None
    def stop(signum, frame):
        raise InterruptedError('Master received signal ' + str(signum))
    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    try:
        print(now(), 'INPUT_CHECK', flush=True)
        rows = verify_setup()
        print(now(), 'VALID_PROJECTS', len(rows), 'LE5000', sum(int(r['sample_count']) <= 5000 for r in rows),
              'GT5000', sum(int(r['sample_count']) > 5000 for r in rows), flush=True)
        set_stage('STAGE1_RUNNING')
        command = ['/usr/bin/parallel', '--plain', '--will-cite', '--jobs', '8', '--colsep', '\t',
                   '--skip-first-line', '--joblog', str(ROOT / 'logs/master' / ('parallel_' + run_id + '.tsv')),
                   str(ROOT / 'scripts/run_one_project.sh'), '{1}', '{2}', '{3}', '{4}',
                   '::::', str(ROOT / 'manifests/projects_le5000.tsv')]
        child = subprocess.Popen(command, cwd=str(ROOT), start_new_session=True)
        parallel_rc = child.wait()
        print(now(), 'STAGE1_PARALLEL_EXIT', parallel_rc, '(project failures do not stop stage 2)', flush=True)
        # wait() above is a hard barrier: no >5000 project starts before every stage-1 job returns.
        merge(rows)
        for row in rows:
            if int(row['sample_count']) <= 5000:
                s = status(row)
                if s['status'] == 'RUNNING' and (project_alive(s.get('worker_pid', ''), row['project']) or project_alive(s.get('r_pid', ''), row['project'])):
                    raise RuntimeError('Live stage-1 process prevents stage-2 start: ' + row['project'])
        (ROOT / 'status/STAGE1_RUNNING').unlink()
        set_stage('STAGE1_FINISHED')
        set_stage('STAGE2_RUNNING')
        for row in read_tsv(ROOT / 'manifests/projects_gt5000.tsv'):
            child = subprocess.Popen([str(ROOT / 'scripts/run_one_project.sh'), row['project'],
                        row['abundance_path'], row['species_path'], row['sample_count']],
                        cwd=str(ROOT), start_new_session=True)
            rc = child.wait()
            print(now(), 'STAGE2_PROJECT_EXIT', row['project'], rc, flush=True)
        merge(rows, finalize=True)
        (ROOT / 'status/STAGE2_RUNNING').unlink()
        set_stage('STAGE2_FINISHED')
        postprocess(rows)
        print(now(), 'PRODUCTION_FINISHED', flush=True)
        return 0
    except BaseException:
        if child is not None and child.poll() is None:
            os.killpg(child.pid, signal.SIGTERM)
            child.wait()
        detail = traceback.format_exc()
        atomic_text(ROOT / 'status/MASTER_FAILED.flag', now() + '\n' + detail)
        atomic_text(ROOT / 'status/current_stage.txt', 'MASTER_FAILED\n')
        print(detail, file=sys.stderr, flush=True)
        return 1
    finally:
        lock.close()

def main():
    mode = sys.argv[1]
    if mode == 'worker':
        def stop(signum, frame):
            raise InterruptedError('Worker received signal ' + str(signum))
        signal.signal(signal.SIGTERM, stop)
        return worker(*sys.argv[2:])
    if mode == 'master':
        return master()
    if mode == 'status':
        check_status()
        return 0
    raise ValueError('Unknown mode')

if __name__ == '__main__':
    sys.exit(main())
