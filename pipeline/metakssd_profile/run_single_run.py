#!/usr/bin/env python3
"""Pinned ERR14036375 example: official MetaKSSD and unchanged MACK-DB R Profile.

No downloads, installations, batch jobs, scientific reimplementation or SRA
execution. Relative data paths resolve from the private configuration file.
"""
import argparse
import array
import csv
import ctypes
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import struct
import subprocess
import sys
import tempfile
import time
import uuid

ROOT = Path(__file__).resolve().parents[2]
EXAMPLE = ROOT/'examples/real_fastq/ERR14036375'
PROFILE = ROOT/'pipeline/metakssd_profile/profile.R'
TIMEOUT = 600


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path = Path(path)
    temporary = path.with_name('.'+path.name+'.writing-'+uuid.uuid4().hex)
    try:
        with temporary.open('x') as f:
            json.dump(value, f, indent=2)
            f.write('\n')
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def pins():
    return json.loads((EXAMPLE/'reference_pins.json').read_text())


def resolve_file(value, base):
    path = Path(value).expanduser()
    return (base/path).resolve() if not path.is_absolute() else path.resolve()


def executable(value, base):
    if '/' in value:
        path = resolve_file(value, base)
    else:
        found = shutil.which(value)
        require(found is not None, 'Executable not found on PATH: '+value)
        path = Path(found).resolve()
    require(path.is_file() and os.access(path, os.X_OK), 'Not executable: '+str(path))
    return path


def settings(config):
    config = Path(config).resolve()
    cfg = json.loads(config.read_text())
    required = {'run_accession', 'fastq_r1', 'fastq_r2', 'metakssd_bin',
                'shuffle', 'markerdb_dir', 'rscript', 'scratch_dir'}
    require(set(cfg) in [required, required | {'metakssd_sha256'}],
            'Config must contain only documented path/run fields and optional metakssd_sha256')
    p = pins()
    require(cfg['run_accession'] == p['run_accession'], 'This example accepts only the fixed ERR14036375 Run')
    require(all(isinstance(v, str) and v for v in cfg.values()), 'Config values must be nonempty strings')
    result = {k: resolve_file(cfg[k], config.parent) for k in
              ['fastq_r1', 'fastq_r2', 'shuffle', 'markerdb_dir', 'scratch_dir']}
    for k in ['metakssd_bin', 'rscript']:
        result[k] = executable(cfg[k], config.parent)
    result.update(run_accession=cfg['run_accession'], config_path=config,
                  metakssd_sha256=cfg.get('metakssd_sha256', p['metakssd_sha256']))
    require(re.fullmatch('[0-9a-f]{64}', result['metakssd_sha256']), 'Invalid executable SHA256 pin')
    return result


class Header(ctypes.LittleEndianStructure):
    _fields_ = [('shuf_id', ctypes.c_uint), ('koc', ctypes.c_bool),
                ('kmerlen', ctypes.c_int), ('dim_rd_len', ctypes.c_int),
                ('comp_num', ctypes.c_int), ('infile_num', ctypes.c_int),
                ('all_ctx_ct', ctypes.c_ulonglong)]


def header(path):
    require(sys.byteorder == 'little' and ctypes.sizeof(Header) == 32,
            'Validated Sketch header layout requires little-endian Linux and 32-byte co_dstat_t')
    with Path(path).open('rb') as f:
        value = Header.from_buffer_copy(f.read(32))
    return {k: getattr(value, k) for k, _ in Header._fields_}


def probe(argv):
    q = subprocess.run(list(map(str, argv)), capture_output=True, text=True, timeout=30,
                       env=dict(os.environ, LC_ALL='C'))
    require(q.returncode == 0, 'Dependency probe failed: '+str(argv[0])+'; '+q.stderr[-1000:])
    return q.stdout.strip()


def preflight(config):
    require(sys.platform.startswith('linux'), 'This execution recipe is validated for Linux')
    c, p = settings(config), pins()
    require(Path('/usr/bin/time').is_file(), '/usr/bin/time is required')
    require(shutil.which('gzip') is not None and shutil.which('zcat') is not None, 'gzip and zcat are required')
    require(c['scratch_dir'].is_dir() and os.access(c['scratch_dir'], os.W_OK), 'Scratch directory is not writable')
    require(re.fullmatch(r'[A-Za-z0-9_./-]+', str(c['scratch_dir'])),
            'Use a short scratch path without whitespace or shell metacharacters (MetaKSSD internal popen)')
    require(len(str(c['scratch_dir'])) < 100, 'Scratch root too long for MetaKSSD internal PATHLEN=256')
    require(shutil.disk_usage(c['scratch_dir']).free >= 3*1024**3, 'At least 3GiB free scratch space required')
    checks = []
    for role, key in [('R1', 'fastq_r1'), ('R2', 'fastq_r2')]:
        f, spec = c[key], p['fastq'][role]
        require(f.is_file(), 'Missing '+role+' FASTQ: '+str(f))
        require(f.stat().st_size == spec['bytes'], role+' FASTQ byte count mismatch')
        require(digest(f, 'md5') == spec['md5'], role+' FASTQ MD5 mismatch')
        require(digest(f) == spec['sha256'], role+' FASTQ SHA256 mismatch')
        checks.append({'role': role, 'path': str(f), 'bytes': f.stat().st_size, 'sha256': digest(f)})
    require(c['fastq_r1'] != c['fastq_r2'], 'The two mates must be different files')
    for f, expected in [(c['shuffle'], p['shuffle_sha256']), (PROFILE, p['profile_core_sha256']),
                        (c['metakssd_bin'], c['metakssd_sha256'])] + [
                            (c['markerdb_dir']/name, value) for name, value in p['markerdb_sha256'].items()]:
        require(f.is_file() and digest(f) == expected, 'Pinned file missing or SHA256 mismatch: '+str(f))
        checks.append({'role': 'tool_reference_or_core', 'path': str(f), 'bytes': f.stat().st_size, 'sha256': expected})
    require(len(str(c['shuffle'])) < 200 and len(str(c['markerdb_dir'])) < 200,
            'Native MetaKSSD reference/shuffle paths must fit PATHLEN=256; use shorter configured paths')
    version = probe([c['metakssd_bin'], '--version'])
    require(version == 'metakssd version 2.21', 'MetaKSSD version must be exactly2.21')
    reference = header(c['markerdb_dir']/'cofiles.stat')
    require(reference['shuf_id'] == p['shuffle_id'] and reference['kmerlen'] == 22 and
            reference['dim_rd_len'] == 6, 'MarkerDB header is incompatible')
    r_version = probe([c['rscript'], '-e',
        'cat(R.version.string,"\\n");for(p in c("tidyverse","optparse")){stopifnot(requireNamespace(p,quietly=TRUE));cat(p,as.character(packageVersion(p)),"\\n")}'])
    return c, {'status': 'PREFLIGHT_PASS', 'run_accession': p['run_accession'], 'files': checks,
               'config_sha256': digest(c['config_path']), 'metakssd_version': version,
               'metakssd_build_matches_recorded_sha': c['metakssd_sha256'] == p['metakssd_sha256'],
               'R_environment': r_version, 'reference_header': reference,
               'parameters': {'threads': 1, 'memory_G': 2, 'sharekmer': 12, 'abundance_mode': True,
                              'shuffle_id': p['shuffle_id'], 'stage_limit_seconds': TIMEOUT},
               'network_requests': 0}


def validate_pair(first, second, expected, limit=TIMEOUT):
    start = time.monotonic()
    counts, bases, longest = [0, 0], [0, 0], [0, 0]
    def record(f):
        h = f.readline()
        if not h:
            return None
        seq, plus, qual = f.readline(), f.readline(), f.readline()
        require(h.startswith(b'@') and plus.startswith(b'+') and qual.endswith(b'\n'), 'Malformed/truncated FASTQ record')
        seq, qual = seq.rstrip(b'\r\n'), qual.rstrip(b'\r\n')
        require(seq and len(seq) == len(qual), 'Sequence/quality length mismatch')
        name = h.split()[0][1:]
        if name.endswith((b'/1', b'/2')):
            name = name[:-2]
        return name, len(seq)
    with gzip.open(first, 'rb') as a, gzip.open(second, 'rb') as b:
        while True:
            x, y = record(a), record(b)
            if x is None and y is None:
                break
            require(x is not None and y is not None and x[0] == y[0], 'Paired record IDs/counts do not match; no automatic reordering')
            for i, r in enumerate((x, y)):
                counts[i] += 1
                bases[i] += r[1]
                longest[i] = max(longest[i], r[1])
            if counts[0] % 65536 == 0:
                require(time.monotonic()-start < limit, 'FASTQ validation exceeded its bounded runtime')
    require(sum(counts) == expected['read_count'] and sum(bases) == expected['base_count'], 'FASTQ read/base totals differ from frozen metadata')
    require(max(longest) < 4094, 'Read length exceeds the validated MetaKSSD short-read buffer contract')
    return {'status': 'PASS', 'mate_records': counts, 'mate_bases': bases,
            'max_read_lengths': longest, 'wall_seconds': time.monotonic()-start,
            'gzip_trailers_validated': True, 'paired_record_ids_match': True}


def run_command(logs, label, argv, stdout_path=None, timeout=TIMEOUT):
    logs = Path(logs)
    stdout_path = Path(stdout_path) if stdout_path else logs/(label+'.stdout')
    stderr_path, time_path = logs/(label+'.stderr'), logs/(label+'.time.txt')
    row = {'stage': label, 'argv': list(map(str, argv)), 'stdout': str(stdout_path),
           'stderr': str(stderr_path), 'time_v_log': str(time_path), 'limit_seconds': timeout,
           'started_UTC': datetime.now(timezone.utc).isoformat()}
    timed = ['/usr/bin/time', '-v', '-o', str(time_path)] + row['argv']
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', LC_ALL='C')
    with (logs/'commands.jsonl').open('a') as f:
        f.write(json.dumps(dict(row, env={k: env[k] for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'LC_ALL']}))+'\n')
    start = time.monotonic()
    with stdout_path.open('xb') as out, stderr_path.open('xb') as err:
        process = subprocess.Popen(timed, stdout=out, stderr=err, env=env, start_new_session=True)
        try:
            process.wait(timeout=timeout)
        except (subprocess.TimeoutExpired, KeyboardInterrupt):
            # Signal only the session created immediately above, including native decompression children.
            try:
                os.killpg(process.pid, signal.SIGTERM)
                process.wait(timeout=5)
            except ProcessLookupError:
                pass
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            row.update(status='TIMEOUT_OR_INTERRUPTED', wall_seconds=time.monotonic()-start, exit_code=process.returncode,
                       ended_UTC=datetime.now(timezone.utc).isoformat())
            write_json(logs/(label+'.status.json'), row)
            raise RuntimeError(label+' timed out/interrupted; no retry; inspect isolated logs')
    row.update(status='COMMAND_COMPLETE', wall_seconds=time.monotonic()-start, exit_code=process.returncode,
               ended_UTC=datetime.now(timezone.utc).isoformat())
    write_json(logs/(label+'.status.json'), row)
    require(process.returncode == 0, label+' failed; no retry; inspect '+str(stderr_path))
    return row


def compare_outputs(sketch, query, profile):
    p = pins()
    h = header(sketch/'cofiles.stat')
    require(all(h[k] == p['sketch_header'][k] for k in p['sketch_header']), 'New Sketch header differs')
    keys, counts = array.array('I'), array.array('H')
    keys.frombytes((sketch/'combco.0').read_bytes())
    counts.frombytes((sketch/'combco.0.a').read_bytes())
    require(len(keys) == len(counts) == h['all_ctx_ct'] and len(set(keys)) == len(keys), 'Invalid Sketch key/count layout')
    canonical = hashlib.sha256()
    for key, count in sorted(zip(keys, counts)):
        canonical.update(struct.pack('<IH', key, count))
    require(canonical.hexdigest() == p['sorted_key_count_sha256'], 'Sketch key/count scientific content differs; no parameter search')
    actual = list(csv.reader(query.open(), delimiter='\t'))
    expected = list(csv.reader((EXAMPLE/'expected_query_without_input_path.tsv').open(), delimiter='\t'))
    require(actual and all(len(r) == 7 for r in actual) and len({r[1] for r in actual}) == len(actual), 'Query shape/taxon IDs invalid')
    q = {r[1]: r[2:] for r in actual}
    e = {r[0]: r[1:] for r in expected}
    require(q.keys() == e.keys(), 'Query taxon names differ')
    require(all(math.isfinite(float(x)) and float(x) == float(y)
                for name in e for x, y in zip(q[name], e[name])), 'Query scientific values differ')
    require(digest(profile) == p['expected_profile_sha256'], 'Species Profile is not byte-identical to the accepted reference')
    return {'status': 'PASS', 'sketch_key_count_exact_match': True, 'sketch_keys': len(keys),
            'sorted_key_count_sha256': canonical.hexdigest(), 'query_taxa': len(q),
            'query_all_numeric_fields_equal': True, 'query_input_path_exempt': True,
            'profile_byte_identical': True, 'profile_max_abs_diff': 0,
            'scope': 'One real Run; serialized Sketch order/path/padding are not required to be byte-identical'}


def safe_cleanup(workspace, scratch_root, token):
    workspace, scratch_root = Path(workspace), Path(scratch_root).resolve()
    require(not workspace.is_symlink() and workspace.parent.resolve() == scratch_root and
            workspace.name.startswith('mackdb-upstream-'), 'Refusing cleanup outside the owned scratch directory')
    require((workspace/'.owner').read_text() == token, 'Scratch ownership marker differs; retain for inspection')
    shutil.rmtree(workspace)  # Only the uniquely created, ownership-checked directory; never an input path.


def execute(config, output, keep_scratch=False):
    output = Path(output).expanduser().resolve()
    require(not output.exists(), 'Output directory already exists; no overwrite or automatic rerun')
    c, check = preflight(config)
    if output.is_relative_to(ROOT):
        require(output.is_relative_to(ROOT/'work') or output.is_relative_to(ROOT/'tests/results'),
                'Within the repository, put generated outputs only under ignored work/ or tests/results/')
    output.parent.mkdir(parents=True, exist_ok=True)
    require(shutil.disk_usage(output.parent).free >= 3*1024**3, 'At least 3GiB free output space required')
    output.mkdir()
    logs = output/'logs'; logs.mkdir()
    state = {'status': 'RUNNING', 'run_accession': c['run_accession'], 'preflight': check, 'stages': [],
             'started_UTC': datetime.now(timezone.utc).isoformat()}
    token = uuid.uuid4().hex
    workspace = Path(tempfile.mkdtemp(prefix='mackdb-upstream-', dir=c['scratch_dir']))
    (workspace/'.owner').write_text(token)
    state['owned_scratch'] = str(workspace)
    status = output/'RUN_STATUS.json'
    write_json(status, state)
    def terminate(signum, frame):
        raise KeyboardInterrupt('SIGTERM: retain evidence and terminate only owned children')
    previous_sigterm = signal.signal(signal.SIGTERM, terminate)
    try:
        state['pair_validation'] = validate_pair(c['fastq_r1'], c['fastq_r2'], pins())
        write_json(status, state)
        fq = workspace/(c['run_accession']+'.fastq')
        state['stages'].append(run_command(logs, 'decompress', [shutil.which('gzip'), '-cd', c['fastq_r1'], c['fastq_r2']], fq))
        state['generated_read_collection'] = {'bytes': fq.stat().st_size, 'sha256': digest(fq),
                                             'operation': 'Concatenate complete FASTQ records; no trimming, filtering or joined read sequences'}
        sketch = workspace/(c['run_accession']+'_sketch')
        require(len(str(fq))+10 < 256 and len(str(sketch))+20 < 256, 'Native MetaKSSD PATHLEN exceeded')
        state['stages'].append(run_command(logs, 'sketch', [c['metakssd_bin'], 'dist', '-L', c['shuffle'], '-A', '-p', '1', '-m', '2', '-o', sketch, fq]))
        staging = output/'.staging'; staging.mkdir()
        query_dir = staging/'query'; query_dir.mkdir()
        query = query_dir/(c['run_accession']+'_all.txt')
        state['stages'].append(run_command(logs, 'query', [c['metakssd_bin'], 'composite', '-r', c['markerdb_dir'], '-q', sketch], query))
        runlist = output/'SINGLE_RUN.txt'; runlist.write_text(c['run_accession']+'\n')
        state['stages'].append(run_command(logs, 'profile', [c['rscript'], PROFILE, '-l', runlist, '-i', query_dir, '-s', '12', '-o', staging/'profile']))
        profile = staging/'profile'/(c['run_accession']+'_profile.txt')
        state['comparison'] = compare_outputs(sketch, query, profile)
        require(all(digest(Path(r['path'])) == r['sha256'] for r in check['files']), 'Input/tool/reference changed during execution')
        shutil.copytree(sketch, staging/'sketch')
        for name in ['sketch', 'query', 'profile']:
            os.replace(staging/name, output/name)
        staging.rmdir()
        if keep_scratch:
            state['scratch_disposition'] = 'RETAINED_BY_EXPLICIT_OPTION'
        else:
            safe_cleanup(workspace, c['scratch_dir'], token)
            state['scratch_disposition'] = 'OWNED_GENERATED_SCRATCH_REMOVED_AFTER_SUCCESS'
        state['status'] = 'PASS_REAL_FASTQ_TO_PROFILE_REFERENCE_MATCH'
        state['ended_UTC'] = datetime.now(timezone.utc).isoformat()
        state['metakssd_build_matches_recorded_sha'] = check['metakssd_build_matches_recorded_sha']
        write_json(status, state)
        with (output/'OUTPUT_SHA256.tsv').open('x') as f:
            w = csv.writer(f, delimiter='\t'); w.writerow(['path', 'bytes', 'sha256'])
            for file in sorted(output.rglob('*')):
                if file.is_file() and file != output/'OUTPUT_SHA256.tsv':
                    w.writerow([file.relative_to(output).as_posix(), file.stat().st_size, digest(file)])
        return state
    except BaseException as exc:
        state.update(status='FAILED_OR_INTERRUPTED_RETAINED_FOR_INSPECTION', error=str(exc),
                     ended_UTC=datetime.now(timezone.utc).isoformat(),
                     scratch_disposition='RETAINED_NO_AUTOMATIC_CLEANUP_OR_RETRY')
        write_json(status, state)
        raise
    finally:
        signal.signal(signal.SIGTERM, previous_sigterm)


def verify_output(output):
    output = Path(output).resolve()
    state = json.loads((output/'RUN_STATUS.json').read_text())
    require(state['status'] == 'PASS_REAL_FASTQ_TO_PROFILE_REFERENCE_MATCH', 'Run has no final scientific success evidence')
    for row in csv.DictReader((output/'OUTPUT_SHA256.tsv').open(), delimiter='\t'):
        path = (output/row['path']).resolve()
        require(path.is_relative_to(output) and path.is_file(), 'Manifest path escapes output or is missing')
        require(path.stat().st_size == int(row['bytes']) and digest(path) == row['sha256'], 'Output integrity mismatch: '+row['path'])
    require(compare_outputs(output/'sketch', output/'query/ERR14036375_all.txt', output/'profile/ERR14036375_profile.txt')['status'] == 'PASS', 'Independent readback failed')
    return {'status': 'PASS_OUTPUT_INTEGRITY_AND_REFERENCE_MATCH', 'run_accession': 'ERR14036375'}


def sra_plan(sra, fastq_dump):
    # Deliberately a plan, with no SRA download, model execution, cleanup, or PASS claim.
    source = Path(sra).expanduser().resolve()
    return {'status': 'NOT_VERIFIED_SRA_STREAM_PLAN_ONLY', 'execution_supported': False,
            'local_sra_file': str(source), 'local_sra_exists': source.is_file(),
            'historical_pipecmd_tokens': [fastq_dump, '--skip-technical', '--split-spot', '-Z'],
            'historical_pipecmd_string': shlex.join([fastq_dump, '--skip-technical', '--split-spot', '-Z']),
            'historical_input_form': 'A local SRA path passed as the positional input after --pipecmd; old RUN could also be an accession, which must not trigger unintended retrieval here',
            'aws_download_command': None,
            'boundary': 'Full AWS download/cleanup wrapper is not reconstructed from log messages; only recorded source commands may be documented. No SRA or AWS operation is performed.'}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest='mode', required=True)
    sub.add_parser('data-info', help='Print the fixed one-Run URLs/checksums; no network')
    pp = sub.add_parser('preflight', help='Verify paths, versions, inputs, cores and references; no science execution')
    pp.add_argument('--config', required=True)
    rp = sub.add_parser('run', help='Run this single real example in a new isolated output directory')
    rp.add_argument('--config', required=True); rp.add_argument('--output', required=True)
    rp.add_argument('--keep-scratch', action='store_true')
    vp = sub.add_parser('verify', help='Read-only independent output/hash/reference validation')
    vp.add_argument('--output', required=True)
    sp = sub.add_parser('plan-sra', help='Inspect the historical SRA input form; never executes it')
    sp.add_argument('--sra-file', required=True); sp.add_argument('--fastq-dump', default='fastq-dump')
    args = ap.parse_args()
    if args.mode == 'data-info': result = pins()
    elif args.mode == 'preflight': result = preflight(args.config)[1]
    elif args.mode == 'run': result = execute(args.config, args.output, args.keep_scratch)
    elif args.mode == 'verify': result = verify_output(args.output)
    else: result = sra_plan(args.sra_file, args.fastq_dump)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, RuntimeError, subprocess.TimeoutExpired, EOFError, KeyError) as exc:
        print('ERROR: '+str(exc), file=sys.stderr)
        sys.exit(1)
