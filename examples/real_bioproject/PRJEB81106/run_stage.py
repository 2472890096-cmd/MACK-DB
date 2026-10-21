#!/usr/bin/env python3
"""Bounded IO/command orchestration for the archived PRJEB81106 matrix example.
All scientific calculations are delegated to existing accepted pipeline code.
No downloads, installations, model training, prediction or automatic retries.
"""
from pathlib import Path
import argparse, csv, hashlib, json, math, os, shutil, signal, subprocess, sys, time

ROOT = Path(__file__).resolve().parents[3]
PROJECT = 'PRJEB81106'


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', required=True, choices=['validate', 'count', 'fastspar',
                        'fastspar-post', 'fastcclasso', 'epi', 'zipi', 'dki-input', 'dki-replay'])
    parser.add_argument('--level', choices=['species', 'genus'], required=True)
    parser.add_argument('--inputs', type=Path, required=True,
                        help='User-supplied private JSON of paths, SHA256 and source class')
    parser.add_argument('--output', type=Path, required=True, help='New isolated directory')
    parser.add_argument('--rscript', default=os.environ.get('MACKDB_RSCRIPT', 'Rscript'))
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--timeout-seconds', type=int, default=600)
    parser.add_argument('--fastspar-bin-dir', type=Path)
    parser.add_argument('--fastcclasso-core', type=Path)
    parser.add_argument('--taxonomy', type=Path)
    parser.add_argument('--allow-archived-genus-replay', action='store_true',
                        help='Replay historical genus inputs only; GTDB suffix aggregation remains under review.')
    args = parser.parse_args()
    if not 1 <= args.timeout_seconds <= 600 or not 1 <= args.threads <= 8:
        parser.error('Timeout must be1..600s; threads must be1..8. Scientific parameters are fixed.')
    cfg = json.loads(args.inputs.read_text())
    if cfg['project'] != PROJECT:
        parser.error('This archived example is only PRJEB81106.')
    allowed = {r['role']:r for r in tsv(Path(__file__).with_name('INPUT_CHECKSUMS.tsv'))}
    files = {}
    for role, item in cfg['files'].items():
        if role not in allowed or item['sha256'] != allowed[role]['sha256']:
            parser.error('Input role/hash differs from the archived example manifest: '+role)
        path = Path(item['path']).expanduser()
        if not path.is_absolute():
            path = args.inputs.resolve().parent/path
        path = path.resolve()
        if not path.is_file() or path.stat().st_size != int(allowed[role]['bytes']) or sha(path) != item['sha256']:
            parser.error('Input missing/hash mismatch: '+role)
        files[role] = path
    level = args.level
    if level == 'genus' and args.stage != 'validate' and not args.allow_archived_genus_replay:
        parser.error('Genus GTDB aggregation review is unresolved; use the explicit archived-replay scope only.')
    required = [level+'_abundance', level+'_mapping', 'runs']
    if args.stage in {'fastspar-post', 'zipi'}:
        required += [level+'_fastspar_correlation', level+'_fastspar_pvalues']
    if args.stage == 'dki-replay':
        required += [level+'_dki_'+k for k in ['Ptrain','Ptest','species_id','sample_id','qtrn','qtst']]
    missing = [k for k in required if k not in files]
    if missing:
        parser.error('Missing input roles for this stage: '+', '.join(missing))
    abundance = files[level+'_abundance']
    mapping = files[level+'_mapping']
    runs = files['runs'].read_text().split()
    with abundance.open() as handle:
        rows = list(csv.reader(handle, delimiter='\t'))
    with mapping.open() as handle:
        names = list(csv.reader(handle, delimiter='\t'))
    taxa = 53 if level == 'species' else 38
    if len(runs) != len(set(runs)) or len(runs) != 33:
        parser.error('Expected33 unique archived runs.')
    if len(rows)-1 != taxa or rows[0][1:] != list(map(str, range(33))):
        parser.error('Archived matrix shape/column IDs mismatch; no automatic transpose.')
    if [r[0] for r in rows[1:]] != [r[0] for r in names]:
        parser.error('Taxon row IDs do not match mapping.')
    if [r[0] for r in names] != list(map(str, range(taxa))) or len({r[1] for r in names}) != taxa:
        parser.error('Mapping must have exact unique0-based IDs and names.')
    if any(len(r) != 34 or any(not math.isfinite(float(x)) or float(x) < 0 for x in r[1:]) for r in rows[1:]):
        parser.error('Matrix contains invalid values or row widths.')
    rscript = shutil.which(args.rscript)
    if args.stage != 'validate' and not rscript:
        parser.error('Rscript unavailable; no installation attempted.')
    out = args.output.expanduser().resolve()
    if ROOT in out.parents and ROOT/'work' not in out.parents:
        parser.error('Inside this checkout, outputs must be under ignored work/; otherwise use an external directory.')
    if any(p == out or out in p.parents for p in files.values()):
        parser.error('Output directory must not contain any input file.')
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    deadline = started+args.timeout_seconds
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               NUMEXPR_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg')
    commands = []
    def command(argv):
        remaining = deadline-time.monotonic()
        if remaining <= 0:
            raise TimeoutError('Stage deadline reached; no automatic retry.')
        argv = list(map(str, argv))
        commands.append(argv)
        with (out/'execution.log').open('a') as log:
            log.write('COMMAND '+json.dumps(argv)+'\n'); log.flush()
            proc = subprocess.Popen(argv, cwd=out, env=env, stdout=log,
                                    stderr=subprocess.STDOUT, start_new_session=True)
            try:
                code = proc.wait(timeout=remaining)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    proc.wait(timeout=3)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL); proc.wait()
                raise TimeoutError('Stopped complete process group at stage deadline; do not auto-retry.')
            if code:
                raise RuntimeError('Command exit '+str(code)+'; inspect execution.log')
    state = dict(project=PROJECT, stage=args.stage, level=level, scientific_parameter_set='FORMAL_UNCHANGED',
                 timeout_seconds=args.timeout_seconds, threads=args.threads, status='RUNNING',
                 trained_model=False, generated_neural_predictions=False,
                 input_sources={k:cfg['files'][k].get('source_class','USER_SUPPLIED') for k in cfg['files']})
    state['input_scope'] = ('ARCHIVED_GENUS_REPLAY_GTDB_SUFFIX_GROUPING_NOT_ACCEPTED_CANONICAL'
                            if level == 'genus' else 'ARCHIVED_SPECIES_MATRIX_REPRODUCTION')
    (out/'INPUTS_USED.json').write_text(json.dumps(cfg, indent=2)+'\n')
    try:
        local = out/'inputs'
        for dirname in ['abundance', 'species', 'genus']:
            (local/dirname).mkdir(parents=True)
        shutil.copy2(abundance, local/'abundance'/f'{PROJECT}_abundance.matrix')
        shutil.copy2(mapping, local/level/f'{PROJECT}_{level}.matrix')
        projectlist = out/'project.txt'; projectlist.write_text(PROJECT+'\n')
        envdir = out/'environment'; envdir.mkdir(); (envdir/'human_gut.list').write_text(PROJECT+'\n')
        def count():
            target = out/f'{PROJECT}_count.matrix'
            command([rscript, ROOT/'pipeline/abundance/count_sparcc.R', '-i', abundance,
                     '-n', '1000000', '-o', target])
            return target
        def post(corr, pv):
            converted = out/'network_csv'
            command([sys.executable, ROOT/'tools/network_matrix_io.py', '--corr', corr,
                     '--pvalue', pv, '--mapping', mapping, '--output', converted])
            command([rscript, ROOT/'pipeline/associations/fastspar/sparcc_corr_pvalue.R',
                     '-c', converted/'correlation.csv', '-p', converted/'pvalues.csv', '-s', mapping,
                     '--corr_num', '.3', '--pvalue_num', '.0001', '-o', out/'selected'])
            command([rscript, ROOT/'pipeline/associations/fastspar/corr_inter.R',
                     '-f', out/'selected.csv', '-n', out/PROJECT])
        if args.stage == 'count':
            count()
        elif args.stage == 'fastspar':
            binaries = args.fastspar_bin_dir
            if not binaries:
                raise FileNotFoundError('--fastspar-bin-dir required; no guessed tool path.')
            for name, var in [('fastspar','FASTSPAR_BIN'), ('fastspar_bootstrap','FASTSPAR_BOOTSTRAP_BIN'),
                              ('fastspar_pvalues','FASTSPAR_PVALUES_BIN')]:
                binary = (binaries/name).resolve()
                if not os.access(binary, os.X_OK):
                    raise FileNotFoundError(binary)
                env[var] = str(binary)
            env['MACKDB_THREADS'] = str(args.threads)
            table = count()
            command(['bash', ROOT/'pipeline/associations/fastspar/run_fastspar.sh', table, out/PROJECT])
            post(out/f'{PROJECT}_correlation.tsv', out/f'{PROJECT}_pvalues.tsv')
        elif args.stage == 'fastspar-post':
            state['scope'] = 'POSTPROCESS_SAVED_NETWORK_NOT_NEW_INFERENCE'
            post(files[level+'_fastspar_correlation'], files[level+'_fastspar_pvalues'])
        elif args.stage == 'fastcclasso':
            if not args.fastcclasso_core or not args.fastcclasso_core.is_file():
                raise FileNotFoundError('--fastcclasso-core required; external core not bundled.')
            expected = '41004c99db932346c651447ef561d7965fddee0dc3252aad63ac6681b57a8638'
            if sha(args.fastcclasso_core) != expected:
                raise ValueError('Official fastCCLasso core differs from accepted hash.')
            runner = 'run_fastcclasso_parallel.R' if level == 'species' else 'run_fastcclasso_parallel_genus.R'
            seed = '20260804' if level == 'species' else '20260806'
            command([rscript, ROOT/'pipeline/associations/fastcclasso'/runner, '--project', PROJECT,
                     '--abundance_dir', local/'abundance', '--species_dir', local/level,
                     '--output', out/'fit', '--fast_core', args.fastcclasso_core,
                     '--nboot', '100', '--kmax', '20', '--pvalue', '.0001', '--cor_cutoff', '.1',
                     '--workers', '1', '--seed', seed, '--write_full_matrices', 'T', '--write_pair_csv', 'T'])
            command([sys.executable, ROOT/'pipeline/associations/fastcclasso/refilter_fastcclasso_cor0p3.py',
                     '--level', level, '--work-root', out/'fit'])
            status = tsv(out/'fit/project_status'/f'{PROJECT}.tsv')[0]
            if status['status'] != 'pass':
                raise ValueError('FastCCLasso project did not report pass: '+str(status))
        elif args.stage == 'epi':
            command([sys.executable, ROOT/'tools/epi_project.py', '--level', level,
                     '--project', PROJECT, '--abundance', abundance, '--mapping', mapping,
                     '--run-root', out/'epi', '--rscript', rscript])
            row = tsv(out/'epi/status/project_status.tsv')[0]
            if row['status'] != 'SUCCESS' or row['exit_code'] != '0':
                raise ValueError('EPI project status is not successful.')
        elif args.stage == 'zipi':
            if not args.taxonomy or not args.taxonomy.is_file():
                raise FileNotFoundError('--taxonomy required; no guessed reference.')
            command([rscript, ROOT/'pipeline/zipi/FastSpar_signed_spinglass_network.R', '-p', PROJECT,
                     '-l', level, '-C', files[level+'_fastspar_correlation'], '-P', files[level+'_fastspar_pvalues'],
                     '--species-base', local, '--genus-base', local, '--taxonomy-file', args.taxonomy,
                     '-f', '.3', '--spins', '25', '--gamma', '1', '--gamma-minus', '1', '--seed', '12345',
                     '-a', out/'plots', '-b', out/'network_html', '-c', out/'chord_html',
                     '-g', out/'gexf', '-d', out/'nodes'])
        elif args.stage == 'dki-input':
            # Genus file-name adaptation only: the generic generator expects _species.matrix.
            if level == 'genus':
                shutil.copy2(mapping, local/'species'/f'{PROJECT}_species.matrix')
            command([sys.executable, ROOT/'pipeline/dki/generate_dki_input_from_abundance.py',
                     '--project', PROJECT, '--abundance-dir', local/'abundance',
                     '--species-dir', local/'species', '--out-base', out/'dki'])
            state['scope'] = 'PERTURBATION_INPUT_GENERATION_ONLY_NO_MODEL'
        elif args.stage == 'dki-replay':
            data = out/'data_base'/PROJECT/'data'; data.mkdir(parents=True)
            pred = out/'pred_base'/PROJECT/'results'; pred.mkdir(parents=True)
            for name in ['Ptrain', 'Ptest', 'species_id', 'sample_id']:
                shutil.copy2(files[f'{level}_dki_{name}'], data/(name+'.csv'))
            for name in ['qtrn', 'qtst']:
                shutil.copy2(files[f'{level}_dki_{name}'], pred/(name+'.csv'))
            command([rscript, ROOT/'pipeline/dki/calculate_structural_keystoneness_rk4_fast.R',
                     '--project-list', projectlist, '--data-base', out/'data_base',
                     '--pred-base', out/'pred_base', '--threads', '1'])
            row = tsv(out/'pred_base/reports/ks_rk4_fast_status.tsv')[0]
            if row['status'] != 'DONE':
                raise ValueError('Scorer did not complete: '+str(row))
            argv = [sys.executable, ROOT/'pipeline/dki'/f'map_keystoneness_{level}_names.py',
                    '--project-list', projectlist, '--pred-base', out/'pred_base', '--status-file', out/'mapping_status.tsv']
            if level == 'species':
                argv += ['--dki-base', out/'data_base', '--species-dir', local/'species']
            else:
                argv += ['--genus-dir', local/'genus']
            command(argv)
            command([sys.executable, ROOT/'pipeline/dki/prepare_dki_keys_for_mysql.py', '--level', level,
                     '--'+level+'-root', out/'pred_base', '--env-dir', envdir,
                     '--output-dir', out/'candidates', '--sd-mode', 'sample'])
            state['scope'] = 'SCORING_REPLAY_OF_SAVED_PREDICTIONS_NO_TRAINING_OR_PREDICTION'
        state['status'] = 'COMMANDS_COMPLETE_OUTPUT_VALIDATION_STILL_REQUIRED'
    except Exception as exc:
        state['status'] = ('TIMEOUT' if isinstance(exc, TimeoutError) else
                           'BLOCKED' if isinstance(exc, FileNotFoundError) else 'FAILED')
        state['error'] = type(exc).__name__+': '+str(exc)
        raise
    finally:
        state.update(elapsed_seconds=time.monotonic()-started, commands=commands)
        (out/'RUN_STATUS.json').write_text(json.dumps(state, indent=2)+'\n')
    print(json.dumps(state, indent=2))


if __name__ == '__main__':
    main()
