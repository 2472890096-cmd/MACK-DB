#!/usr/bin/env python3
"""Portable complete-project matrix example; delegates unchanged scientific code.

Default invocation checks prerequisites only. --steps all explicitly includes
new DKI training/prediction. No downloads, installations, saved-prediction
substitution, parameter tuning, production writes or automatic retries.
"""
from pathlib import Path
import argparse
import csv
import fcntl
import hashlib
import json
import math
import os
import platform
import shutil
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
BUNDLE = HERE / 'bundled'
PROJECT = 'PRJEB81106'
STEPS = ['prepare', 'fastspar', 'fastcclasso', 'zipi', 'epi',
         'dki-input', 'dki-train', 'dki-score']
REQUIRES = {'prepare': [], 'fastspar': ['prepare'], 'fastcclasso': ['prepare'],
            'zipi': ['fastspar'], 'epi': ['prepare'], 'dki-input': ['prepare'],
            'dki-train': ['dki-input'], 'dki-score': ['dki-train']}
FASTCC_SHA = '41004c99db932346c651447ef561d7965fddee0dc3252aad63ac6681b57a8638'
CONVERTER_SHA = 'dcc828b84119f8a6ac842bec27a2a4e675c2e99a485b655d3edac9eab7ebd790'


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def atomic_json(path, data):
    path = Path(path)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def plan_steps(requested):
    wanted = set()
    def visit(stage):
        if stage in wanted:
            return
        if stage not in REQUIRES:
            raise ValueError('Unknown workflow step: ' + stage)
        for dep in REQUIRES[stage]:
            visit(dep)
        wanted.add(stage)
    for stage in requested:
        if stage == 'all':
            for item in STEPS:
                visit(item)
        else:
            visit(stage)
    return [s for s in STEPS if s in wanted]


def checked_bundle():
    manifest = BUNDLE / 'SHA256SUMS.tsv'
    for row in read_tsv(manifest):
        p = (BUNDLE / row['path']).resolve()
        if BUNDLE.resolve() not in p.parents or not p.is_file():
            raise ValueError('Invalid/missing bundled path: ' + row['path'])
        if p.stat().st_size != int(row['bytes']) or digest(p) != row['sha256']:
            raise ValueError('Bundled input checksum mismatch: ' + row['path'])
    runs = (BUNDLE / 'metadata/runs.txt').read_text().split()
    if len(runs) != len(set(runs)) or len(runs) != 33:
        raise ValueError('Expected exactly the archived33 unique Runs')
    sm = read_tsv(BUNDLE / 'metadata/sample_mapping.tsv')
    if [r['run'] for r in sm] != runs or [r['sample_index'] for r in sm] != list(map(str, range(33))):
        raise ValueError('Run order and matrix-column mapping disagree')
    for level, taxa in [('species', 53), ('genus', 38)]:
        matrix, mapping = input_rows(level)
        if len(matrix) != taxa + 1 or len(mapping) != taxa:
            raise ValueError('Input matrix/mapping dimension mismatch: ' + level)
        ids = list(map(str, range(taxa)))
        if matrix[0][1:] != list(map(str, range(33))) or [r[0] for r in matrix[1:]] != ids or [r[0] for r in mapping] != ids:
            raise ValueError('Ordered row/column identifiers disagree: ' + level)
        if any(len(r) != 2 for r in mapping) or len({r[1] for r in mapping}) != taxa:
            raise ValueError('Invalid/duplicate mapping names: ' + level)
        if any(len(r) != 34 or any(not math.isfinite(float(v)) or float(v) < 0 for v in r[1:]) for r in matrix[1:]):
            raise ValueError('Invalid abundance values or width: ' + level)
    return digest(manifest)


def input_rows(level):
    with (BUNDLE / 'inputs' / level / 'abundance.matrix').open() as f:
        matrix = list(csv.reader(f, delimiter='\t'))
    with (BUNDLE / 'inputs' / level / 'taxa.tsv').open() as f:
        mapping = list(csv.reader(f, delimiter='\t'))
    return matrix, mapping


def checked_scientific_sources():
    rows = read_tsv(ROOT / 'docs/RELEASE_SCRIPT_MANIFEST.tsv')
    for row in rows:
        if digest(ROOT / row['path']) != row['sha256']:
            raise ValueError('Registered source modified: ' + row['path'])
    converter = ROOT / 'pipeline/abundance/SpeciesAbundanceMatrix_TO_GenusAbundanceMatrix.py'
    if digest(converter) != CONVERTER_SHA:
        raise ValueError('Author-selected original genus converter changed')
    return digest(ROOT / 'docs/RELEASE_SCRIPT_MANIFEST.tsv')


def executable(value):
    found = shutil.which(str(value))
    if not found:
        raise FileNotFoundError('Executable unavailable: ' + str(value))
    return str(Path(found).resolve())


def tool_config(args):
    cfg = {}
    base = ROOT
    if args.config:
        cfg = json.loads(args.config.read_text())
        base = args.config.resolve().parent
    def value(key, default):
        item = str(cfg.get(key, default))
        return str((base / item).resolve()) if '/' in item and not Path(item).is_absolute() else item
    tools = dict(python=args.python or value('python', os.environ.get('MACKDB_PYTHON', sys.executable)),
                 rscript=args.rscript or value('rscript', os.environ.get('MACKDB_RSCRIPT', 'Rscript')),
                 dki_python=args.dki_python or value('dki_python', os.environ.get('MACKDB_DKI_PYTHON', sys.executable)),
                 fastspar=value('fastspar', os.environ.get('FASTSPAR_BIN', 'fastspar')),
                 fastspar_bootstrap=value('fastspar_bootstrap', os.environ.get('FASTSPAR_BOOTSTRAP_BIN', 'fastspar_bootstrap')),
                 fastspar_pvalues=value('fastspar_pvalues', os.environ.get('FASTSPAR_PVALUES_BIN', 'fastspar_pvalues')),
                 parallel=value('parallel', os.environ.get('PARALLEL_BIN', 'parallel')),
                 fastcclasso_core=str(args.fastcclasso_core.resolve()) if args.fastcclasso_core else value('fastcclasso_core', os.environ.get('FASTCCLASSO_CORE', str(ROOT / 'external/fastCCLasso/RCode/fastCCLasso.R'))))
    if args.fastspar_bin_dir:
        for name in ['fastspar', 'fastspar_bootstrap', 'fastspar_pvalues']:
            tools[name] = str((args.fastspar_bin_dir / name).resolve())
    # Resolve available programs consistently even when their stage is not yet
    # requested. Adding a later stage must not change a symlink's resume identity.
    # Missing programs remain stage-specific preflight errors, not global ones.
    for name in tools:
        if name != 'fastcclasso_core':
            found = shutil.which(tools[name])
            if found:
                tools[name] = str(Path(found).resolve())
    return tools


def probe(argv):
    result = subprocess.run(list(map(str, argv)), capture_output=True, text=True, timeout=60)
    if result.returncode:
        raise RuntimeError('Dependency probe failed: ' + json.dumps(list(map(str, argv))) + '\n' + result.stderr[-2000:])
    return (result.stdout + result.stderr).strip()


def preflight(tools, stages):
    if platform.system() != 'Linux':
        raise RuntimeError('This workflow supports Linux; other operating systems are not validated')
    facts = {'platform': platform.platform(), 'steps': stages, 'installed_nothing': True}
    tools['python'] = executable(tools['python'])
    facts['python'] = probe([tools['python'], '-c', 'import sys,numpy,pandas;assert sys.version_info>=(3,11);print(sys.version);print("numpy",numpy.__version__,"pandas",pandas.__version__)'])
    rpkgs = set()
    if any(s in stages for s in ['fastspar', 'fastcclasso', 'zipi']):
        rpkgs |= {'optparse', 'tidyverse'}
    if 'fastspar' in stages:
        rpkgs.add('VennDiagram')
        for name in ['fastspar', 'fastspar_bootstrap', 'fastspar_pvalues']:
            tools[name] = executable(tools[name])
            version = probe([tools[name], '--version'])
            if '1.0.0' not in version:
                raise ValueError(name + ' must be the documented1.0.0 build')
            facts[name] = version
        tools['parallel'] = executable(tools['parallel'])
        facts['parallel'] = probe([tools['parallel'], '--version']).splitlines()[0]
    if 'fastcclasso' in stages:
        path = Path(tools['fastcclasso_core'])
        if not path.is_file() or digest(path) != FASTCC_SHA:
            raise ValueError('Obtain the pinned external fastCCLasso core; see environments/portable-matrix/README.md')
        facts['fastcclasso_sha256'] = FASTCC_SHA
    if 'epi' in stages:
        rpkgs.add('vegan')
    if 'zipi' in stages:
        rpkgs |= {'microeco', 'igraph', 'htmlwidgets', 'networkD3', 'gridExtra', 'ggrepel', 'RColorBrewer', 'chorddiag'}
    if 'dki-score' in stages:
        rpkgs.add('dplyr')
    if rpkgs:
        tools['rscript'] = executable(tools['rscript'])
        code = 'p<-c(' + ','.join(json.dumps(p) for p in sorted(rpkgs)) + ');for(x in p){stopifnot(requireNamespace(x,quietly=TRUE));cat(x,as.character(packageVersion(x)),"\\n")};cat("R",as.character(getRversion()),"\\n")'
        facts['R_packages'] = probe([tools['rscript'], '--vanilla', '-e', code])
    if 'dki-train' in stages:
        tools['dki_python'] = executable(tools['dki_python'])
        facts['DKI_runtime'] = probe([tools['dki_python'], '-c', 'import torch,torchdiffeq,numpy;print("torch",torch.__version__,"torchdiffeq",torchdiffeq.__version__,"numpy",numpy.__version__)'])
    return facts


def artifact_records(folder):
    return [dict(path=str(p.relative_to(folder)), bytes=p.stat().st_size, sha256=digest(p))
            for p in sorted(folder.rglob('*')) if p.is_file() and p.name not in {'SUCCESS.json', 'FAILED.json'}]


def verify_stage(folder):
    record = json.loads((folder / 'SUCCESS.json').read_text())
    if record['status'] != 'COMPLETE_VALIDATED':
        raise ValueError('Stage has no validated completion record: ' + str(folder))
    for item in record['artifacts']:
        p = folder / item['path']
        if not p.is_file() or p.stat().st_size != item['bytes'] or digest(p) != item['sha256']:
            raise ValueError('Completed stage artifact changed/missing: ' + str(p))
    return record


def finite_csv(path, rows, cols, normalized=False):
    with Path(path).open() as f:
        data = list(csv.reader(f))
    if len(data) != rows or any(len(r) != cols for r in data):
        raise ValueError('Prediction/input CSV dimensions invalid: ' + str(path))
    values = [[float(v) for v in r] for r in data]
    if any(not math.isfinite(v) or v < 0 for r in values for v in r):
        raise ValueError('Invalid negative/nonfinite prediction: ' + str(path))
    if normalized and any(abs(sum(r) - 1) > 1e-5 for r in values):
        raise ValueError('Prediction rows are not normalized: ' + str(path))


def numerical_comparison(actual, reference):
    with Path(actual).open() as f:
        a = list(csv.reader(f, delimiter='\t'))
    with Path(reference).open() as f:
        b = list(csv.reader(f, delimiter='\t'))
    if len(a) != len(b) or any(len(x) != len(y) for x, y in zip(a, b)):
        return {'same_shape': False, 'reference_is_not_a_new_output_input': True}
    finite_pairs = []
    same_missing = True
    for x, y in zip(a[1:], b[1:]):
        for u, v in zip(x[1:], y[1:]):
            try:
                u, v = float(u), float(v)
            except ValueError:
                same_missing &= u in {'NA', 'NaN', ''} and v in {'NA', 'NaN', ''}
                continue
            if math.isfinite(u) and math.isfinite(v):
                finite_pairs.append(abs(u-v))
            else:
                same_missing &= not math.isfinite(u) and not math.isfinite(v)
    return dict(same_shape=True, same_missing_pattern=bool(same_missing),
                max_abs_difference=max(finite_pairs, default=0),
                reference_is_not_a_new_output_input=True)


class Workflow:
    def __init__(self, args, tools, facts, stages, bundle_sha, source_sha):
        self.args, self.tools, self.facts, self.stages = args, tools, facts, stages
        self.out = args.output.expanduser().resolve()
        if ROOT in self.out.parents and ROOT / 'work' not in self.out.parents:
            raise ValueError('Checkout outputs must be under work/; inputs/source must remain unchanged')
        self.identity = dict(project=PROJECT, levels=args.levels, bundle_sha256=bundle_sha,
                             scientific_source_manifest_sha256=source_sha, tools=tools,
                             formal_parameters='UNCHANGED_ACCEPTED_SCIENTIFIC_PROGRAMS')
        existed = self.out.exists()
        if existed and not args.resume:
            raise FileExistsError('Use a new output directory, or --resume for unchanged completed stages')
        self.out.mkdir(parents=True, exist_ok=True)
        self.lock = (self.out / '.workflow.lock').open('a')
        fcntl.flock(self.lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        pin = self.out / 'RUN_IDENTITY.json'
        if existed:
            if not pin.is_file() or json.loads(pin.read_text()) != self.identity:
                raise ValueError('Resume identity differs; do not mix inputs, code or tool paths')
        else:
            atomic_json(pin, self.identity)
        atomic_json(self.out / 'PREFLIGHT.json', facts)
        self.env = os.environ.copy()
        self.env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                        NUMEXPR_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg')

    def stage_path(self, level, stage):
        return self.out / level / stage

    def command(self, argv, folder):
        remain = self.deadline-time.monotonic()
        if remain <= 0:
            raise TimeoutError('Configured stage deadline reached')
        argv = list(map(str, argv))
        with (folder / 'commands.jsonl').open('a') as f:
            f.write(json.dumps({'argv': argv, 'utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}) + '\n')
        with (folder / 'execution.log').open('a') as log:
            proc = subprocess.Popen(argv, cwd=folder, env=self.env, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            try:
                code = proc.wait(timeout=remain)
            except (subprocess.TimeoutExpired, KeyboardInterrupt):
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
                raise TimeoutError('Child process group stopped; no retry or parameter reduction')
            if code:
                raise RuntimeError('Scientific program exited ' + str(code) + '; inspect ' + str(folder / 'execution.log'))

    def local_inputs(self, level):
        return self.stage_path(level, 'prepare') / 'inputs'

    def execute(self, level, stage):
        folder = self.stage_path(level, stage)
        if (folder / 'SUCCESS.json').is_file():
            old = verify_stage(folder)
            print(json.dumps(dict(level=level, stage=stage, status='REUSED_VALIDATED', results=old['results'])), flush=True)
            return old
        if folder.exists():
            raise FileExistsError('Incomplete prior stage retained; use a fresh output root: ' + str(folder))
        folder.mkdir(parents=True)
        self.deadline = time.monotonic() + self.args.timeout_seconds
        started = time.monotonic()
        print(json.dumps(dict(level=level, stage=stage, status='STARTED')), flush=True)
        results = {}
        try:
            local = self.local_inputs(level)
            mapping = local / level / f'{PROJECT}_{level}.matrix'
            abundance = local / 'abundance' / f'{PROJECT}_abundance.matrix'
            projectlist = folder / 'project.txt'
            projectlist.write_text(PROJECT + '\n')
            if stage == 'prepare':
                target = folder / 'inputs'
                for name in ['abundance', 'species', 'genus']:
                    (target / name).mkdir(parents=True)
                if level == 'species':
                    shutil.copyfile(BUNDLE / 'inputs/species/abundance.matrix', abundance)
                    shutil.copyfile(BUNDLE / 'inputs/species/taxa.tsv', mapping)
                else:
                    legacy = folder / 'converter_input'
                    legacy.mkdir()
                    shutil.copyfile(BUNDLE / 'inputs/species/abundance.matrix', legacy / f'{PROJECT}_profile.txt')
                    shutil.copyfile(BUNDLE / 'inputs/species/taxa.tsv', legacy / f'{PROJECT}_species.txt')
                    self.command([self.tools['python'], ROOT / 'pipeline/abundance/SpeciesAbundanceMatrix_TO_GenusAbundanceMatrix.py', '-i', projectlist, '-d', legacy, '-o', target], folder)
                    produced = target / 'abundance' / f'{PROJECT}_abundance.txt'
                    names = target / 'genus' / f'{PROJECT}_genus.txt'
                    shutil.copyfile(produced, abundance)
                    shutil.copyfile(names, mapping)
                    # Numerical equality is required; newline/float text representation is not science.
                    actual = [r for r in csv.reader(abundance.open(), delimiter='\t')]
                    expected, names_expected = input_rows('genus')
                    names_actual = list(csv.reader(mapping.open(), delimiter='\t'))
                    if names_actual != names_expected or actual[0] != expected[0] or len(actual) != len(expected):
                        raise ValueError('Original genus converter labels/order differ from archive')
                    if any(a[0] != b[0] or len(a) != len(b) or any(float(u) != float(v) for u, v in zip(a[1:], b[1:])) for a, b in zip(actual[1:], expected[1:])):
                        raise ValueError('Original genus converter does not reproduce every archived cell')
                    results['all_1254_genus_cells_exact'] = True
                results.update(samples=33, taxa=53 if level == 'species' else 38)
            elif stage == 'fastspar':
                for name, var in [('fastspar', 'FASTSPAR_BIN'), ('fastspar_bootstrap', 'FASTSPAR_BOOTSTRAP_BIN'), ('fastspar_pvalues', 'FASTSPAR_PVALUES_BIN'), ('parallel', 'PARALLEL_BIN')]:
                    self.env[var] = self.tools[name]
                self.env['MACKDB_THREADS'] = str(self.args.threads)
                count = folder / f'{PROJECT}_count.matrix'
                self.command([self.tools['rscript'], ROOT / 'pipeline/abundance/count_sparcc.R', '-i', abundance, '-n', '1000000', '-o', count], folder)
                # The unchanged GNU parallel shell entry embeds the prefix in
                # command text. Use a simple relative prefix inside this owned
                # stage directory so checkout/output paths may contain spaces.
                # The scientific arguments and resulting physical files agree.
                self.command(['bash', ROOT / 'pipeline/associations/fastspar/run_fastspar.sh', count, PROJECT], folder)
                csvdir = folder / 'network_csv'
                corr, pv = [folder / f'{PROJECT}_{suffix}.tsv' for suffix in ['correlation', 'pvalues']]
                self.command([self.tools['python'], ROOT / 'tools/network_matrix_io.py', '--corr', corr, '--pvalue', pv, '--mapping', mapping, '--output', csvdir], folder)
                self.command([self.tools['rscript'], ROOT / 'pipeline/associations/fastspar/sparcc_corr_pvalue.R', '-c', csvdir / 'correlation.csv', '-p', csvdir / 'pvalues.csv', '-s', mapping, '--corr_num', '0.3', '--pvalue_num', '0.0001', '-o', folder / 'selected'], folder)
                self.command([self.tools['rscript'], ROOT / 'pipeline/associations/fastspar/corr_inter.R', '-f', folder / 'selected.csv', '-n', folder / PROJECT], folder)
                results['selected_edges'] = len(list(csv.DictReader((folder / 'selected.csv').open())))
                results['reference_comparison'] = {name: numerical_comparison(folder / f'{PROJECT}_{name}.tsv', BUNDLE / 'expected' / level / f'fastspar_{name}.tsv') for name in ['correlation', 'covariance', 'pvalues']}
            elif stage == 'fastcclasso':
                runner = 'run_fastcclasso_parallel' + ('_genus' if level == 'genus' else '') + '.R'
                self.command([self.tools['rscript'], ROOT / 'pipeline/associations/fastcclasso' / runner, '--project', PROJECT, '--abundance_dir', local / 'abundance', '--species_dir', local / level, '--output', folder / 'fit', '--fast_core', self.tools['fastcclasso_core'], '--nboot', '100', '--kmax', '20', '--pvalue', '0.0001', '--cor_cutoff', '0.1', '--workers', '1', '--seed', '20260804' if level == 'species' else '20260806', '--write_full_matrices', 'T', '--write_pair_csv', 'T'], folder)
                self.command([self.tools['python'], ROOT / 'pipeline/associations/fastcclasso/refilter_fastcclasso_cor0p3.py', '--level', level, '--work-root', folder / 'fit'], folder)
                record = read_tsv(folder / 'fit/project_status' / f'{PROJECT}.tsv')[0]
                if record['status'] != 'pass':
                    raise ValueError('fastCCLasso did not report project pass')
                results['project_status'] = record
            elif stage == 'zipi':
                network = self.stage_path(level, 'fastspar')
                verify_stage(network)
                self.command([self.tools['rscript'], ROOT / 'pipeline/zipi/FastSpar_signed_spinglass_network.R', '-p', PROJECT, '-l', level, '-C', network / f'{PROJECT}_correlation.tsv', '-P', network / f'{PROJECT}_pvalues.tsv', '--species-base', local, '--genus-base', local, '--taxonomy-file', BUNDLE / 'reference' / f'{level}_taxonomy.tsv', '-f', '0.3', '--spins', '25', '--gamma', '1', '--gamma-minus', '1', '--seed', '12345', '-a', folder / 'plots', '-b', folder / 'network_html', '-c', folder / 'chord_html', '-g', folder / 'gexf', '-d', folder / 'nodes'], folder)
                if '节点角色计算失败' in (folder / 'execution.log').read_text():
                    raise ValueError('ZiPi role calculation failed; fallback basic nodes are not complete scientific output')
                nodepath = folder / 'nodes' / f'{PROJECT}_node.txt'
                with nodepath.open() as f:
                    reader = csv.DictReader(f, delimiter='\t')
                    if not {'name', 'z', 'p', 'taxa_roles'} <= set(reader.fieldnames or []):
                        raise ValueError('ZiPi core node fields are missing')
                    nodes = list(reader)
                if len({r['name'] for r in nodes}) != len(nodes):
                    raise ValueError('Duplicate ZiPi node names')
                for r in nodes:
                    role = r['taxa_roles'].strip().lower()
                    if r['z'] in {'NA', 'NaN', ''} or r['p'] in {'NA', 'NaN', ''}:
                        if role not in {'na', 'nan', ''}:
                            raise ValueError('Role declared for undefined Zi/Pi')
                        continue
                    z, p = float(r['z']), float(r['p'])
                    if not math.isfinite(z) or not math.isfinite(p) or not -1e-12 <= p <= 1 + 1e-12:
                        raise ValueError('Invalid finite Zi/Pi value')
                    expected = ('network hub' if p > .62 else 'module hub') if z > 2.5 else ('connector' if p > .62 else 'peripheral node')
                    if role.rstrip('s') != expected:
                        raise ValueError('Saved ZiPi role contradicts saved Zi/Pi')
                envdir = folder / 'environment'
                envdir.mkdir()
                (envdir / 'human_gut.list').write_text(PROJECT + '\n')
                self.command([self.tools['python'], ROOT / 'pipeline/zipi/extract_signed_zipi_keys_by_environment.py', '--level', level, '--env-dir', envdir, '--' + level + '-node-dir', folder / 'nodes', '--' + level + '-out-dir', folder / 'keys', '--write-empty'], folder)
                self.command([self.tools['python'], ROOT / 'pipeline/zipi/prepare_signed_keys_for_mysql.py', '--level', level, '--' + level + '-dir', folder / 'keys', '--output-dir', folder / 'candidates'], folder)
                summary = read_tsv(folder / 'keys/_project_extraction_summary.tsv')[0]
                if summary['status'] not in {'OK', 'EMPTY_WRITTEN'}:
                    raise ValueError('ZiPi extraction was not successful: ' + str(summary))
                results.update(nodes=len(nodes), extracted_candidates=int(summary['key_rows']), legal_empty_network=not bool(nodes))
            elif stage == 'epi':
                self.command([self.tools['python'], ROOT / 'tools/epi_project.py', '--level', level, '--project', PROJECT, '--abundance', abundance, '--mapping', mapping, '--run-root', folder / 'epi', '--rscript', self.tools['rscript']], folder)
                row = read_tsv(folder / 'epi/status/project_status.tsv')[0]
                if row['status'] != 'SUCCESS' or row['exit_code'] != '0':
                    raise ValueError('EPI lacks complete success evidence')
                results['project_status'] = row
                for metric in ['D1', 'D2', 'Q']:
                    path = folder / 'epi/keys' / f'{PROJECT}_EPI-{metric}.key'
                    if not path.is_file():
                        raise ValueError('EPI missing metric: ' + metric)
            elif stage == 'dki-input':
                mapdir = folder / 'mapping'
                mapdir.mkdir()
                shutil.copyfile(mapping, mapdir / f'{PROJECT}_species.matrix')
                self.command([self.tools['python'], ROOT / 'pipeline/dki/generate_dki_input_from_abundance.py', '--project', PROJECT, '--abundance-dir', local / 'abundance', '--species-dir', mapdir, '--out-base', folder / 'data_base'], folder)
                data = folder / 'data_base' / PROJECT / 'data'
                for name in ['Ptrain', 'Ptest', 'species_id', 'sample_id']:
                    if not (data / (name + '.csv')).is_file():
                        raise ValueError('Missing DKI perturbation file: ' + name)
                results['generated_perturbations_from_matrix'] = True
            elif stage == 'dki-train':
                data = self.stage_path(level, 'dki-input') / 'data_base' / PROJECT / 'data'
                pred = folder / 'pred_base' / PROJECT / 'results'
                self.command([self.tools['dki_python'], ROOT / 'pipeline/dki/DKI_formula_fixed.py', '--train-file', data / 'Ptrain.csv', '--test-file', data / 'Ptest.csv', '--output-dir', pred, '--project', PROJECT, '--max-epochs', '1000', '--early-stop-patience', '80', '--lr', '0.005', '--batch-size', '32', '--activation', 'relu', '--bottleneck', 'auto', '--normalize-z', 'true', '--time-mode', 'endpoint', '--ode-method', 'rk4', '--ode-step-size', '0.05', '--device', 'cpu', '--seed', '20260629'], folder)
                metrics = json.loads((pred / 'model_metrics.json').read_text())
                ntaxa = 53 if level == 'species' else 38
                with (data / 'Ptest.csv').open() as f:
                    cases = len(next(csv.reader(f)))
                finite_csv(pred / 'qtrn.csv', 33, ntaxa, normalized=True)
                finite_csv(pred / 'qtst.csv', cases, ntaxa, normalized=True)
                if metrics['project'] != PROJECT or metrics['N_species'] != ntaxa or metrics['seed'] != 20260629 or metrics['ode_method'] != 'rk4' or metrics['ode_step_size'] != .05:
                    raise ValueError('DKI model metrics contradict production parameters')
                results.update(trained_model=True, generated_neural_predictions=True, model_metrics=metrics,
                               historical_prediction_files_used=False)
            elif stage == 'dki-score':
                verify_stage(self.stage_path(level, 'dki-train'))
                data_base = self.stage_path(level, 'dki-input') / 'data_base'
                pred_base = folder / 'pred_base'
                pred = pred_base / PROJECT / 'results'
                pred.mkdir(parents=True)
                fresh = self.stage_path(level, 'dki-train') / 'pred_base' / PROJECT / 'results'
                for name in ['qtrn.csv', 'qtst.csv']:
                    shutil.copyfile(fresh / name, pred / name)
                self.command([self.tools['rscript'], ROOT / 'pipeline/dki/calculate_structural_keystoneness_rk4_fast.R', '--project-list', projectlist, '--data-base', data_base, '--pred-base', pred_base, '--threads', '1'], folder)
                status = read_tsv(pred_base / 'reports/ks_rk4_fast_status.tsv')[0]
                if status['status'] != 'DONE':
                    raise ValueError('DKI scorer did not finish: ' + str(status))
                argv = [self.tools['python'], ROOT / 'pipeline/dki' / f'map_keystoneness_{level}_names.py', '--project-list', projectlist, '--pred-base', pred_base, '--status-file', folder / 'mapping_status.tsv']
                argv += ['--dki-base', data_base, '--species-dir', local / 'species'] if level == 'species' else ['--genus-dir', local / 'genus']
                self.command(argv, folder)
                mapped = read_tsv(folder / 'mapping_status.tsv')[0]
                if mapped['status'] != 'DONE' or int(mapped['unmapped_n']) != 0:
                    raise ValueError('DKI taxon-name mapping incomplete: ' + str(mapped))
                envdir = folder / 'environment'
                envdir.mkdir()
                (envdir / 'human_gut.list').write_text(PROJECT + '\n')
                self.command([self.tools['python'], ROOT / 'pipeline/dki/prepare_dki_keys_for_mysql.py', '--level', level, '--' + level + '-root', pred_base, '--env-dir', envdir, '--output-dir', folder / 'candidates', '--sd-mode', 'sample'], folder)
                selection = read_tsv(folder / 'candidates' / f'key_{level}_dki_project_summary.tsv')[0]
                if selection['status'] not in {'OK', 'NO_CANDIDATES'} or int(selection['N']) != (53 if level == 'species' else 38):
                    raise ValueError('DKI candidate selection incomplete: ' + str(selection))
                results.update(scorer_status=status, selection_status=selection,
                               fresh_model_predictions=True, historical_predictions_used=False)
            record = dict(project=PROJECT, level=level, stage=stage, status='COMPLETE_VALIDATED',
                          elapsed_seconds=time.monotonic()-started, results=results,
                          artifact_validation='Core exits/status + documented schema; reference equality separate',
                          artifacts=artifact_records(folder))
            atomic_json(folder / 'SUCCESS.json', record)
            print(json.dumps(dict(level=level, stage=stage, status='COMPLETE_VALIDATED', elapsed_seconds=record['elapsed_seconds'], results=results)), flush=True)
            return record
        except BaseException as exc:
            atomic_json(folder / 'FAILED.json', dict(project=PROJECT, stage=stage, level=level,
                        status='TIMEOUT' if isinstance(exc, TimeoutError) else 'FAILED',
                        error=type(exc).__name__ + ': ' + str(exc), elapsed_seconds=time.monotonic()-started))
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--steps', nargs='+', choices=['preflight', 'all'] + STEPS, default=['preflight'])
    parser.add_argument('--levels', nargs='+', choices=['species', 'genus'], default=['species', 'genus'])
    parser.add_argument('--preflight', action='store_true', help='Only check prerequisites for --steps')
    parser.add_argument('--output', type=Path, default=ROOT / 'work/PRJEB81106_portable')
    parser.add_argument('--config', type=Path)
    parser.add_argument('--python')
    parser.add_argument('--rscript')
    parser.add_argument('--dki-python')
    parser.add_argument('--fastspar-bin-dir', type=Path)
    parser.add_argument('--fastcclasso-core', type=Path)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--timeout-seconds', type=float, default=600, help='Finite timeout per stage; increase explicitly on slower servers; no retries')
    parser.add_argument('--resume', action='store_true', help='Reuse only checksum-validated completed stages with identical inputs/code/tool paths')
    args = parser.parse_args()
    if not 1 <= args.threads <= 64 or not math.isfinite(args.timeout_seconds) or args.timeout_seconds <= 0:
        parser.error('Positive finite stage timeout and1..64threads are required')
    if len(args.levels) != len(set(args.levels)):
        parser.error('Duplicate levels')
    checking = args.preflight or args.steps == ['preflight']
    if 'preflight' in args.steps and args.steps != ['preflight']:
        parser.error('Use --preflight --steps ... for a subset of dependency checks')
    stages = plan_steps(['all'] if args.steps == ['preflight'] else args.steps)
    bundle_sha = checked_bundle()
    source_sha = checked_scientific_sources()
    tools = tool_config(args)
    facts = preflight(tools, stages)
    if checking:
        print(json.dumps(dict(status='PREFLIGHT_PASS', bundle_sha256=bundle_sha, **facts), indent=2))
        return
    flow = Workflow(args, tools, facts, stages, bundle_sha, source_sha)
    done = []
    try:
        for level in args.levels:
            for stage in stages:
                done.append(flow.execute(level, stage))
        atomic_json(flow.out / 'WORKFLOW_STATUS.json', dict(project=PROJECT, status='REQUESTED_STEPS_COMPLETE_VALIDATED',
                    levels=args.levels, stages=stages, input='BUNDLED_COMPLETE33RUN_MATRIX',
                    fresh_DKI_training_requested='dki-train' in stages, downloads=False,
                    stages_complete=[{'level': r['level'], 'stage': r['stage']} for r in done]))
    except BaseException as exc:
        atomic_json(flow.out / 'WORKFLOW_STATUS.json', dict(project=PROJECT, status='FAILED', error=str(exc),
                    completed=[{'level': r['level'], 'stage': r['stage']} for r in done]))
        raise


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        print(type(exc).__name__ + ': ' + str(exc), file=sys.stderr)
        sys.exit(1)
