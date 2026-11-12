#!/usr/bin/env python3
"""Recalculate paper summaries from archived simulations and support inputs.

No random graphs, network inference or source exclusion is run here. Scientific
summary/rank/test functions are compiled directly from the retained repository
sources; their module initializers and historical main functions are not run.
"""
from __future__ import annotations

import argparse
import ast
import collections
import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
LEGACY = ROOT / 'analysis/external_validation/shared_protocol/chen_wallen_null_models_final_v1.py'
SUPPORT = ROOT / 'analysis/external_validation/literature_networks/stage4d_support_rank/scripts/04d_literature_recurrence_support_rank.py'
SOURCE_SHA256 = {
    LEGACY: 'dcda62e6a6981cb31b177d3afb4efe358cb21c8df61c518b6c7f10d673278bc5',
    SUPPORT: 'e5326f7ba3d5c1ded4519f2da3b8bd279f3e257904680b0693067b5da9b48f67',
}
NETWORKS = ['Zhou_MS', 'Zhou_HHC', 'Guo_AggControl', 'Guo_AggMECFS',
            'Lacruz_MHNO', 'Lacruz_MHO', 'Lacruz_MUNO', 'Lacruz_MUO', 'Li_IBS']
INPUTS: dict[Path, str] = {}
QC: list[dict] = []
OUTPUT: Path
csv.field_size_limit(10000000)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def register(path, role):
    path = Path(path).resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    INPUTS[path] = role
    return path


def read(path):
    with Path(path).open(newline='', encoding='utf-8') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write(path, rows, fields=None):
    path, rows = Path(path).resolve(), list(rows)
    if OUTPUT not in path.parents:
        raise ValueError('Output escapes new output directory')
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def check(name, ok, detail=''):
    QC.append(dict(check=name, status='PASS' if ok else 'FAIL', detail=str(detail)))
    if not ok:
        raise ValueError(name + ': ' + str(detail))


def load_functions(path, names, namespace):
    """Compile named existing functions only; never execute module top-level IO."""
    register(path, 'unchanged_scientific_function_source')
    check('scientific_source_SHA256_' + path.stem, sha(path) == SOURCE_SHA256[path])
    tree = ast.parse(path.read_text())
    functions = {node.name: node for node in tree.body if isinstance(node, ast.FunctionDef)}
    selected = [functions[name] for name in names]
    module = ast.fix_missing_locations(ast.Module(body=selected, type_ignores=[]))
    exec(compile(module, str(path), 'exec'), namespace)
    write(OUTPUT / ('reused_' + path.stem + '.tsv'), [
        dict(function=name, repository_path=str(path.relative_to(ROOT)),
             source_sha256=sha(path), start_line=functions[name].lineno,
             end_line=functions[name].end_lineno,
             function_sha256=hashlib.sha256(ast.get_source_segment(path.read_text(), functions[name]).encode()).hexdigest())
        for name in names])


def close(a, b):
    return (math.isnan(a) and math.isnan(b)) or a == b or math.isclose(a, b, rel_tol=1e-11, abs_tol=1e-12)


def compare_rows(actual, expected, keys, label):
    left = {tuple(str(r[k]) for k in keys): r for r in actual}
    right = {tuple(str(r[k]) for k in keys): r for r in expected}
    check(label + '_unique_keys', len(left) == len(actual) and len(right) == len(expected))
    check(label + '_same_keys', left.keys() == right.keys())
    diffs = []
    numeric = 0
    for key in left:
        for column in left[key].keys() & right[key].keys():
            a, b = str(left[key][column]), str(right[key][column])
            if a == b:
                continue
            try:
                x, y = float(a), float(b)
            except ValueError:
                check(label + '_text_' + column, False, key)
            check(label + '_number_' + column, close(x, y), (key, x, y))
            numeric += 1
            diffs.append(dict(table=label, key=' | '.join(key), column=column,
                              recalculated=a, frozen=b, absolute_difference=abs(x-y)))
    return diffs, numeric


def nulls(args):
    ns = dict(pd=pd, np=np)
    load_functions(LEGACY, ['summarize'], ns)
    lit, shared = args.literature_root.resolve(), args.chen_wallen_root.resolve()
    stage = lit / 'stage3c_null_model'
    summary9 = read(register(stage / 'summary/03c_null_model_summary.tsv', 'frozen_nine_network_summary'))
    summary2 = read(register(shared / '05_null_model_summary.tsv', 'frozen_Chen_Wallen_summary'))
    config = json.loads(register(stage / 'qc/03c_run_config.json', 'nine_network_seed_protocol').read_text())
    check('nine_network_config', config['iterations'] == 1000 and config['master_seed'] == 20260917
          and config['successful_swaps_per_edge'] == 50 and config['network_order'] == NETWORKS)
    summary = []
    comparisons = []
    all_seeds = set()
    seed_keys = set()
    for network in ['Chen', 'Wallen'] + NETWORKS:
        for null in 'AB':
            if network in ['Chen', 'Wallen']:
                frozen = next(r for r in summary2 if r['study'] == network and r['null_model'] == 'Null ' + null)
                prefix = {('Chen', 'A'): '01', ('Chen', 'B'): '02', ('Wallen', 'A'): '03', ('Wallen', 'B'): '04'}[network, null]
                path = shared / f'{prefix}_{network}_Null{null}_replicates.tsv'
                E, observed, n = (int(frozen[k]) for k in ['observed_edge_count', 'observed_recovered_count', 'external_species_count'])
                seed_base = 20260910 + 100000 * (1 + 2 * ['Chen', 'Wallen'].index(network) + (null == 'B'))
                iteration = 'replicate'
            else:
                frozen9 = next(r for r in summary9 if r['network'] == network and r['null_type'] == 'Null ' + null)
                path = stage / ('null_a' if null == 'A' else 'null_b') / f'{network}_null{null}_iterations.tsv'
                E, observed, n = (int(frozen9[k]) for k in ['edge_count', 'observed_recovered', 'node_count'])
                seed_base = 20260917 + 100000 * (2 * NETWORKS.index(network) + (1 if null == 'A' else 2))
                iteration = 'iteration'
            rows = read(register(path, 'stored_1000_randomized_graph_outcomes_' + network + '_' + null))
            check(network + null + '_1000_unique_iterations', len(rows) == len({r[iteration] for r in rows}) == 1000
                  and {int(r[iteration]) for r in rows} == set(range(1, 1001)))
            for r in rows:
                seed = int(r['seed'])
                check(network + null + '_seed', seed == seed_base + int(r[iteration]))
                # Chen/Wallen and the later nine-network extension have distinct
                # masters seven integers apart. Uniqueness is guaranteed within
                # each recorded run family, not between the two historic jobs.
                family = 'Chen_Wallen' if network in ['Chen', 'Wallen'] else 'nine_networks'
                check(network + null + '_unique_seed_within_family', (family, seed) not in seed_keys)
                seed_keys.add((family, seed))
                all_seeds.add(seed)
                count = int(r['recovered_count'])
                check(network + null + '_valid_count_rate', 0 <= count <= E and close(float(r['recovery_rate']), count/E))
            study = SimpleNamespace(name=network, nodes=range(n), E=E,
                                    observed_recovered=observed, observed_rate=observed/E)
            value = ns['summarize'](study, 'Null ' + null, pd.DataFrame(rows), 1000)
            if network in ['Chen', 'Wallen']:
                comparisons.extend(compare_rows([value], [frozen], ['study', 'null_model'], network + null)[0])
            else:
                mapping = dict(null_mean_recovered_count='null_mean_recovered', null_sd_recovered_count='null_sd_recovered',
                               null_mean_recovery_rate='null_mean_recovery_rate', null_sd_recovery_rate='null_sd_recovery_rate',
                               **{'null_2.5pct_recovery_rate': 'null_q025_recovery_rate',
                                  'null_97.5pct_recovery_rate': 'null_q975_recovery_rate'},
                               observed_to_null_mean_enrichment='recovery_enrichment', z_score='z_score',
                               null_replicates_ge_observed='n_ge_observed', empirical_one_sided_P='empirical_p')
                check(network + null + '_frozen_summary', all(close(float(value[a]), float(frozen9[b])) for a, b in mapping.items()))
            obs, mean = value['observed_recovery_rate'], value['null_mean_recovery_rate']
            low, high = value['null_2.5pct_recovery_rate'], value['null_97.5pct_recovery_rate']
            value.update(network=network, null_type='Null ' + null,
                         excess_recovery_rate=obs-mean, excess_q025=obs-high, excess_q975=obs-low,
                         excess_recovery_percentage_points=100*(obs-mean),
                         excess_q025_percentage_points=100*(obs-high), excess_q975_percentage_points=100*(obs-low),
                         empirical_p=value['empirical_one_sided_P'])
            summary.append(value)
    plotted = read(register(lit / 'figures/Figure5_external_randomization/Figure5_excess_values.tsv', 'frozen_S1_22_display_records'))
    diffs, _ = compare_rows(summary, plotted, ['network', 'null_type'], 'S1_22_plot_records')
    comparisons.extend(diffs)
    write(OUTPUT / 'NULL_REPLAY_22_SUMMARIES.tsv', summary)
    write(OUTPUT / 'NULL_FLOAT_SERIALIZATION_DIFFERENCES.tsv', comparisons,
          ['table', 'key', 'column', 'recalculated', 'frozen', 'absolute_difference'])
    check('22000_iterations_with_family_unique_seeds', len(seed_keys) == 22000)
    (OUTPUT / 'SEED_SCOPE.json').write_text(json.dumps(dict(
        iteration_records=22000, family_seed_keys=len(seed_keys),
        globally_distinct_integer_seeds=len(all_seeds),
        note='Seed uniqueness is checked within each historical run family. Integer seeds overlap across the two different masters; no random simulations were changed or rerun.'), indent=2)+'\n')


def support(args):
    base = args.literature_root.resolve()
    before = base / 'stage4b_recurrence_recovery'
    old = base / 'stage4d_support_rank'
    ns = dict(Path=Path, collections=collections, csv=csv, json=json, math=math,
              np=np, stats=stats, multipletests=multipletests, read=read,
              write=write, check=check, log=print, OUT=OUTPUT,
              OUTCOMES=['source_excluded_support_count', 'support_rank_percentile', 'source_excluded_support_fraction'],
              GROUPS=['1', '2', '3', '>=4'],
              DECODING=register(base / 'stage3a/qc/03a_database_label_decoding.tsv', 'exact_GTDB_R226_label_inverse'),
              DB=register(args.rank_reference, 'entire_653071_pair_baseline_before_source_exclusion'),
              DB_HASH='4b67935940a765d28518a2dd612587f347c0e47e985323910fcdde51950a7834')
    load_functions(SUPPORT, ['fmt', 'yn', 'ids', 'group', 'database', 'pair_data', 'describe', 'summaries', 'statistics'], ns)
    check('baseline_database_SHA256', sha(ns['DB']) == ns['DB_HASH'])
    populations = {}
    for analysis, expected in [('primary', [5114, 749, 144, 55, 5]), ('sensitivity', [5487, 771, 152, 59, 5])]:
        path = before / f'pair_level/04b_{analysis}_pair_recovery.tsv'
        rows = [r for r in read(register(path, 'frozen_source_excluded_pair_recovery_' + analysis))
                if r['combined_source_excluded_recovered'] == 'Yes']
        check(analysis + '_recurrence_population', collections.Counter(int(r['literature_count']) for r in rows) == dict(enumerate(expected, 1)))
        check(analysis + '_unique_pair_keys', len(rows) == len({r['pair_key'] for r in rows}))
        populations[analysis] = sorted(rows, key=lambda r: r['pair_key'])
    tested = read(register(old / 'project_eligibility/04d_project_tested_species_sets.tsv', '859_frozen_actual_tested_species_sets'))
    check('859_unique_tested_projects', len(tested) == len({r['project_id'] for r in tested}) == 859)
    species_projects = collections.defaultdict(set)
    for row in tested:
        names = json.loads(row['eligible_species_set'])
        check('tested_set_unique_' + row['project_id'], len(names) == len(set(names)) == int(row['eligible_species_count']))
        for name in names:
            species_projects[name].add(row['project_id'])
    check('tested_union_12862', len(species_projects) == 12862)
    matches, ranks, _ = ns['database'](populations)
    differences = []
    for analysis, population in populations.items():
        rows, _ = ns['pair_data'](analysis, population, species_projects, matches, ranks)
        ns['summaries'](analysis, rows)
        ns['statistics'](analysis, rows)
        pattern = f'04d_{analysis}_*.tsv'
        candidates = [p for folder in ['data', 'summary', 'statistics'] for p in (OUTPUT/folder).glob(pattern)]
        for path in candidates:
            reference = old / path.relative_to(OUTPUT)
            if not reference.is_file():
                continue
            register(reference, 'frozen_' + path.name)
            actual, frozen = read(path), read(reference)
            if 'pair_key' in actual[0]:
                keys = ['pair_key']
            elif 'comparison' in actual[0]:
                keys = ['comparison'] + (['outcome'] if 'outcome' in actual[0] else [])
            elif 'outcome' in actual[0]:
                keys = ['outcome']
            else:
                keys = ['literature_count' if 'literature_count' in actual[0] else 'recurrence_group']
            differences.extend(compare_rows(actual, frozen, keys, path.name)[0])
    write(OUTPUT / 'SUPPORT_FLOAT_SERIALIZATION_DIFFERENCES.tsv', differences,
          ['table', 'key', 'column', 'recalculated', 'frozen', 'absolute_difference'])


def main():
    global OUTPUT
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='mode', required=True)
    for name in ['nulls', 'support']:
        command = sub.add_parser(name)
        command.add_argument('--literature-root', type=Path, required=True)
        command.add_argument('--output', type=Path, required=True)
        if name == 'nulls':
            command.add_argument('--chen-wallen-root', type=Path, required=True)
        else:
            command.add_argument('--rank-reference', type=Path, required=True)
    args = parser.parse_args()
    OUTPUT = args.output.resolve()
    if OUTPUT.exists():
        parser.error('--output must be a new directory; existing results are never overwritten')
    OUTPUT.mkdir(parents=True)
    status = 'FAIL'
    try:
        globals()[args.mode](args)
        status = 'PASS'
    finally:
        write(OUTPUT / 'CHECKS.tsv', QC, ['check', 'status', 'detail'])
        write(OUTPUT / 'INPUT_SHA256.tsv', [dict(path=str(p), role=role, size_bytes=p.stat().st_size, sha256=sha(p))
              for p, role in sorted(INPUTS.items())], ['path', 'role', 'size_bytes', 'sha256'])
        (OUTPUT / 'RESULT.json').write_text(json.dumps(dict(status=status, mode=args.mode,
              simulation_executed=False, original_results_modified=False,
              numpy=np.__version__, pandas=pd.__version__, python=sys.version), indent=2)+'\n')
    print(json.dumps(dict(status=status, mode=args.mode, output=str(OUTPUT))))


if __name__ == '__main__':
    main()
