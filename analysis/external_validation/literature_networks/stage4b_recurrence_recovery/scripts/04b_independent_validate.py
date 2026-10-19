#!/usr/bin/env python3
"""Independent read-back verification of Stage 4B files and summary tables.

Does not import the construction script or use its in-memory counts. Reloads
the original database, frozen recurrence tables, source membership, and new
pair-level output. Uses tuple filtering and Counter aggregation independently.
"""

# v0.1: repository-local IO configuration; scientific functions below are unchanged.
import sys as _path_sys
from pathlib import Path as _RepoPath
_repo_root = next(p for p in _RepoPath(__file__).resolve().parents if (p / "lib/mackdb_paths.py").is_file())
_path_sys.path.insert(0, str(_repo_root / "lib"))
from mackdb_paths import resolve_path
import collections
import csv
import json
import math
from pathlib import Path

BASE = Path(str(resolve_path('literature')))
OUT = BASE / 'stage4b_recurrence_recovery'
A = BASE / 'stage4a_recurrence'
DB = Path(str(resolve_path('environment_network')))
PAPERS = ['Chen', 'Wallen', 'Zhou', 'Guo', 'Lacruz', 'Li']
VECTORS = ['proj_ID', 'proj_corr', 'proj_pv', 'proj_sam']
CHECKS = []


def read(p):
    with Path(p).open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def check(name, ok, detail=''):
    CHECKS.append(dict(check=name, status='PASS' if ok else 'FAIL', details=str(detail)))
    if not ok:
        raise AssertionError(name + ': ' + str(detail))


def split(s):
    return s.split(';') if s else []


def eq_rate(observed, numerator, denominator):
    if denominator == 0:
        return observed == 'NA'
    try:
        return math.isclose(float(observed), numerator / denominator, rel_tol=0, abs_tol=1e-14)
    except ValueError:
        return False


def main():
    pair_tables = {a: read(OUT / ('pair_level/04b_' + a + '_pair_recovery.tsv')) for a in ['primary', 'sensitivity']}
    wanted = {r['pair_key'] for rows in pair_tables.values() for r in rows}
    decoder = {r['database_label']: r['gtdb_r226_species'] for r in read(BASE / 'stage3a/qc/03a_database_label_decoding.tsv')}
    original_db = {}
    with DB.open(newline='', encoding='utf-8') as f:
        for raw in csv.DictReader(f, delimiter='\t'):
            endpoints = sorted([decoder[raw['species_a']], decoder[raw['species_b']]])
            key = '\t'.join(endpoints)
            if key in wanted:
                if key in original_db:
                    raise AssertionError('Duplicate exact database pair')
                original_db[key] = raw
    config = {r['network']: set(split(r['reported_source_projects'])) for r in read(OUT / 'source_sets/04b_network_source_config.tsv')}
    universe = {r['project'] for r in read(OUT / 'qc/04b_actual_project_support_universe.tsv')}
    for analysis, rows in pair_tables.items():
        frozen = {r['pair_key']: r for r in read(A / ('membership/04a_' + analysis + '_unique_pair_recurrence.tsv'))}
        member = {(r['pair_key'], r['paper']): r['subnetwork_sources'].split('|') for r in read(A / ('membership/04a_' + analysis + '_pair_paper_membership.tsv'))}
        source_table = {r['pair_key']: r for r in read(OUT / ('source_sets/04b_' + analysis + '_pair_source_sets.tsv'))}
        check(analysis + '_same_unique_frozen_pair_universe', len(rows) == len(frozen) == len({r['pair_key'] for r in rows}) and {r['pair_key'] for r in rows} == set(frozen) == set(source_table))
        counts, baselines, recovered, losses = [collections.Counter() for _ in range(4)]
        combo_n, combo_b, combo_r, combo_l = [collections.Counter() for _ in range(4)]
        errors = collections.Counter()
        for r in rows:
            key = r['pair_key']
            old, sr = frozen[key], source_table[key]
            errors['all_frozen_master_fields'] += not all(r[k] == v for k, v in old.items())
            errors['canonical_TAB_key'] += key != r['canonical_species_a'] + '\t' + r['canonical_species_b'] or not r['canonical_species_a'] < r['canonical_species_b']
            per_paper = {}
            provenance = {}
            for p in PAPERS:
                nets = member.get((key, p), [])
                provenance[p] = nets
                expected = set()
                for n in nets:
                    expected.update(config[n])
                per_paper[p] = expected
                errors['paper_sources'] += set(split(sr[p + '_sources'])) != expected
                errors['paper_overlapping_sources'] += set(split(sr[p + '_overlapping_sources'])) != expected & universe
                errors['paper_subnetwork_sources'] += sr[p + '_subnetwork_sources'] != '|'.join(nets)
            reported = set.union(*per_paper.values())
            removal = reported & universe
            errors['MHO_policy'] += analysis == 'primary' and 'Lacruz_MHO' in provenance['Lacruz']
            errors['provenance_JSON'] += json.loads(r['subnetwork_sources_by_paper']) != provenance
            errors['union_source_sets'] += set(split(r['source_projects_to_remove'])) != removal or set(split(sr['overlapping_source_projects_union'])) != removal or set(split(sr['reported_source_projects_union'])) != reported
            errors['unique_source_count'] += len(split(r['source_projects_to_remove'])) != len(removal) or int(r['source_project_remove_count']) != len(removal)
            raw = original_db.get(key)
            baseline_flag = raw is not None
            original_tuples = list(zip(*(split(raw[c]) for c in VECTORS))) if raw else []
            kept = [t for t in original_tuples if t[0] not in removal]
            deleted = [t for t in original_tuples if t[0] in removal]
            for i, c in enumerate(VECTORS):
                errors['baseline_vectors_match_original_database'] += r['baseline_' + c] != (raw[c] if raw else '')
                errors['remaining_vectors_from_independent_tuple_filter'] += split(r['remaining_' + c]) != [t[i] for t in kept]
            errors['remaining_count'] += int(r['remaining_proj_count']) != len(kept)
            errors['baseline_count'] += int(r['baseline_proj_count']) != len(original_tuples)
            errors['actual_removed_support'] += int(r['removed_source_support_count']) != len(deleted) or split(r['removed_source_support_projects']) != [t[0] for t in deleted]
            errors['removed_indices'] += r['removed_support_indices_zero_based'] != ';'.join(str(i) for i, t in enumerate(original_tuples) if t[0] in removal)
            recovery_flag = len(kept) >= 2
            loss_flag = baseline_flag and not recovery_flag
            errors['baseline_outcome'] += r['baseline_mackdb_present'] != ('Yes' if baseline_flag else 'No')
            errors['combined_outcome'] += r['combined_source_excluded_recovered'] != ('Yes' if recovery_flag else 'No')
            errors['loss_outcome'] += r['loss_due_to_combined_source_exclusion'] != ('Yes' if loss_flag else 'No')
            errors['optional_mean_policy'] += r['remaining_aver_corr'] != 'NA' or r['remaining_aver_pv'] != 'NA'
            level, combo = int(r['literature_count']), r['supporting_papers']
            errors['paper_unit_count'] += level != sum(int(r[p]) for p in PAPERS) or combo.split('|') != [p for p in PAPERS if r[p] == '1']
            counts[level] += 1
            baselines[level] += baseline_flag
            recovered[level] += recovery_flag
            losses[level] += loss_flag
            combo_n[combo] += 1
            combo_b[combo] += baseline_flag
            combo_r[combo] += recovery_flag
            combo_l[combo] += loss_flag
        for name in sorted(errors):
            check(analysis + '_all_rows_' + name, errors[name] == 0, errors[name])
        frozen_distribution = read(A / ('summary/04a_' + analysis + '_literature_count_distribution.tsv'))
        check(analysis + '_literature_distribution_independent_readback', {int(r['literature_count']): int(r['unique_pair_count']) for r in frozen_distribution} == {k: counts[k] for k in range(1, 7)})
        exact = read(OUT / ('summary/04b_' + analysis + '_recovery_by_literature_count.tsv'))
        check(analysis + '_exact_bins_include_zero_six', [int(r['literature_count']) for r in exact] == list(range(1, 7)))
        for r in exact:
            k = int(r['literature_count'])
            check(analysis + '_exact_' + str(k) + '_integer_counts', [int(r[c]) for c in ['external_pair_count', 'baseline_mackdb_present_count', 'combined_source_excluded_recovered_count', 'source_exclusion_loss_count']] == [counts[k], baselines[k], recovered[k], losses[k]])
            check(analysis + '_exact_' + str(k) + '_rates', eq_rate(r['baseline_recovery_rate'], baselines[k], counts[k]) and eq_rate(r['combined_source_excluded_recovery_rate'], recovered[k], counts[k]) and eq_rate(r['source_exclusion_loss_fraction_of_baseline_recovered'], losses[k], baselines[k]))
        cum = read(OUT / ('summary/04b_' + analysis + '_cumulative_recurrence_recovery.tsv'))
        check(analysis + '_cumulative_thresholds_1_to_6', [int(r['minimum_literature_count']) for r in cum] == list(range(1, 7)))
        for r in cum:
            k = int(r['minimum_literature_count'])
            n, b, rec, loss = [sum(counter[x] for x in range(k, 7)) for counter in [counts, baselines, recovered, losses]]
            check(analysis + '_cumulative_' + str(k) + '_integer_counts', [int(r[c]) for c in ['external_pair_count', 'baseline_mackdb_present_count', 'recovered_count', 'source_exclusion_loss_count']] == [n, b, rec, loss])
            check(analysis + '_cumulative_' + str(k) + '_rates', eq_rate(r['recovery_rate'], rec, n) and eq_rate(r['baseline_recovery_rate'], b, n))
        combinations = read(OUT / ('summary/04b_' + analysis + '_support_combination_recovery.tsv'))
        check(analysis + '_exact_combination_universe', len(combinations) == len(combo_n) and {r['supporting_papers'] for r in combinations} == set(combo_n))
        for r in combinations:
            c = r['supporting_papers']
            check(analysis + '_combination_' + c, int(r['pair_count']) == combo_n[c] and int(r['baseline_mackdb_present_count']) == combo_b[c] and int(r['recovered_count']) == combo_r[c] and int(r['source_exclusion_loss_count']) == combo_l[c] and eq_rate(r['recovery_rate'], combo_r[c], combo_n[c]) and eq_rate(r['baseline_recovery_rate'], combo_b[c], combo_n[c]))
        lost_rows = read(OUT / ('pair_level/04b_' + analysis + '_pairs_lost_after_source_exclusion.tsv'))
        check(analysis + '_lost_pair_audit_complete', {r['pair_key'] for r in lost_rows} == {r['pair_key'] for r in rows if r['loss_due_to_combined_source_exclusion'] == 'Yes'} and len(lost_rows) == sum(losses.values()))
    print('Independent validation PASS: {} checks; original DB tuples, source provenance, frozen identity and all descriptive summaries.'.format(len(CHECKS)))


if __name__ == '__main__':
    try:
        main()
    finally:
        path = OUT / 'qc/04b_independent_summary_validation.tsv'
        with path.open('w', encoding='utf-8', newline='') as f:
            w = csv.DictWriter(f, fieldnames=['check', 'status', 'details'], delimiter='\t', lineterminator='\n')
            w.writeheader()
            w.writerows(CHECKS)
