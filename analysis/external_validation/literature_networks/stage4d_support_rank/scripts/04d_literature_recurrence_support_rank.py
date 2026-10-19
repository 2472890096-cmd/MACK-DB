#!/home/wangmengyang/anaconda3/bin/python
"""Stage 4D. Read frozen Stage 4B outcomes and Stage 3A.1 tested sets only.

Run: /home/wangmengyang/anaconda3/bin/python -s -B <this file>
Writes exclusively stage4d_support_rank/ and the named Stage 4D report.
"""

# v0.1: repository-local IO configuration; scientific functions below are unchanged.
import sys as _path_sys
from pathlib import Path as _RepoPath
_repo_root = next(p for p in _RepoPath(__file__).resolve().parents if (p / "lib/mackdb_paths.py").is_file())
_path_sys.path.insert(0, str(_repo_root / "lib"))
from mackdb_paths import resolve_path
import os
for _var in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[_var] = '1'
import collections
import csv
import hashlib
import json
import math
import platform
import subprocess
import sys
import traceback
from pathlib import Path
import numpy as np
import scipy
from scipy import stats
import statsmodels
from statsmodels.stats.multitest import multipletests

sys.dont_write_bytecode = True
csv.field_size_limit(10000000)
BASE = Path(str(resolve_path('literature')))
OUT = BASE / 'stage4d_support_rank'
REPORT = BASE / 'reports/04d_literature_recurrence_support_rank.md'
B = BASE / 'stage4b_recurrence_recovery'
C = BASE / 'stage4c_recurrence_statistics'
A1 = BASE / 'stage3a1_zhou_miss'
DB = Path(str(resolve_path('environment_network')))
DB_HASH = '4b67935940a765d28518a2dd612587f347c0e47e985323910fcdde51950a7834'
PROJECTS = Path(str(resolve_path('environment', 'human_gut.list')))
MISSING = DB.parent / 'human_gut_links_missing_projects.list'
AUDIT = A1 / 'qc/03a1_project_input_audit.tsv'
DECODING = BASE / 'stage3a/qc/03a_database_label_decoding.tsv'
SCRIPT = Path(__file__).resolve()
VALIDATOR = SCRIPT.with_name('04d_independent_recalculation.py')
EXPECTED = {'primary': [5114, 749, 144, 55, 5], 'sensitivity': [5487, 771, 152, 59, 5]}
GROUPS = ['1', '2', '3', '>=4']
OUTCOMES = ['source_excluded_support_count', 'support_rank_percentile', 'source_excluded_support_fraction']
QC, LOG, OUTPUTS, INPUT_ROLES = [], [], set(), {}


def read(path):
    with Path(path).open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write(path, rows, fields=None):
    path, rows = Path(path), list(rows)
    assert OUT in path.parents, 'Unauthorized write: ' + str(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    OUTPUTS.add(path)


def textfile(path, text):
    assert OUT in path.parents or path == REPORT
    path.write_text(text, encoding='utf-8')
    OUTPUTS.add(path)


def check(name, ok, detail=''):
    QC.append(dict(check=name, status='PASS' if ok else 'FAIL', details=str(detail)))
    if not ok:
        raise RuntimeError(name + ': ' + str(detail))


def log(text):
    LOG.append(text)
    print(text, flush=True)


def fmt(x):
    return format(float(x), '.17g')


def yn(x):
    return 'Yes' if x else 'No'


def ids(text):
    return set() if text in {'', 'NA'} else set(text.split(';'))


def group(k):
    return str(k) if int(k) < 4 else '>=4'


def snapshot(paths):
    return {p: (p.stat().st_size, sha(p)) for p in sorted(paths)}


def prepare():
    for d in ['data', 'project_eligibility', 'summary', 'statistics', 'qc', 'scripts', 'logs']:
        (OUT / d).mkdir(parents=True, exist_ok=True)
    prior = {}
    for p in [C / 'qc/04c_input_integrity.tsv', C / 'qc/04c_formal_output_manifest.tsv']:
        INPUT_ROLES[p] = 'frozen_integrity_evidence'
        for r in read(p):
            prior[Path(r.get('input_file') or r['absolute_path'])] = (r.get('sha256_after') or r['sha256'], str(p))
    for p, role in [(DB, 'complete_FastSpar_rank_reference'), (PROJECTS, 'formal_project_list'),
                    (MISSING, 'formal_unavailable_projects'), (AUDIT, 'frozen_actual_tested_species_sets'),
                    (DECODING, 'frozen_exact_database_label_decoding'),
                    (A1 / 'qc/03a1_confirmed_rules.tsv', 'frozen_eligibility_rule'),
                    (A1 / 'qc/03a1_input_manifest.tsv', 'frozen_species_file_hashes'),
                    (B / 'qc/04b_actual_project_support_universe.tsv', 'frozen_support_universe'),
                    (SCRIPT, 'implementation'), (VALIDATOR, 'independent_implementation')]:
        INPUT_ROLES[p] = role
    for a in EXPECTED:
        INPUT_ROLES[B / f'pair_level/04b_{a}_pair_recovery.tsv'] = 'formal_recovery_recurrence_support'
        INPUT_ROLES[B / f'summary/04b_{a}_recovery_by_literature_count.tsv'] = 'frozen_recovered_count_audit'
    species_prior = {Path(r['input_file']): r for r in read(A1 / 'qc/03a1_input_manifest.tsv')
                     if r['role'] == 'project_species' and r['exists'] == 'True'}
    for r in read(AUDIT):
        if r['input_status'] == 'complete':
            p = Path(r['species_file'])
            INPUT_ROLES[p] = 'actual_species_input_verified_against_frozen_tested_set'
            expected = species_prior[p]['sha256_after']
            check('species_file_prior_hash_consistent:' + r['project_id'], p not in prior or prior[p][0] == expected)
            prior[p] = (expected, str(A1 / 'qc/03a1_input_manifest.tsv'))
    protected = set(prior) | set(INPUT_ROLES)
    protected.update(p for p in BASE.rglob('*') if p.is_file() and OUT not in p.parents and p != REPORT)
    check('protected_inputs_exist', all(p.is_file() for p in protected))
    log('Hashing frozen inputs, prior outputs and 859 small species input files.')
    before = snapshot(protected)
    for p in sorted(prior):
        check('prior_frozen_hash:' + str(p), before[p][1] == prior[p][0])
    check('official_database_SHA256', before[DB][1] == DB_HASH)
    write(OUT / 'qc/04d_input_manifest.tsv', [dict(input_file=str(p), role=INPUT_ROLES[p],
          size_bytes=before[p][0], sha256=before[p][1]) for p in sorted(INPUT_ROLES)])
    versions = [dict(software='python', version=platform.python_version(), path=sys.executable)]
    for mod in [np, scipy, statsmodels]:
        versions.append(dict(software=mod.__name__, version=mod.__version__, path=mod.__file__))
    write(OUT / 'qc/04d_software_versions.tsv', versions)
    return protected, before, prior


def eligibility(before):
    listed = PROJECTS.read_text().splitlines()
    missing = set(MISSING.read_text().splitlines())
    audit = read(AUDIT)
    check('project_list_875_unique', len(listed) == len(set(listed)) == 875)
    check('Stage3A1_project_audit_complete', len(audit) == 875 and {r['project_id'] for r in audit} == set(listed))
    check('unavailable_16', len(missing) == 16 and missing <= set(listed))
    sets, out, set_rows = {}, [], []
    for r in sorted(audit, key=lambda r: r['project_id']):
        p = r['project_id']
        available = p not in missing
        check('eligibility_availability:' + p, available == (r['input_status'] == 'complete')
              == (r['included_in_formal_f2_project_universe'] == 'True'))
        if available:
            names = json.loads(r['eligible_species_set'])
            check('frozen_tested_set:' + p, len(names) == len(set(names)) == int(r['eligible_species_count'])
                  and r['matrix_labels_align'] == 'True')
            # Verify the already established serialization, without taxonomy lookup.
            # Eligibility itself is taken directly from the frozen Stage 3A.1 set.
            sp = [line.split('\t') for line in Path(r['species_file']).read_text().splitlines()]
            check('actual_species_rows:' + p, len(sp) == len(names) and all(len(x) == 2 for x in sp)
                  and [x[0] for x in sp] == [str(i) for i in range(len(sp))])
            check('actual_serialization_matches_frozen_set:' + p,
                  {x[1].split('_', 1)[1] for x in sp} == {s.replace(' ', '_') for s in names}
                  and len({s.replace(' ', '_') for s in names}) == len(names)
                  and all(x[1].split('_', 1)[0].isdigit() for x in sp))
            sets[p] = frozenset(names)
            set_rows.append(dict(project_id=p, eligible_species_count=len(names),
                eligible_species_set=json.dumps(sorted(names), ensure_ascii=False, separators=(',', ':')),
                source=str(AUDIT), source_sha256=before[AUDIT][1], species_file=r['species_file'],
                species_file_sha256=before[Path(r['species_file'])][1]))
        else:
            check('unavailable_NA_not_absent:' + p, r['eligible_species_set'] == 'NA' and r['matrix_labels_align'] == 'NA')
        out.append(dict(project_id=p, listed_in_human_gut='Yes', included_in_denominator_universe=yn(available),
            input_status=r['input_status'], eligible_species_count=len(sets[p]) if available else 'NA',
            tested_labels_read=yn(available) if available else 'NA',
            frozen_matrix_labels_align=r['matrix_labels_align'],
            actual_species_file_hash_matches_stage3a1='Yes' if available else 'NA',
            actual_species_rows_match_frozen_set='Yes' if available else 'NA',
            eligibility_source=str(AUDIT), species_file=r['species_file'],
            unavailable_treatment='exclude_entirely_not_absent_or_unsupported' if not available else 'actual_tested_set'))
    universe = set().union(*sets.values())
    check('all_859_tested_sets_read', len(sets) == 859 and set(sets) == set(listed) - missing)
    check('tested_species_union_12862_not_database_6210', len(universe) == 12862)
    actual = {r['project'] for r in read(B / 'qc/04b_actual_project_support_universe.tsv')}
    check('formal_859_equal_frozen_database_support_universe', actual == set(sets))
    write(OUT / 'project_eligibility/04d_project_tested_species_sets.tsv', set_rows)
    write(OUT / 'qc/04d_project_eligibility_audit.tsv', out)
    species_projects = collections.defaultdict(set)
    for p, names in sets.items():
        for s in names:
            species_projects[s].add(p)
    log('Reused all 859 frozen actual tested sets (12862 species); verified 859 species files; excluded 16 unavailable projects.')
    return species_projects


def population():
    populations = {}
    for a, expected in EXPECTED.items():
        rows = read(B / f'pair_level/04b_{a}_pair_recovery.tsv')
        check(a + '_input_unique_pairs', len(rows) == len({r['pair_key'] for r in rows}))
        check(a + '_formal_recovery_Yes_No', all(r['combined_source_excluded_recovered'] in {'Yes', 'No'} for r in rows))
        recovered = [r for r in rows if r['combined_source_excluded_recovered'] == 'Yes']
        counts = collections.Counter(int(r['literature_count']) for r in recovered)
        check(a + '_recovered_counts_match_requested', counts == dict(enumerate(expected, 1)), dict(counts))
        frozen = read(B / f'summary/04b_{a}_recovery_by_literature_count.tsv')
        check(a + '_recovered_counts_match_Stage4B', [int(r['combined_source_excluded_recovered_count']) for r in frozen] == expected + [0])
        check(a + '_canonical_keys', all(r['canonical_species_a'] < r['canonical_species_b'] and
              r['pair_key'] == r['canonical_species_a'] + '\t' + r['canonical_species_b'] for r in recovered))
        populations[a] = sorted(recovered, key=lambda r: r['pair_key'])
    log('Frozen recovered counts verified: primary 6067; sensitivity 6474. No unrecovered pair included.')
    return populations


def database(populations):
    decoder = {r['database_label']: r['gtdb_r226_species'] for r in read(DECODING)}
    targets = {r['pair_key'] for rows in populations.values() for r in rows}
    counts, keys, matches = [], set(), {}
    with DB.open(newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            a, b = sorted([decoder[row['species_a']], decoder[row['species_b']]])
            key = a + '\t' + b
            if a == b or key in keys:
                check('database_no_self_loop_or_duplicate', False, key)
            keys.add(key)
            n = int(row['proj_count'])
            counts.append(n)
            if key in targets:
                matches[key] = (n, row['proj_ID'])
    check('database_653071_unique_no_self_loops', len(counts) == len(keys) == 653071)
    check('all_recovered_exactly_matched_database', set(matches) == targets)
    v = np.array(counts, dtype=int)
    ranks = stats.rankdata(-v, method='average')
    by_support = {}
    for n, rank in zip(v, ranks):
        if int(n) in by_support and by_support[int(n)] != float(rank):
            check('equal_support_equal_rank', False)
        by_support[int(n)] = float(rank)
    ordered = sorted(by_support)
    check('support_rank_monotonic', all(by_support[a] > by_support[b] for a, b in zip(ordered, ordered[1:])))
    counts_by_support = collections.Counter(counts)
    rank_rows = [dict(baseline_support_count=n, database_pair_count=counts_by_support[n],
        fastspar_support_rank=fmt(by_support[n]), support_rank_percentile=fmt(100*(len(v)-by_support[n])/(len(v)-1)),
        **{f'top_{t}_percent': yn(100*(len(v)-by_support[n])/(len(v)-1) >= 100-t) for t in [1, 5, 10, 25]}) for n in ordered]
    write(OUT / 'data/04d_database_support_rank_distribution.tsv', rank_rows)
    summary = dict(N_DB=len(v), unique_canonical_pairs=len(keys), self_loops=0, duplicates=0,
        minimum_support_count=int(v.min()), maximum_support_count=int(v.max()), median_support_count=fmt(np.median(v)),
        Q1_support_count=fmt(np.quantile(v, .25)), Q3_support_count=fmt(np.quantile(v, .75)), mean_support_count=fmt(v.mean()),
        support_rank_monotonic='Yes', tie_method='average',
        percentile_formula='100*(N_DB-rank)/(N_DB-1)', rank_reference='entire_baseline_FastSpar_database',
        original_database_sha256=DB_HASH)
    for t in [1, 5, 10, 25]:
        n = sum(r['database_pair_count'] for r in rank_rows if r[f'top_{t}_percent'] == 'Yes')
        summary[f'top_{t}_database_count'] = n
        summary[f'top_{t}_database_fraction'] = fmt(n / len(v))
    write(OUT / 'qc/04d_database_rank_audit.tsv', [summary])
    log('Database reference: {} pairs, support min/max/median = {}/{}/{}.'.format(len(v), v.min(), v.max(), np.median(v)))
    return matches, by_support, summary


def pair_data(analysis, inputs, species_projects, matches, ranks):
    rows, evidence, fraction_qc = [], [], []
    for r in inputs:
        a, b, key = r['canonical_species_a'], r['canonical_species_b'], r['pair_key']
        baseline = int(r['baseline_proj_count'])
        remaining = int(r['remaining_proj_count'])  # Direct frozen value. No new support removal.
        base_ids, remaining_ids, source = ids(r['baseline_proj_ID']), ids(r['remaining_proj_ID']), ids(r['source_projects_to_remove'])
        exact = matches[key]
        joint = species_projects[a] & species_projects[b]
        removed_eligible = joint & source
        excluded_joint = joint - source  # Authorized denominator-only exclusion.
        n0, n1 = len(joint), len(excluded_joint)
        tests = dict(baseline_support_matches_DB=baseline == exact[0], baseline_support_IDs_match_DB=base_ids == ids(exact[1]),
            frozen_support_vectors_unique=baseline == len(base_ids) == len(r['baseline_proj_ID'].split(';'))
                and remaining == len(remaining_ids) == len(r['remaining_proj_ID'].split(';')),
            frozen_source_set_count=int(r['source_project_remove_count']) == len(source),
            frozen_recovered_baseline_present=r['baseline_mackdb_present'] == 'Yes' and remaining >= 2,
            all_baseline_support_projects_eligible=base_ids <= joint,
            all_remaining_support_projects_eligible=remaining_ids <= excluded_joint,
            baseline_support_le_joint=baseline <= n0, remaining_support_le_joint=remaining <= n1,
            source_excluded_denominator_ge2=n1 >= 2)
        if not all(tests.values()):
            fraction_qc.append(dict(analysis_set=analysis, pair_key=key, status='FAIL', **tests))
            write(OUT / 'qc/04d_support_fraction_qc.tsv', fraction_qc)
            check(analysis + '_pair_support_eligibility', False, (key, tests))
        corr = [float(x) for x in r['remaining_proj_corr'].split(';')]
        check(analysis + '_correlation_vector:' + key, len(corr) == remaining and all(math.isfinite(x) for x in corr))
        pos, neg, zero = sum(x > 0 for x in corr), sum(x < 0 for x in corr), sum(x == 0 for x in corr)
        rank = ranks[baseline]
        percentile = 100*(653071-rank)/653070
        row = {c: r[c] for c in ['canonical_species_a', 'canonical_species_b', 'pair_key', 'literature_count', 'supporting_papers']}
        row.update(analysis_set=analysis, recurrence_group=group(r['literature_count']), baseline_support_count=baseline,
            source_excluded_support_count=remaining, fastspar_support_rank=fmt(rank), support_rank_percentile=fmt(percentile),
            **{f'top_{t}_percent': yn(percentile >= 100-t) for t in [1, 5, 10, 25]},
            baseline_jointly_testable_projects=n0, source_excluded_jointly_testable_projects=n1,
            baseline_support_fraction=fmt(baseline/n0) if n0 else 'NA',
            source_excluded_support_fraction=fmt(remaining/n1) if n1 else 'NA',
            source_projects_to_remove=r['source_projects_to_remove'], source_project_remove_count=r['source_project_remove_count'],
            remaining_positive_project_count=pos, remaining_negative_project_count=neg, remaining_zero_project_count=zero,
            dominant_sign_fraction=fmt(max(pos, neg)/remaining), baseline_proj_ID=r['baseline_proj_ID'],
            remaining_proj_ID=r['remaining_proj_ID'], combined_source_excluded_recovered=r['combined_source_excluded_recovered'])
        rows.append(row)
        evidence.append(dict(analysis_set=analysis, pair_key=key, baseline_jointly_testable_projects=n0,
            source_excluded_jointly_testable_projects=n1, baseline_eligible_project_ID=';'.join(sorted(joint)),
            source_projects_to_remove=r['source_projects_to_remove'], removed_eligible_project_ID=';'.join(sorted(removed_eligible)),
            source_excluded_eligible_project_ID=';'.join(sorted(excluded_joint))))
        fraction_qc.append(dict(analysis_set=analysis, pair_key=key, status='PASS', **tests))
    write(OUT / f'data/04d_{analysis}_recovered_pair_support_rank.tsv', rows)
    write(OUT / f'project_eligibility/04d_{analysis}_pair_project_eligibility.tsv', evidence)
    check(analysis + '_all_support_fractions_in_range', all(0 <= float(r[c]) <= 1 for r in rows
          for c in ['baseline_support_fraction', 'source_excluded_support_fraction']))
    return rows, fraction_qc


def describe(values):
    v = np.array(values, dtype=float)
    return dict(median=fmt(np.median(v)), Q1=fmt(np.quantile(v, .25)), Q3=fmt(np.quantile(v, .75)),
                mean=fmt(v.mean()), min=fmt(v.min()), max=fmt(v.max()))


def summaries(analysis, rows):
    exact, grouped, paper = [], [], []
    for mode, levels, result in [('exact', list(range(1, 6)), exact), ('grouped', GROUPS, grouped)]:
        for k in levels:
            subset = [r for r in rows if (int(r['literature_count']) == k if mode == 'exact' else r['recurrence_group'] == k)]
            rec = dict(analysis_set=analysis, **{('literature_count' if mode == 'exact' else 'recurrence_group'): k}, N=len(subset))
            for outcome in ['baseline_support_count'] + OUTCOMES:
                for stat, value in describe([float(r[outcome]) for r in subset]).items():
                    rec[outcome + '_' + stat] = value
            for t in [1, 5, 10, 25]:
                n = sum(r[f'top_{t}_percent'] == 'Yes' for r in subset)
                rec[f'top_{t}_percent_count'], rec[f'top_{t}_percent_fraction'] = n, fmt(n/len(subset))
            result.append(rec)
            if mode == 'exact':
                paper.append(dict(analysis_set=analysis, literature_count=k, N_recovered_pairs=len(subset),
                    median_source_excluded_support=rec['source_excluded_support_count_median'],
                    IQR_source_excluded_support=rec['source_excluded_support_count_Q1'] + ';' + rec['source_excluded_support_count_Q3'],
                    median_rank_percentile=rec['support_rank_percentile_median'],
                    IQR_rank_percentile=rec['support_rank_percentile_Q1'] + ';' + rec['support_rank_percentile_Q3'],
                    median_support_fraction=rec['source_excluded_support_fraction_median'],
                    IQR_support_fraction=rec['source_excluded_support_fraction_Q1'] + ';' + rec['source_excluded_support_fraction_Q3'],
                    top10_fraction=rec['top_10_percent_fraction'], IQR_representation='Q1;Q3 (linear interpolation)'))
    for suffix, prefixes in [('support', ['baseline_support_count_', 'source_excluded_support_count_']),
                             ('rank', ['support_rank_percentile_', 'top_']),
                             ('support_fraction', ['source_excluded_support_fraction_'])]:
        fields = ['analysis_set', 'literature_count', 'N'] + [c for c in exact[0] if any(c.startswith(p) for p in prefixes)]
        write(OUT / f'summary/04d_{analysis}_{suffix}_by_literature_count.tsv', [{c: r[c] for c in fields} for r in exact])
    write(OUT / f'summary/04d_{analysis}_grouped_summary.tsv', grouped)
    fields = ['analysis_set', 'pair_key', 'literature_count', 'recurrence_group'] + OUTCOMES
    write(OUT / f'summary/04d_{analysis}_plot_data.tsv', [{c: r[c] for c in fields} for r in rows])
    fields = ['analysis_set', 'literature_count', 'N'] + [c for c in exact[0] if c.startswith('top_')]
    write(OUT / f'summary/04d_{analysis}_top_rank_fraction_plot_data.tsv', [{c: r[c] for c in fields} for r in exact])
    return dict(exact=exact, grouped=grouped, paper=paper)


def statistics(analysis, rows):
    spearman, kruskal, mann, fisher, long = [], [], [], [], []
    x = np.array([int(r['literature_count']) for r in rows])
    g = np.array([r['recurrence_group'] for r in rows])
    for outcome in OUTCOMES:
        v = np.array([float(r[outcome]) for r in rows])
        check(analysis + '_finite_outcome:' + outcome, np.isfinite(v).all())
        sp = stats.spearmanr(x, v, alternative='two-sided')
        spearman.append(dict(analysis_set=analysis, outcome=outcome, N=len(rows), rho=fmt(sp.statistic),
            two_sided_P=fmt(sp.pvalue), recurrence='exact_1_to_5', P_method='asymptotic_Student_t_N_minus_2_with_average_tie_ranks'))
        samples = [v[g == k] for k in GROUPS]
        kw = stats.kruskal(*samples)
        kruskal.append(dict(analysis_set=analysis, outcome=outcome, N=len(rows), H=fmt(kw.statistic), df=3,
            P=fmt(kw.pvalue), groups='1;2;3;>=4', tie_correction='Yes'))
        family = []
        for k, target in zip(GROUPS[1:], samples[1:]):
            reference = samples[0]
            mw = stats.mannwhitneyu(target, reference, alternative='two-sided', method='asymptotic', use_continuity=True)
            family.append(dict(analysis_set=analysis, outcome=outcome, comparison=k + '_vs_1',
                group_A=k, group_B='1', N_A=len(target), N_B=len(reference), N=len(target)+len(reference),
                median_group_A=fmt(np.median(target)), median_group_B=fmt(np.median(reference)),
                U=fmt(mw.statistic), raw_P=fmt(mw.pvalue), BH_adjusted_P='',
                BH_family=analysis + '_' + outcome + '_three_comparisons', BH_family_size=3,
                U_A_over_nA_nB=fmt(mw.statistic/(len(target)*len(reference))),
                direction='higher' if mw.statistic > len(target)*len(reference)/2 else 'lower' if mw.statistic < len(target)*len(reference)/2 else 'equal',
                method='asymptotic_tie_corrected', continuity_correction='0.5', alternative='two-sided'))
        for r, adjusted in zip(family, multipletests([float(r['raw_P']) for r in family], method='fdr_bh')[1]):
            r['BH_adjusted_P'] = fmt(adjusted)
        mann.extend(family)
    top = np.array([r['top_10_percent'] == 'Yes' for r in rows])
    for k in GROUPS[1:]:
        a, b = int(top[g == k].sum()), int((~top[g == k]).sum())
        c, d = int(top[g == '1'].sum()), int((~top[g == '1']).sum())
        table = [[a, b], [c, d]]
        ft = stats.fisher_exact(table, alternative='two-sided')
        conditional = stats.contingency.odds_ratio(table, kind='conditional')
        ci = conditional.confidence_interval(confidence_level=.95, alternative='two-sided')
        fisher.append(dict(analysis_set=analysis, comparison=k + '_vs_1', N=a+b+c+d,
            target_top10=a, target_not_top10=b, reference_top10=c, reference_not_top10=d,
            target_top10_fraction=fmt(a/(a+b)), reference_top10_fraction=fmt(c/(c+d)),
            odds_ratio=fmt(conditional.statistic), OR_95CI_lower=fmt(ci.low), OR_95CI_upper=fmt(ci.high),
            OR_estimator='conditional_maximum_likelihood', CI_method='exact_conditional_equal_tailed_95_percent',
            sample_cross_product_OR=fmt(ft.statistic), raw_P=fmt(ft.pvalue), BH_adjusted_P='',
            BH_family=analysis + '_three_Top10_Fisher_comparisons', BH_family_size=3, alternative='two-sided'))
    for r, adjusted in zip(fisher, multipletests([float(r['raw_P']) for r in fisher], method='fdr_bh')[1]):
        r['BH_adjusted_P'] = fmt(adjusted)
    for records, test, effect, pcol in [(spearman, 'Spearman', 'rho', 'two_sided_P'),
        (kruskal, 'Kruskal-Wallis', 'H', 'P'), (mann, 'Mann-Whitney U', 'U', 'raw_P'),
        (fisher, 'Top10 Fisher exact', 'odds_ratio', 'raw_P')]:
        for r in records:
            check(analysis + '_valid_statistics:' + test + ':' + r.get('outcome', '') + ':' + r.get('comparison', ''),
                  0 <= float(r[pcol]) <= 1 and not math.isnan(float(r[effect])))
            long.append(dict(analysis_set=analysis, outcome=r.get('outcome', 'top_10_percent'), test=test,
                comparison=r.get('comparison', 'all_recovered_pairs'), effect_statistic=r[effect],
                effect_name=effect, P=r[pcol], adjusted_P_if_applicable=r.get('BH_adjusted_P', 'NA'), N=r['N'],
                OR_95CI_lower=r.get('OR_95CI_lower', 'NA'), OR_95CI_upper=r.get('OR_95CI_upper', 'NA')))
    for suffix, records in [('spearman', spearman), ('kruskal_wallis', kruskal), ('mann_whitney', mann), ('top10_fisher', fisher)]:
        write(OUT / f'statistics/04d_{analysis}_{suffix}.tsv', records)
    return dict(spearman=spearman, kruskal=kruskal, mann=mann, fisher=fisher, long=long)


def comparison(results):
    output = []
    for i, outcome in enumerate(OUTCOMES):
        p, s = results['primary'], results['sensitivity']
        pr, sr = float(p['spearman'][i]['rho']), float(s['spearman'][i]['rho'])
        pm = [float(r[outcome + '_median']) for r in p['grouped']]
        sm = [float(r[outcome + '_median']) for r in s['grouped']]
        pp, ss = [r for r in p['mann'] if r['outcome'] == outcome], [r for r in s['mann'] if r['outcome'] == outcome]
        same_trend = np.sign(pr) == np.sign(sr)
        same_median_order = np.array_equal(np.sign(np.diff(pm)), np.sign(np.diff(sm)))
        same_pairs = all(a['direction'] == b['direction'] for a, b in zip(pp, ss))
        output.append(dict(outcome=outcome, primary_rho=fmt(pr), sensitivity_rho=fmt(sr),
            primary_grouped_medians=';'.join(map(fmt, pm)), sensitivity_grouped_medians=';'.join(map(fmt, sm)),
            primary_exact_medians=';'.join(r[outcome + '_median'] for r in p['exact']),
            sensitivity_exact_medians=';'.join(r[outcome + '_median'] for r in s['exact']),
            Spearman_direction_consistent=yn(same_trend), grouped_median_order_consistent=yn(same_median_order),
            pairwise_direction_consistent=yn(same_pairs), direction_consistent=yn(same_trend and same_median_order and same_pairs),
            primary_all_three_pairwise_higher_BH_lt_0_05=yn(all(r['direction'] == 'higher' and float(r['BH_adjusted_P']) < .05 for r in pp)),
            sensitivity_all_three_pairwise_higher_BH_lt_0_05=yn(all(r['direction'] == 'higher' and float(r['BH_adjusted_P']) < .05 for r in ss))))
    write(OUT / 'summary/04d_primary_sensitivity_comparison.tsv', output)
    return output


METHODS = '''# Stage 4D methods

Population: unique Stage 4B pairs with combined_source_excluded_recovered == Yes.
literature_count, supporting_papers, source sets and remaining support are copied,
never recomputed. Unrecovered pairs have no assigned rank or support outcome here.
The primary set excludes exploratory Lacruz_MHO; sensitivity includes it.

Rank: scipy.stats.rankdata(-all_653071_proj_counts, method="average").
Percentile: 100*(653071-rank)/653070. Tied extremes need not attain 0/100;
threshold groups need not contain exactly 1/5/10/25% of the full database.
Top categories use percentiles >=99/95/90/75, without external-subset reranking.
Baseline rank is monotone in baseline count; source-excluded count uses the
frozen remaining_proj_count and need not have exactly the same ranks.

Eligibility: reuse the full exact species sets in the frozen Stage 3A.1 audit,
whose rows and FastSpar axes were aligned in that stage. Additionally verify all
859 original small species files against their frozen hashes and serialized
species sets. Do not rerun matrices, filtering, FastSpar, taxonomy or thresholds.
Unavailable 16 projects remain NA and are excluded entirely. For each pair,
joint=set(projects_testing_a) intersect set(projects_testing_b). The denominator
after source removal is len(joint - exact_frozen_Stage4B_source_projects_to_remove).
Only the denominator is newly computed; neither support removal nor recovery
is repeated. Each fraction is its support count / respective denominator.
All baseline/remaining support projects must belong to the respective eligible set.

Summaries: median and linear-interpolated Q1/Q3 (NumPy method="linear", default).
IQR columns use Q1;Q3, not Q3-Q1. Exact count 1..5 for Spearman; 1,2,3,>=4
for Kruskal-Wallis, two-sided Mann-Whitney and Top10 Fisher comparisons.
Spearman uses average ranks and the standard two-sided Student-t approximation
with df=N-2. Kruskal-Wallis uses tie correction and chi-square df=3 upper tail.
Mann-Whitney reports U for target A (2/3/>=4) versus reference B (1); asymptotic
normal approximation with tie correction and continuity correction 0.5.
BH adjusts only the three requested contrasts per outcome per analysis set;
the three Top10 Fisher contrasts form a separate BH family per analysis set.
No BH adjustment is assigned to the separately specified omnibus/trend tests.

Top10 table: [[target_Yes,target_No],[reference_Yes,reference_No]]. Main OR is
conditional maximum likelihood, and 95% CI inverts equal-tailed noncentral
hypergeometric tests; no pseudo-count is added. The sample cross-product OR is
also recorded separately. Two-sided Fisher P sums fixed-margin table probabilities
no greater than that of the observed table. Exact conditional OR/CI remain
well-defined at zero cells, allowing boundary 0/infinity as appropriate.

Independent validation: separate script imports no analysis functions. It reads
Stage 4B, original database and original frozen tested-set audit anew; computes
database ranks from frequency cumulative counts, denominators with per-project
bit masks, summaries with hand-coded quantiles, and tests from rank/variance
formulas. Fisher uses enumerated log combinatorial weights; exact OR and CI use
independent exponential tilting and root finding. BH uses sorted reverse minima.

These are unadjusted association tests with unique pairs as observational units.
Pairs can share species and studies; the requested tests do not correct network
dependence. Support fraction normalizes tested opportunity, not confounding.
No causal effect or direct ecological interaction is inferred. Exact count 5
has only five recovered pairs and supports descriptive reporting only.

Primary documentation:
- https://docs.scipy.org/doc/scipy-1.15.3/reference/generated/scipy.stats.mannwhitneyu.html
- https://docs.scipy.org/doc/scipy-1.15.3/reference/generated/scipy.stats.contingency.odds_ratio.html
'''


def table(fields, rows):
    def value(x):
        if isinstance(x, (float, np.floating)):
            return format(x, '.6g')
        return str(x).replace('|', '\\|')
    return '\n'.join(['| ' + ' | '.join(fields) + ' |', '| ' + ' | '.join(['---']*len(fields)) + ' |'] +
                     ['| ' + ' | '.join(value(r[c]) for c in fields) + ' |' for r in rows])


def report(results, db, compare, status):
    lines = ['# Stage 4D — Literature recurrence, BioProject support and database rank', '',
        'overall_status = ' + status['overall_status'], '',
        'Stage 4B 回答 pair 是否进入 MACK-DB；Stage 4D 仅在已经 combined-source-excluded recovered 的 pair 中，'
        '检验跨文献重复性与支持强度、全库排名及共同可检测机会归一化支持比例的关联。未恢复 pair 不赋零支持或末位排名。', '',
        'Primary 排除 exploratory Lacruz_MHO；Sensitivity 包含 MHO。正式 Stage 4B literature_count、supporting_papers、'
        'source set、remaining_proj_count、recovery 状态全部直接复用，未重新计算。', '',
        '## 定义与项目分母', '',
        '仅 FastSpar。baseline_support_count 由原始 human_gut.info 的 proj_count exact canonical lookup 取得，'
        '逐条等于 Stage 4B baseline_proj_count。主要支持指标 source_excluded_support_count 直接来自 Stage 4B remaining_proj_count。', '',
        '完整数据库按 proj_count 降序，ties 使用 average rank；percentile = 100 × (653071 − rank) / 653070。'
        'percentile 越高越靠前，Top 1/5/10/25% 分别用 ≥99/95/90/75。并列可能使实际入选比例偏离名义比例。'
        '排名是 baseline 全库排名，支持数主要指标是 source-excluded count，两者不保证 Spearman rho 完全相同。', '',
        '沿用 Stage 3A.1 已确认实际 FastSpar tested sets：875 listed，859 analyzable，16 unavailable。'
        '成功读取全部 859 个冻结完整 tested sets（并集 12862 species），并核验 859 个原 species 文件的冻结 SHA256、行数和序列化物种集。'
        '未用数据库 6210 species 代替 tested universe，未重新设置 prevalence 阈值；FastSpar 两轴对齐沿用冻结审计。', '',
        'baseline_jointly_testable_projects 为 859 个项目中同时测试两物种的项目数；source-excluded denominator '
        '再排除该 pair 在 Stage 4B 的 source_projects_to_remove。此步只计算 denominator，不重新删除支持项目或变更 recovery。'
        '16 unavailable 均排除，未视为 absent 或 unsupported。baseline/source-excluded support fraction 分别是对应支持数除以对应分母。', '',
        '## 当前数据库参考', '', table(list(db), [db]), '',
        '实际最大 support count 来自当前数据库，未使用旧网页 slider 的 557。', '',
        '## Primary / Sensitivity 描述与统计', '']
    for a, result in results.items():
        lines += ['### ' + a, '', 'Recovered N = ' + str(sum(EXPECTED[a])) + '。IQR 表示 Q1–Q3；分位数采用线性插值。', '',
            table(['literature_count', 'N_recovered_pairs', 'median_source_excluded_support', 'IQR_source_excluded_support',
                   'median_rank_percentile', 'median_support_fraction', 'top10_fraction'], result['paper']), '',
            'Spearman 使用 exact literature_count 1–5；P 为双侧。', '',
            table(['outcome', 'rho', 'two_sided_P', 'N'], result['spearman']), '',
            'Kruskal-Wallis：分组 1、2、3、>=4，ties 校正，df=3。', '',
            table(['outcome', 'H', 'df', 'P'], result['kruskal']), '',
            'Mann-Whitney：双侧、ties 校正、0.5 continuity correction；每种 outcome 的三个比较分别 BH。'
            'U 对应 target group A；direction 依据 U/(nA*nB) 与 0.5 的关系，而非把分布检验解释成纯中位数检验。', '',
            table(['outcome', 'comparison', 'median_group_A', 'median_group_B', 'U', 'raw_P', 'BH_adjusted_P', 'direction'], result['mann']), '',
            'Grouped Top10 比例：', '',
            table(['recurrence_group', 'N', 'top_10_percent_count', 'top_10_percent_fraction'], result['grouped']), '',
            'Top10 Fisher 双侧检验；OR 为 conditional MLE，CI 为 exact conditional 95% 区间，三个 P 单独 BH。', '',
            table(['comparison', 'odds_ratio', 'OR_95CI_lower', 'OR_95CI_upper', 'raw_P', 'BH_adjusted_P'], result['fisher']), '']
    lines += ['## 方向一致性与解释边界', '', table(list(compare[0]), compare), '',
        'direction_consistent 要求 Spearman 符号、grouped median 相邻顺序和 pairwise U 方向均一致；完整 exact medians 同时保留。'
        '这是 recovered pairs 内的关联分析，不进行因果解释；support fraction 是检测机会归一化，并未消除所有混杂。'
        'pair 可能共享物种或研究，当前指定检验不额外校正这种依赖。literature_count=5 仅五条，不能单独作强统计结论。'
        'sign consistency 仅作为描述字段保存，未检验。', '',
        '## QC 与独立复算', '',
        f"所有 {status['independent_checks']} 项独立复算通过，{status['protected_file_count']} 个受保护文件运行前后 SHA256 不变。"
        '原数据库正式 SHA256：`' + DB_HASH + '`。支持数超过 jointly testable count 的异常为 0；所有 recovered 分母 ≥2。'
        '全部 source-independent 支持项目均位于对应 source-excluded eligible set。独立复算覆盖输入人数/分组、支持数、全库 rank、'
        'percentile、逐项目分母、fraction、分位数、Spearman/P、Kruskal-Wallis/P、Mann-Whitney/P、BH 及 Fisher/OR/CI。', '',
        '统计方法及官方文档链接：`stage4d_support_rank/qc/04d_statistical_methods.md`。'
        '可重复命令：`/home/wangmengyang/anaconda3/bin/python -s -B stage4d_support_rank/scripts/04d_literature_recurrence_support_rank.py`。', '',
        '## 正式输出', '',
        '完整绝对路径、字节数和 SHA256 见 `stage4d_support_rank/qc/04d_formal_output_manifest.tsv`（清单自身不递归列入）。', '']
    textfile(REPORT, '\n'.join(lines))


def main():
    protected, before, prior = prepare()
    populations = population()
    species_projects = eligibility(before)
    matches, ranks, db = database(populations)
    results, all_fraction_qc = {}, []
    for a, inputs in populations.items():
        rows, fraction_qc = pair_data(a, inputs, species_projects, matches, ranks)
        all_fraction_qc.extend(fraction_qc)
        results[a] = {**summaries(a, rows), **statistics(a, rows)}
        log(a + ': support, rank, opportunity denominator, summaries and requested tests completed.')
    write(OUT / 'qc/04d_support_fraction_qc.tsv', all_fraction_qc)
    write(OUT / 'summary/04d_support_rank_paper_table.tsv', [r for v in results.values() for r in v['paper']])
    write(OUT / 'summary/04d_support_rank_statistics.tsv', [r for v in results.values() for r in v['long']])
    compare = comparison(results)
    textfile(OUT / 'qc/04d_statistical_methods.md', METHODS)
    log('Running independent formula/data-path validation.')
    completed = subprocess.run([sys.executable, '-s', '-B', str(VALIDATOR)], text=True, capture_output=True)
    if completed.stdout.strip():
        log(completed.stdout.strip())
    check('independent_validation_process', completed.returncode == 0, completed.stderr.strip())
    independent_path = OUT / 'qc/04d_independent_recalculation.tsv'
    OUTPUTS.add(independent_path)
    independent = read(independent_path)
    check('all_independent_recalculations_PASS', bool(independent) and all(r['status'] == 'PASS' for r in independent))
    log('Rechecking all frozen hashes after analysis.')
    after = snapshot(protected)
    integrity = [dict(input_file=str(p), scope=INPUT_ROLES.get(p, 'frozen_upstream_result_or_prior_input'),
        size_bytes=before[p][0], sha256_before=before[p][1], sha256_after=after[p][1],
        prior_recorded_sha256=prior[p][0] if p in prior else 'NA',
        prior_hash_evidence=prior[p][1] if p in prior else 'NA',
        matches_prior_recorded_hash=yn(before[p][1] == prior[p][0]) if p in prior else 'NA',
        unchanged=yn(before[p] == after[p])) for p in sorted(protected)]
    write(OUT / 'qc/04d_input_integrity.tsv', integrity)
    check('all_frozen_inputs_unchanged', before == after)
    check('original_database_final_hash', after[DB][1] == DB_HASH)
    status = dict(overall_status='PASS', primary_recovered=6067, sensitivity_recovered=6474,
        listed_projects=875, tested_labels_read_projects=859, unavailable_projects_excluded=16,
        support_exceeding_eligible_anomalies=0, independent_checks=len(independent), protected_file_count=len(protected),
        all_frozen_hashes_unchanged=True, direction_consistent=all(r['direction_consistent'] == 'Yes' for r in compare),
        database=db)
    textfile(OUT / 'qc/04d_run_status.json', json.dumps(status, ensure_ascii=False, indent=2) + '\n')
    report(results, db, compare, status)
    write(OUT / 'qc/04d_statistical_qc.tsv', QC)
    log('overall_status=PASS. Stage 4D complete; no downstream analysis executed.')
    textfile(OUT / 'logs/04d_literature_recurrence_support_rank.log', '\n'.join(LOG) + '\n')
    OUTPUTS.update([SCRIPT, VALIDATOR])
    manifest = OUT / 'qc/04d_formal_output_manifest.tsv'
    write(manifest, [dict(absolute_path=str(p), size_bytes=p.stat().st_size, sha256=sha(p),
          role='report' if p == REPORT else p.parent.name) for p in sorted(OUTPUTS)])
    actual = {p for p in OUT.rglob('*') if p.is_file()} | {REPORT}
    assert actual == OUTPUTS, 'Unmanifested or missing Stage 4D output'


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        OUT.joinpath('qc').mkdir(parents=True, exist_ok=True)
        write(OUT / 'qc/04d_statistical_qc.tsv', QC, ['check', 'status', 'details'])
        textfile(OUT / 'qc/04d_run_status.json', json.dumps(dict(overall_status='FAIL', error=str(exc)), ensure_ascii=False, indent=2) + '\n')
        textfile(OUT / 'logs/04d_literature_recurrence_support_rank.log', '\n'.join(LOG) + '\n' + traceback.format_exc())
        textfile(REPORT, '# Stage 4D\n\noverall_status = FAIL\n\n' + str(exc) + '\n')
        raise
