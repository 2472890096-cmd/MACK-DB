#!/home/wangmengyang/anaconda3/bin/python
"""Independent Stage 4D audit. Does not import the production implementation."""

# v0.1: repository-local IO configuration; scientific functions below are unchanged.
import sys as _path_sys
from pathlib import Path as _RepoPath
_repo_root = next(p for p in _RepoPath(__file__).resolve().parents if (p / "lib/mackdb_paths.py").is_file())
_path_sys.path.insert(0, str(_repo_root / "lib"))
from mackdb_paths import resolve_path
import os
for _var in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[_var] = '1'
import csv
import collections
import json
import math
import sys
from pathlib import Path
import numpy as np
from scipy import special, optimize

sys.dont_write_bytecode = True
csv.field_size_limit(10000000)
BASE = Path(str(resolve_path('literature')))
OUT = BASE / 'stage4d_support_rank'
B = BASE / 'stage4b_recurrence_recovery'
DB = Path(str(resolve_path('environment_network')))
EXPECTED = {'primary': [5114, 749, 144, 55, 5], 'sensitivity': [5487, 771, 152, 59, 5]}
GROUPS = ['1', '2', '3', '>=4']
OUTCOMES = ['source_excluded_support_count', 'support_rank_percentile', 'source_excluded_support_fraction']
CHECKS = []


def read(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def verify(name, ok, detail=''):
    CHECKS.append(dict(check=name, status='PASS' if ok else 'FAIL', detail=str(detail)))


def close(name, expected, observed, rel=2e-7, absolute=1e-10):
    a, b = float(expected), float(observed)
    ok = (a == b or (math.isfinite(a) and math.isfinite(b) and math.isclose(a, b, rel_tol=rel, abs_tol=absolute)))
    verify(name, ok, 'independent={:.17g}; formal={:.17g}'.format(a, b))


def pclose(name, expected, observed):
    close(name, expected, observed, rel=2e-6, absolute=1e-300)


def rank(values):
    counts = collections.Counter(values)
    before, ranks = 0, {}
    for value in sorted(counts):
        n = counts[value]
        ranks[value] = before + (n+1)/2
        before += n
    return np.array([ranks[x] for x in values], dtype=float)


def quantile(values, probability):
    v = sorted(values)
    index = (len(v)-1)*probability
    lo, hi = math.floor(index), math.ceil(index)
    return v[lo] + (v[hi]-v[lo])*(index-lo)


def describe(v):
    return dict(median=quantile(v, .5), Q1=quantile(v, .25), Q3=quantile(v, .75),
                mean=math.fsum(v)/len(v), min=min(v), max=max(v))


def bh(values):
    n = len(values)
    order = sorted(range(n), key=lambda i: values[i])
    out, running = [0.]*n, 1.
    for i in reversed(range(n)):
        j = order[i]
        running = min(running, n*values[j]/(i+1))
        out[j] = running
    return out


def manual_spearman(x, y):
    u, v = rank(x), rank(y)
    u -= u.mean()
    v -= v.mean()
    rho = float(np.dot(u, v)/math.sqrt(np.dot(u, u)*np.dot(v, v)))
    return rho, float(special.betainc((len(x)-2)/2, .5, max(0., 1-rho*rho)))


def manual_kw(samples):
    all_values = [x for sample in samples for x in sample]
    ranks = rank(all_values)
    n, index, total = len(all_values), 0, 0.
    for sample in samples:
        total += math.fsum(ranks[index:index+len(sample)])**2 / len(sample)
        index += len(sample)
    tie = sum(t**3-t for t in collections.Counter(all_values).values())
    correction = 1-tie/(n**3-n)
    h = (12*total/(n*(n+1))-3*(n+1))/correction
    return h, float(special.gammaincc(1.5, h/2))


def manual_mw(a, b):
    values = a + b
    n1, n2, n = len(a), len(b), len(values)
    u = math.fsum(rank(values)[:n1])-n1*(n1+1)/2
    tie = sum(t**3-t for t in collections.Counter(values).values())
    variance = n1*n2/12*((n+1)-tie/(n*(n-1)))
    z = (abs(u-n1*n2/2)-.5)/math.sqrt(variance)
    return u, min(1., math.erfc(z/math.sqrt(2)))


def exact_table(a, b, c, d):
    """Enumerate fixed margins, with independent exponential tilting for OR/CI."""
    total, row, col = a+b+c+d, a+b, a+c
    lo, hi = max(0, row+col-total), min(row, col)
    x = np.arange(lo, hi+1, dtype=float)
    def logchoose(n, k):
        return special.gammaln(n+1)-special.gammaln(k+1)-special.gammaln(n-k+1)
    logw = logchoose(col, x)+logchoose(total-col, row-x)
    logprob = logw-special.logsumexp(logw)
    p = min(1., float(np.exp(special.logsumexp(logprob[logprob <= logprob[a-lo]+1e-10]))))
    def probability(eta):
        weights = logw + eta*x
        return np.exp(weights-special.logsumexp(weights))
    def root(fn):
        return math.exp(optimize.brentq(fn, -700., 700., xtol=1e-12))
    estimate = 0. if a == lo else math.inf if a == hi else root(lambda e: float(np.dot(x, probability(e)))-a)
    lower = 0. if a == lo else root(lambda e: float(probability(e)[x >= a].sum())-.025)
    upper = math.inf if a == hi else root(lambda e: float(probability(e)[x <= a].sum())-.025)
    return p, estimate, lower, upper


def audit_inputs():
    inputs = {a: {r['pair_key']: r for r in read(B / f'pair_level/04b_{a}_pair_recovery.tsv')
                   if r['combined_source_excluded_recovered'] == 'Yes'} for a in EXPECTED}
    targets = set().union(*(set(rows) for rows in inputs.values()))
    decoder = {r['database_label']: r['gtdb_r226_species'] for r in read(BASE / 'stage3a/qc/03a_database_label_decoding.tsv')}
    frequencies, match, unique, selfloops = collections.Counter(), {}, set(), 0
    with DB.open(newline='') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            a, b = sorted([decoder[r['species_a']], decoder[r['species_b']]])
            pair = a+'\t'+b
            unique.add(pair)
            selfloops += a == b
            n = int(r['proj_count'])
            frequencies[n] += 1
            if pair in targets:
                match[pair] = n, r['proj_ID']
    total = sum(frequencies.values())
    verify('independent_database_count_no_duplicates_no_selfloops', total == len(unique) == 653071 and selfloops == 0)
    verify('all_recovered_pairs_exact_database_match', set(match) == targets)
    ranks, percentiles, below = {}, {}, 0
    for n in sorted(frequencies):
        count = frequencies[n]
        ranks[n] = total-below-(count-1)/2
        percentiles[n] = 100*(below+(count-1)/2)/(total-1)
        below += count
    for r in read(OUT / 'data/04d_database_support_rank_distribution.tsv'):
        n = int(r['baseline_support_count'])
        close('database_rank_for_support_' + str(n), ranks[n], r['fastspar_support_rank'])
        close('database_percentile_for_support_' + str(n), percentiles[n], r['support_rank_percentile'])
        verify('database_frequency_' + str(n), frequencies[n] == int(r['database_pair_count']))
    all_support = [n for n, count in frequencies.items() for _ in range(count)]
    dbsummary = read(OUT / 'qc/04d_database_rank_audit.tsv')[0]
    for c, v in describe(all_support).items():
        key = {'min': 'minimum', 'max': 'maximum'}.get(c, c)+'_support_count'
        close('database_summary_' + c, v, dbsummary[key])
    audit = read(BASE / 'stage3a1_zhou_miss/qc/03a1_project_input_audit.tsv')
    available = sorted([r for r in audit if r['included_in_formal_f2_project_universe'] == 'True'], key=lambda r: r['project_id'])
    verify('independent_875_859_16', len(audit) == 875 and len(available) == 859 and sum(r['eligible_species_set'] == 'NA' for r in audit) == 16)
    masks, project_bits = collections.defaultdict(int), {}
    for i, r in enumerate(available):
        bit = 1 << i
        project_bits[r['project_id']] = bit
        for s in json.loads(r['eligible_species_set']):
            masks[s] |= bit
    for a, frozen in inputs.items():
        output = read(OUT / f'data/04d_{a}_recovered_pair_support_rank.tsv')
        counts = collections.Counter(int(r['literature_count']) for r in frozen.values())
        verify(a+'_independent_recovered_groups', counts == dict(enumerate(EXPECTED[a], 1)))
        verify(a+'_output_unique_complete_recovered_only', len(output) == len(frozen) == sum(EXPECTED[a]) and {r['pair_key'] for r in output} == set(frozen))
        eligibility_rows = {r['pair_key']: r for r in read(OUT / f'project_eligibility/04d_{a}_pair_project_eligibility.tsv')}
        mismatches = collections.Counter()
        for r in output:
            key, old = r['pair_key'], frozen[r['pair_key']]
            for c in ['canonical_species_a', 'canonical_species_b', 'pair_key', 'literature_count', 'supporting_papers',
                      'source_projects_to_remove', 'source_project_remove_count', 'baseline_proj_ID', 'remaining_proj_ID', 'combined_source_excluded_recovered']:
                mismatches[c] += r[c] != old[c]
            base, left = int(old['baseline_proj_count']), int(old['remaining_proj_count'])
            source = [] if old['source_projects_to_remove'] in {'', 'NA'} else old['source_projects_to_remove'].split(';')
            joint = masks[old['canonical_species_a']] & masks[old['canonical_species_b']]
            source_mask = 0
            for p in source:
                source_mask |= project_bits.get(p, 0)
            excluded = joint & ~source_mask
            n0, n1 = joint.bit_count(), excluded.bit_count()
            expected = dict(baseline_support_count=base, source_excluded_support_count=left,
                fastspar_support_rank=ranks[base], support_rank_percentile=percentiles[base],
                baseline_jointly_testable_projects=n0, source_excluded_jointly_testable_projects=n1,
                baseline_support_fraction=base/n0, source_excluded_support_fraction=left/n1)
            for c, value in expected.items():
                mismatches[c] += not math.isclose(float(r[c]), value, rel_tol=1e-12, abs_tol=1e-12)
            mismatches['baseline_original_DB'] += base != match[key][0] or old['baseline_proj_ID'] != match[key][1]
            mismatches['vector_counts'] += base != len(set(old['baseline_proj_ID'].split(';'))) or left != len(set(old['remaining_proj_ID'].split(';')))
            mismatches['fraction_bounds'] += not (2 <= left <= n1 <= n0 <= 859 and base <= n0)
            for field, mask in [('baseline_eligible_project_ID', joint), ('source_excluded_eligible_project_ID', excluded),
                                ('removed_eligible_project_ID', joint & source_mask)]:
                wanted = {p for p, bit in project_bits.items() if mask & bit}
                actual = set(eligibility_rows[key][field].split(';')) - {''}
                mismatches[field] += wanted != actual
            mismatches['group'] += r['recurrence_group'] != (str(old['literature_count']) if int(old['literature_count']) < 4 else '>=4')
            for t in [1, 5, 10, 25]:
                mismatches['top_'+str(t)] += (r[f'top_{t}_percent'] == 'Yes') != (percentiles[base] >= 100-t)
            corr = [float(x) for x in old['remaining_proj_corr'].split(';')]
            pos, neg, zero = sum(x > 0 for x in corr), sum(x < 0 for x in corr), sum(x == 0 for x in corr)
            mismatches['descriptive_sign_counts'] += (pos, neg, zero) != tuple(int(r[c]) for c in ['remaining_positive_project_count', 'remaining_negative_project_count', 'remaining_zero_project_count'])
            mismatches['dominant_sign_fraction'] += not math.isclose(float(r['dominant_sign_fraction']), max(pos, neg)/left, rel_tol=1e-12)
        for field, count in mismatches.items():
            verify(a+'_all_pairs_'+field, count == 0, 'mismatches='+str(count))
    return inputs


def audit_summaries(a, rows):
    exact = {k: [r for r in rows if int(r['literature_count']) == k] for k in range(1, 6)}
    grouped = {k: [r for r in rows if r['recurrence_group'] == k] for k in GROUPS}
    for suffix in ['support', 'rank', 'support_fraction']:
        for r in read(OUT / f'summary/04d_{a}_{suffix}_by_literature_count.tsv'):
            subset = exact[int(r['literature_count'])]
            verify(a+'_summary_N_'+suffix+'_'+r['literature_count'], len(subset) == int(r['N']))
            for outcome in ['baseline_support_count'] + OUTCOMES:
                for stat, value in describe([float(t[outcome]) for t in subset]).items():
                    col = outcome+'_'+stat
                    if col in r:
                        close(a+'_'+suffix+'_'+r['literature_count']+'_'+col, value, r[col])
            for t in [1, 5, 10, 25]:
                col = f'top_{t}_percent_count'
                if col in r:
                    n = sum(x[f'top_{t}_percent'] == 'Yes' for x in subset)
                    close(a+'_exact_'+r['literature_count']+'_'+col, n, r[col])
                    close(a+'_exact_'+r['literature_count']+'_top_fraction_'+str(t), n/len(subset), r[f'top_{t}_percent_fraction'])
    for r in read(OUT / f'summary/04d_{a}_grouped_summary.tsv'):
        subset = grouped[r['recurrence_group']]
        verify(a+'_grouped_N_'+r['recurrence_group'], len(subset) == int(r['N']))
        for outcome in ['baseline_support_count'] + OUTCOMES:
            for stat, value in describe([float(t[outcome]) for t in subset]).items():
                close(a+'_grouped_'+r['recurrence_group']+'_'+outcome+'_'+stat, value, r[outcome+'_'+stat])
        for t in [1, 5, 10, 25]:
            n = sum(x[f'top_{t}_percent'] == 'Yes' for x in subset)
            close(a+'_grouped_'+r['recurrence_group']+'_top_fraction_'+str(t), n/len(subset), r[f'top_{t}_percent_fraction'])


def audit_statistics(a, rows):
    spearman = {r['outcome']: r for r in read(OUT / f'statistics/04d_{a}_spearman.tsv')}
    kw = {r['outcome']: r for r in read(OUT / f'statistics/04d_{a}_kruskal_wallis.tsv')}
    mw = read(OUT / f'statistics/04d_{a}_mann_whitney.tsv')
    for outcome in OUTCOMES:
        x = [int(r['literature_count']) for r in rows]
        y = [float(r[outcome]) for r in rows]
        rho, p = manual_spearman(x, y)
        close(a+'_'+outcome+'_Spearman_rho', rho, spearman[outcome]['rho'])
        pclose(a+'_'+outcome+'_Spearman_P', p, spearman[outcome]['two_sided_P'])
        samples = [[float(r[outcome]) for r in rows if r['recurrence_group'] == k] for k in GROUPS]
        h, p = manual_kw(samples)
        close(a+'_'+outcome+'_Kruskal_H', h, kw[outcome]['H'])
        pclose(a+'_'+outcome+'_Kruskal_P', p, kw[outcome]['P'])
        family = [r for r in mw if r['outcome'] == outcome]
        independent_p = []
        for k, sample, r in zip(GROUPS[1:], samples[1:], family):
            verify(a+'_'+outcome+'_contrast_'+k, r['comparison'] == k+'_vs_1')
            u, p = manual_mw(sample, samples[0])
            independent_p.append(p)
            close(a+'_'+outcome+'_MW_U_'+k, u, r['U'])
            pclose(a+'_'+outcome+'_MW_P_'+k, p, r['raw_P'])
            close(a+'_'+outcome+'_MW_median_A_'+k, quantile(sample, .5), r['median_group_A'])
            close(a+'_'+outcome+'_MW_median_B_'+k, quantile(samples[0], .5), r['median_group_B'])
        for r, p in zip(family, bh(independent_p)):
            pclose(a+'_'+outcome+'_MW_BH_'+r['comparison'], p, r['BH_adjusted_P'])
    fisher = read(OUT / f'statistics/04d_{a}_top10_fisher.tsv')
    independent_p = []
    for k, r in zip(GROUPS[1:], fisher):
        target = [s for s in rows if s['recurrence_group'] == k]
        reference = [s for s in rows if s['recurrence_group'] == '1']
        aa, cc = sum(s['top_10_percent'] == 'Yes' for s in target), sum(s['top_10_percent'] == 'Yes' for s in reference)
        bb, dd = len(target)-aa, len(reference)-cc
        verify(a+'_Top10_cells_'+k, [aa, bb, cc, dd] == [int(r[c]) for c in ['target_top10', 'target_not_top10', 'reference_top10', 'reference_not_top10']])
        p, estimate, low, high = exact_table(aa, bb, cc, dd)
        independent_p.append(p)
        pclose(a+'_Fisher_P_'+k, p, r['raw_P'])
        close(a+'_conditional_OR_'+k, estimate, r['odds_ratio'])
        close(a+'_exact_CI_lower_'+k, low, r['OR_95CI_lower'])
        close(a+'_exact_CI_upper_'+k, high, r['OR_95CI_upper'])
    for r, p in zip(fisher, bh(independent_p)):
        pclose(a+'_Fisher_BH_'+r['comparison'], p, r['BH_adjusted_P'])


def main():
    audit_inputs()
    for a in EXPECTED:
        rows = read(OUT / f'data/04d_{a}_recovered_pair_support_rank.tsv')
        audit_summaries(a, rows)
        audit_statistics(a, rows)
    path = OUT / 'qc/04d_independent_recalculation.tsv'
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['check', 'status', 'detail'], delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(CHECKS)
    failures = [r for r in CHECKS if r['status'] != 'PASS']
    print('Independent validation: {} checks, {} failures.'.format(len(CHECKS), len(failures)), flush=True)
    for r in failures:
        print(json.dumps(r, ensure_ascii=False))
    if failures:
        sys.exit(1)


if __name__ == '__main__':
    main()
