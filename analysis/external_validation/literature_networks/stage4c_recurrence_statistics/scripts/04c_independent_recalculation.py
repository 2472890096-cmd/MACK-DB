#!/home/wangmengyang/anaconda3/bin/python
"""Independent formulas/read-back checks; never import the main statistics script."""

# v0.1: repository-local IO configuration; scientific functions below are unchanged.
import sys as _path_sys
from pathlib import Path as _RepoPath
_repo_root = next(p for p in _RepoPath(__file__).resolve().parents if (p / "lib/mackdb_paths.py").is_file())
_path_sys.path.insert(0, str(_repo_root / "lib"))
from mackdb_paths import resolve_path
import os
for _name in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[_name] = '1'
import collections
import csv
import math
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from statsmodels.stats.contingency_tables import Table

BASE = Path(str(resolve_path('literature')))
OUT = BASE / 'stage4c_recurrence_statistics'
PREVIOUS = BASE / 'stage4b_recurrence_recovery'
CHECKS = []


def read(path):
    with Path(path).open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def check(name, ok, actual='', expected='', tolerance='exact'):
    CHECKS.append(dict(check=name, formal_value=str(actual), independent_value=str(expected),
                       tolerance=tolerance, status='PASS' if ok else 'FAIL'))
    if not ok:
        raise AssertionError(name + ': ' + str(actual) + ' != ' + str(expected))


def compare(name, actual, expected, probability=False):
    a, e = float(actual), float(expected)
    rtol, atol = (2e-7, 0.) if probability else (2e-9, 2e-9)
    check(name, math.isfinite(a) and math.isfinite(e) and math.isclose(a, e, rel_tol=rtol, abs_tol=atol),
          format(a, '.17g'), format(e, '.17g'), 'relative={};absolute={}'.format(rtol, atol))


def logistic(t):
    if t >= 0:
        return 1/(1+math.exp(-t))
    e = math.exp(t)
    return e/(1+e)


def normal_two_sided(z):
    return math.erfc(abs(z)/math.sqrt(2))


def chi_square_tail(statistic, df):
    if df == 1:
        return math.erfc(math.sqrt(statistic/2))
    if df == 3:
        return math.erfc(math.sqrt(statistic/2)) + math.sqrt(2*statistic/math.pi)*math.exp(-statistic/2)
    raise ValueError('Only requested df=1 and df=3 are supported')


def continuous_roots(n, success):
    """Nested scalar score roots, independent of pair-level Newton optimization."""
    scores = list(range(1, 6))
    total_success = sum(success)

    def intercept_for_slope(slope):
        def intercept_score(intercept):
            return math.fsum(ni*logistic(intercept+slope*x) for ni, x in zip(n, scores))-total_success
        return brentq(intercept_score, -200., 200., xtol=1e-13, rtol=1e-14)

    def slope_score(slope):
        intercept = intercept_for_slope(slope)
        return math.fsum(x*(r-ni*logistic(intercept+slope*x)) for x, ni, r in zip(scores, n, success))

    slope = brentq(slope_score, -20., 20., xtol=1e-13, rtol=1e-14)
    intercept = intercept_for_slope(slope)
    probs = [logistic(intercept+slope*x) for x in scores]
    weights = [ni*p*(1-p) for ni, p in zip(n, probs)]
    a = math.fsum(weights)
    b = math.fsum(w*x for w, x in zip(weights, scores))
    c = math.fsum(w*x*x for w, x in zip(weights, scores))
    det = a*c-b*b
    covariance = np.array([[c, -b], [-b, a]])/det
    ll = math.fsum(r*math.log(p)+(ni-r)*math.log1p(-p) for ni, r, p in zip(n, success, probs))
    return [intercept, slope], covariance, ll


def exact_fisher(n_target, success_target, n_reference, success_reference):
    """Exact integer hypergeometric mass ordering, without scipy Fisher/hypergeom."""
    population = n_target+n_reference
    all_success = success_target+success_reference
    low, high = max(0, n_target-(population-all_success)), min(n_target, all_success)
    observed_weight = math.comb(all_success, success_target)*math.comb(population-all_success, n_target-success_target)
    tail_weight = 0
    for x in range(low, high+1):
        weight = math.comb(all_success, x)*math.comb(population-all_success, n_target-x)
        if weight <= observed_weight:
            tail_weight += weight
    probability = tail_weight/math.comb(population, n_target)
    odds_ratio = success_target*(n_reference-success_reference)/((n_target-success_target)*success_reference)
    return odds_ratio, probability


def main():
    results = {}
    for analysis in ['primary', 'sensitivity']:
        original_rows = read(PREVIOUS / ('pair_level/04b_' + analysis + '_pair_recovery.tsv'))
        original = {r['pair_key']: r for r in original_rows}
        data = read(OUT / ('data/04c_' + analysis + '_model_input.tsv'))
        check(analysis + '_unique_frozen_input_identity', len(data) == len(original_rows) == len(original) == len({r['pair_key'] for r in data}) and {r['pair_key'] for r in data} == set(original))
        check(analysis + '_unchanged_outcome_and_paper_count', all(
            r['literature_count'] == original[r['pair_key']]['literature_count'] and
            r['supporting_papers'] == original[r['pair_key']]['supporting_papers'] and
            r['combined_source_excluded_recovered'] == original[r['pair_key']]['combined_source_excluded_recovered'] and
            int(r['recovered_binary']) == int(original[r['pair_key']]['combined_source_excluded_recovered'] == 'Yes') and
            r['recurrence_group'] == ('>=4' if int(r['literature_count']) >= 4 else r['literature_count']) for r in data))
        pair_count = collections.Counter(int(r['literature_count']) for r in original_rows)
        success_count = collections.Counter(int(r['literature_count']) for r in original_rows if r['combined_source_excluded_recovered'] == 'Yes')
        n, success = [pair_count[k] for k in range(1, 6)], [success_count[k] for k in range(1, 6)]
        total, all_success = sum(n), sum(success)
        check(analysis + '_levels_1_to_5_without_6', set(pair_count) == {1, 2, 3, 4, 5})
        prior_summary = read(PREVIOUS / ('summary/04b_' + analysis + '_recovery_by_literature_count.tsv'))
        for r in prior_summary:
            k = int(r['literature_count'])
            check(analysis + '_Stage4B_counts_' + str(k), int(r['external_pair_count']) == pair_count[k] and int(r['combined_source_excluded_recovered_count']) == success_count[k])
        observed = read(OUT / ('summary/04c_' + analysis + '_recovery_table.tsv'))
        check(analysis + '_exact_recovery_rows', len(observed) == 5)
        for r in observed:
            k = int(r['literature_count'])
            check(analysis + '_observed_counts_' + str(k), int(r['pair_count']) == pair_count[k] and int(r['recovered']) == success_count[k] and int(r['not_recovered']) == pair_count[k]-success_count[k])
            compare(analysis + '_observed_rate_' + str(k), r['recovery_rate'], success_count[k]/pair_count[k])
        group_n, group_y = n[:3]+[n[3]+n[4]], success[:3]+[success[3]+success[4]]
        grouped = read(OUT / ('summary/04c_' + analysis + '_grouped_recovery_table.tsv'))
        check(analysis + '_grouped_counts', len(grouped) == 4 and [int(r['pair_count']) for r in grouped] == group_n and [int(r['recovered']) for r in grouped] == group_y and [r['recurrence_group'] for r in grouped] == ['1', '2', '3', '>=4'])
        beta, covariance, ll = continuous_roots(n, success)
        se = [math.sqrt(covariance[j, j]) for j in range(2)]
        p0 = all_success/total
        ll0 = all_success*math.log(p0)+(total-all_success)*math.log1p(-p0)
        lr = 2*(ll-ll0)
        continuous = read(OUT / ('models/04c_' + analysis + '_continuous_logistic.tsv'))[0]
        expected = dict(beta0=beta[0], SE_beta0=se[0], beta1=beta[1], SE_beta1=se[1],
            OR_per_additional_paper=math.exp(beta[1]), OR_95CI_lower=math.exp(beta[1]-1.96*se[1]),
            OR_95CI_upper=math.exp(beta[1]+1.96*se[1]), percent_change_in_odds=math.expm1(beta[1])*100,
            Wald_z=beta[1]/se[1], two_sided_Wald_P=normal_two_sided(beta[1]/se[1]), logLik_model0=ll0,
            logLik_model1=ll, LR_statistic=lr, LR_P=chi_square_tail(lr, 1), AIC=-2*ll+4,
            BIC=-2*ll+2*math.log(total), McFadden_pseudo_R2=1-ll/ll0)
        for field, value in expected.items():
            compare(analysis + '_continuous_' + field, continuous[field], value, probability=field.endswith('_P'))
        check(analysis + '_continuous_fitted_N_and_df', int(continuous['N']) == total and continuous['LR_df'] == '1')
        group_log_odds = [math.log(r/(ni-r)) for ni, r in zip(group_n, group_y)]
        group_variance = [1/r+1/(ni-r) for ni, r in zip(group_n, group_y)]
        cat_beta = [group_log_odds[0]] + [a-group_log_odds[0] for a in group_log_odds[1:]]
        cat_cov = np.empty((4, 4))
        cat_cov[0, 0] = group_variance[0]
        cat_cov[0, 1:] = cat_cov[1:, 0] = -group_variance[0]
        cat_cov[1:, 1:] = group_variance[0]
        for j in range(1, 4):
            cat_cov[j, j] += group_variance[j]
        cat_ll = math.fsum(r*math.log(r/ni)+(ni-r)*math.log((ni-r)/ni) for ni, r in zip(group_n, group_y))
        cat_lr = 2*(cat_ll-ll0)
        categorical = read(OUT / ('models/04c_' + analysis + '_categorical_logistic.tsv'))
        check(analysis + '_three_categorical_contrasts', [r['contrast'] for r in categorical] == ['2_vs_1', '3_vs_1', 'ge4_vs_1'])
        for j, r in enumerate(categorical, 1):
            sigma = math.sqrt(cat_cov[j, j])
            values = dict(beta0=cat_beta[0], SE_beta0=math.sqrt(cat_cov[0, 0]), beta=cat_beta[j], SE=sigma,
                OR=math.exp(cat_beta[j]), CI_lower=math.exp(cat_beta[j]-1.96*sigma), CI_upper=math.exp(cat_beta[j]+1.96*sigma),
                Wald_z=cat_beta[j]/sigma, two_sided_Wald_P=normal_two_sided(cat_beta[j]/sigma),
                logLik_model0=ll0, logLik_categorical=cat_ll, categorical_overall_LR_statistic=cat_lr,
                categorical_overall_LR_P=chi_square_tail(cat_lr, 3), AIC=-2*cat_ll+8,
                BIC=-2*cat_ll+4*math.log(total), McFadden_pseudo_R2=1-cat_ll/ll0)
            for field, v in values.items():
                compare(analysis + '_' + r['contrast'] + '_' + field, r[field], v, probability=field.endswith('_P'))
            check(analysis + '_' + r['contrast'] + '_N_df_and_group_counts', int(r['N']) == total and r['categorical_overall_LR_df'] == '3' and int(r['target_total']) == group_n[j] and int(r['target_recovered']) == group_y[j])
        for r in read(OUT / ('models/04c_' + analysis + '_coefficient_covariance.tsv')):
            terms = ['intercept', 'literature_count'] if r['model'] == 'continuous' else ['intercept', 'group_2', 'group_3', 'group_ge4']
            cov = covariance if r['model'] == 'continuous' else cat_cov
            compare(analysis + '_' + r['model'] + '_covariance_' + r['term_a'] + '_' + r['term_b'], r['covariance'], cov[terms.index(r['term_a']), terms.index(r['term_b'])])
        fisher = read(OUT / ('models/04c_' + analysis + '_fisher_pairwise.tsv'))
        fisher_probs = []
        for j, r in enumerate(fisher, 1):
            odds, prob = exact_fisher(group_n[j], group_y[j], group_n[0], group_y[0])
            fisher_probs.append(prob)
            compare(analysis + '_Fisher_' + r['contrast'] + '_odds', r['odds_ratio'], odds)
            compare(analysis + '_Fisher_' + r['contrast'] + '_exact_integer_P', r['raw_P'], prob, probability=True)
        order = sorted(range(3), key=lambda j: fisher_probs[j])
        adjusted = [0.]*3
        running_min = 1.
        for rank in range(3, 0, -1):
            j = order[rank-1]
            running_min = min(running_min, fisher_probs[j]*3/rank)
            adjusted[j] = running_min
        for j, r in enumerate(fisher):
            compare(analysis + '_BH_' + r['contrast'], r['BH_adjusted_P'], adjusted[j], probability=True)
        table = np.column_stack((np.array(n)-np.array(success), success)).astype(float)
        original_table = table.copy()
        ordinal = Table(table, shift_zeros=False).test_ordinal_association(row_scores=np.arange(1., 6.), col_scores=np.array([0., 1.]))
        check(analysis + '_CA_zero_cells_unchanged', np.array_equal(table, original_table))
        z_independent = ordinal.zscore*math.sqrt(total/(total-1))
        trend = read(OUT / ('models/04c_' + analysis + '_cochran_armitage.tsv'))[0]
        compare(analysis + '_CA_Z_after_explicit_variance_conversion', trend['Z'], z_independent)
        compare(analysis + '_CA_P_erfc', trend['two_sided_P'], normal_two_sided(z_independent), probability=True)
        check(analysis + '_CA_variance_convention_audit', True, ordinal.zscore, z_independent, 'fixed-margin library Z multiplied by sqrt(N/(N-1)); no data change')
        predictions = read(OUT / ('summary/04c_' + analysis + '_observed_predicted_recovery.tsv'))
        for r in predictions:
            k = int(r['literature_count'])
            eta = beta[0]+beta[1]*k
            variance = covariance[0, 0]+2*k*covariance[0, 1]+k*k*covariance[1, 1]
            sigma = math.sqrt(variance)
            for field, value in [('predicted_probability', logistic(eta)), ('predicted_95CI_lower', logistic(eta-1.96*sigma)), ('predicted_95CI_upper', logistic(eta+1.96*sigma))]:
                compare(analysis + '_prediction_' + str(k) + '_' + field, r[field], value)
        results[analysis] = expected
    comparison = read(OUT / 'summary/04c_primary_sensitivity_comparison.tsv')[0]
    relative = results['sensitivity']['OR_per_additional_paper']/results['primary']['OR_per_additional_paper']-1
    compare('primary_sensitivity_OR_relative_difference', comparison['OR_relative_difference'], relative)
    check('positive_direction_consistency', comparison['direction_consistent'] == ('Yes' if results['primary']['beta1'] > 0 and results['sensitivity']['beta1'] > 0 else 'No'))
    diagnostics = read(OUT / 'models/04c_model_diagnostics.tsv')
    check('six_models_converged_no_complete_or_quasi_separation', len(diagnostics) == 6 and all(r['converged'] == 'True' and r['complete_separation'] == 'False' and r['quasi_separation'] == 'False' and int(r['mixed_outcome_design_rank']) == int(r['parameter_count']) for r in diagnostics))
    print('Independent recalculation PASS: {} checks; nested score roots, closed-form group odds, exact-integer Fisher, manual BH, and independently converted CA.'.format(len(CHECKS)))


if __name__ == '__main__':
    try:
        main()
    finally:
        with (OUT / 'qc/04c_independent_recalculation.tsv').open('w', newline='', encoding='utf-8') as f:
            w = csv.DictWriter(f, fieldnames=['check', 'formal_value', 'independent_value', 'tolerance', 'status'], delimiter='\t', lineterminator='\n')
            w.writeheader()
            w.writerows(CHECKS)
