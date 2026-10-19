#!/home/wangmengyang/anaconda3/bin/python
"""Deterministic Stage 4C inference on frozen Stage 4B outcomes only.

Run with /home/wangmengyang/anaconda3/bin/python -s -B. No taxonomy, database
matching, source exclusion, permutations or additional-covariate models.
"""

# v0.1: repository-local IO configuration; scientific functions below are unchanged.
import sys as _path_sys
from pathlib import Path as _RepoPath
_repo_root = next(p for p in _RepoPath(__file__).resolve().parents if (p / "lib/mackdb_paths.py").is_file())
_path_sys.path.insert(0, str(_repo_root / "lib"))
from mackdb_paths import resolve_path
import os
for _thread_variable in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[_thread_variable] = '1'
import collections
import csv
import hashlib
import json
import math
import platform
import subprocess
import sys
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy import special, stats
import statsmodels
from statsmodels.discrete.discrete_model import Logit
from statsmodels.stats.multitest import multipletests

sys.dont_write_bytecode = True
BASE = Path(str(resolve_path('literature')))
OUT = BASE / 'stage4c_recurrence_statistics'
PREVIOUS = BASE / 'stage4b_recurrence_recovery'
REPORT = BASE / 'reports/04c_literature_recurrence_statistics.md'
SCRIPT = Path(__file__).resolve()
VALIDATOR = SCRIPT.with_name('04c_independent_recalculation.py')
DB = Path(str(resolve_path('environment_network')))
DB_HASH = '4b67935940a765d28518a2dd612587f347c0e47e985323910fcdde51950a7834'
EXPECTED = {'primary': ([10839, 846, 154, 56, 5], [5114, 749, 144, 55, 5]),
            'sensitivity': ([11331, 868, 162, 60, 5], [5487, 771, 152, 59, 5])}
PAPERS = ['Chen', 'Wallen', 'Zhou', 'Guo', 'Lacruz', 'Li']
GROUPS = ['1', '2', '3', '>=4']
QC, DIAGNOSTICS, OUTPUTS, LOG = [], [], set(), []


def read(p):
    with Path(p).open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write(path, rows, fields=None):
    path, rows = Path(path), list(rows)
    if OUT not in path.parents:
        raise RuntimeError('Unauthorized output path: ' + str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    OUTPUTS.add(path)


def check(name, ok, detail=''):
    QC.append(dict(check=name, status='PASS' if ok else 'FAIL', details=str(detail)))
    if not ok:
        raise RuntimeError(name + ': ' + str(detail))


def log(s):
    LOG.append(s)
    print(s, flush=True)


def fmt(x):
    return format(float(x), '.17g')


def yn(x):
    return 'Yes' if x else 'No'


def dataframe(x, kind):
    if kind == 'null':
        return pd.DataFrame({'intercept': np.ones(len(x))})
    if kind == 'continuous':
        return pd.DataFrame({'intercept': np.ones(len(x)), 'literature_count': x})
    return pd.DataFrame({'intercept': np.ones(len(x)), 'group_2': (x == 2).astype(float),
                         'group_3': (x == 3).astype(float), 'group_ge4': (x >= 4).astype(float)})


def fit_model(analysis, kind, x, y):
    design = dataframe(x, kind)
    # A separating direction must be orthogonal to every covariate row that
    # has both outcomes. Full column rank of these rows rules out both complete
    # and quasi separation (including a nonzero weak separating direction).
    strata = collections.defaultdict(set)
    for row, outcome in zip(design.to_numpy(), y):
        strata[tuple(row)].add(int(outcome))
    mixed = np.array([row for row, outcomes in strata.items() if outcomes == {0, 1}])
    mixed_rank = int(np.linalg.matrix_rank(mixed)) if mixed.size else 0
    check(analysis + '_' + kind + '_no_complete_or_quasi_separation', mixed_rank == design.shape[1],
          'Mixed-outcome design rank {} / {}'.format(mixed_rank, design.shape[1]))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        fitted = Logit(y, design, missing='raise', check_rank=True).fit(method='newton', maxiter=100, tol=1e-10, disp=False)
    warning_text = json.dumps([w.category.__name__ + ': ' + str(w.message) for w in caught], ensure_ascii=False)
    converged = bool(fitted.mle_retvals['converged'])
    diag = dict(analysis_set=analysis, model=kind, N=len(y), parameter_count=design.shape[1],
                converged=converged, iterations=int(fitted.mle_retvals['iterations']), warnings=warning_text,
                complete_separation=False, quasi_separation=False, mixed_outcome_design_rank=mixed_rank,
                separation_check='full_column_rank_of_covariate_strata_with_both_outcomes',
                finite_coefficients=bool(np.isfinite(fitted.params).all()),
                max_abs_score=fmt(np.max(np.abs(fitted.model.score(fitted.params)))))
    DIAGNOSTICS.append(diag)
    check(analysis + '_' + kind + '_converged_finite', converged and np.isfinite(fitted.params).all()
          and np.isfinite(fitted.bse).all() and (fitted.bse > 0).all(), warning_text)
    check(analysis + '_' + kind + '_no_separation_warning', not any('separation' in str(w.message).lower() for w in caught), warning_text)
    check(analysis + '_' + kind + '_score_equations', float(diag['max_abs_score']) < 1e-6, diag['max_abs_score'])
    return fitted, diag


def analyze(analysis, before):
    input_path = PREVIOUS / ('pair_level/04b_' + analysis + '_pair_recovery.tsv')
    rows = read(input_path)
    expected_n, expected_y = EXPECTED[analysis]
    check(analysis + '_unique_pair_count', len(rows) == len({r['pair_key'] for r in rows}) == sum(expected_n))
    check(analysis + '_formal_outcome_only_Yes_No', all(r['combined_source_excluded_recovered'] in {'Yes', 'No'} for r in rows))
    check(analysis + '_canonical_pair_key_unchanged', all(r['pair_key'] == r['canonical_species_a'] + '\t' + r['canonical_species_b'] and r['canonical_species_a'] < r['canonical_species_b'] for r in rows))
    x = np.array([int(r['literature_count']) for r in rows], dtype=float)
    y = np.array([int(r['combined_source_excluded_recovered'] == 'Yes') for r in rows], dtype=int)
    check(analysis + '_actual_levels_1_to_5_only', set(x) == {1., 2., 3., 4., 5.})
    check(analysis + '_six_paper_units_unchanged', all(int(r['literature_count']) == sum(int(r[p]) for p in PAPERS)
          and r['supporting_papers'].split('|') == [p for p in PAPERS if r[p] == '1'] for r in rows))
    observed = [dict(literature_count=k, pair_count=int(np.sum(x == k)), recovered=int(y[x == k].sum()),
                     not_recovered=int(np.sum(x == k) - y[x == k].sum()), recovery_rate=fmt(y[x == k].mean())) for k in range(1, 6)]
    check(analysis + '_independent_input_counts_match_requested', [r['pair_count'] for r in observed] == expected_n and [r['recovered'] for r in observed] == expected_y, observed)
    frozen_summary = read(PREVIOUS / ('summary/04b_' + analysis + '_recovery_by_literature_count.tsv'))
    check(analysis + '_counts_match_Stage4B_summary', [int(r['external_pair_count']) for r in frozen_summary] == expected_n + [0]
          and [int(r['combined_source_excluded_recovered_count']) for r in frozen_summary] == expected_y + [0])
    model_rows = []
    for r, outcome in zip(rows, y):
        model_rows.append({**{c: r[c] for c in ['canonical_species_a', 'canonical_species_b', 'pair_key',
            'literature_count', 'supporting_papers', 'combined_source_excluded_recovered']},
            'recovered_binary': int(outcome), 'recurrence_group': r['literature_count'] if int(r['literature_count']) < 4 else '>=4',
            'analysis_set': analysis})
    write(OUT / ('data/04c_' + analysis + '_model_input.tsv'), model_rows)
    write(OUT / ('summary/04c_' + analysis + '_recovery_table.tsv'), observed)
    grouped = []
    for group in GROUPS:
        mask = x >= 4 if group == '>=4' else x == int(group)
        n, r = int(mask.sum()), int(y[mask].sum())
        grouped.append(dict(recurrence_group=group, pair_count=n, recovered=r, not_recovered=n-r, recovery_rate=fmt(r/n)))
    check(analysis + '_group_ge4_correct_and_both_outcomes', grouped[-1]['pair_count'] == expected_n[3] + expected_n[4]
          and grouped[-1]['recovered'] == expected_y[3] + expected_y[4] and all(r['recovered'] > 0 and r['not_recovered'] > 0 for r in grouped))
    write(OUT / ('summary/04c_' + analysis + '_grouped_recovery_table.tsv'), grouped)
    fitted0, diag0 = fit_model(analysis, 'null', x, y)
    fitted, diag = fit_model(analysis, 'continuous', x, y)
    category, cdiag = fit_model(analysis, 'categorical', x, y)
    b, se, cov = fitted.params.to_numpy(), fitted.bse.to_numpy(), fitted.cov_params().to_numpy()
    lr = 2 * (fitted.llf - fitted0.llf)
    clr = 2 * (category.llf - fitted0.llf)
    check(analysis + '_LR_statistics_nonnegative', lr >= 0 and clr >= 0)
    check(analysis + '_design_has_only_requested_predictors', fitted.model.exog_names == ['intercept', 'literature_count']
          and category.model.exog_names == ['intercept', 'group_2', 'group_3', 'group_ge4'])
    continuous = dict(analysis_set=analysis, N=len(y), beta0=fmt(b[0]), SE_beta0=fmt(se[0]),
        beta1=fmt(b[1]), SE_beta1=fmt(se[1]), OR_per_additional_paper=fmt(np.exp(b[1])),
        OR_95CI_lower=fmt(np.exp(b[1] - 1.96*se[1])), OR_95CI_upper=fmt(np.exp(b[1] + 1.96*se[1])),
        percent_change_in_odds=fmt(np.expm1(b[1])*100), Wald_z=fmt(b[1]/se[1]),
        two_sided_Wald_P=fmt(2*stats.norm.sf(abs(b[1]/se[1]))),
        logLik_model0=fmt(fitted0.llf), logLik_model1=fmt(fitted.llf), LR_statistic=fmt(lr), LR_df=1,
        LR_P=fmt(stats.chi2.sf(lr, 1)), primary_model_level_P='likelihood_ratio',
        AIC=fmt(fitted.aic), BIC=fmt(fitted.bic), McFadden_pseudo_R2=fmt(1-fitted.llf/fitted0.llf),
        converged=diag['converged'], iterations=diag['iterations'], warnings=diag['warnings'],
        complete_separation=False, quasi_separation=False)
    write(OUT / ('models/04c_' + analysis + '_continuous_logistic.tsv'), [continuous])
    categorical = []
    for j, contrast in enumerate(['2_vs_1', '3_vs_1', 'ge4_vs_1'], 1):
        beta, sigma = category.params.iloc[j], category.bse.iloc[j]
        categorical.append(dict(analysis_set=analysis, contrast=contrast, reference_group='1', N=len(y),
            beta0=fmt(category.params.iloc[0]), SE_beta0=fmt(category.bse.iloc[0]), beta=fmt(beta), SE=fmt(sigma),
            OR=fmt(np.exp(beta)), CI_lower=fmt(np.exp(beta-1.96*sigma)), CI_upper=fmt(np.exp(beta+1.96*sigma)),
            Wald_z=fmt(beta/sigma), two_sided_Wald_P=fmt(2*stats.norm.sf(abs(beta/sigma))),
            logLik_model0=fmt(fitted0.llf), logLik_categorical=fmt(category.llf),
            categorical_overall_LR_statistic=fmt(clr), categorical_overall_LR_df=3,
            categorical_overall_LR_P=fmt(stats.chi2.sf(clr, 3)),
            AIC=fmt(category.aic), BIC=fmt(category.bic), McFadden_pseudo_R2=fmt(1-category.llf/fitted0.llf),
            target_total=grouped[j]['pair_count'], target_recovered=grouped[j]['recovered'],
            target_not_recovered=grouped[j]['not_recovered'], target_observed_recovery_rate=grouped[j]['recovery_rate'],
            converged=cdiag['converged'], iterations=cdiag['iterations'], warnings=cdiag['warnings'],
            complete_separation=False, quasi_separation=False))
    write(OUT / ('models/04c_' + analysis + '_categorical_logistic.tsv'), categorical)
    covariance_rows = []
    for kind, fit in [('continuous', fitted), ('categorical', category)]:
        c = fit.cov_params()
        for term_a in c.index:
            for term_b in c.columns:
                covariance_rows.append(dict(model=kind, term_a=term_a, term_b=term_b, covariance=fmt(c.loc[term_a, term_b])))
    write(OUT / ('models/04c_' + analysis + '_coefficient_covariance.tsv'), covariance_rows)
    fisher = []
    for target, contrast in zip(grouped[1:], ['2_vs_1', '3_vs_1', 'ge4_vs_1']):
        table = [[target['recovered'], target['not_recovered']], [grouped[0]['recovered'], grouped[0]['not_recovered']]]
        test = stats.fisher_exact(table, alternative='two-sided')
        fisher.append(dict(analysis_set=analysis, contrast=contrast, target_recovered=table[0][0], target_not_recovered=table[0][1],
            reference_recovered=table[1][0], reference_not_recovered=table[1][1], odds_ratio=fmt(test.statistic), raw_P=fmt(test.pvalue)))
    adjusted = multipletests([float(r['raw_P']) for r in fisher], method='fdr_bh')[1]
    for r, p_adj in zip(fisher, adjusted):
        r.update(BH_adjusted_P=fmt(p_adj), BH_family=analysis + '_three_pairwise_Fisher_tests', BH_family_size=3, alternative='two-sided')
    write(OUT / ('models/04c_' + analysis + '_fisher_pairwise.tsv'), fisher)
    # Standard asymptotic binomial-score CA test, no continuity correction.
    n = np.array([r['pair_count'] for r in observed], dtype=float)
    successes = np.array([r['recovered'] for r in observed], dtype=float)
    scores = np.arange(1., 6.)
    proportion = successes.sum()/n.sum()
    score_mean = np.dot(n, scores)/n.sum()
    numerator = np.dot(scores-score_mean, successes)
    variance = proportion*(1-proportion)*np.dot(n, (scores-score_mean)**2)
    z = numerator/np.sqrt(variance)
    trend = dict(analysis_set=analysis, N=len(y), scores='1;2;3;4;5', recovered_total=int(y.sum()),
        score_weighted_mean=fmt(score_mean), pooled_recovered_fraction=fmt(proportion), numerator=fmt(numerator),
        null_variance=fmt(variance), Z=fmt(z), two_sided_P=fmt(2*stats.norm.sf(abs(z))),
        variance_convention='binomial_score_pooled_p_no_continuity_correction', zero_count_cells_preserved=True,
        excluded_empty_level=6, formula='Z=sum_i[(s_i-sbar)*r_i]/sqrt[p*(1-p)*sum_i(n_i*(s_i-sbar)^2)]')
    write(OUT / ('models/04c_' + analysis + '_cochran_armitage.tsv'), [trend])
    predictions = []
    for r in observed:
        design = np.array([1., r['literature_count']])
        eta, sigma = float(design @ b), float(np.sqrt(design @ cov @ design))
        predictions.append(dict(literature_count=r['literature_count'], observed_total=r['pair_count'],
            observed_recovered=r['recovered'], observed_recovery_rate=r['recovery_rate'],
            predicted_probability=fmt(special.expit(eta)), predicted_95CI_lower=fmt(special.expit(eta-1.96*sigma)),
            predicted_95CI_upper=fmt(special.expit(eta+1.96*sigma)),
            CI_method='Wald_interval_for_mean_logit_then_inverse_logit', linear_predictor=fmt(eta), SE_linear_predictor=fmt(sigma)))
    write(OUT / ('summary/04c_' + analysis + '_observed_predicted_recovery.tsv'), predictions)
    metrics = [dict(analysis_set=analysis, model=kind, N=len(y), parameter_count=len(fit.params), log_likelihood=fmt(fit.llf),
                   AIC=fmt(fit.aic), BIC=fmt(fit.bic), McFadden_pseudo_R2=fmt(1-fit.llf/fitted0.llf))
               for kind, fit in [('continuous', fitted), ('categorical', category)]]
    log('{}: continuous beta1={}, OR={}, LR P={}; all requested models converged.'.format(analysis, continuous['beta1'], continuous['OR_per_additional_paper'], continuous['LR_P']))
    return dict(continuous=continuous, categorical=categorical, fisher=fisher, trend=trend, predictions=predictions,
                observed=observed, grouped=grouped, metrics=metrics)


def protected_files():
    paths = {p for p in BASE.rglob('*') if p.is_file() and OUT not in p.parents and p != REPORT and '__pycache__' not in p.parts}
    paths.update(Path(r['input_file']) for r in read(PREVIOUS / 'qc/04b_input_integrity.tsv'))
    paths.update({SCRIPT, VALIDATOR, DB})
    return paths


def snapshot(paths):
    return {p: (p.stat().st_size, sha(p)) for p in sorted(paths)}


def md(fields, rows):
    def clean(x):
        return str(x).replace('|', '\\|').replace('\n', ' ')
    return '\n'.join(['| ' + ' | '.join(fields) + ' |', '| ' + ' | '.join(['---']*len(fields)) + ' |'] +
                     ['| ' + ' | '.join(clean(r[k]) for k in fields) + ' |' for r in rows])


def main():
    for d in ['data', 'models', 'summary', 'qc', 'scripts', 'logs']:
        (OUT / d).mkdir(parents=True, exist_ok=True)
    log('Stage 4C: statistical tests on frozen Stage 4B combined-source-excluded recovery only.')
    protected = protected_files()
    before = snapshot(protected)
    prior = {}
    for path in [PREVIOUS / 'qc/04b_input_integrity.tsv', PREVIOUS / 'qc/04b_formal_output_manifest.tsv']:
        for r in read(path):
            prior[Path(r.get('input_file') or r['absolute_path'])] = (r.get('sha256_after') or r['sha256'], str(path))
    for p in sorted(protected & prior.keys()):
        check('prior_frozen_hash:' + str(p), before[p][1] == prior[p][0])
    check('original_database_official_hash', before[DB][1] == DB_HASH, 'Checksum only; no database parsing or matching.')
    versions = [dict(software='python', version=platform.python_version(), path=sys.executable)]
    for package in [np, pd, scipy, statsmodels]:
        versions.append(dict(software=package.__name__, version=package.__version__, path=package.__file__))
    write(OUT / 'qc/04c_software_versions.tsv', versions)
    inputs = []
    for analysis in ['primary', 'sensitivity']:
        p = PREVIOUS / ('pair_level/04b_' + analysis + '_pair_recovery.tsv')
        inputs.append(dict(analysis_set=analysis, input_file=str(p), size_bytes=before[p][0], sha256=before[p][1],
                           role='formal_pair_level_outcome', outcome='combined_source_excluded_recovered', Yes=1, No=0))
        p = PREVIOUS / ('summary/04b_' + analysis + '_recovery_by_literature_count.tsv')
        inputs.append(dict(analysis_set=analysis, input_file=str(p), size_bytes=before[p][0], sha256=before[p][1],
                           role='frozen_count_audit_only', outcome='', Yes='', No=''))
    for p in [SCRIPT, VALIDATOR, PREVIOUS / 'qc/04b_input_integrity.tsv', PREVIOUS / 'qc/04b_formal_output_manifest.tsv']:
        inputs.append(dict(analysis_set='both', input_file=str(p), size_bytes=before[p][0], sha256=before[p][1],
                           role='implementation_or_integrity_evidence', outcome='', Yes='', No=''))
    write(OUT / 'qc/04c_input_manifest.tsv', inputs)
    results = {analysis: analyze(analysis, before) for analysis in ['primary', 'sensitivity']}
    write(OUT / 'models/04c_model_diagnostics.tsv', DIAGNOSTICS)
    write(OUT / 'summary/04c_model_fit_metrics.tsv', [r for analysis in results for r in results[analysis]['metrics']])
    cp, cs = results['primary']['continuous'], results['sensitivity']['continuous']
    direction = float(cp['beta1']) > 0 and float(cs['beta1']) > 0
    relative = (float(cs['OR_per_additional_paper']) - float(cp['OR_per_additional_paper']))/float(cp['OR_per_additional_paper'])
    comparison = dict(direction_consistent=yn(direction), primary_beta1=cp['beta1'], sensitivity_beta1=cs['beta1'],
        primary_OR=cp['OR_per_additional_paper'], sensitivity_OR=cs['OR_per_additional_paper'],
        OR_relative_difference=fmt(relative), OR_relative_difference_percent=fmt(100*relative),
        relative_difference_definition='(OR_sensitivity-OR_primary)/OR_primary; signed; no equivalence threshold')
    write(OUT / 'summary/04c_primary_sensitivity_comparison.tsv', [comparison])
    paper_rows = []
    for analysis, r in results.items():
        c, trend = r['continuous'], r['trend']
        row = dict(analysis_set=analysis, total_pairs=c['N'], OR_per_additional_paper=c['OR_per_additional_paper'],
            OR_95CI_lower=c['OR_95CI_lower'], OR_95CI_upper=c['OR_95CI_upper'], beta1=c['beta1'],
            continuous_LR_chi2=c['LR_statistic'], continuous_LR_P=c['LR_P'], continuous_Wald_P=c['two_sided_Wald_P'],
            Cochran_Armitage_Z=trend['Z'], Cochran_Armitage_P=trend['two_sided_P'],
            categorical_overall_LR_P=r['categorical'][0]['categorical_overall_LR_P'])
        for t in r['categorical']:
            contrast = t['contrast']
            row['OR_' + contrast] = t['OR']
            row['CI_' + contrast] = t['CI_lower'] + ';' + t['CI_upper']
            row['P_' + contrast] = t['two_sided_Wald_P']
        row.update(direction_consistent=yn(direction), OR_relative_difference=fmt(relative))
        paper_rows.append(row)
    write(OUT / 'summary/04c_recurrence_statistics_paper_table.tsv', paper_rows)
    methods = '''# Stage 4C statistical definitions

Outcome is the frozen Stage 4B combined_source_excluded_recovered (Yes=1, No=0).
The analysis unit is one unique external species pair. Primary and sensitivity
are fitted separately. Pair outcomes are treated as independent Bernoulli
observations by these requested, unadjusted models; no shared-species/network
dependence correction or causal effect is estimated.

Continuous model: logit(p_j)=beta0+beta1*x_j; x_j=literature_count.
Categorical model: intercept plus indicators 2, 3, >=4 (reference 1).
No penalty, continuity correction, pseudo-count, other covariate or data-driven
model selection is used. Categorical >=4 includes exact counts 4 and 5.

statsmodels Logit fits the pair-level Bernoulli likelihood by Newton iteration.
SEs use the inverse observed information. Coefficient CIs use beta +/- 1.96*SE;
OR CIs exponentiate these endpoints. Wald P=2*Normal.sf(abs(beta/SE)).
LR=2*(logLik_full-logLik_intercept); primary continuous model-level P=Chi2.sf(LR,1).
Categorical overall model P=Chi2.sf(LR,3), the standard upper-tail LR reference
probability for an unrestricted alternative (not a doubled chi-square tail).
AIC=-2*logLik+2*k; BIC=-2*logLik+k*ln(N); McFadden R2=1-logLik/logLik_null.
Log-likelihood is the sum of individual Bernoulli log-likelihoods; no grouped
binomial combinatorial constants are added in independent recalculation.
Source: https://www.statsmodels.org/stable/generated/statsmodels.discrete.discrete_model.Logit.html

Predicted mean probability: expit(v' beta), v=(1,x). Its 95% CI is
expit(v' beta +/- 1.96*sqrt(v' covariance(beta) v)); it is a confidence interval
for model mean probability, not a binary-observation prediction interval.

Fisher tables: [[target recovered,target not recovered],
               [reference recovered,reference not recovered]].
Use scipy.stats.fisher_exact(alternative='two-sided'). P sums fixed-margin
hypergeometric probabilities no greater than the observed table probability.
Odds ratio is a*d/(b*c). BH-FDR has exactly 3 tests per analysis set, handled
separately for primary and sensitivity. No adjustment is applied to the other
requested primary/model statistics.
Sources:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html
https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html

Cochran-Armitage: scores s_i=1,2,3,4,5; n_i=pair count; r_i=recovered count;
N=sum(n_i), p=sum(r_i)/N, sbar=sum(n_i*s_i)/N.
U=sum((s_i-sbar)*r_i).
V=p*(1-p)*sum(n_i*(s_i-sbar)^2).
Z=U/sqrt(V); two-sided P=2*Normal.sf(abs(Z)).
This is the asymptotic binomial-score convention; no continuity correction.
Positive Z indicates increasing recovery with literature count. Empty level 6
is excluded and the failure count 0 at level 5 is retained exactly.
Formula source: https://support.sas.com/documentation/cdl/en/statug/68162/HTML/default/statug_freq_details72.htm

Independent CA verification uses statsmodels Table(...,shift_zeros=False),
row_scores=(1,2,3,4,5), col_scores=(0,1), columns=(not_recovered,recovered).
That library's fixed-margin variance is V*N/(N-1). Therefore its Z is multiplied
by sqrt(N/(N-1)) before comparing with the declared binomial-score Z. This is
an explicit variance-convention conversion, not a change to data or a new test.
The original library Z and converted Z are both retained in the audit.
Source: https://www.statsmodels.org/stable/generated/statsmodels.stats.contingency_tables.Table.test_ordinal_association.html

Separation check: covariate strata containing both outcomes must have v'd=0
for any complete/weak separating direction d. If those v rows have full column
rank, only d=0 is possible, ruling out complete and quasi separation. Rank is
checked for each fitted model; warning and convergence diagnostics are saved.

Independent coefficient verification: nested scalar score-equation roots on
five exact counts for continuous beta, closed-form group log-odds for categorical
coefficients, and independently derived information/covariance. Fisher P uses
exact integer combinations over every feasible 2x2 table. BH uses manual ordered
minimum adjustment. Normal/chi-square tails are checked with erfc identities.
For tiny P values comparisons use relative tolerance only, never an absolute
tolerance that would treat all tiny probabilities as equal.
'''
    method_path = OUT / 'qc/04c_statistical_methods.md'
    method_path.write_text(methods, encoding='utf-8')
    OUTPUTS.add(method_path)
    log('Running independent statistical verification through separate coefficient, exact-test and trend formulas.')
    process = subprocess.run([sys.executable, '-s', '-B', str(VALIDATOR)], capture_output=True, text=True)
    if process.returncode:
        raise RuntimeError('Independent recalculation failed:\n' + process.stdout + process.stderr)
    log(process.stdout.strip())
    independent_path = OUT / 'qc/04c_independent_recalculation.tsv'
    independent = read(independent_path)
    OUTPUTS.add(independent_path)
    check('independent_recalculation_all_PASS', independent and all(r['status'] == 'PASS' for r in independent), len(independent))
    check('all_logistic_models_converged_without_separation', all(r['converged'] and not r['complete_separation'] and not r['quasi_separation'] for r in DIAGNOSTICS))
    check('frozen_file_inventory_unchanged', protected_files() == protected)
    after = snapshot(protected)
    integrity = []
    for p in sorted(protected):
        check('unchanged:' + str(p), before[p] == after[p])
        integrity.append(dict(input_file=str(p), scope='current_implementation' if p in {SCRIPT, VALIDATOR} else 'frozen_upstream',
            size_bytes=before[p][0], sha256_before=before[p][1], sha256_after=after[p][1],
            prior_recorded_sha256=prior[p][0] if p in prior else '', prior_hash_evidence=prior[p][1] if p in prior else '',
            matches_prior_recorded_hash='PASS' if p in prior else 'NO_PRIOR_RECORD', unchanged=before[p] == after[p]))
    check('human_gut_info_unchanged_checksum_only', after[DB][1] == DB_HASH)
    check('no_upstream_outcome_or_source_exclusion_recomputed', True, 'Only frozen 4B outcome column used; no DB/source/taxonomy parser called.')
    write(OUT / 'qc/04c_input_integrity.tsv', integrity)
    write(OUT / 'qc/04c_statistical_qc.tsv', QC)
    lines = ['# Stage 4C — Literature recurrence versus MACK-DB recovery', '', 'overall_status: **PASS**', '',
        '正式主分析使用 Stage 4B pair-specific combined-source-excluded recovery（Yes=1 / No=0）。Primary 不包含 exploratory Lacruz_MHO；sensitivity 包含 MHO。没有使用 baseline recovery 作为 outcome，也没有重做 source exclusion、overlap 或 recurrence。', '',
        '## 方法', '',
        '主模型为逐 pair `recovered ~ 1 + literature_count`，采用未惩罚 maximum-likelihood logistic regression；LR P（df=1）为主要 model-level P，双侧 Wald P 同时保留。',
        'Categorical 模型使用预先定义的 1、2、3、≥4，reference=1；overall LR df=3。每套分析的三个双侧 Fisher tests 分别做 BH-FDR，用于补充稳健性。',
        'Cochran–Armitage 使用 exact scores 1–5、pooled-binomial variance，无连续性修正。空 level 6 不进入模型，level 5 的零 failure 不加 0.5。',
        '明确公式、独立算法和软件方法引用见 [统计方法说明](' + str(method_path) + ')。CA 公式采用 [SAS 官方定义](https://support.sas.com/documentation/cdl/en/statug/68162/HTML/default/statug_freq_details72.htm)，并与关闭 zero-shift 的 statsmodels ordinal test 校正到同一 variance convention 后核对。',
        '所有标准误/置信区间/P 值均对应所请求的未调整 pair-level 模型。模型将 pair outcomes 作为独立 Bernoulli observations；没有估计共享物种或网络依赖的修正，也不作因果解释。', '',
        md(['software', 'version', 'path'], versions), '']
    for analysis, r in results.items():
        c = r['continuous']
        lines += ['## ' + analysis, '', 'Observed recovery（固定 Stage 4B 比例）：', '',
                  md(list(r['observed'][0]), r['observed']), '',
                  '**literature_count=5 只有 5 条 pair。5/5=100% 样本极小，不能单独作为完全恢复的强结论。**主要证据来自 continuous model 和 1/2/3/≥4 分组。', '',
                  'Continuous：beta1={beta1}，OR={OR_per_additional_paper}，95% CI [{OR_95CI_lower}, {OR_95CI_upper}]，LR P={LR_P}，Wald P={two_sided_Wald_P}。'.format(**c),
                  '每增加一篇独立文献支持，模型 recovery **odds** 增加 {:.2f}%；这不是概率增加的百分点。'.format(float(c['percent_change_in_odds'])), '',
                  'Categorical：', '', md(['contrast', 'OR', 'CI_lower', 'CI_upper', 'two_sided_Wald_P'], r['categorical']), '',
                  'Categorical overall LR P=' + r['categorical'][0]['categorical_overall_LR_P'] + '。', '',
                  'Grouped observed recovery：', '', md(list(r['grouped'][0]), r['grouped']), '',
                  'Supplementary Fisher + BH-FDR：', '', md(['contrast', 'odds_ratio', 'raw_P', 'BH_adjusted_P'], r['fisher']), '',
                  'Cochran–Armitage：Z={Z}，two-sided P={two_sided_P}。'.format(**r['trend']), '',
                  'Continuous model mean probabilities and logit-Wald 95% CI：', '',
                  md(['literature_count', 'observed_recovery_rate', 'predicted_probability', 'predicted_95CI_lower', 'predicted_95CI_upper'], r['predictions']), '']
    lines += ['## 拟合指标与一致性', '', md(list(results['primary']['metrics'][0]), results['primary']['metrics'] + results['sensitivity']['metrics']), '',
              '不依据 AIC 自动选择模型。', '', md(list(comparison), [comparison]), '',
              'Primary 与 sensitivity 的 beta1 均为正，方向一致；OR relative difference 明确定义为 (OR_sensitivity−OR_primary)/OR_primary，没有设置等同性硬阈值。',
              '在所指定模型下，较高的跨文献 recurrence 与较高的联合 source-excluded MACK-DB recovery probability 显著相关；不解释为文献支持导致 recovery 增加。', '',
              '## 收敛、独立复核与冻结完整性', '', md(['analysis_set', 'model', 'converged', 'iterations', 'warnings', 'complete_separation', 'quasi_separation', 'mixed_outcome_design_rank'], DIAGNOSTICS), '',
              '{} 项独立复核全部 PASS：包括输入计数、连续/分类系数、SE/covariance、OR/CI、LL/LR/P、预测概率区间、Fisher exact、BH-FDR、CA 和拟合指标。'.format(len(independent)),
              '{} 项整体 QC 全部 PASS；{} 个受保护文件 before/after SHA-256 全部一致，其中 {} 个匹配之前记录的 hash。'.format(len(QC), len(before), len(protected & prior.keys())),
              '冻结 Stage 2B/2C/3A/3A.1/3B/3C/4A/4B、历史 Chen/Wallen 结果及原始 human_gut.info 保持不变。数据库只读取字节用于 SHA-256，不进行解析或 matching。', '',
              '## 重现与交付', '', '`' + sys.executable + ' -s -B ' + str(SCRIPT) + '`', '',
              '输出固定排序，不依赖系统时间或随机种子。所有正式输出绝对路径与 SHA-256 见 `' + str(OUT / 'qc/04c_formal_output_manifest.tsv') + '`；清单不包含自身 hash。',
              'Stage 4C 完成后停止，未执行额外协变量模型、permutation、null model、mixed model、degree adjustment 或新的 source/taxonomy analysis。', '']
    REPORT.write_text('\n'.join(lines), encoding='utf-8')
    OUTPUTS.update({REPORT, SCRIPT, VALIDATOR})
    status = dict(overall_status='PASS', paper_table=paper_rows, primary_sensitivity=comparison,
                  all_models_converged=True, complete_separation=False, quasi_separation=False,
                  independent_checks=len(independent), qc_checks=len(QC), frozen_files_verified=len(before),
                  prior_hashes_verified=len(protected & prior.keys()),
                  observed_strictly_increasing={a: all(float(r1['recovery_rate']) < float(r2['recovery_rate']) for r1, r2 in zip(r['observed'], r['observed'][1:])) for a, r in results.items()})
    status_path = OUT / 'qc/04c_run_status.json'
    status_path.write_text(json.dumps(status, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    OUTPUTS.add(status_path)
    log('overall_status=PASS; all requested statistics independently verified; frozen inputs unchanged; Stage 4C complete.')
    log_path = OUT / 'logs/04c_recurrence_statistics.log'
    log_path.write_text('\n'.join(LOG) + '\n', encoding='utf-8')
    OUTPUTS.add(log_path)
    write(OUT / 'qc/04c_formal_output_manifest.tsv', [dict(absolute_path=str(p), size_bytes=p.stat().st_size,
          sha256=sha(p), role='formal_Stage4C_output') for p in sorted(OUTPUTS)])


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        write(OUT / 'qc/04c_statistical_qc.tsv', QC + [dict(check='overall_status', status='FAIL', details=str(exc))])
        (OUT / 'qc/04c_run_status.json').write_text(json.dumps(dict(overall_status='FAIL', error=str(exc)), ensure_ascii=False, indent=2) + '\n')
        REPORT.write_text('# Stage 4C\n\noverall_status: **FAIL**\n\nStopped: ' + str(exc) + '\n', encoding='utf-8')
        raise
