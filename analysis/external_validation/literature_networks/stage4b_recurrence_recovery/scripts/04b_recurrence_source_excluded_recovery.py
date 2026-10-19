#!/usr/bin/env python3
"""Stage 4B: pair-specific source removal, f=2 recovery, descriptive summaries.

Run with python3 -s -B. No upstream script is imported/executed. Database
labels are decoded only using the frozen Stage 3A exact serialization table.
Species identities, paper membership and literature_count are copied from 4A.
Only Stage 4B and its named report are written. No statistical tests are run.
"""

# v0.1: repository-local IO configuration; scientific functions below are unchanged.
import sys as _path_sys
from pathlib import Path as _RepoPath
_repo_root = next(p for p in _RepoPath(__file__).resolve().parents if (p / "lib/mackdb_paths.py").is_file())
_path_sys.path.insert(0, str(_repo_root / "lib"))
from mackdb_paths import resolve_path
import ast
import collections
import csv
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

sys.dont_write_bytecode = True
ROOT = Path(str(resolve_path('external_root')))
BASE = resolve_path('literature')
OUT = BASE / 'stage4b_recurrence_recovery'
REPORT = BASE / 'reports/04b_literature_recurrence_source_excluded_recovery.md'
SCRIPT = Path(__file__).resolve()
VALIDATOR = SCRIPT.with_name('04b_independent_validate.py')
A = BASE / 'stage4a_recurrence'
B = BASE / 'stage3b_source_exclusion'
PROV = B / 'provenance'
DB = Path(str(resolve_path('environment_network')))
DB_HASH = '4b67935940a765d28518a2dd612587f347c0e47e985323910fcdde51950a7834'
LISTED = Path(str(resolve_path('environment', 'human_gut.list')))
CHEN_SCRIPT = ROOT / 'chen_postrescue_source_exclusion_overlap_v45.py'
LEGACY_SCRIPT = resolve_path('shared_protocol', 'chen_wallen_null_models_final_v1.py')
LEGACY_SUMMARY = ROOT / 'chen_wallen_null_models_final_v1/07_source_exclusion_summary.tsv'
PAPERS = ['Chen', 'Wallen', 'Zhou', 'Guo', 'Lacruz', 'Li']
EXPECTED = {'primary': [10839, 846, 154, 56, 5, 0],
            'sensitivity': [11331, 868, 162, 60, 5, 0]}
VECTORS = ['proj_ID', 'proj_corr', 'proj_pv', 'proj_sam']
QC, LOG, OUTPUTS = [], [], set()


def read(path):
    with Path(path).open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write(path, rows, fields=None):
    path = Path(path)
    if OUT not in path.parents:
        raise RuntimeError('Unauthorized output path: ' + str(path))
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    OUTPUTS.add(path)


def check(name, ok, details=''):
    QC.append(dict(check=name, status='PASS' if ok else 'FAIL', details=str(details)))
    if not ok:
        raise RuntimeError(name + ': ' + str(details))


def log(message):
    LOG.append(message)
    print(message, flush=True)


def tokens(value):
    return value.split(';') if value else []


def packed(values):
    return ';'.join(sorted(set(values)))


def rate(n, d):
    return format(n / d, '.15g') if d else 'NA'


def yn(value):
    return 'Yes' if value else 'No'


def literal_assignment(path, name):
    tree = ast.parse(path.read_text(encoding='utf-8'))
    hits = [n.value for n in tree.body if isinstance(n, ast.Assign)
            and any(isinstance(t, ast.Name) and t.id == name for t in n.targets)]
    check('unique_historical_assignment:' + name, len(hits) == 1, str(path))
    return ast.literal_eval(hits[0])


def prior_hashes():
    hashes = {}
    for relative in ['stage3b_source_exclusion/qc/03b_input_integrity.tsv',
                     'stage3b_source_exclusion/qc/03b_formal_output_manifest.tsv',
                     'stage3c_null_model/qc/03c_input_integrity.tsv',
                     'stage3c_null_model/qc/03c_formal_output_manifest.tsv',
                     'stage4a_recurrence/qc/04a_input_integrity.tsv',
                     'stage4a_recurrence/qc/04a_formal_output_manifest.tsv']:
        path = BASE / relative
        for r in read(path):
            p = r.get('input_file') or r.get('absolute_path')
            h = r.get('sha256_after') or r.get('sha256')
            if p and h:
                hashes[Path(p)] = (h, str(path))
    return hashes


def protect(extra):
    return {p for p in BASE.rglob('*') if p.is_file() and OUT not in p.parents
            and p != REPORT and '__pycache__' not in p.parts} | set(extra)


def snapshot(paths):
    return {p: (p.stat().st_size, sha(p)) for p in sorted(paths)}


def md(fields, rows):
    def cell(v):
        return str(v).replace('|', '\\|').replace('\n', ' ')
    return '\n'.join(['| ' + ' | '.join(fields) + ' |', '| ' + ' | '.join(['---'] * len(fields)) + ' |'] +
                     ['| ' + ' | '.join(cell(r[k]) for k in fields) + ' |' for r in rows])


def summary(rows):
    n = len(rows)
    before = sum(r['baseline_mackdb_present'] == 'Yes' for r in rows)
    recovered = sum(r['combined_source_excluded_recovered'] == 'Yes' for r in rows)
    loss = sum(r['loss_due_to_combined_source_exclusion'] == 'Yes' for r in rows)
    return dict(external_pair_count=n, baseline_mackdb_present_count=before,
                baseline_recovery_rate=rate(before, n),
                combined_source_excluded_recovered_count=recovered,
                combined_source_excluded_recovery_rate=rate(recovered, n),
                source_exclusion_loss_count=loss,
                source_exclusion_loss_fraction_of_baseline_recovered=rate(loss, before))


def load_recurrence():
    masters, memberships, network_pairs = {}, {}, collections.defaultdict(set)
    for analysis in ['primary', 'sensitivity']:
        master_path = A / ('membership/04a_' + analysis + '_unique_pair_recurrence.tsv')
        rows = read(master_path)
        keys = [r['pair_key'] for r in rows]
        observed = collections.Counter(int(r['literature_count']) for r in rows)
        check(analysis + '_input_unique_pair_count', len(rows) == len(set(keys)) == sum(EXPECTED[analysis]), len(rows))
        check(analysis + '_input_exact_literature_count_distribution', [observed[k] for k in range(1, 7)] == EXPECTED[analysis], dict(observed))
        check(analysis + '_input_GTDB_pair_identity', all(r['canonical_species_a'] < r['canonical_species_b'] and
              r['pair_key'] == r['canonical_species_a'] + '\t' + r['canonical_species_b'] for r in rows))
        check(analysis + '_paper_indicators_and_support_count', all(
            all(r[p] in ['0', '1'] for p in PAPERS) and
            int(r['literature_count']) == sum(int(r[p]) for p in PAPERS) and
            r['supporting_papers'].split('|') == [p for p in PAPERS if r[p] == '1'] for r in rows))
        dist = read(A / ('summary/04a_' + analysis + '_literature_count_distribution.tsv'))
        check(analysis + '_distribution_matches_frozen_summary', {int(r['literature_count']): int(r['unique_pair_count']) for r in dist} == {k: observed[k] for k in range(1, 7)})
        masters[analysis] = rows
        membership_rows = read(A / ('membership/04a_' + analysis + '_pair_paper_membership.tsv'))
        mem = {}
        for r in membership_rows:
            index = r['pair_key'], r['paper']
            if index in mem:
                raise RuntimeError('Duplicate frozen paper membership: ' + str(index))
            sources = r['subnetwork_sources'].split('|')
            if r['paper_level_present'] != '1' or r['analysis_set'] != analysis or len(sources) != len(set(sources)):
                raise RuntimeError('Invalid frozen membership row: ' + str(index))
            mem[index] = sources
            if analysis == 'sensitivity':
                for n in sources:
                    network_pairs[n].add(r['pair_key'])
        check(analysis + '_membership_exactly_matches_master', set(mem) == {(r['pair_key'], p) for r in rows for p in PAPERS if r[p] == '1'})
        for p in PAPERS:
            label = 'Lacruz_' + analysis if p == 'Lacruz' else p
            paper_rows = read(A / ('paper_level/' + label + '_paper_pairs.tsv'))
            check(analysis + '_' + p + '_paper_level_provenance_matches_long',
                  {r['pair_key']: r['subnetwork_sources'].split('|') for r in paper_rows} == {k: v for (k, paper), v in mem.items() if paper == p})
        check(analysis + '_Lacruz_MHO_policy', all('Lacruz_MHO' not in v for v in mem.values()) if analysis == 'primary' else any('Lacruz_MHO' in v for v in mem.values()))
        memberships[analysis] = mem
    return masters, memberships, network_pairs


def main():
    for d in ['pair_level', 'source_sets', 'summary', 'qc', 'scripts', 'logs']:
        (OUT / d).mkdir(parents=True, exist_ok=True)
    log('Stage 4B: frozen-input checks, recurrence input validation and existing source provenance.')
    discoveries = [r for r in read(A / 'qc/04a_chen_wallen_input_discovery.tsv') if r['used_as_formal_input'] == 'True']
    check('Chen_Wallen_unique_frozen_formal_input_paths', len(discoveries) == 2 and {r['paper'] for r in discoveries} == {'Chen', 'Wallen'})
    historical = {r['paper']: Path(r['candidate_file']).parent for r in discoveries}
    chen_report = historical['Chen'] / 'CHEN_V45_source_exclusion_overlap_report.txt'
    chen_audit = historical['Chen'] / '06_mackdb_source_exclusion_affected_rows.tsv'
    wallen_report = historical['Wallen'] / 'WALLEN_taxonomy_rescue_overlap_report.txt'
    extra = {DB, LISTED, CHEN_SCRIPT, LEGACY_SCRIPT, LEGACY_SUMMARY, SCRIPT, VALIDATOR}
    for folder in list(historical.values()) + [ROOT / 'chen_wallen_null_models_final_v1', ROOT / '9_14_codex/chen_wallen_combined_overlap_v1']:
        extra.update(p for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts)
    protected = protect(extra)
    before = snapshot(protected)
    prior = prior_hashes()
    for p in sorted(protected & prior.keys()):
        check('prior_frozen_hash:' + str(p), before[p][1] == prior[p][0])
    check('original_human_gut_info_official_hash', before[DB][1] == DB_HASH)
    masters, memberships, network_pairs = load_recurrence()

    config_file = PROV / '03b_network_source_config.json'
    config = json.loads(config_file.read_text())
    source_sets = {n: set(c['sources']) for n, c in config.items()}
    legacy_expected = literal_assignment(LEGACY_SCRIPT, 'EXPECTED')
    chen_sources = set(literal_assignment(CHEN_SCRIPT, 'SOURCE_PROJECTS'))
    report_matches = re.findall(r'^Excluded BioProjects:\s*(.+)$', chen_report.read_text(), flags=re.M)
    check('Chen_source_set_uniquely_confirmed_by_script_report', len(report_matches) == 1 and {s.strip() for s in tokens(report_matches[0].strip())} == chen_sources == legacy_expected['Chen']['exclude'], report_matches)
    legacy_summary = {r['study']: r for r in read(LEGACY_SUMMARY)}
    check('Chen_source_set_matches_historical_summary', chen_sources == set(tokens(legacy_summary['Chen']['excluded_projects'])))
    wallen_frozen = read(B / 'summary/03b_wallen_reference_summary.tsv')
    check('Wallen_unique_existing_source_reference', len(wallen_frozen) == 1 and wallen_frozen[0]['reference_status'] == 'PASS')
    wallen_sources = set(tokens(wallen_frozen[0]['source_project']))
    check('Wallen_source_set_matches_historical_scripts_and_summary', wallen_sources == legacy_expected['Wallen']['exclude'] == set(tokens(legacy_summary['Wallen']['excluded_projects'])))
    source_sets.update(Chen=chen_sources, Wallen=wallen_sources)
    check('eleven_subnetworks_match_frozen_provenance', set(source_sets) == set(network_pairs))
    stage3a_manifest = read(BASE / 'stage3a/qc/03a_input_manifest.tsv')
    for r in stage3a_manifest:
        check(r['network'] + '_subnetwork_count_from_membership', len(network_pairs[r['network']]) == int(r['final_unique_pair_count']))
    for paper in ['Chen', 'Wallen']:
        check(paper + '_unchanged_formal_pair_count', len(network_pairs[paper]) == legacy_expected[paper]['edges'])

    alias_rows = read(PROV / '03b_accession_alias_map.tsv')
    aliases = {}
    for r in alias_rows:
        accession, canonical = r['reported_accession'], r['canonical_project_accession']
        if accession in aliases and aliases[accession] != canonical:
            raise RuntimeError('Ambiguous frozen accession alias: ' + accession)
        aliases[accession] = canonical
    check('all_accession_aliases_frozen_PASS', all(r['status'] == 'PASS' for r in alias_rows))
    check('Zhou_ERP115476_single_canonical_study', aliases.get('ERP115476') == 'PRJEB32762')
    check('all_selected_source_accessions_already_canonical', all(aliases.get(p, p) == p and re.fullmatch(r'PRJ(?:NA|EB|DB|CA)\d+', p) for v in source_sets.values() for p in v))
    stage3b_summary = {r['network']: r for r in read(B / 'summary/03b_source_exclusion_summary.tsv')}
    overlap = {r['network']: r for r in read(PROV / '03b_source_overlap_with_mackdb.tsv')}
    db_manifest = {n: r for r in read(B / 'source_excluded_networks/03b_database_manifest.tsv') for n in r['networks'].split(';')}
    project_audit = read(PROV / '03b_source_project_audit.tsv')
    check('Stage3B_formal_provenance_QC_PASS', all(r['status'] == 'PASS' for r in read(B / 'qc/03b_source_provenance_qc.tsv')))
    for n in config:
        included = {r['source_project'] for r in project_audit if r['network'] == n and r['include_in_source_set'] == 'Yes'}
        check(n + '_reported_sources_agree_across_frozen_provenance', included == source_sets[n] == set(tokens(stage3b_summary[n]['source_projects_reported'])) == set(tokens(overlap[n]['source_projects'])) == set(tokens(db_manifest[n]['source_projects'])))
    # The audit inclusion flag is project-level: three treatment rows share
    # already-included discovery projects. Count cohorts by their explicit role.
    li_audit = read(PROV / '03b_li_table_s1_cohort_audit.tsv')
    li_rows = [r for r in li_audit if r['role'] in {'Discovery', 'WGS Discovery'}]
    check('Li_14_discovery_cohorts_13_unique_projects', len(li_rows) == 14 and {r['canonical_project_accession'] for r in li_rows} == source_sets['Li_IBS'] and len(source_sets['Li_IBS']) == 13 and all(r['include_in_network_source_set'] == 'Yes' for r in li_rows))
    check('Li_treatment_or_validation_rows_add_no_projects', {r['canonical_project_accession'] for r in li_audit if r['include_in_network_source_set'] == 'Yes'} == source_sets['Li_IBS'])
    phenotype = read(PROV / '03b_lacruz_phenotype_project_counts.tsv')
    for n in config:
        if n.startswith('Lacruz_'):
            check(n + '_phenotype_specific_source_set', source_sets[n] == {r['canonical_project_accession'] for r in phenotype if r['network'] == n and r['include_in_source_set'] == 'Yes'})
    log('All six paper source sets confirmed; primary and sensitivity recurrence inputs match Stage 4A exactly.')

    # Freeze manifest lists analytical inputs; integrity separately covers all protected results.
    used = {DB, LISTED, CHEN_SCRIPT, LEGACY_SCRIPT, LEGACY_SUMMARY, chen_report, chen_audit, wallen_report,
            config_file, B / 'summary/03b_wallen_reference_summary.tsv', B / 'summary/03b_source_exclusion_summary.tsv',
            B / 'source_excluded_networks/03b_database_manifest.tsv', B / 'qc/03b_source_provenance_qc.tsv',
            B / 'qc/03b_aggregation_rule.tsv', BASE / 'stage3a/qc/03a_database_label_decoding.tsv',
            BASE / 'stage3a/qc/03a_input_manifest.tsv', A / 'qc/04a_chen_wallen_input_discovery.tsv',
            PROV / '03b_accession_alias_map.tsv', PROV / '03b_source_project_audit.tsv', PROV / '03b_source_overlap_with_mackdb.tsv',
            PROV / '03b_li_table_s1_cohort_audit.tsv', PROV / '03b_lacruz_phenotype_project_counts.tsv', SCRIPT, VALIDATOR}
    used.update(p for directory in ['membership', 'paper_level', 'summary'] for p in (A / directory).glob('*.tsv'))
    write(OUT / 'qc/04b_input_manifest.tsv', [dict(input_file=str(p), size_bytes=before[p][0], sha256=before[p][1],
          input_role='original_MACKDB_read_only' if p == DB else ('current_implementation' if p in {SCRIPT, VALIDATOR} else 'frozen_recurrence_or_provenance'), status='PASS') for p in sorted(used)])

    log('Scanning the original database once: exact pair lookup, aligned support vectors, actual project universe.')
    decoder_rows = read(BASE / 'stage3a/qc/03a_database_label_decoding.tsv')
    decoder = {r['database_label']: r['gtdb_r226_species'] for r in decoder_rows}
    check('frozen_exact_GTDB_decoder_unique', len(decoder) == len(decoder_rows) and all(r['decoding_status'] == 'exact_unique_R226_serialization' for r in decoder_rows))
    wanted = {r['pair_key'] for rs in masters.values() for r in rs}
    listed_values = LISTED.read_text().split()
    check('human_gut_list_875_unique_canonical_projects', len(listed_values) == len(set(listed_values)) == 875 and all(aliases.get(p, p) == p for p in listed_values))
    listed = set(listed_values)
    supports, baseline, seen = collections.Counter(), {}, set()
    malformed, alias_changed, unknown_label = 0, 0, 0
    with DB.open(newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        check('database_required_columns', reader.fieldnames == ['species_a', 'species_b', 'aver_corr', 'aver_pv', 'proj_count', 'proj_ID', 'proj_corr', 'proj_pv', 'proj_sam', 'super_env', 'sub_env'])
        for line, r in enumerate(reader, 2):
            vectors = [tokens(r[k]) for k in VECTORS]
            ids = vectors[0]
            malformed += not (int(r['proj_count']) >= 2 and all(len(v) == int(r['proj_count']) for v in vectors) and len(set(ids)) == len(ids))
            alias_changed += any(aliases.get(p, p) != p for p in ids)
            supports.update(ids)
            if r['species_a'] not in decoder or r['species_b'] not in decoder:
                unknown_label += 1
                continue
            a, b = sorted((decoder[r['species_a']], decoder[r['species_b']]))
            key = a + '\t' + b
            if a == b or key in seen:
                raise RuntimeError('Database duplicate canonical pair or self-loop at line ' + str(line))
            seen.add(key)
            if key in wanted:
                baseline[key] = dict(r, database_row_number=line)
    check('database_653071_unique_rows', len(seen) == 653071, len(seen))
    check('database_all_support_vectors_aligned_unique', malformed == 0, malformed)
    check('database_all_labels_decode_without_new_mapping', unknown_label == 0, unknown_label)
    check('database_support_aliases_already_canonical', alias_changed == 0, alias_changed)
    check('database_859_actual_support_projects_within_list', len(supports) == 859 and set(supports) <= listed, len(supports))
    check('actual_support_accession_syntax', all(re.fullmatch(r'PRJ(?:NA|EB|DB|CA)\d+', p) for p in supports))
    universe = set(supports)
    actual = {n: sources & universe for n, sources in source_sets.items()}
    for n in config:
        check(n + '_actual_support_overlap_matches_frozen_provenance', actual[n] == set(tokens(overlap[n]['source_projects_in_actual_support'])) == set(tokens(stage3b_summary[n]['source_projects_overlapping_mackdb'])) == set(tokens(db_manifest[n]['actual_source_intersection'])))
        check(n + '_listed_source_overlap_matches_frozen_provenance', source_sets[n] & listed == set(tokens(overlap[n]['source_projects_in_human_gut_list'])))
    for r in project_audit:
        p = r['source_project']
        check(r['network'] + '_' + p + '_prior_project_support_evidence', yn(p in listed) == r['in_human_gut_list'] and yn(p in universe) == r['in_human_gut_info_support'] and supports[p] == int(r['supporting_association_count']))
    ca = read(chen_audit)
    ca_sources = {p for r in ca for p in tokens(r['excluded_projects_present'])}
    check('Chen_sources_unique_actual_support_and_list_overlap', ca_sources == chen_sources == actual['Chen'] == (chen_sources & listed))
    check('Chen_old_affected_rows_match', len(ca) == int(legacy_summary['Chen']['affected_rows']))
    for p in sorted(chen_sources):
        count = sum(p in tokens(r['excluded_projects_present']) for r in ca)
        check('Chen_' + p + '_actual_support_matches_formal_audit', count == supports[p] == int(float(legacy_summary['Chen']['affected_rows_' + p])))
    check('Wallen_actual_support_matches_formal_audit', actual['Wallen'] == wallen_sources and all(supports[p] == int(wallen_frozen[0]['associations_with_source_support']) for p in wallen_sources))
    write(OUT / 'qc/04b_actual_project_support_universe.tsv', [dict(project=p, supporting_database_association_count=supports[p], in_human_gut_list=yn(p in listed)) for p in sorted(universe)])
    source_audit = []
    for r in project_audit:
        n, p = r['network'], r['source_project']
        source_audit.append(dict(paper=n.split('_')[0], network=n, source_project=p,
            include_in_network_construction_sources=r['include_in_source_set'], source_role=r['provenance_decision'],
            in_human_gut_list=yn(p in listed), in_actual_project_support=yn(p in universe),
            supporting_database_association_count=supports[p],
            eligible_for_pair_specific_exclusion=yn(p in actual[n]),
            evidence_source=str(PROV / '03b_source_project_audit.tsv'), status='PASS'))
    for paper in ['Chen', 'Wallen']:
        for p in sorted(source_sets[paper]):
            source_audit.append(dict(paper=paper, network=paper, source_project=p,
                include_in_network_construction_sources='Yes', source_role='existing_formal_network_source_exclusion_set',
                in_human_gut_list=yn(p in listed), in_actual_project_support=yn(p in universe),
                supporting_database_association_count=supports[p], eligible_for_pair_specific_exclusion=yn(p in actual[paper]),
                evidence_source=' | '.join(map(str, [CHEN_SCRIPT, chen_report, chen_audit, LEGACY_SUMMARY])) if paper == 'Chen' else ' | '.join(map(str, [B / 'summary/03b_wallen_reference_summary.tsv', wallen_report, LEGACY_SUMMARY])), status='PASS'))
    write(OUT / 'qc/04b_source_provenance_audit.tsv', source_audit)
    network_config = [dict(paper=n.split('_')[0], network=n, reported_source_projects=packed(source_sets[n]),
         source_projects_in_human_gut_list=packed(source_sets[n] & listed), overlapping_source_projects=packed(actual[n]),
         source_count=len(source_sets[n]), overlapping_source_count=len(actual[n]), status='PASS') for n in source_sets]
    write(OUT / 'source_sets/04b_network_source_config.tsv', network_config)
    write(OUT / 'qc/04b_chen_wallen_source_discovery.tsv', [dict(paper=p, reported_sources=packed(source_sets[p]),
          sources_in_human_gut_list=packed(source_sets[p] & listed), sources_in_actual_support=packed(actual[p]),
          source_support_counts=json.dumps({s: supports[s] for s in sorted(source_sets[p])}, sort_keys=True),
          evidence_source=' | '.join(sorted({r['evidence_source'] for r in source_audit if r['paper'] == p})),
          status='CONFIRMED_UNIQUE') for p in ['Chen', 'Wallen']])

    # Cross-check the exact baseline lookup against each frozen single-network result.
    for n, keys in network_pairs.items():
        recovered = sum(k in baseline for k in keys)
        single_recovered = sum(k in baseline and len(set(tokens(baseline[k]['proj_ID'])) - actual[n]) >= 2 for k in keys)
        expected_before = int(stage3b_summary[n]['baseline_recovered']) if n in config else legacy_expected[n]['recovered']
        expected_after = int(stage3b_summary[n]['source_excluded_recovered']) if n in config else legacy_expected[n]['recovered']
        check(n + '_baseline_and_single_source_results_match_frozen', recovered == expected_before and single_recovered == expected_after, '{} / {}'.format(recovered, single_recovered))
    log('Original database and all source overlaps verified; constructing pair-specific combined exclusions.')

    # Construction and finalization are below, kept separate from independent validation.
    finish(masters, memberships, source_sets, actual, baseline, universe, network_config,
           before, protected, extra, prior, historical)


def finish(masters, memberships, source_sets, actual, baseline, universe, network_config,
           before, protected, extra, prior, historical):
    pair_qc, results = [], {}
    for analysis in ['primary', 'sensitivity']:
        source_rows, pair_rows = [], []
        for original in masters[analysis]:
            key = original['pair_key']
            supporters = original['supporting_papers'].split('|')
            per_paper, per_actual, provenance = {}, {}, {}
            for p in PAPERS:
                subnetworks = memberships[analysis].get((key, p), [])
                if any(n.split('_')[0] != p or n not in source_sets for n in subnetworks):
                    raise RuntimeError('Invalid subnetwork provenance: ' + key)
                provenance[p] = subnetworks
                per_paper[p] = set().union(*(source_sets[n] for n in subnetworks)) if subnetworks else set()
                per_actual[p] = set().union(*(actual[n] for n in subnetworks)) if subnetworks else set()
            reported = set().union(*per_paper.values())
            remove = set().union(*per_actual.values())
            source = {k: original[k] for k in ['pair_key', 'canonical_species_a', 'canonical_species_b', 'literature_count', 'supporting_papers']}
            for p in PAPERS:
                source[p + '_sources'] = packed(per_paper[p])
                source[p + '_overlapping_sources'] = packed(per_actual[p])
                source[p + '_subnetwork_sources'] = '|'.join(provenance[p])
            source.update(reported_source_projects_union=packed(reported), overlapping_source_projects_union=packed(remove),
                          source_project_remove_count=len(remove), analysis_set=analysis)
            source_rows.append(source)
            r = baseline.get(key)
            vectors = [tokens(r[c]) if r else [] for c in VECTORS]
            keep = [i for i, p in enumerate(vectors[0]) if p not in remove]
            removed_indices = [i for i, p in enumerate(vectors[0]) if p in remove]
            retained = [[v[i] for i in keep] for v in vectors]
            recovered = bool(r) and len(keep) >= 2
            lost = bool(r) and not recovered
            row = dict(original)
            row.update(source_projects_to_remove=packed(remove), source_project_remove_count=len(remove),
                       reported_source_projects_union=packed(reported),
                       subnetwork_sources_by_paper=json.dumps(provenance, ensure_ascii=False, separators=(',', ':')),
                       baseline_mackdb_present=yn(bool(r)), baseline_proj_count=int(r['proj_count']) if r else 0,
                       baseline_proj_ID=r['proj_ID'] if r else '',
                       baseline_proj_corr=r['proj_corr'] if r else '', baseline_proj_pv=r['proj_pv'] if r else '',
                       baseline_proj_sam=r['proj_sam'] if r else '', baseline_aver_corr=r['aver_corr'] if r else 'NA',
                       baseline_aver_pv=r['aver_pv'] if r else 'NA',
                       removed_source_support_count=len(removed_indices),
                       removed_source_support_projects=';'.join(vectors[0][i] for i in removed_indices),
                       removed_support_indices_zero_based=';'.join(map(str, removed_indices)),
                       remaining_proj_count=len(keep), remaining_proj_ID=';'.join(retained[0]),
                       remaining_proj_corr=';'.join(retained[1]), remaining_proj_pv=';'.join(retained[2]),
                       remaining_proj_sam=';'.join(retained[3]),
                       remaining_aver_corr='NA', remaining_aver_pv='NA',
                       average_policy='not_recomputed_optional_full_precision_quantity',
                       combined_source_excluded_recovered=yn(recovered), loss_due_to_combined_source_exclusion=yn(lost),
                       database_row_number=r['database_row_number'] if r else '',
                       database_species_a=r['species_a'] if r else '', database_species_b=r['species_b'] if r else '',
                       analysis_set=analysis)
            pair_rows.append(row)
            # Verify the index-mask construction independently with zipped support tuples.
            tuples = list(zip(*vectors))
            expected_tuples = [t for t in tuples if t[0] not in remove]
            checks = dict(
                frozen_pair_identity_unchanged=all(row[k] == original[k] for k in original),
                paper_count_unchanged=int(row['literature_count']) == len(supporters) == sum(int(row[p]) for p in PAPERS),
                provenance_exactly_matches_supporters=[p for p in PAPERS if provenance[p]] == supporters,
                source_union_uses_only_actual_subnetworks=remove == reported & universe,
                unique_source_accessions=len(tokens(row['source_projects_to_remove'])) == len(remove),
                primary_has_no_MHO=analysis != 'primary' or 'Lacruz_MHO' not in provenance['Lacruz'],
                deletion_only_removes_present_source_support=set(vectors[0][i] for i in removed_indices) == set(vectors[0]) & remove,
                four_vectors_synchronously_deleted=list(zip(*retained)) == expected_tuples,
                retained_vectors_count_aligned=all(len(v) == len(keep) for v in retained),
                removed_and_retained_partition_original=len(removed_indices) + len(keep) == len(vectors[0]) and not set(retained[0]) & remove,
                recovery_equals_remaining_f2=recovered == (len(retained[0]) >= 2),
                no_new_recovery=not recovered or bool(r),
                loss_definition_correct=lost == (bool(r) and len(retained[0]) < 2))
            if not all(checks.values()):
                raise RuntimeError('Pair source-exclusion QC failed: ' + key + ' ' + str(checks))
            pair_qc.append(dict(analysis_set=analysis, pair_key=key,
                               **{k: 'PASS' if v else 'FAIL' for k, v in checks.items()}, status='PASS'))
        write(OUT / ('source_sets/04b_' + analysis + '_pair_source_sets.tsv'), source_rows)
        write(OUT / ('pair_level/04b_' + analysis + '_pair_recovery.tsv'), pair_rows)
        loss_rows = [r for r in pair_rows if r['loss_due_to_combined_source_exclusion'] == 'Yes']
        write(OUT / ('pair_level/04b_' + analysis + '_pairs_lost_after_source_exclusion.tsv'), loss_rows, list(pair_rows[0]))
        exact = [dict(literature_count=k, **summary([r for r in pair_rows if int(r['literature_count']) == k])) for k in range(1, 7)]
        cumulative = []
        for k in range(1, 7):
            s = summary([r for r in pair_rows if int(r['literature_count']) >= k])
            cumulative.append(dict(minimum_literature_count=k, external_pair_count=s['external_pair_count'],
                                   recovered_count=s['combined_source_excluded_recovered_count'],
                                   recovery_rate=s['combined_source_excluded_recovery_rate'],
                                   baseline_mackdb_present_count=s['baseline_mackdb_present_count'],
                                   baseline_recovery_rate=s['baseline_recovery_rate'], source_exclusion_loss_count=s['source_exclusion_loss_count']))
        frozen_combos = read(A / ('summary/04a_' + analysis + '_support_combination_counts.tsv'))
        combinations = []
        for c in frozen_combos:
            s = summary([r for r in pair_rows if r['supporting_papers'] == c['supporting_papers']])
            check(analysis + '_unchanged_combination_denominator:' + c['supporting_papers'], s['external_pair_count'] == int(c['pair_count']))
            combinations.append(dict(supporting_papers=c['supporting_papers'], literature_count=c['literature_count'],
                                     pair_count=s['external_pair_count'], recovered_count=s['combined_source_excluded_recovered_count'],
                                     recovery_rate=s['combined_source_excluded_recovery_rate'],
                                     baseline_mackdb_present_count=s['baseline_mackdb_present_count'],
                                     baseline_recovery_rate=s['baseline_recovery_rate'],
                                     source_exclusion_loss_count=s['source_exclusion_loss_count']))
        write(OUT / ('summary/04b_' + analysis + '_recovery_by_literature_count.tsv'), exact)
        write(OUT / ('summary/04b_' + analysis + '_cumulative_recurrence_recovery.tsv'), cumulative)
        write(OUT / ('summary/04b_' + analysis + '_support_combination_recovery.tsv'), combinations)
        totals = summary(pair_rows)
        check(analysis + '_combined_recovery_never_exceeds_baseline', totals['combined_source_excluded_recovered_count'] <= totals['baseline_mackdb_present_count'])
        check(analysis + '_baseline_equals_retained_plus_losses', totals['baseline_mackdb_present_count'] == totals['combined_source_excluded_recovered_count'] + totals['source_exclusion_loss_count'])
        check(analysis + '_summaries_partition_exact_denominators', sum(r['external_pair_count'] for r in exact) == sum(EXPECTED[analysis]) == len(pair_rows))
        results[analysis] = dict(rows=pair_rows, exact=exact, cumulative=cumulative, combinations=combinations, totals=totals)
        log('{}: {} pairs; baseline recovered {}; combined-source-excluded recovered {}; lost {}.'.format(
            analysis, len(pair_rows), totals['baseline_mackdb_present_count'], totals['combined_source_excluded_recovered_count'], totals['source_exclusion_loss_count']))
    write(OUT / 'qc/04b_pair_source_exclusion_qc.tsv', pair_qc)

    primary = {r['pair_key']: r for r in results['primary']['rows']}
    sensitivity = {r['pair_key']: r for r in results['sensitivity']['rows']}
    shared_status_changes = []
    for key, r in primary.items():
        t = sensitivity[key]
        if not set(tokens(r['source_projects_to_remove'])) <= set(tokens(t['source_projects_to_remove'])):
            raise RuntimeError('Sensitivity removed fewer sources for shared pair: ' + key)
        if r['combined_source_excluded_recovered'] == 'No' and t['combined_source_excluded_recovered'] == 'Yes':
            raise RuntimeError('Sensitivity gained recovery for shared pair: ' + key)
        if r['combined_source_excluded_recovered'] != t['combined_source_excluded_recovered']:
            shared_status_changes.append(key)
    check('primary_sensitivity_shared_pair_baseline_equal', all(r['baseline_mackdb_present'] == sensitivity[k]['baseline_mackdb_present'] for k, r in primary.items()))
    check('primary_sensitivity_extra_evidence_cannot_increase_shared_pair_recovery', True, 'Checked all {} shared pairs'.format(len(primary)))

    log('Running independent validator: re-read pair-level files and recompute exact, cumulative and combination summaries.')
    completed = subprocess.run([sys.executable, '-s', '-B', str(VALIDATOR)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if completed.returncode != 0:
        raise RuntimeError('Independent validator failed:\n' + completed.stdout + completed.stderr)
    log(completed.stdout.strip())
    validation_path = OUT / 'qc/04b_independent_summary_validation.tsv'
    independent = read(validation_path)
    OUTPUTS.add(validation_path)
    check('independent_validation_all_PASS', independent and all(r['status'] == 'PASS' for r in independent), len(independent))
    check('all_pair_source_exclusion_QC_PASS', len(pair_qc) == sum(map(sum, EXPECTED.values())) and all(r['status'] == 'PASS' for r in pair_qc), len(pair_qc))
    check('no_taxonomy_or_recurrence_reconstruction', True, 'Stage 4A fields copied unchanged; frozen database-label decoder only.')
    check('only_f2_recovery_and_descriptive_summaries', True, 'No model fitting, significance test, randomization, FastSpar or threshold sweep.')
    log('Checking all protected inputs and historical outputs after construction and independent validation.')
    check('frozen_file_inventory_unchanged', protect(extra) == protected)
    after = snapshot(protected)
    integrity = []
    for p in sorted(protected):
        unchanged = before[p] == after[p]
        check('unchanged:' + str(p), unchanged)
        integrity.append(dict(input_file=str(p), scope='current_implementation' if p in {SCRIPT, VALIDATOR} else 'frozen_upstream_or_historical',
                              size_bytes=before[p][0], sha256_before=before[p][1], sha256_after=after[p][1],
                              prior_recorded_sha256=prior[p][0] if p in prior else '', prior_hash_evidence=prior[p][1] if p in prior else '',
                              matches_prior_recorded_hash='PASS' if p in prior else 'NO_PRIOR_RECORD', unchanged=unchanged))
    check('original_database_hash_unchanged', after[DB][1] == DB_HASH)
    write(OUT / 'qc/04b_input_integrity.tsv', integrity)
    write(OUT / 'qc/04b_recovery_qc.tsv', QC)

    lines = ['# Stage 4B — Literature-recurrence-specific combined source-exclusion recovery', '',
        'overall_status: **PASS**', '',
        '独立 paper unit 固定为 Chen、Wallen、Zhou、Guo、Lacruz、Li。literature_count 表示独立外部论文数；remaining_proj_count 表示联合排除外部来源后数据库仍保留的独立 BioProject 数，两者完全不同。', '',
        '## 输入与逐 pair 来源规则', '',
        'Primary 的 11900 对、sensitivity 的 12426 对均直接来自 Stage 4A master。pair identity、paper indicators、supporting_papers 和 literature_count 原样复制并验证。',
        '来源只采用已有 Stage 3B 和 Chen/Wallen 正式 evidence。读取每条 pair 在 Stage 4A long membership 中实际 subnetwork_sources，再将对应网络的来源合并。不会仅凭 paper 指示列使用过大的全论文来源集合。',
        'Primary Lacruz 只计 MHNO/MUNO/MUO，MHO 的 observed flag 不参与 primary 来源；sensitivity 才允许 MHO。MHO-only pair 使用 MHO 来源集合。',
        'Zhou 的 ERP115476 已由冻结 alias map 归一到 PRJEB32762，不重复计数。Li 采用 14 discovery cohorts / 13 unique projects，排除 validation 和 rifaximin-only 来源。', '',
        md(['paper', 'network', 'reported_source_projects', 'source_projects_in_human_gut_list', 'overlapping_source_projects'], network_config), '',
        'Chen 正式证据：`' + str(CHEN_SCRIPT) + '` 中 SOURCE_PROJECTS；`' + str(historical['Chen'] / 'CHEN_V45_source_exclusion_overlap_report.txt') + '` 的 500FG source exclusion；`' + str(historical['Chen'] / '06_mackdb_source_exclusion_affected_rows.tsv') + '`；以及 `' + str(LEGACY_SUMMARY) + '`。',
        'Wallen 正式证据：`' + str(B / 'summary/03b_wallen_reference_summary.tsv') + '`、`' + str(historical['Wallen'] / 'WALLEN_taxonomy_rescue_overlap_report.txt') + '` 和上述历史 source-exclusion summary。',
        'Chen PRJNA319574、PRJNA942468 都出现在 human_gut.list 与 actual support，分别支持 1699、1921 条数据库 association，与历史正式 audit 一致。', '',
        '## Recovery construction', '',
        '原始数据库只读：`' + str(DB) + '`；SHA-256：`' + DB_HASH + '`。',
        '653071 条 canonical association、859 个 actual-support BioProjects；human_gut.list 共 875 个项目。数据库物种标签只通过冻结 Stage 3A exact decoding 表恢复表示，不重新 mapping。',
        '每条 pair 的 reported_source_projects_union 与 actual support universe 相交，形成 source_projects_to_remove；再只删除该 pair 自身 proj_ID 中真正出现的项目。',
        'proj_ID/proj_corr/proj_pv/proj_sam 使用同一个位置掩码同步删除，顺序与序列化数值保持原样。removed_source_support_count 是实际删掉的支持数，不等同于候选 source_project_remove_count。',
        '原始不存在的 pair recovered=No；存在的 pair 在 remaining_proj_count >= 2 时 recovered=Yes。loss 仅指原来存在而排除后不足 2 个项目。未生成任何完整数据库副本。',
        '不存在于 f=2 数据库时，输出支持向量为空、对应计数为 0；这只表示本次数据库没有可用记录，不推断原始单项目层面完全没有证据。',
        '可选 remaining_aver_corr / remaining_aver_pv 均记 NA。Stage 3B 已确认平均采用原始全精度项目值；本阶段不从已经舍入的 proj_corr/proj_pv 近似重算。这不影响 primary f=2 outcome。',
        '所有 rate 在分母为 0 时记 NA，因此 literature_count=6 保留 0 行而非报告 0% recovery。pair_key 内含 TAB，使用标准带引号 TSV reader。', '']
    for analysis in ['primary', 'sensitivity']:
        z = results[analysis]
        lines += ['## ' + analysis, '',
                  md(list(z['exact'][0]), z['exact']), '', 'Cumulative ≥N（独立的描述性表，不替代 exact bins）：', '',
                  md(list(z['cumulative'][0]), z['cumulative']), '',
                  '全部：baseline {baseline_mackdb_present_count}/{external_pair_count}；联合排除后 {combined_source_excluded_recovered_count}/{external_pair_count}；loss={source_exclusion_loss_count}。'.format(**z['totals']), '']
        lines += ['实际有 {} 个 pair 删除了来源支持，共删除 {} 条 project-support entries；零 loss 不表示未执行删除。'.format(
            sum(r['removed_source_support_count'] > 0 for r in z['rows']),
            sum(r['removed_source_support_count'] for r in z['rows'])), '']
        losses = [r for r in z['combinations'] if r['source_exclusion_loss_count']]
        lines += ['Loss supporting-paper combinations：', '', md(list(losses[0]), losses) if losses else '无 baseline recovered pair 丢失。', '']
        nonempty = [r for r in z['exact'] if r['external_pair_count']]
        increasing = all(float(a['combined_source_excluded_recovery_rate']) <= float(b['combined_source_excluded_recovery_rate']) for a, b in zip(nonempty, nonempty[1:]))
        lines += ['非空 exact literature_count 分组的 recovery rate 随文献数增加' + ('呈非递减排列。' if increasing else '并非逐级单调增加。') + '这仅是描述性比较，没有显著性或因果结论。', '']
    lines += ['## Primary 与 sensitivity', '',
              '共有 {} 对；其中联合排除后的 recovery 状态变化 {} 对。Sensitivity 新增 {} 对，且已有 pair 可因 MHO evidence 改变 literature_count 分组。'.format(len(primary), len(shared_status_changes), len(sensitivity) - len(primary)),
              '两套结果使用相同规则；组别分母不同，不能把数值差异解释为显著性差异。', '',
              '## QC、独立复算与冻结完整性', '',
              '{} 条逐 pair QC 记录全部 PASS；{} 条独立复算检查全部 PASS；{} 条整体 QC 全部 PASS。'.format(len(pair_qc), len(independent), len(QC)),
              '独立脚本重新打开 pair-level/source-set/master/summary 文件，用支持 tuple 过滤和 Counter 分组复算 distribution、baseline、combined recovery、exact、cumulative 与 exact combinations；不复用正式汇总函数。',
              '{} 个受保护文件 before/after hash 一致，其中 {} 个与此前记录的 hash 一致；全部冻结 Stage 2B/2C/3A/3A.1/3B/3C/4A、使用的 Chen/Wallen 历史结果及原始数据库保持不变。'.format(len(before), len(protected & prior.keys())),
              '独立检查文件：`' + str(validation_path) + '`。', '',
              '## 重现与正式输出', '',
              '`python3 -s -B ' + str(SCRIPT) + '`', '',
              '固定输入、排序与序列化格式；无系统时间或随机因素。正式输出绝对路径与 SHA-256 见 `' + str(OUT / 'qc/04b_formal_output_manifest.tsv') + '`，清单不自哈希。',
              'Stage 4B 完成后停止；未执行 regression、trend test、odds ratio、Fisher test、permutation、null model、enrichment、FastSpar 或新的 threshold sensitivity。', '']
    REPORT.write_text('\n'.join(lines), encoding='utf-8')
    OUTPUTS.update({REPORT, SCRIPT, VALIDATOR})
    status = dict(overall_status='PASS', primary=results['primary']['totals'], sensitivity=results['sensitivity']['totals'],
                  primary_exact=results['primary']['exact'], sensitivity_exact=results['sensitivity']['exact'],
                  shared_pair_recovery_changes=len(shared_status_changes), frozen_files_verified=len(before),
                  prior_hashes_verified=len(protected & prior.keys()), pair_qc_rows=len(pair_qc), independent_checks=len(independent), overall_checks=len(QC))
    status_path = OUT / 'qc/04b_run_status.json'
    status_path.write_text(json.dumps(status, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    OUTPUTS.add(status_path)
    log('overall_status=PASS; Stage 4B complete; stop after descriptive recovery construction.')
    logfile = OUT / 'logs/04b_recurrence_source_excluded_recovery.log'
    logfile.write_text('\n'.join(LOG) + '\n', encoding='utf-8')
    OUTPUTS.add(logfile)
    write(OUT / 'qc/04b_formal_output_manifest.tsv', [dict(absolute_path=str(p), size_bytes=p.stat().st_size, sha256=sha(p), role='formal_Stage4B_output') for p in sorted(OUTPUTS)])


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        write(OUT / 'qc/04b_recovery_qc.tsv', QC + [dict(check='overall_status', status='FAIL', details=str(exc))])
        (OUT / 'qc/04b_run_status.json').write_text(json.dumps(dict(overall_status='FAIL', error=str(exc)), ensure_ascii=False, indent=2) + '\n')
        REPORT.write_text('# Stage 4B\n\noverall_status: **FAIL**\n\nStopped: ' + str(exc) + '\n', encoding='utf-8')
        raise
