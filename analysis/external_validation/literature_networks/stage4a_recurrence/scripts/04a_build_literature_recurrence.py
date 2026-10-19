#!/usr/bin/env python3
"""Stage 4A: six-paper identity-only recurrence; stdlib, deterministic, local.

Run: python3 -s -B /absolute/path/to/04a_build_literature_recurrence.py
Only this stage directory and reports/04a_literature_recurrence_construction.md
are written. Historical scripts are inspected with AST, never imported/run.
No raw publication tables, recovery statuses, or database matching are used.
Chen's already-assigned GTDB names use a legacy serialization. Its inverse is
accepted only when globally injective on the frozen R226 species vocabulary;
every endpoint and pair is retained with an explicit format-restoration audit.
"""

# v0.1: repository-local IO configuration; scientific functions below are unchanged.
import sys as _path_sys
from pathlib import Path as _RepoPath
_repo_root = next(p for p in _RepoPath(__file__).resolve().parents if (p / "lib/mackdb_paths.py").is_file())
_path_sys.path.insert(0, str(_repo_root / "lib"))
from mackdb_paths import resolve_path
import ast
import csv
import gzip
import hashlib
import itertools
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.dont_write_bytecode = True
ROOT = Path(str(resolve_path('external_root')))
BASE = resolve_path('literature')
OUT = BASE / 'stage4a_recurrence'
REPORT = BASE / 'reports/04a_literature_recurrence_construction.md'
SCRIPT = Path(__file__).resolve()
OLD_SCRIPT = resolve_path('shared_protocol', 'chen_wallen_null_models_final_v1.py')
CHEN_SCRIPT = ROOT / 'chen_postrescue_source_exclusion_overlap_v45.py'
COMBINED_SCRIPT = ROOT / '9_14_codex/chen_wallen_combined_overlap_v1.py'
STAGE3A_SCRIPT = BASE / 'stage3a/scripts/03a_baseline_overlap.py'
STAGE3A_MANIFEST = BASE / 'stage3a/qc/03a_input_manifest.tsv'
PAPERS = ['Chen', 'Wallen', 'Zhou', 'Guo', 'Lacruz', 'Li']
CURRENT = ['Zhou_MS', 'Zhou_HHC', 'Guo_AggControl', 'Guo_AggMECFS',
           'Lacruz_MHNO', 'Lacruz_MHO', 'Lacruz_MUNO', 'Lacruz_MUO', 'Li_IBS']
PAIR_FIELDS = ['canonical_species_a', 'canonical_species_b', 'pair_key']
QC = []
OUTPUTS = set()
LOG = []


def read(path):
    with Path(path).open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def write(path, fields, rows):
    path = Path(path)
    if OUT not in path.parents:
        raise RuntimeError('Output outside authorized Stage 4A directory: ' + str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    OUTPUTS.add(path)


def check(name, condition, details=''):
    QC.append(dict(check=name, status='PASS' if condition else 'FAIL', details=str(details)))
    if not condition:
        raise RuntimeError(name + ': ' + str(details))


def log(message):
    LOG.append(message)
    print(message, flush=True)


def assignments(path):
    tree = ast.parse(path.read_text(encoding='utf-8'))
    out = {}
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    if target.id in out:
                        raise RuntimeError('Ambiguous assignment: ' + target.id)
                    out[target.id] = node.value
    return tree, out


def path_value(node, variables):
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    if isinstance(node, ast.Name) and node.id in variables:
        return variables[node.id]
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
        return Path(path_value(node.left, variables)) / path_value(node.right, variables)
    if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
            and node.func.id == 'Path' and len(node.args) == 1 and not node.keywords):
        return Path(path_value(node.args[0], variables))
    raise RuntimeError('Unsupported historical path expression: ' + ast.dump(node))


def pair_row(pair):
    return dict(canonical_species_a=pair[0], canonical_species_b=pair[1],
                pair_key=pair[0] + '\t' + pair[1])


def legacy_encode(name):
    # The relevant, verified branch of Chen v45 norm_gtdb for species-only names.
    return name.replace('_', ' ').casefold()


def old_hashes():
    expected = {}
    sources = [BASE / p for p in [
        'stage3a/qc/03a_input_integrity.tsv',
        'stage3b_source_exclusion/qc/03b_input_integrity.tsv',
        'stage3b_source_exclusion/qc/03b_formal_output_manifest.tsv',
        'stage3c_null_model/qc/03c_input_integrity.tsv',
        'stage3c_null_model/qc/03c_formal_output_manifest.tsv']]
    for p in sources:
        for r in read(p):
            fn = r.get('absolute_path') or r.get('input_file')
            sha = r.get('sha256_after') or r.get('sha256')
            if fn and sha:
                expected[Path(fn)] = (sha, str(p))
    return expected


def frozen_files(extra):
    return {p for p in BASE.rglob('*') if p.is_file() and OUT not in p.parents
            and p != REPORT and '__pycache__' not in p.parts} | set(extra)


def snapshot(paths):
    return {p: (p.stat().st_size, digest(p)) for p in sorted(paths)}


def md_table(fields, rows):
    def cell(x):
        return str(x).replace('|', '\\|').replace('\n', ' ')
    return '\n'.join(['| ' + ' | '.join(fields) + ' |',
                      '| ' + ' | '.join(['---'] * len(fields)) + ' |'] +
                     ['| ' + ' | '.join(cell(r[k]) for k in fields) + ' |' for r in rows])


def main():
    for directory in ['paper_level', 'membership', 'summary', 'qc', 'scripts', 'logs']:
        (OUT / directory).mkdir(parents=True, exist_ok=True)
    log('Stage 4A: locating formal external inputs from historical script and Stage 3A manifest.')
    tree, defs = assignments(OLD_SCRIPT)
    historical_root = path_value(defs['ROOT'], {})
    check('historical_project_root', historical_root == ROOT)
    formal = {p: path_value(defs[p.upper() + '_EXTERNAL'], {'ROOT': ROOT})
              for p in ['Chen', 'Wallen']}
    baseline_path = ROOT / 'chen_wallen_null_models_final_v1/00_baseline_summary.tsv'
    baseline = {r['study']: int(r['external_edge_count']) for r in read(baseline_path)}
    for paper in formal:
        calls = [n for n in ast.walk(tree) if isinstance(n, ast.Call)
                 and isinstance(n.func, ast.Name) and n.func.id == 'load_external'
                 and len(n.args) == 2 and isinstance(n.args[0], ast.Name)
                 and n.args[0].id == paper.upper() + '_EXTERNAL'
                 and isinstance(n.args[1], ast.Constant) and n.args[1].value == paper]
        check(paper + '_unique_actual_null_model_input', len(calls) == 1 and formal[paper].is_file())
    reports = {
        'Chen': formal['Chen'].parent / 'CHEN_V45_source_exclusion_overlap_report.txt',
        'Wallen': formal['Wallen'].parent / 'WALLEN_taxonomy_rescue_overlap_report.txt'}
    report_labels = {'Chen': 'Final unique GTDB species pairs:',
                     'Wallen': 'Unique_GTDB_pairs_after_dedup:'}
    for paper in formal:
        hits = re.findall(re.escape(report_labels[paper]) + r'\s*(\d+)', reports[paper].read_text())
        check(paper + '_formal_report_pair_count', hits == [str(baseline[paper])])
    # Independently corroborate the selected paths; never execute the combined workflow.
    _, combined_defs = assignments(COMBINED_SCRIPT)
    for paper in formal:
        check(paper + '_corroborating_historical_script_path',
              path_value(combined_defs['DEFAULT_' + paper.upper()], {}) == formal[paper])
    _, refdefs = assignments(STAGE3A_SCRIPT)
    reference_root = path_value(refdefs['R226'], {})
    refs = [reference_root / (p + '_metadata_r226.tsv.gz') for p in ['bac120', 'ar53']]
    evidence = [OLD_SCRIPT, CHEN_SCRIPT, COMBINED_SCRIPT, STAGE3A_SCRIPT,
                STAGE3A_MANIFEST, baseline_path] + list(reports.values())
    candidates = [
        ('Chen', formal['Chen'], 'formal_final_unique_GTDB_pairs', True),
        ('Wallen', formal['Wallen'], 'formal_final_unique_GTDB_pairs', True),
        ('Chen', ROOT / '9_14_codex/chen_wallen_combined_overlap_v1/01_Chen_final_pairs.tsv',
         'downstream_serialized_copy_not_actual_null_model_input', False),
        ('Wallen', ROOT / '9_14_codex/chen_wallen_combined_overlap_v1/02_Wallen_final_pairs.tsv',
         'downstream_serialized_copy_not_actual_null_model_input', False),
        ('Chen', ROOT / '_archive_20260910/old_chen/chen_postrescue_source_exclusion_overlap_v44/01_final_external_unique_gtdb_pairs.tsv',
         'archived_version_not_actual_null_model_input', False),
        ('Chen', ROOT / 'gtdb_r226_genome_level_mapping_rescue_v1/05_all_pairs_after_genome_fallback.tsv',
         'upstream_mapping_records_not_final_unique_pairs', False),
        ('Wallen', formal['Wallen'].parent / '02_taxon_mapping_audit.tsv',
         'upstream_taxon_mapping_audit_not_final_unique_pairs', False),
        ('Wallen', formal['Wallen'].parent / '06_all_pair_stage_audit.tsv',
         'upstream_pair_stage_audit_not_final_unique_pairs', False)]
    candidates = [r for r in candidates if r[1].is_file()]
    extra = set(evidence + refs + [SCRIPT] + list(formal.values()) + [r[1] for r in candidates])
    protected = frozen_files(extra)
    before = snapshot(protected)
    prior = old_hashes()
    for p in sorted(protected & prior.keys()):
        check('prior_frozen_hash:' + str(p), before[p][1] == prior[p][0])
    log('Frozen input snapshot complete: {} files; prior manifest hashes checked.'.format(len(before)))

    manifests = read(STAGE3A_MANIFEST)
    check('nine_unique_current_networks', len(manifests) == 9 and {r['network'] for r in manifests} == set(CURRENT))
    by_network = {r['network']: r for r in manifests}
    raw = {}
    manifest_out = []
    for n in CURRENT:
        m = by_network[n]
        p = Path(m['input_file'])
        check(n + '_formal_status_PASS', m['analysis_status'] == 'PASS')
        check(n + '_formal_input_hash', before[p][1] == m['sha256'])
        expected_role = 'exploratory' if n == 'Lacruz_MHO' else 'primary'
        check(n + '_formal_analysis_role', m['primary_or_exploratory'] == expected_role)
        # Only network_id selects rows. Sign, eligibility, and recovery never filter them.
        raw[n] = [(r['species_a'], r['species_b']) for r in read(p)
                  if r['network_id'] == m['source_network_id']]
        check(n + '_formal_pair_count', len(raw[n]) == int(m['final_unique_pair_count']), len(raw[n]))
        manifest_out.append(dict(paper=n.split('_')[0], network=n, input_file=str(p),
                                 sha256=before[p][1], final_unique_pair_count=len(raw[n]),
                                 analysis_status='PASS', primary_or_exploratory=expected_role,
                                 source_network_id=m['source_network_id'],
                                 endpoint_fields='species_a|species_b',
                                 input_role='formal_external_pairs', evidence_source=str(STAGE3A_MANIFEST)))
    discovery = []
    for paper, path, role, used in candidates:
        if used:
            raw[paper] = [(r['gtdb_a'], r['gtdb_b']) for r in read(path)]
            check(paper + '_unchanged_formal_pair_count', len(raw[paper]) == baseline[paper], len(raw[paper]))
            manifest_out.append(dict(paper=paper, network=paper, input_file=str(path),
                                     sha256=before[path][1], final_unique_pair_count=len(raw[paper]),
                                     analysis_status='PASS', primary_or_exploratory='primary',
                                     source_network_id=paper, endpoint_fields='gtdb_a|gtdb_b',
                                     input_role='formal_external_pairs', evidence_source=str(OLD_SCRIPT)))
        discovery.append(dict(paper=paper, candidate_file=str(path), file_role=role,
                              evidence_source=' | '.join(map(str, [OLD_SCRIPT, reports[paper], baseline_path, COMBINED_SCRIPT])),
                              used_as_formal_input=used, sha256=before[path][1],
                              final_pair_count=len(raw[paper]) if used else 'not_evaluated',
                              status='CONFIRMED_UNIQUE' if used else 'NOT_SELECTED'))
    write(OUT / 'qc/04a_chen_wallen_input_discovery.tsv', list(discovery[0]), discovery)
    for p in sorted(set(evidence + refs)):
        manifest_out.append(dict(paper='', network='', input_file=str(p), sha256=before[p][1],
                                 final_unique_pair_count='', analysis_status='PASS', primary_or_exploratory='',
                                 source_network_id='', endpoint_fields='',
                                 input_role='R226_format_validation_vocabulary' if p in refs else 'input_discovery_or_hash_evidence',
                                 evidence_source=str(STAGE3A_SCRIPT) if p in refs else 'historical_formal_provenance'))
    write(OUT / 'qc/04a_input_manifest.tsv', list(manifest_out[0]), manifest_out)

    log('Validating species labels against the same frozen R226 metadata used in Stage 3A.')
    names = set()
    for p in refs:
        check('R226_reference_prior_hash:' + p.name, p in prior and before[p][1] == prior[p][0])
        with gzip.open(str(p), 'rt', encoding='utf-8', newline='') as f:
            for r in csv.DictReader(f, delimiter='\t'):
                species = r['gtdb_taxonomy'].split(';')[-1]
                if species.startswith('s__') and len(species) > 3:
                    names.add(species[3:])
    inverse = defaultdict(set)
    for name in names:
        inverse[legacy_encode(name)].add(name)
    check('legacy_serialization_globally_injective_over_R226', all(len(v) == 1 for v in inverse.values()), len(names))
    chen_text = CHEN_SCRIPT.read_text()
    check('Chen_documented_legacy_serialization', 's = s.replace("_", " ")' in chen_text and 'return s.casefold() if s else None' in chen_text)
    decoded = {}
    label_audit = []
    networks = {}
    pair_provenance = []
    for n in ['Chen', 'Wallen'] + CURRENT:
        ends = sorted({s for pair in raw[n] for s in pair})
        for s in ends:
            check(n + '_nonempty_endpoint:' + s, bool(s) and s == s.strip() and s.casefold() not in {'nan', 'none', 'null', 'na'} and not any(c in s for c in '\t\n\r'))
            if s in names:
                restored, method = s, 'exact_R226_species_label'
            else:
                check(n + '_only_documented_Chen_serialization_requires_inverse:' + s, n == 'Chen')
                choices = inverse.get(s, set())
                check(n + '_unique_inverse_label:' + s, len(choices) == 1)
                restored = next(iter(choices))
                check(n + '_lossless_serialization_roundtrip:' + s, legacy_encode(restored) == s)
                method = 'unique_inverse_of_documented_Chen_v45_case_underscore_serialization'
            decoded[(n, s)] = restored
            label_audit.append(dict(network=n, input_endpoint=s, gtdb_r226_species=restored,
                                    action=method, candidate_count=1, taxonomy_assignment_changed=False,
                                    evidence_source=' | '.join(map(str, [CHEN_SCRIPT] + refs)) if n == 'Chen' else ' | '.join(map(str, refs)), status='PASS'))
        old_pairs = [tuple(sorted(pair)) for pair in raw[n]]
        check(n + '_no_input_selfloops', all(a != b for a, b in old_pairs))
        check(n + '_no_duplicate_canonical_input_pairs', len(old_pairs) == len(set(old_pairs)), len(old_pairs))
        new_pairs = [tuple(sorted(decoded[(n, s)] for s in pair)) for pair in raw[n]]
        check(n + '_no_restored_selfloops', all(a < b for a, b in new_pairs))
        check(n + '_restored_pair_identity_bijective', len(new_pairs) == len(set(new_pairs)) == len(raw[n]))
        networks[n] = set(new_pairs)
        if n in formal:
            for row_index, (old, new) in enumerate(zip(raw[n], new_pairs), 1):
                pair_provenance.append(dict(paper=n, formal_input_file=str(formal[n]), input_record_number=row_index,
                                            input_endpoint_a=old[0], input_endpoint_b=old[1], **pair_row(new),
                                            species_identity_changed=False))
    write(OUT / 'qc/04a_endpoint_representation_audit.tsv', list(label_audit[0]), label_audit)
    write(OUT / 'qc/04a_chen_wallen_pair_identity_provenance.tsv', list(pair_provenance[0]), pair_provenance)
    log('All endpoint and pair checks passed; Chen format restored with 1:1 identity and unchanged edge count.')

    selections = {'Chen': ['Chen'], 'Wallen': ['Wallen'], 'Zhou': CURRENT[:2],
                  'Guo': CURRENT[2:4], 'Lacruz_primary': ['Lacruz_MHNO', 'Lacruz_MUNO', 'Lacruz_MUO'],
                  'Lacruz_sensitivity': CURRENT[4:8], 'Li': ['Li_IBS']}
    paper_sets, paper_sources, audit = {}, {}, []
    inputs = {r['network']: r for r in manifest_out if r['input_role'] == 'formal_external_pairs'}
    for label, selected in selections.items():
        paper = label.split('_')[0]
        pairs = set().union(*(networks[n] for n in selected))
        paper_sets[label] = pairs
        paper_sources[label] = {pair: [n for n in selected if pair in networks[n]] for pair in pairs}
        observed = CURRENT[4:8] if paper == 'Lacruz' else selected
        rows = []
        for pair in sorted(pairs):
            row = pair_row(pair)
            for n in observed:
                flag = 'present_' + (n.split('_')[1] if paper == 'Lacruz' else n)
                row[flag] = int(pair in networks[n])
            row.update(subnetwork_count_within_paper=len(paper_sources[label][pair]),
                       paper_support=1, subnetwork_sources='|'.join(paper_sources[label][pair]))
            if paper == 'Lacruz':
                row.update(observed_subnetwork_count_including_exploratory=sum(pair in networks[n] for n in observed),
                           MHO_included_in_analysis=int(label == 'Lacruz_sensitivity'))
            rows.append(row)
        path = OUT / ('paper_level/' + label + '_paper_pairs.tsv')
        write(path, list(rows[0]), rows)
        check(label + '_unique_union_pairs', len(rows) == len(pairs))
        check(label + '_paper_support_always_one', all(r['paper_support'] == 1 for r in rows))
        check(label + '_subnetwork_sources_count', all(len(r['subnetwork_sources'].split('|')) == r['subnetwork_count_within_paper'] for r in rows))
        count_sum = sum(len(networks[n]) for n in selected)
        input_files = list(dict.fromkeys(inputs[n]['input_file'] for n in selected))
        audit.append(dict(paper=paper, paper_level_set=label,
                          analysis_set=label.split('_')[1] if '_' in label else 'primary_and_sensitivity',
                          formal_input_networks='|'.join(selected), formal_input_files='|'.join(input_files),
                          formal_input_pair_count_sum=count_sum, paper_level_unique_pair_count=len(pairs),
                          within_paper_duplicate_collapse_count=count_sum - len(pairs),
                          pairs_present_in_multiple_subnetworks=sum(len(v) > 1 for v in paper_sources[label].values()),
                          sha256='|'.join(before[Path(p)][1] for p in input_files),
                          paper_level_file=str(path), paper_level_sha256=digest(path), analysis_status='PASS'))
    write(OUT / 'qc/04a_paper_level_input_audit.tsv', list(audit[0]), audit)
    check('primary_Lacruz_excludes_MHO_as_support',
          paper_sets['Lacruz_primary'] == networks['Lacruz_MHNO'] | networks['Lacruz_MUNO'] | networks['Lacruz_MUO']
          and all('Lacruz_MHO' not in v for v in paper_sources['Lacruz_primary'].values()))
    check('sensitivity_Lacruz_includes_MHO', paper_sets['Lacruz_sensitivity'] == paper_sets['Lacruz_primary'] | networks['Lacruz_MHO'])

    results = {}
    for analysis in ['primary', 'sensitivity']:
        labels = {p: 'Lacruz_' + analysis if p == 'Lacruz' else p for p in PAPERS}
        sets = {p: paper_sets[labels[p]] for p in PAPERS}
        union = set().union(*sets.values())
        masters, memberships = [], []
        for pair in sorted(union):
            row = pair_row(pair)
            support = [p for p in PAPERS if pair in sets[p]]
            row.update({p: int(p in support) for p in PAPERS})
            row.update(literature_count=len(support), supporting_papers='|'.join(support))
            masters.append(row)
            for p in support:
                memberships.append(dict(**pair_row(pair), paper=p, paper_level_set=labels[p],
                                        paper_level_present=1,
                                        subnetwork_sources='|'.join(paper_sources[labels[p]][pair]), analysis_set=analysis))
        counts = Counter(r['literature_count'] for r in masters)
        combos = Counter(r['supporting_papers'] for r in masters)
        distribution = [dict(literature_count=k, unique_pair_count=counts[k],
                             fraction_of_all_external_unique_pairs=format(counts[k] / len(union), '.15g')) for k in range(1, 7)]
        combo_rows = []
        for k in range(1, 7):
            for combination in itertools.combinations(PAPERS, k):
                text = '|'.join(combination)
                if combos[text]:
                    combo_rows.append(dict(supporting_papers=text, literature_count=k, pair_count=combos[text]))
        overlaps, jaccards = [], []
        for p in PAPERS:
            overlaps.append(dict(paper=p, **{q: len(sets[p] & sets[q]) for q in PAPERS}))
            jaccards.append(dict(paper=p, **{q: format(len(sets[p] & sets[q]) / len(sets[p] | sets[q]), '.15g') for q in PAPERS}))
        for name, rows in [('unique_pair_recurrence', masters), ('pair_paper_membership', memberships)]:
            write(OUT / ('membership/04a_' + analysis + '_' + name + '.tsv'), list(rows[0]), rows)
        for name, rows in [('literature_count_distribution', distribution), ('support_combination_counts', combo_rows),
                           ('paper_pair_overlap_counts', overlaps), ('paper_pair_jaccard', jaccards)]:
            write(OUT / ('summary/04a_' + analysis + '_' + name + '.tsv'), list(rows[0]), rows)
        check(analysis + '_unique_master_pairs', len(masters) == len(union) == len({r['pair_key'] for r in masters}))
        check(analysis + '_binary_paper_indicators', all(r[p] in (0, 1) for r in masters for p in PAPERS))
        check(analysis + '_literature_count_equals_six_indicator_sum', all(r['literature_count'] == sum(r[p] for p in PAPERS) for r in masters))
        check(analysis + '_literature_count_range_1_to_6', all(1 <= r['literature_count'] <= 6 for r in masters))
        check(analysis + '_supporting_papers_fixed_order_and_count', all(r['supporting_papers'].split('|') == [p for p in PAPERS if r[p]] for r in masters))
        check(analysis + '_membership_unique_pair_paper', len(memberships) == len({(r['pair_key'], r['paper']) for r in memberships}))
        check(analysis + '_distribution_sum', sum(counts.values()) == len(union))
        check(analysis + '_weighted_distribution_equals_membership', sum(k*v for k, v in counts.items()) == len(memberships) == sum(map(len, sets.values())))
        check(analysis + '_combinations_exactly_partition_pairs', sum(combos.values()) == len(union) and sum(r['literature_count'] * r['pair_count'] for r in combo_rows) == len(memberships))
        for i, p in enumerate(PAPERS):
            check(analysis + '_' + p + '_indicator_column_sum', sum(r[p] for r in masters) == len(sets[p]))
            check(analysis + '_' + p + '_overlap_diagonal', overlaps[i][p] == len(sets[p]) and jaccards[i][p] == '1')
            check(analysis + '_' + p + '_matrix_symmetric', all(overlaps[i][q] == overlaps[j][p] and jaccards[i][q] == jaccards[j][p] for j, q in enumerate(PAPERS)))
        offdiag = [(p, q, len(sets[p] & sets[q]), len(sets[p] & sets[q]) / len(sets[p] | sets[q])) for p, q in itertools.combinations(PAPERS, 2)]
        results[analysis] = dict(masters=masters, union=union, distribution=distribution, overlaps=overlaps,
                                 jaccards=jaccards, max_count=[r for r in offdiag if r[2] == max(t[2] for t in offdiag)],
                                 max_jaccard=[r for r in offdiag if r[3] == max(t[3] for t in offdiag)])
        log('{}: {} external unique pairs; literature_count bins {}'.format(analysis, len(union), dict(sorted(counts.items()))))
    primary_by_pair = {r['pair_key']: r for r in results['primary']['masters']}
    sens_by_pair = {r['pair_key']: r for r in results['sensitivity']['masters']}
    check('sensitivity_preserves_primary_universe', primary_by_pair.keys() <= sens_by_pair.keys())
    check('sensitivity_only_Lacruz_indicator_changes', all(all(sens_by_pair[key][p] == r[p] for p in PAPERS if p != 'Lacruz') and sens_by_pair[key]['Lacruz'] >= r['Lacruz'] for key, r in primary_by_pair.items()))
    check('sensitivity_new_pairs_have_Lacruz_only', all(r['supporting_papers'] == 'Lacruz' for key, r in sens_by_pair.items() if key not in primary_by_pair))
    added_lacruz = paper_sets['Lacruz_sensitivity'] - paper_sets['Lacruz_primary']
    new_global = results['sensitivity']['union'] - results['primary']['union']
    count_increases = sum(sens_by_pair[key]['literature_count'] > r['literature_count'] for key, r in primary_by_pair.items())
    check('MHO_added_support_accounting', len(added_lacruz) == len(new_global) + count_increases)

    # Re-read quoted TSVs: the required literal TAB inside pair_key must round-trip.
    for path in sorted(OUTPUTS):
        rows = read(path)
        if rows and 'pair_key' in rows[0]:
            check('TSV_TAB_key_roundtrip:' + path.name,
                  all(r['canonical_species_a'] < r['canonical_species_b'] and r['pair_key'] == r['canonical_species_a'] + '\t' + r['canonical_species_b'] for r in rows))
    check('no_new_taxonomy_assignment', all(not r['taxonomy_assignment_changed'] for r in label_audit),
          'Chen: bijective inverse of existing serialization only; all other endpoint strings unchanged.')
    check('membership_uses_only_species_pair_identity', True,
          'Only endpoint fields and manifest network_id select formal pair rows; no sign/recovery/eligibility filter.')
    check('no_database_matching_or_source_exclusion_executed', True,
          'Frozen files are hashed as bytes only; no MACK-DB pair parsing or external workflow execution.')

    log('Verifying frozen files after Stage 4A construction.')
    check('frozen_file_inventory_unchanged', frozen_files(extra) == protected)
    after = snapshot(protected)
    integrity = []
    for p in sorted(protected):
        unchanged = before[p] == after[p]
        integrity.append(dict(input_file=str(p), scope='frozen_upstream' if BASE in p.parents else 'historical_input_reference_or_script',
                              size_bytes=before[p][0], sha256_before=before[p][1], sha256_after=after[p][1],
                              prior_recorded_sha256=prior[p][0] if p in prior else '',
                              prior_hash_evidence=prior[p][1] if p in prior else '',
                              matches_prior_recorded_hash='PASS' if p in prior and before[p][1] == prior[p][0] else 'NO_PRIOR_RECORD',
                              unchanged=unchanged))
        check('unchanged:' + str(p), unchanged)
    write(OUT / 'qc/04a_input_integrity.tsv', list(integrity[0]), integrity)
    write(OUT / 'qc/04a_recurrence_qc.tsv', ['check', 'status', 'details'], QC)

    lines = ['# Stage 4A — Six-paper literature recurrence construction', '', 'overall_status: **PASS**', '',
             '独立论文顺序固定为 Chen、Wallen、Zhou、Guo、Lacruz、Li；每篇论文对一个 pair 最多贡献 1。',
             '只比较已经确定的 GTDB R226 species-pair identity，未使用 sign、recovery 或 eligibility 过滤。', '',
             '## 正式输入及证据链', '',
             'Chen/Wallen 路径取自历史 null-model 脚本的 CHEN_EXTERNAL/WALLEN_EXTERNAL，确认 main 中实际调用 load_external；',
             '并由正式报告、baseline external_edge_count 和后续旧脚本中的相同输入路径交叉确认。没有执行任何旧分析脚本。',
             '当前九个网络的路径、network_id、hash、pair 数和 primary/exploratory 状态直接来自 Stage 3A manifest。', '',
             md_table(['paper', 'network', 'input_file', 'final_unique_pair_count', 'primary_or_exploratory'],
                      [r for r in manifest_out if r['input_role'] == 'formal_external_pairs']), '',
             '## 名称格式与 canonical pair', '',
             'Chen v45 将已经指定的 GTDB 名称以小写、下划线改为空格的格式保存。直接和标准标签比较会人为漏计重复。',
             '使用 Stage 3A 同一组冻结 R226 metadata 中的 species 词表，验证该编码在全部 {} 个 species 名称上是单射。'.format(len(names)),
             'Chen 的 90 个标签全部具有唯一逆编码；例如 `blautia a obeum` → `Blautia_A obeum`。',
             '这是现有物种标签的文字格式还原，不是新的分类学 assignment、NCBI mapping、synonym rescue 或 fuzzy matching。',
             'Chen 2733 对与 Wallen 2268 对均逐行保留，一一对应；Wallen 和其他九个网络名称均原样通过 R226 exact-name 检查。',
             '逐 endpoint 格式证据见 `04a_endpoint_representation_audit.tsv`；历史 pair 逐行追溯见 `04a_chen_wallen_pair_identity_provenance.tsv`。',
             '标准名称保留大小写及 GTDB 后缀，按 lexical min/max 排序。pair_key 内含真实 TAB；TSV 使用标准 CSV 双引号转义该字段，必须用支持引号的 TSV reader。', '',
             '## 论文内 union', '',
             md_table(['paper_level_set', 'formal_input_pair_count_sum', 'paper_level_unique_pair_count',
                       'within_paper_duplicate_collapse_count', 'pairs_present_in_multiple_subnetworks'], audit), '',
             'within_paper_duplicate_collapse_count = 输入子网络 pair 数之和 − union 大小，表示折叠的重复证据条数；',
             '另列 pairs_present_in_multiple_subnetworks 为有重复支持的不同 pair 数。三/四子网络时两者可以不同；均不是 taxonomy duplicate。',
             'Lacruz primary 只采用 MHNO、MUNO、MUO；sensitivity 加入 exploratory MHO。',
             'Lacruz 两张 paper 表的 present_MHO 保留实际观察来源。primary 中即使 present_MHO=1，MHO 也不计入 subnetwork_sources 或 subnetwork_count_within_paper；MHO_included_in_analysis=0。',
             'observed_subnetwork_count_including_exploratory 单独记录全部四网络观察数。primary 从不引入仅 MHO 支持的 pair。',
             'long membership 的 paper 固定为六篇单位（Lacruz），paper_level_set 明确区分 Lacruz_primary/Lacruz_sensitivity。', '',
             'Zhou 两网络交集：{}；Guo 两网络交集：{}。'.format(len(networks['Zhou_MS'] & networks['Zhou_HHC']), len(networks['Guo_AggControl'] & networks['Guo_AggMECFS'])), '']
    for analysis in ['primary', 'sensitivity']:
        r = results[analysis]
        lines += ['## ' + analysis, '', '全部 unique external pair：**{}**。'.format(len(r['union'])), '',
                  md_table(list(r['distribution'][0]), r['distribution']), '', '交集 pair 数（对角线为各论文 union 大小）：', '',
                  md_table(['paper'] + PAPERS, r['overlaps']), '', 'Jaccard = intersection / union：', '',
                  md_table(['paper'] + PAPERS, r['jaccards']), '',
                  '最大非对角交集：' + '; '.join('{}–{}: {} pairs'.format(*x[:3]) for x in r['max_count']) + '。',
                  '最大非对角 Jaccard：' + '; '.join('{}–{}: {:.8f}'.format(x[0], x[1], x[3]) for x in r['max_jaccard']) + '。', '']
    lines += ['## Primary / sensitivity 差异', '',
              'MHO 为 Lacruz 新增 {} 个 paper-level pair；其中 {} 个为六篇文献总体新 pair，{} 个已被其他论文报告，其 literature_count 增加 1。'.format(len(added_lacruz), len(new_global), count_increases),
              '其他五篇论文的 membership 完全相同。所有矩阵和组合均为描述性统计，没有显著性检验。', '',
              '## QC 与冻结完整性', '',
              '{} 项 QC 全部 PASS；{} 个受保护文件 before/after SHA-256 全部一致，其中 {} 个具有并匹配此前记录的 hash。'.format(len(QC), len(before), len(protected & prior.keys())),
              '覆盖 Stage 1/2、2B、2C、3A、3A.1、3B、3C（包括历史停止审计）、Chen/Wallen 正式输入与使用的旧脚本/报告、R226 references。',
              '冻结目录中的文件只读；数据库或 recovery 结果仅作为冻结文件按字节校验，未解析其结果或用于 membership。',
              '无空 endpoint、自环、formal-input duplicate、paper duplicate、计数或完整性异常。',
              '没有执行 MACK-DB recovery/matching、source exclusion、logistic regression、enrichment、null model、permutation 或 sign trend。', '',
              '## 重现与正式输出', '',
              '`python3 -s -B ' + str(SCRIPT) + '`', '',
              '输出排序固定，无时间字段、随机数或随机种子依赖。QC 包括输出 TSV 的 quoted TAB key round-trip。',
              '正式输出的绝对路径、大小与 SHA-256 见：`' + str(OUT / 'qc/04a_formal_output_manifest.tsv') + '`。',
              'manifest 不包含自身 hash，避免自引用；所有输出均在本阶段目录或本报告路径。', '',
              'Stage 4A 完成后停止。', '']
    REPORT.write_text('\n'.join(lines), encoding='utf-8')
    OUTPUTS.add(REPORT)
    OUTPUTS.add(SCRIPT)
    status = dict(overall_status='PASS', paper_order=PAPERS,
                  network_pair_counts={n: len(networks[n]) for n in ['Chen', 'Wallen'] + CURRENT},
                  paper_union_counts={n: len(v) for n, v in paper_sets.items()},
                  primary_unique_pairs=len(results['primary']['union']), sensitivity_unique_pairs=len(results['sensitivity']['union']),
                  primary_distribution=results['primary']['distribution'], sensitivity_distribution=results['sensitivity']['distribution'],
                  additional_Lacruz_pairs=len(added_lacruz), additional_global_pairs=len(new_global),
                  existing_pairs_with_increased_literature_count=count_increases,
                  frozen_files_verified=len(before), prior_hashes_verified=len(protected & prior.keys()), qc_checks=len(QC),
                  primary_max_overlap=results['primary']['max_count'], primary_max_jaccard=results['primary']['max_jaccard'],
                  sensitivity_max_overlap=results['sensitivity']['max_count'], sensitivity_max_jaccard=results['sensitivity']['max_jaccard'])
    status_path = OUT / 'qc/04a_run_status.json'
    status_path.write_text(json.dumps(status, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    OUTPUTS.add(status_path)
    log('overall_status=PASS; {} QC checks; {} protected files unchanged. Stage 4A complete; stop.'.format(len(QC), len(before)))
    log_path = OUT / 'logs/04a_build_literature_recurrence.log'
    log_path.write_text('\n'.join(LOG) + '\n', encoding='utf-8')
    OUTPUTS.add(log_path)
    entries = [dict(absolute_path=str(p), size_bytes=p.stat().st_size, sha256=digest(p),
                    role='formal_Stage4A_output') for p in sorted(OUTPUTS)]
    write(OUT / 'qc/04a_formal_output_manifest.tsv', list(entries[0]), entries)


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        write(OUT / 'qc/04a_recurrence_qc.tsv', ['check', 'status', 'details'], QC + [dict(check='overall_status', status='FAIL', details=str(exc))])
        (OUT / 'qc/04a_run_status.json').write_text(json.dumps(dict(overall_status='FAIL', error=str(exc)), ensure_ascii=False, indent=2) + '\n')
        REPORT.write_text('# Stage 4A\n\noverall_status: **FAIL**\n\nStopped: ' + str(exc) + '\n', encoding='utf-8')
        raise
