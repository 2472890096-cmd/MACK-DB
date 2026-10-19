#!/usr/bin/env python3
"""Offline Stage 3A exact-pair baseline overlap; no taxonomy/source rescue.

Run: python3 -s stage3a/scripts/03a_baseline_overlap.py
All earlier-stage files and the full human_gut.info are read-only. Database
numeric-ID/underscore serialization is inverted by exact lookup in the R226
species-name vocabulary, retaining genus and species suffixes and case. No
NCBI mapping, synonyms, fuzzy matches, support exclusions or reaggregation.
"""

# v0.1: repository-local IO configuration; scientific functions below are unchanged.
import sys as _path_sys
from pathlib import Path as _RepoPath
_repo_root = next(p for p in _RepoPath(__file__).resolve().parents if (p / "lib/mackdb_paths.py").is_file())
_path_sys.path.insert(0, str(_repo_root / "lib"))
from mackdb_paths import resolve_path
import csv
import gzip
import hashlib
import json
import logging
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.dont_write_bytecode = True
BASE = Path(str(resolve_path('literature')))
OUT = BASE / 'stage3a'
REPORT = BASE / 'reports/03a_baseline_overlap_qc.md'
DB = Path(str(resolve_path('environment_network')))
R226 = Path(str(resolve_path('gtdb_r226')))
REFS = [R226 / (name + '_metadata_r226.tsv.gz') for name in ['bac120', 'ar53']]
PAIR2B = BASE / 'stage2b/04_all_external_edges_gtdb_r226_stage2b.tsv'
AUDIT2B = BASE / 'stage2b/04_edge_mapping_audit_stage2b.tsv'
LI = BASE / 'stage2c/li2024'
PAIR2C = LI / 'taxonomy/li2024_gtdb_r226_pairs.tsv'
NETWORKS = [('Zhou_MS', 'Zhou2022_MS'), ('Zhou_HHC', 'Zhou2022_HHC'),
            ('Guo_AggControl', 'Guo2023_AggControl'), ('Guo_AggMECFS', 'Guo2023_AggMECFS'),
            ('Lacruz_MHNO', 'Lacruz2026_MHNO'), ('Lacruz_MHO', 'Lacruz2026_MHO'),
            ('Lacruz_MUNO', 'Lacruz2026_MUNO'), ('Lacruz_MUO', 'Lacruz2026_MUO'),
            ('Li_IBS', 'Li2024_published')]
NA = 'NA'
QC, INVENTORY, MANIFEST, SUMMARY = [], [], [], []
STATE = {}
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')
LOG = logging.getLogger('stage3a')


def j(x):
    return json.dumps(x, ensure_ascii=False, sort_keys=True, separators=(',', ':'))


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    with path.open(encoding='utf-8', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        return reader.fieldnames, list(reader)


def output_path(path):
    assert OUT in path.resolve().parents, 'Write outside Stage 3A: ' + str(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def write(path, rows, fields=None):
    output_path(path)
    fields = fields or list(dict.fromkeys(k for r in rows for k in r))
    assert fields, path
    encoded = [{k: NA if r.get(k) is None else str(r[k]) for k in fields} for r in rows]
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(encoded)
    assert read(path) == (fields, encoded), 'TSV round-trip: ' + str(path)


def check(name, condition, detail='', network='ALL'):
    QC.append(dict(network=network, check=name, status='PASS' if condition else 'FAIL', detail=detail))
    if not condition:
        raise AssertionError(name + ': ' + str(detail))


def canonical(a, b):
    # Exact, case-sensitive lexical sorting. No name normalization is applied.
    return (a, b) if a <= b else (b, a)


def pair_text(pair):
    assert all(' || ' not in x for x in pair)
    return pair[0] + ' || ' + pair[1]


def sign(x):
    value = float(x)
    assert math.isfinite(value), 'Nonfinite association value'
    return 'positive' if value > 0 else 'negative' if value < 0 else 'zero'


def frozen_paths():
    return {p for name in ['parsed', 'mapping', 'scripts', 'stage2b', 'stage2c', 'raw', 'reference', 'reports']
            for p in (BASE / name).rglob('*') if p.is_file() and p != REPORT}


def snapshot():
    previous = read(LI / 'qc/li2024_input_integrity.tsv')[1]
    check('stage2c_recorded_integrity_pass', all(r['unchanged'] == 'True' and r['matches_prior_manifest'] == 'True' for r in previous))
    expected = {Path(r['absolute_path']): r['sha256_after'] for r in previous}
    paths = frozen_paths() | set(expected) | set(REFS) | {DB}
    for path in sorted(paths):
        digest = sha(path)
        INVENTORY.append(dict(absolute_path=str(path), input_scope='human_gut_database' if path == DB else 'frozen_stage1_2_2b_2c_or_reference',
                              size_bytes=path.stat().st_size, sha256_before=digest,
                              prior_recorded_sha256=expected.get(path, NA),
                              matches_prior_recorded_hash=path not in expected or expected[path] == digest,
                              sha256_after=NA, unchanged='PENDING'))
    write(OUT / 'qc/03a_input_integrity.tsv', INVENTORY)
    check('frozen_inputs_match_previous_hashes', all(r['matches_prior_recorded_hash'] for r in INVENTORY))
    STATE['hashes'] = {r['absolute_path']: r['sha256_before'] for r in INVENTORY}
    STATE['frozen_paths_before'] = {str(p) for p in frozen_paths()}
    LOG.info('Snapshot: %d frozen files/reference inputs plus database', len(INVENTORY) - 1)


def load_external():
    report2b = (BASE / 'reports/02b_taxonomy_rescue_qc.md').read_text(encoding='utf-8')
    report2c = (BASE / 'reports/02c_li2024_taxonomy_qc.md').read_text(encoding='utf-8')
    check('stage2b_report_overall_PASS', bool(re.search(r'overall_status\s*=\s*\*\*PASS\*\*', report2b)))
    check('stage2c_report_overall_PASS', bool(re.search(r'overall_status\s*=\s*PASS', report2c)))
    qc2c = read(LI / 'qc/li2024_qc_summary.tsv')[1]
    check('stage2c_machine_QC_PASS', all(r['status'] == 'PASS' for r in qc2c) and
          any(r['check'] == 'overall_status' and r['status'] == 'PASS' for r in qc2c))
    qc2b = {r['network_id']: r for r in read(BASE / 'reports/02b_taxonomy_rescue_qc.tsv')[1]}
    gain2b = {r['network']: r for r in read(BASE / 'stage2b/stage2_vs_stage2b_network_summary.tsv')[1]}
    check('stage2b_machine_QC_eight_networks_balanced', len(qc2b) == 8 and all(r['edge_fate_balance_ok'] == 'True' for r in qc2b.values()))
    rows2b = read(PAIR2B)[1]
    rows2c = read(PAIR2C)[1]
    raw2b = read(AUDIT2B)[1]
    raw_by_id = {r['raw_record_id']: r for r in raw2b}
    check('frozen_stage2b_audit_ids_unique', len(raw_by_id) == len(raw2b))
    raw_counts = Counter(r['network_id'] for r in raw2b)
    published_li = read(LI / 'processed/li2024_published_edges_raw.tsv')[1]
    li_by_id = {r['source_row_id']: r for r in published_li}
    li_fates = {r['primary_fate']: int(r['count']) for r in read(LI / 'qc/li2024_edge_fate_summary.tsv')[1]}
    check('Li_final_pair_denominator_599', len(rows2c) == li_fates['kept_unique_gtdb_pair'] == 599)
    check('Li_published_raw_count_matches_frozen_QC', len(published_li) == sum(li_fates.values()) == 742 and
          any(r['check'] == 'published_edges_exact_742' and r['detail'] == '742' for r in qc2c))
    check('Li_report_pair_count_matches_file', bool(re.search(r'final unique pairs = {}\b'.format(len(rows2c)), report2c)))
    observed_networks = {r['network_id'] for r in rows2b + rows2c}
    check('exactly_nine_expected_networks', observed_networks == {old for _, old in NETWORKS})
    external = {}
    for network, old_network in NETWORKS:
        li = network == 'Li_IBS'
        path = PAIR2C if li else PAIR2B
        rows = [r for r in (rows2c if li else rows2b) if r['network_id'] == old_network]
        expected_count = li_fates['kept_unique_gtdb_pair'] if li else int(qc2b[old_network]['final_unique_gtdb_pairs'])
        raw_count = len(published_li) if li else int(qc2b[old_network]['raw_edges'])
        label = 'exploratory' if network == 'Lacruz_MHO' else 'primary'
        check('final_pair_counts_match_latest_taxonomy_QC', len(rows) == expected_count and
              (li or expected_count == int(gain2b[old_network]['Stage2B_final_pairs'])), str(len(rows)), network)
        check('published_raw_counts_match_frozen_audit', li or raw_counts[old_network] == raw_count, raw_count, network)
        check('primary_exploratory_unchanged', {r['primary_status'] for r in rows} == {label} and
              (li or qc2b[old_network]['primary_status'] == gain2b[old_network]['primary_status'] == label), label, network)
        check('frozen_final_pairs_unique_canonical_valid', len({canonical(r['species_a'], r['species_b']) for r in rows}) == len(rows) and
              all(r['species_a'] not in ('', NA) and r['species_b'] not in ('', NA) and r['species_a'] < r['species_b'] and
                  pair_text(canonical(r['species_a'], r['species_b'])) == r['pair_key'] for r in rows), network=network)
        items = []
        for r in rows:
            member_ids = json.loads(r['raw_record_ids'])
            members = [li_by_id[x] if li else raw_by_id[x] for x in member_ids]
            assert len(member_ids) == int(r['raw_record_count']) and len(set(member_ids)) == len(member_ids)
            assert member_ids[0] == r['representative_record_id']
            values = [x['cor'] if li else x['association_value'] for x in members]
            signs = {sign(x) for x in values}
            # All current final pairs have a consistent signed association.
            # Never choose one sign to override a frozen collapsed conflict.
            assert (r['sign_status'] == 'consistent' and len(signs) == 1 and r['sign'] == next(iter(signs))) or (
                    r['sign_status'] == 'conflict' and len(signs) > 1 and r['sign'] == NA)
            if not li:
                assert all(x['network_id'] == old_network and x['pair_key'] == r['pair_key'] for x in members)
                assert members[0]['primary_fate'] == 'kept_unique_gtdb_pair'
            else:
                assert r['cor'] == values[0] and r['external_sign'] == r['sign']
            association_type = 'SparCC_correlation' if li else members[0]['association_type']
            items.append(dict(network=network, source_network_id=old_network, primary_or_exploratory=label,
                              external_species_a=r['species_a'], external_species_b=r['species_b'], external_pair_key=r['pair_key'],
                              external_cor=values[0] if r['sign_status'] == 'consistent' else NA, external_sign=r['sign'],
                              external_association_type=association_type, external_sign_status=r['sign_status'],
                              external_cor_policy='frozen first-record association value; no weight aggregation',
                              external_representative_cor=values[0], external_member_values=j(values),
                              external_raw_record_count=len(member_ids), external_representative_record_id=member_ids[0],
                              external_raw_record_ids=r['raw_record_ids'], external_input_file=str(path)))
        check('external_signs_match_frozen_member_values', True, 'Representative value; collapsed sign preserved', network)
        check('frozen_sign_conflict_counts_preserved', li or sum(r['sign_status'] == 'conflict' for r in rows) == int(qc2b[old_network]['sign_conflict_pairs']), network=network)
        external[network] = items
        MANIFEST.append(dict(network=network, source_network_id=old_network, input_file=str(path), sha256=STATE['hashes'][str(path)],
                             final_unique_pair_count=len(rows), analysis_status='PASS', primary_or_exploratory=label,
                             published_raw_edges=raw_count, taxonomy_pair_coverage=len(rows) / raw_count,
                             taxonomy_QC_report=str(BASE / 'reports' / ('02c_li2024_taxonomy_qc.md' if li else '02b_taxonomy_rescue_qc.md')),
                             taxonomy_QC_overall_status='PASS'))
    write(OUT / 'qc/03a_input_manifest.tsv', MANIFEST)
    LOG.info('External manifest written BEFORE overlap: %s', {r['network']: r['final_unique_pair_count'] for r in MANIFEST})
    return external


def vocabulary():
    names = defaultdict(set)
    sources = defaultdict(set)
    for path in REFS:
        with gzip.open(path, 'rt', encoding='utf-8') as f:
            columns = next(f).rstrip('\n').split('\t')
            index = columns.index('gtdb_taxonomy')
            for line in f:
                label = line.rstrip('\n').split('\t')[index].split(';')[-1]
                assert label.startswith('s__')
                species = label[3:]
                if species:
                    payload = species.replace(' ', '_')
                    names[payload].add(species)
                    sources[payload].add(str(path))
    check('R226_serialization_is_injective', all(len(v) == 1 for v in names.values()), len(names))
    return {k: next(iter(v)) for k, v in names.items()}, sources


def audit_database():
    names, sources = vocabulary()
    decoded, label_rows, seen, id_to_name = {}, [], {}, {}
    counts, missing, duplicate_examples = Counter(), Counter(), []
    pc_min, pc_max = None, None
    canonical_path = OUT / 'database/human_gut_canonical_keys.tsv'
    key_fields = ['source_row_id', 'database_species_a_raw', 'database_species_b_raw',
                  'decoded_endpoint1', 'decoded_endpoint2', 'canonical_species_a', 'canonical_species_b', 'canonical_pair_key']
    with DB.open(encoding='utf-8', newline='') as f, output_path(canonical_path).open('w', encoding='utf-8', newline='') as out:
        reader = csv.DictReader(f, delimiter='\t')
        fields = reader.fieldnames
        check('database_required_header', {'species_a', 'species_b', 'aver_corr', 'aver_pv', 'proj_count', 'proj_ID'}.issubset(fields), j(fields))
        writer = csv.DictWriter(out, fieldnames=key_fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for source_row, r in enumerate(reader, 2):
            assert None not in r and all(v is not None for v in r.values()), 'Malformed database row'
            counts['data_records'] += 1
            missing.update(k for k in fields if r[k].strip().lower() in {'', 'na', 'nan', 'null', 'none'})
            endpoints = []
            for field in ['species_a', 'species_b']:
                raw = r[field]
                if raw not in decoded:
                    match = re.fullmatch(r'(\d+)_(.+)', raw)
                    value = names.get(match[2]) if match else None
                    decoded[raw] = value
                    if value is not None:
                        assert match[2] == value.replace(' ', '_'), 'Serialization round-trip mismatch'
                        assert match[1] not in id_to_name or id_to_name[match[1]] == value, 'Database ID has conflicting labels'
                        id_to_name[match[1]] = value
                    label_rows.append(dict(database_label=raw, database_numeric_id=match[1] if match else NA,
                                           encoded_payload=match[2] if match else NA, gtdb_r226_species=value or NA,
                                           decoding_status='exact_unique_R226_serialization' if value else 'FAIL',
                                           evidence_file=';'.join(sorted(sources.get(match[2], []))) if match else NA,
                                           rule='strip numeric ID prefix; exact inverse of R226 species.replace(space, underscore); preserve suffixes/case'))
                endpoints.append(decoded[raw])
            if None in endpoints:
                counts['undecodable_endpoint_records'] += 1
                if any(r[k].strip().lower() in {'', 'na', 'nan', 'null', 'none'} for k in ['species_a', 'species_b']):
                    counts['missing_endpoint_records'] += 1
                continue
            pair = canonical(*endpoints)
            counts['self_loop_records'] += pair[0] == pair[1]
            counts['rows_reordered_for_canonical_species'] += endpoints[0] > endpoints[1]
            if pair in seen:
                counts['duplicate_canonical_records_excess'] += 1
                if len(duplicate_examples) < 20:
                    duplicate_examples.append(dict(first_row=seen[pair], repeated_row=source_row, pair_key=pair_text(pair)))
            else:
                seen[pair] = source_row
            counts['database_' + sign(r['aver_corr']) + '_records'] += 1
            assert math.isfinite(float(r['aver_pv']))
            project_count = int(r['proj_count'])
            pc_min = project_count if pc_min is None else min(pc_min, project_count)
            pc_max = project_count if pc_max is None else max(pc_max, project_count)
            project_fields = [k for k in fields if k.startswith('proj_') and k != 'proj_count']
            counts['project_support_vector_length_mismatch_rows'] += any(len(r[k].split(';')) != project_count for k in project_fields)
            counts['repeated_project_ID_rows'] += len(set(r['proj_ID'].split(';'))) != len(r['proj_ID'].split(';'))
            writer.writerow(dict(source_row_id=source_row, database_species_a_raw=r['species_a'], database_species_b_raw=r['species_b'],
                                 decoded_endpoint1=endpoints[0], decoded_endpoint2=endpoints[1], canonical_species_a=pair[0],
                                 canonical_species_b=pair[1], canonical_pair_key=pair_text(pair)))
    with DB.open('rb') as f:
        total_lines = sum(1 for _ in f)
    audit = dict(input_file=str(DB), sha256=STATE['hashes'][str(DB)], size_bytes=DB.stat().st_size,
                 header=j(fields), total_lines_including_header=total_lines, data_records=counts['data_records'],
                 species_pair_fields='species_a;species_b', correlation_field='aver_corr', pvalue_field='aver_pv',
                 project_count_field='proj_count', project_ID_field='proj_ID',
                 other_project_support_fields=j([k for k in fields if k.startswith('proj_') and k not in ['proj_count', 'proj_ID']]),
                 other_fields=j([k for k in fields if k not in ['species_a', 'species_b', 'aver_corr', 'aver_pv'] and not k.startswith('proj_')]),
                 unique_endpoint_labels=len(decoded), unique_decoded_species=len({v for v in decoded.values() if v is not None}),
                 undecodable_endpoint_labels=sum(v is None for v in decoded.values()), canonical_unique_pairs=len(seen),
                 pair_representation='numeric_ID + underscore-serialized GTDB R226 species; unordered lexical species key derived locally',
                 missing_cells_by_field=j({k: missing[k] for k in fields}), duplicate_examples=j(duplicate_examples),
                 project_count_min=pc_min, project_count_max=pc_max,
                 derived_canonical_file=str(canonical_path), all_project_support_values_preserved=True,
                 source_exclusion_performed=False, database_reaggregation_performed=False)
    for k in ['self_loop_records', 'duplicate_canonical_records_excess', 'missing_endpoint_records',
              'undecodable_endpoint_records', 'rows_reordered_for_canonical_species', 'database_positive_records',
              'database_negative_records', 'database_zero_records', 'project_support_vector_length_mismatch_rows', 'repeated_project_ID_rows']:
        audit[k] = counts[k]
    STATE['database_audit'] = audit
    write(OUT / 'qc/03a_database_input_audit.tsv', [dict(metric=k, value=v) for k, v in audit.items()])
    write(OUT / 'qc/03a_database_label_decoding.tsv', sorted(label_rows, key=lambda r: r['database_label']))
    check('database_duplicate_canonical_pairs_absent_STOP_gate', counts['duplicate_canonical_records_excess'] == 0, j(duplicate_examples))
    check('database_endpoints_present_and_exactly_decoded', counts['missing_endpoint_records'] == counts['undecodable_endpoint_records'] == 0)
    check('database_self_loops_absent', counts['self_loop_records'] == 0)
    check('database_physical_lines_reconcile', total_lines == counts['data_records'] + 1)
    check('database_project_support_vectors_audited', True, j({k: audit[k] for k in ['project_support_vector_length_mismatch_rows', 'repeated_project_ID_rows', 'project_count_min', 'project_count_max']}))
    STATE['database_fields'] = fields
    LOG.info('Database audit passed: %d unique pairs; %d labels; %d rows reordered', len(seen), len(decoded), counts['rows_reordered_for_canonical_species'])
    return decoded, seen


def overlap(external, decoded, database_pairs):
    wanted = {(r['external_species_a'], r['external_species_b']) for rows in external.values() for r in rows}
    matching = {}
    # Audited database is now read as-is. No support/correlation/p-value filter.
    with DB.open(encoding='utf-8', newline='') as f:
        for source_row, raw in enumerate(csv.DictReader(f, delimiter='\t'), 2):
            pair = canonical(decoded[raw['species_a']], decoded[raw['species_b']])
            if pair in wanted:
                assert pair not in matching, 'Duplicate database pair; never collapse'
                matching[pair] = (source_row, raw)
    check('exact_lookup_matches_full_database_membership', set(matching) == wanted.intersection(database_pairs))
    detail_rows = {}
    for manifest in MANIFEST:
        network = manifest['network']
        details = []
        for source in external[network]:
            pair = canonical(source['external_species_a'], source['external_species_b'])
            result = dict(source, recovered='Yes' if pair in matching else 'No', mack_source_row_id=NA)
            result.update({'mack_' + k: NA for k in STATE['database_fields']})
            result.update(mack_sign=NA, sign_concordant=NA)
            if result['recovered'] == 'Yes':
                source_row, raw = matching[pair]
                result.update({'mack_' + k: v for k, v in raw.items()})
                result['mack_source_row_id'] = source_row
                result['mack_sign'] = sign(raw['aver_corr'])
                if result['external_sign'] != NA:
                    result['sign_concordant'] = 'Yes' if result['external_sign'] == result['mack_sign'] else 'No'
            details.append(result)
        write(OUT / 'baseline_overlap' / (network + '_baseline_overlap.tsv'), details)
        detail_rows[network] = details
        recovered = [r for r in details if r['recovered'] == 'Yes']
        n, r = len(details), len(recovered)
        concordant = sum(x['sign_concordant'] == 'Yes' for x in recovered)
        discordant = sum(x['sign_concordant'] == 'No' for x in recovered)
        row = dict(network=network, source_network_id=manifest['source_network_id'], primary_or_exploratory=manifest['primary_or_exploratory'],
                   published_raw_edges=manifest['published_raw_edges'], taxonomy_mapped_unique_pairs=n,
                   taxonomy_pair_coverage=n / manifest['published_raw_edges'], baseline_recovered_pairs=r,
                   baseline_not_recovered_pairs=n-r, baseline_recovery_rate=r/n,
                   primary_recovery_denominator='taxonomy_mapped_unique_pairs',
                   external_positive_pairs=sum(x['external_sign'] == 'positive' for x in details),
                   external_negative_pairs=sum(x['external_sign'] == 'negative' for x in details),
                   external_zero_pairs=sum(x['external_sign'] == 'zero' for x in details),
                   external_unresolved_sign_pairs=sum(x['external_sign'] == NA for x in details),
                   recovered_external_positive=sum(x['external_sign'] == 'positive' for x in recovered),
                   recovered_external_negative=sum(x['external_sign'] == 'negative' for x in recovered),
                   recovered_external_zero=sum(x['external_sign'] == 'zero' for x in recovered),
                   sign_concordant_pairs=concordant, sign_discordant_pairs=discordant,
                   recovered_sign_unassessable_pairs=r-concordant-discordant,
                   sign_concordance_rate_among_recovered=concordant/r if r else NA,
                   sign_concordance_role='descriptive_only; not a recovery criterion',
                   end_to_end_recovery_rate=r/manifest['published_raw_edges'],
                   end_to_end_metric_role='secondary descriptive metric', edge_pair_balance_ok=r+(n-r)==n)
        SUMMARY.append(row)
        check('recovered_plus_not_recovered_equals_final_pairs', row['edge_pair_balance_ok'] and n == manifest['final_unique_pair_count'], str(n), network)
        check('every_external_pair_once', len({x['external_pair_key'] for x in details}) == n, network=network)
        check('every_recovery_exact_and_every_nonrecovery_absent', all((x['recovered'] == 'Yes') ==
              ((x['external_species_a'], x['external_species_b']) in database_pairs) for x in details), network=network)
        check('sign_does_not_filter_recovery', all(x['recovered'] == 'Yes' for x in details if x['sign_concordant'] == 'No') and
              {x['external_pair_key'] for x in recovered} == {pair_text(pair) for pair in
              {(x['external_species_a'], x['external_species_b']) for x in details}.intersection(database_pairs)},
              'Discordant recovered pairs retained: ' + str(discordant), network)
        check('sign_counts_reconcile', row['external_positive_pairs']+row['external_negative_pairs']+row['external_zero_pairs']+
              row['external_unresolved_sign_pairs'] == n and concordant+discordant+row['recovered_sign_unassessable_pairs'] == r, network=network)
        LOG.info('%s: %d/%d recovered (%.4f%%); concordant=%d discordant=%d', network, r, n, 100*r/n, concordant, discordant)
    write(OUT / 'summary/03a_baseline_overlap_summary.tsv', SUMMARY)
    return detail_rows


def independent_validation(details, database_pairs):
    # Validate serialized results using a separate database scan and unordered
    # frozensets. This does not use the overlap() matching cache or its sort key.
    decoded_rows = read(OUT / 'qc/03a_database_label_decoding.tsv')[1]
    lookup = {r['database_label']: r['gtdb_r226_species'] for r in decoded_rows}
    targets = defaultdict(list)
    for network, _ in NETWORKS:
        columns, rows = read(OUT / 'baseline_overlap' / (network + '_baseline_overlap.tsv'))
        assert len(rows) == len(details[network])
        for row in rows:
            targets[frozenset([row['external_species_a'], row['external_species_b']])].append(row)
    observed, recovered_per_network = set(), Counter()
    with DB.open(encoding='utf-8', newline='') as f, (OUT / 'database/human_gut_canonical_keys.tsv').open(encoding='utf-8', newline='') as keyfile:
        reader, keys = csv.DictReader(f, delimiter='\t'), csv.DictReader(keyfile, delimiter='\t')
        count = 0
        for source_row, raw in enumerate(reader, 2):
            derived = next(keys)
            count += 1
            a, b = lookup[raw['species_a']], lookup[raw['species_b']]
            unordered = frozenset([a, b])
            assert len(unordered) == 2
            assert derived['source_row_id'] == str(source_row)
            assert derived['database_species_a_raw'] == raw['species_a'] and derived['database_species_b_raw'] == raw['species_b']
            assert derived['decoded_endpoint1'] == a and derived['decoded_endpoint2'] == b
            assert frozenset([derived['canonical_species_a'], derived['canonical_species_b']]) == unordered
            assert derived['canonical_species_a'] < derived['canonical_species_b']
            assert derived['canonical_pair_key'] == derived['canonical_species_a'] + ' || ' + derived['canonical_species_b']
            if unordered not in targets:
                continue
            assert unordered not in observed
            observed.add(unordered)
            for r in targets[unordered]:
                assert r['recovered'] == 'Yes' and r['mack_source_row_id'] == str(source_row)
                assert all(r['mack_' + k] == v for k, v in raw.items()), 'Database values/support changed'
                assert r['mack_sign'] == sign(raw['aver_corr'])
                expected = NA if r['external_sign'] == NA else 'Yes' if r['external_sign'] == r['mack_sign'] else 'No'
                assert r['sign_concordant'] == expected
                recovered_per_network[r['network']] += 1
        assert next(keys, None) is None and count == len(database_pairs)
    for pair, rows in targets.items():
        for r in rows:
            assert (r['recovered'] == 'Yes') == (pair in observed)
            if pair not in observed:
                assert all(r['mack_' + k] == NA for k in STATE['database_fields'])
                assert r['mack_sign'] == r['sign_concordant'] == NA
    for s in SUMMARY:
        check('independent_database_scan_and_support_preservation', recovered_per_network[s['network']] == s['baseline_recovered_pairs'],
              'All database source columns compared verbatim; frozenset identity independently verified', s['network'])
    check('all_database_canonical_derived_rows_independently_verified', True, count)
    check('nine_networks_complete', len(SUMMARY) == 9 and all(r['edge_pair_balance_ok'] for r in SUMMARY), '9/9')
    check('Lacruz_MHO_exploratory_retained', next(r for r in SUMMARY if r['network'] == 'Lacruz_MHO')['primary_or_exploratory'] == 'exploratory')
    check('exact_identity_only_no_taxonomy_rescue', True, 'Exact R226 serialization decoding + case-sensitive unordered species identity only')
    check('no_source_exclusion_or_support_reaggregation', True, 'All recovered project fields copied verbatim from full human_gut.info')


def integrity():
    for r in INVENTORY:
        p = Path(r['absolute_path'])
        r['sha256_after'] = sha(p) if p.is_file() else 'MISSING'
        r['unchanged'] = r['sha256_before'] == r['sha256_after']
    if INVENTORY:
        write(OUT / 'qc/03a_input_integrity.tsv', INVENTORY)
        check('all_stage1_2_2b_2c_frozen_hashes_unchanged', all(r['unchanged'] for r in INVENTORY if r['input_scope'] != 'human_gut_database'), len(INVENTORY)-1)
        check('human_gut_info_hash_unchanged', all(r['unchanged'] for r in INVENTORY if r['input_scope'] == 'human_gut_database'))
        check('no_new_or_missing_files_in_frozen_directories', {str(p) for p in frozen_paths()} == STATE.get('frozen_paths_before'))


def finish(error):
    status = 'FAIL' if error or any(r['status'] == 'FAIL' for r in QC) else 'PASS'
    QC.append(dict(network='ALL', check='overall_status', status=status, detail=error or 'Stage 3A baseline overlap and QC completed'))
    write(OUT / 'qc/03a_overlap_qc.tsv', QC)
    audit = STATE.get('database_audit', {})
    lines = ['# Stage 3A — Baseline MACK-DB overlap validation', '', 'overall_status = ' + status, '',
             '完整 human_gut f=2 数据库 baseline；primary recovery = exact recovered pairs / taxonomy-mapped final unique GTDB pairs。',
             '不排除任何 source/BioProject，不按相关符号、correlation 大小或 p-value 筛选，不重新聚合数据库。', '',
             '## 数据库只读审计', '', '```json', json.dumps(audit, ensure_ascii=False, indent=2), '```', '',
             '数据库标签为数字 ID 前缀与 GTDB 名称的下划线序列化。只删除前缀，并用本地 R226 species 原始名称词汇表做精确反序列化；'
             '该词汇表的编码无冲突，保留大小写及 genus/species 后缀。这不是 synonym mapping 或 taxonomy rescue。',
             '每个标签的证据见 qc/03a_database_label_decoding.tsv；所有原始行的 canonical key 见 database/human_gut_canonical_keys.tsv。'
             '两端点按完全相同的 lexical min/max 规则组成 pair。发现 duplicate canonical pair 将立即 FAIL，不做 collapse。', '',
             '## External 输入与 association 值', '',
             '自动核验 Stage 2B/2C 报告 PASS、机器 QC、pair 文件计数和 primary/exploratory 标签；'
             '正式匹配前写出 qc/03a_input_manifest.tsv。Stage 2B 使用 rescue 后 final 文件；Li 使用 Stage 2C 的 599 pairs。',
             'Stage 2B pair 表的数值从已冻结 audit 按 representative_record_id 取回；仍使用原定首条代表记录，不重新平均。'
             'external_member_values、external_raw_record_ids 保留全部成员数值与来源；external_sign 保留冻结的 collapsed sign。',
             'external_cor 按请求统一命名，但 Guo 是 FlashWeave association weight、Lacruz 是 SPIEC-EASI conditional-independence weight；'
             '它们不应被解释成 Pearson correlation。external_association_type 明确保留数值类型。Zhou 使用已冻结 covariate-corrected SparCC 值；Li 使用文件 cor。', '',
             '## 各网络结果', '',
             '| network | 标签 | raw edges | final pairs | recovered | not recovered | baseline recovery | end-to-end（secondary） | sign concordant / recovered（descriptive） |',
             '| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for r in SUMMARY:
        sign_rate = '{:.4%}'.format(r['sign_concordance_rate_among_recovered']) if r['sign_concordance_rate_among_recovered'] != NA else NA
        lines.append('| {} | {} | {} | {} | {} | {} | {:.4%} | {:.4%} | {} / {} ({}) |'.format(
            r['network'], r['primary_or_exploratory'], r['published_raw_edges'], r['taxonomy_mapped_unique_pairs'],
            r['baseline_recovered_pairs'], r['baseline_not_recovered_pairs'], r['baseline_recovery_rate'],
            r['end_to_end_recovery_rate'], r['sign_concordant_pairs'], r['baseline_recovered_pairs'], sign_rate))
    lines += ['', 'baseline_recovery_rate 的分母始终是 final unique pairs；taxonomy_pair_coverage 和 secondary end-to-end 的分母是 published raw edges。',
              'sign_concordant/sign_discordant 仅描述 recovered pairs。符号不一致仍 recovered=Yes；未恢复的 mack/sign_concordant 字段记 NA。',
              '完整正负/零符号计数、符号不可评估数量和比例见 summary/03a_baseline_overlap_summary.tsv。Lacruz_MHO 始终 exploratory。', '',
              '## QC 与冻结完整性', '',
              '逐网络 recovered + not_recovered = final unique pairs，并独立重读数据库以 frozenset 无向 identity 核验每条 Yes/No。'
              '所有匹配的数据库列（含 proj_ID/proj_corr/proj_pv/proj_sam 等）逐值核对原文件；canonical 派生文件逐行对账。',
              '已对 {} 个冻结文件/参考输入及 human_gut.info 做运行前后 SHA256 检查；并核对 Stage 2C 记录的既有输入 hash。'.format(len(INVENTORY)-1),
              '结果见 qc/03a_input_integrity.tsv 和 qc/03a_overlap_qc.tsv。',
              '截至本报告仅完成 Stage 3A；没有 source exclusion、recurrence、null model 或 enrichment。']
    if error:
        lines += ['', 'FAIL reason: ' + error]
    lines += ['', '## 全部正式输出文件', '']
    lines += ['- `' + str(p) + '`' for p in sorted(OUT.rglob('*')) if p.is_file()]
    lines += ['- `' + str(REPORT) + '`', '']
    REPORT.write_text('\n'.join(lines), encoding='utf-8')
    print(j(dict(overall_status=status, database_records=audit.get('data_records'), frozen_inputs=len(INVENTORY)-1,
                 networks=[{k: r[k] for k in ['network', 'published_raw_edges', 'taxonomy_mapped_unique_pairs',
                     'baseline_recovered_pairs', 'baseline_recovery_rate', 'end_to_end_recovery_rate',
                     'sign_concordant_pairs', 'sign_discordant_pairs', 'sign_concordance_rate_among_recovered']} for r in SUMMARY])), flush=True)
    return 0 if status == 'PASS' else 1


def main():
    error = ''
    try:
        snapshot()
        external = load_external()
        decoded, db_pairs = audit_database()
        details = overlap(external, decoded, db_pairs)
        independent_validation(details, db_pairs)
    except Exception as exc:
        error = '{}: {}'.format(type(exc).__name__, exc)
        LOG.exception('Stage 3A stopped')
        QC.append(dict(network='ALL', check='execution', status='FAIL', detail=error))
    try:
        integrity()
    except Exception as exc:
        error += '; integrity: ' + str(exc)
    return finish(error)


if __name__ == '__main__':
    sys.exit(main())
