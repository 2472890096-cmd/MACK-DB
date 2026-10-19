#!/usr/bin/env python3
"""Offline Stage 2B; all Stage 1/2 files are read-only and hash-verified.

Run: python3 stage2b/02b_taxonomy_rescue.py
Reuse the frozen Stage-2 Chen-compatible NCBI/official/fallback implementation
without calling its main(), importing the side-effectful original Chen scripts,
or writing bytecode. Stage 2B adds evidenced bracket-genus candidates for only
the top 30 unresolved Zhou taxa, exact-version SGB evidence, and an exact genome
accession R220->R226 bridge. No overlap, source exclusion, or recurrence.
"""

# v0.1: repository-local IO configuration; scientific functions below are unchanged.
import sys as _path_sys
from pathlib import Path as _RepoPath
_repo_root = next(p for p in _RepoPath(__file__).resolve().parents if (p / "lib/mackdb_paths.py").is_file())
_path_sys.path.insert(0, str(_repo_root / "lib"))
from mackdb_paths import resolve_path
import ast
import bz2
import csv
import gzip
import hashlib
import importlib.util
import json
import logging
import os
import re
import statistics
import sys
import tempfile
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

BASE = Path(str(resolve_path('literature')))
ROOT = BASE.parent
OUT = BASE / 'stage2b'
REPORTS = BASE / 'reports'
REF = BASE / 'reference'
VERSION = 'mpa_vJan25_CHOCOPhlAnSGB_202503'
SPECIES = REF / 'metaphlan_vjan25' / (VERSION + '_species.txt.bz2')
SGB_TABLE = REF / 'metaphlan_vjan25' / (VERSION + '_SGB2GTDB.tsv')
WHEEL = REF / 'metaphlan_vjan25/MetaPhlAn-4.2.2-py3-none-any.whl'
R220 = {d: REF / ('gtdb_r220/{}_metadata_r220.tsv.gz'.format(p)) for d, p in [('Bacteria', 'bac120'), ('Archaea', 'ar53')]}
R226_DIR = Path(str(resolve_path('gtdb_r226')))
R226 = {d: R226_DIR / ('{}_metadata_r226.tsv.gz'.format(p)) for d, p in [('Bacteria', 'bac120'), ('Archaea', 'ar53')]}
HASHES = {SPECIES: '1207a258f7b08b5c4398207af85a9c2c2fbc0d750ecaba53043a6dabafe901cb',
          SGB_TABLE: '037947d047f48719359905c63a0ba51ae0e12cbfae15d8fda02c2e82524313c2',
          WHEEL: '389d8abce7eb76e236816cf94d0d0cdff609c7f3bd24a7731ba8274d97cf510b',
          R220['Bacteria']: 'f133b46c54a7b7df3ee09e92d1930d84bb9fa987baa55c817dc509ccde49fd81',
          R220['Archaea']: 'cf2bb46db5e731034c527e243c194720834f593daac87a177bf3033565a91c6c'}
BASELINE_PAIRS = {'Zhou2022_MS': 4554, 'Zhou2022_HHC': 3801, 'Guo2023_AggControl': 335,
                  'Guo2023_AggMECFS': 312, 'Lacruz2026_MHNO': 242, 'Lacruz2026_MHO': 290,
                  'Lacruz2026_MUNO': 241, 'Lacruz2026_MUO': 241}
BASELINE_FATES = {'unresolved_ncbi': 2393, 'unresolved_gtdb': 1306, 'pending_sgb': 1241,
                  'self_loop_after_mapping': 58, 'duplicate_after_mapping': 520, 'kept_unique_gtdb_pair': 10016}
OLD_ACCEPTED = {'mapped_official_unique', 'mapped_official_unique_majority', 'mapped_genome_fallback'}
ACCEPTED = OLD_ACCEPTED | {'mapped_speciesfile', 'mapped_r220_r226_bridge'}
FATES = ['unresolved_ncbi', 'unresolved_gtdb', 'unresolved_sgb', 'pending_missing_reference',
         'self_loop_after_mapping', 'duplicate_after_mapping', 'kept_unique_gtdb_pair']
GROUPS = ['MHNO', 'MHO', 'MUNO', 'MUO']
NA = 'NA'
LOG = logging.getLogger('stage2b')
STATE = {'overall_status': 'RUNNING', 'warnings': [], 'qc': [], 'summary': []}


def check(value, reason):
    assert value, reason


def j(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'))


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    with path.open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        return reader.fieldnames, list(reader)


def write(path, rows, fields=None):
    check(OUT in path.resolve().parents or (path.parent == REPORTS and path.name.startswith('02b_')), 'Write outside Stage 2B: ' + str(path))
    fields = fields or list(dict.fromkeys(k for r in rows for k in r))
    check(fields, 'Missing output schema: ' + str(path))
    encoded = [{f: NA if r.get(f) is None else str(r[f]) for f in fields} for r in rows]
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(encoded)
    columns, observed = read(path)
    check(columns == fields and observed == encoded, 'TSV round-trip mismatch: ' + path.name)


def load_framework():
    sys.dont_write_bytecode = True
    path = BASE / 'scripts/02_harmonize_external_taxonomy.py'
    spec = importlib.util.spec_from_file_location('frozen_stage2_framework', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # Definition-only module; main() is not called.
    module.LOG = LOG
    # Explicitly read the original reference scripts without executing them.
    for path in [module.REF_NAME, module.REF_GENOME]:
        ast.parse(path.read_text(encoding='utf-8'))
    return module


def snapshot(framework):
    paths = []
    for directory in [BASE / 'parsed', BASE / 'mapping', BASE / 'scripts']:
        paths.extend(p for p in directory.iterdir() if p.is_file())
    paths += [p for p in REPORTS.iterdir() if p.is_file() and not p.name.startswith('02b_')]
    paths += list(HASHES) + [framework.REF_NAME, framework.REF_GENOME, framework.OFFICIAL, framework.FORWARD,
                            framework.TAXDUMP / 'names.dmp', framework.TAXDUMP / 'nodes.dmp', framework.TAXDUMP / 'merged.dmp', R226['Bacteria']]
    if R226['Archaea'].is_file():
        paths.append(R226['Archaea'])
    inventory = []
    for path in sorted(set(paths)):
        check(path.is_file(), 'Missing fixed input: ' + str(path))
        digest = sha(path)
        inventory.append(dict(absolute_path=str(path), size_bytes=path.stat().st_size, sha256=digest,
                              expected_sha256=HASHES.get(path, NA), hash_matches_expected=digest == HASHES.get(path, digest)))
    STATE['inventory'] = inventory
    write(REPORTS / '02b_input_inventory.tsv', inventory)
    check(all(r['hash_matches_expected'] for r in inventory), 'A fixed-reference SHA256 does not match the supplied value')
    with zipfile.ZipFile(WHEEL) as archive:
        members = [n for n in archive.namelist() if n.endswith(VERSION + '_SGB2GTDB.tsv')]
        check(len(members) == 1, 'Wheel has no unique exact-version SGB2GTDB member')
        check(hashlib.sha256(archive.read(members[0])).hexdigest() == HASHES[SGB_TABLE], 'SGB2GTDB differs from verified wheel member')
    STATE['wheel_member'] = members[0]
    STATE['r226_archaea_metadata_exists'] = R226['Archaea'].is_file()
    LOG.info('Verified %d protected/reference hashes; all 5 pinned hashes match; wheel member verified; R226 Archaea exists=%s',
             len(inventory), R226['Archaea'].is_file())


def verify_inputs():
    rows = []
    for old in STATE['inventory']:
        path = Path(old['absolute_path'])
        digest = sha(path)
        rows.append(dict(old, sha256_after=digest, unchanged=digest == old['sha256'] and path.stat().st_size == old['size_bytes']))
    write(REPORTS / '02b_input_integrity_qc.tsv', rows)
    check(all(r['unchanged'] for r in rows), 'A frozen Stage 1/2 file or taxonomy reference changed')
    STATE['inputs_unchanged'] = True


def fixed_inputs(framework):
    fields, edges = read(BASE / 'parsed/all_new_external_edges_raw.tsv')
    map_fields, taxa = read(BASE / 'mapping/03_external_taxa_to_gtdb_r226.tsv')
    _, before_pairs = read(BASE / 'mapping/04_all_external_edges_gtdb_r226.tsv')
    _, before_audit = read(BASE / 'mapping/04_edge_mapping_audit.tsv')
    _, unresolved = read(REPORTS / '02_unresolved_taxa.tsv')
    _, old_qc = read(REPORTS / '02_taxonomy_mapping_qc.tsv')
    (REPORTS / '02_taxonomy_mapping_qc.md').read_text(encoding='utf-8')
    check(len(edges) == len(before_audit) == 15534, 'Raw/audit baseline size changed')
    check(Counter(r['primary_status'] for r in edges) == {'primary': 14836, 'exploratory': 698}, 'Baseline primary/exploratory counts changed')
    check(Counter(r['network_id'] for r in edges) == framework.EXPECTED, 'Baseline raw network counts changed')
    check(Counter(r['network_id'] for r in before_pairs) == BASELINE_PAIRS and len(before_pairs) == 10016, 'Stage 2 pair baseline changed')
    check(Counter(r['primary_fate'] for r in before_audit) == BASELINE_FATES, 'Stage 2 fate baseline changed')
    check({r['network_id']: int(r['final_unique_gtdb_pairs']) for r in old_qc} == BASELINE_PAIRS, 'Stage 2 QC baseline changed')
    check({(r['paper_id'], r['taxon']) for r in unresolved} == {(r['paper_id'], r['raw_taxon']) for r in taxa if r['mapping_status'] not in OLD_ACCEPTED},
          'Frozen unresolved taxon inventory disagrees with frozen mapping')
    for i, (raw, old) in enumerate(zip(edges, before_audit), 1):
        check(all(raw[f] == old[f] for f in fields) and old['raw_record_id'] == 'E{:05d}'.format(i), 'Frozen edge/audit row mismatch')
        raw.update(raw_record_id=old['raw_record_id'], input_tsv_row=i+1)
    STATE.update(input_fields=fields, map_fields=map_fields, baseline_qc=old_qc)
    return edges, taxa, before_pairs, before_audit


def zhou_degrees(edges, taxa, stage):
    unresolved = [r for r in taxa if r['paper_id'] == 'Zhou2022' and r['mapping_status'] not in OLD_ACCEPTED]
    check(len(unresolved) == 148, 'Zhou unresolved taxon count must be 148')
    records, neighbors = defaultdict(lambda: defaultdict(set)), defaultdict(lambda: defaultdict(set))
    for row in edges:
        if row['paper_id'] != 'Zhou2022':
            continue
        group = row['condition']
        for a, b in [('taxon1_raw', 'taxon2_raw'), ('taxon2_raw', 'taxon1_raw')]:
            records[row[a]][group].add(row['raw_record_id'])
            neighbors[row[a]][group].add(row[b])
    audit = []
    for old in unresolved:
        raw = old['raw_taxon']
        audit.append(dict(raw_taxon=raw, MS_degree=len(neighbors[raw]['MS']), HHC_degree=len(neighbors[raw]['HHC']),
                          MS_affected_edges=len(records[raw]['MS']), HHC_affected_edges=len(records[raw]['HHC']),
                          total_unique_affected_edges=len(records[raw]['MS'] | records[raw]['HHC']),
                          old_ncbi_status=old['ncbi_stage_status'], old_gtdb_status=old['mapping_status'], old_reason=old['mapping_reason']))
    audit.sort(key=lambda r: (-r['total_unique_affected_edges'], r['raw_taxon']))
    for i, row in enumerate(audit, 1):
        row['degree_rank'] = i
        row['top30_review'] = i <= 30
    write(stage / 'zhou_unresolved_degree_audit.tsv', audit)
    values = sorted(r['total_unique_affected_edges'] for r in audit)
    STATE['zhou_degree_distribution'] = dict(count=len(values), minimum=min(values), maximum=max(values),
                                            median=statistics.median(values), mean=statistics.mean(values),
                                            histogram=j(dict(sorted(Counter(values).items()))))
    STATE['zhou_top30'] = audit[:30]
    return audit[:30]


def pending_manifest(edges, taxa, stage):
    pending = [r for r in taxa if r['paper_id'] == 'Lacruz2026' and r['mapping_status'] == 'pending_SGB_mapping']
    check(len(pending) == 138, 'Initial pending SGB taxa must equal 138')
    counts = defaultdict(Counter)
    affected = set()
    labels = {r['raw_taxon'] for r in pending}
    for edge in edges:
        if edge['paper_id'] == 'Lacruz2026':
            for raw in {edge['taxon1_raw'], edge['taxon2_raw']} & labels:
                counts[raw][edge['condition']] += 1
                affected.add(edge['raw_record_id'])
    manifest = []
    for old in pending:
        parts = dict((p[:1], p[3:]) for p in old['raw_taxon'].split('|'))
        check(parts.get('t') == old['metaphlan_sgb'] and re.fullmatch(r'SGB\d+', parts.get('t', '')), 'Invalid frozen pending SGB ID')
        c = counts[old['raw_taxon']]
        manifest.append(dict(raw_taxon=old['raw_taxon'], species_field=parts['s'], sgb_id=parts['t'],
                             ggb_id=parts['g'] if re.fullmatch(r'GGB\d+', parts.get('g', '')) else NA,
                             metaphlan_genus_label=parts.get('g', NA), domain=parts.get('k', NA), network_count=len(c),
                             **{g+'_edges': c[g] for g in GROUPS},
                             primary_affected_edges=sum(c[g] for g in ['MHNO', 'MUNO', 'MUO']),
                             exploratory_affected_edges=c['MHO'], total_affected_edges=sum(c.values())))
    check(len({r['sgb_id'] for r in manifest}) == 138, 'Unique pending SGB IDs must equal 138')
    check(len(affected) == 1241, 'Pending SGB affected edge count must equal 1241')
    status_counts = Counter(r['primary_status'] for r in edges if r['raw_record_id'] in affected)
    check(status_counts == {'primary': 874, 'exploratory': 367}, 'Pending SGB primary/exploratory baseline changed')
    manifest.sort(key=lambda r: (-r['total_affected_edges'], r['sgb_id']))
    STATE['archaeal_pending_sgb'] = sum(r['domain'] == 'Archaea' for r in manifest)
    STATE['pending_record_ids'] = affected
    write(stage / 'lacruz_pending_sgb_manifest.tsv', manifest)
    LOG.info('Frozen pending SGB=%d, affected records=1241 (874 primary + 367 exploratory); archaeal=%d', len(manifest), STATE['archaeal_pending_sgb'])
    return manifest


def species_schema_and_matches(pending, stage):
    """Inspect the actual two-column schema before interpreting assignments.
    Comma-delimited full k__ lineages are separate species assignments. All
    matching source lines, including placeholder labels, are retained verbatim.
    """
    matches, labels, label_sgbs, columns = defaultdict(list), defaultdict(set), defaultdict(set), Counter()
    id_prefixes, domain_prefixes = Counter(), Counter()
    with bz2.open(SPECIES, 'rt', encoding='utf-8', newline='') as handle:
        for line_no, line in enumerate(handle, 1):
            raw = line.rstrip('\r\n')
            if line_no <= 20:
                LOG.info('species.txt first20 line %d: %s', line_no, raw)
            fields = raw.split('\t')
            columns[len(fields)] += 1
            check(len(fields) == 2 and fields[0] and not any(c.isspace() for c in fields[0]), 'Unexpected actual species-file column schema at line ' + str(line_no))
            prefix = re.match(r'[A-Za-z]+', fields[0])
            id_prefixes[prefix.group() if prefix else 'other'] += 1
            assignments = fields[1].split(',')
            check(all(a.startswith('k__') for a in assignments), 'Unrecognized species-file assignment delimiter')
            for assignment in assignments:
                names = [x[3:] for x in assignment.split('|') if x.startswith('s__')]
                check(len(names) == 1, 'Species-file assignment has no unique s__ field')
                labels[fields[0]].add(names[0])
                domain_prefixes[assignment.split('|', 1)[0]] += 1
                if fields[0].startswith('SGB'):
                    label_sgbs[names[0]].add(fields[0])
            matches[fields[0]].append(dict(line_number=line_no, raw_record=raw, assignments=assignments))
    schema = dict(line_count=line_no, column_count_distribution=dict(columns), delimiter='TAB', header='absent',
                  sgb_field='column 1; cluster ID (SGB and EUK observed); pending taxa joined by exact SGB ID',
                  species_label_field='column 2; s__ token of each comma-separated full lineage',
                  lineage_delimiter='|', assignment_delimiter=',', unique_cluster_ids=len(matches),
                  unique_sgb_ids=sum(s.startswith('SGB') for s in matches), cluster_id_prefix_counts=dict(id_prefixes),
                  assignment_domain_counts=dict(domain_prefixes),
                  species_labels_with_multiple_sgbs=sum(len(v) > 1 for v in label_sgbs.values()),
                  sgbs_with_multiple_species_labels=sum(len(v) > 1 for s, v in labels.items() if s.startswith('SGB')))
    (stage / 'metaphlan_species_file_schema.txt').write_text(json.dumps(schema, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    STATE['species_schema'] = schema
    audit = []
    for row in pending:
        sgb = row['sgb_id']
        assigned = sorted(labels.get(sgb, set()))
        named = [name for name in assigned if name and not re.search(r'(?:^|_)SGB\d+(?:_|$)', name)]
        audit.append(dict(sgb_id=sgb, raw_taxon=row['raw_taxon'], species_file_found=sgb in matches,
                          species_file_assignment_count=sum(len(m['assignments']) for m in matches.get(sgb, [])),
                          species_file_species_labels=j(assigned), species_file_named_species_labels=j(named),
                          species_file_raw_records=j(matches.get(sgb, [])), named_species_assignment_count=len(named)))
    write(stage / 'lacruz_sgb_speciesfile_audit.tsv', audit)
    STATE['speciesfile_counts'] = dict(found_sgb=sum(r['species_file_found'] for r in audit),
                                      not_found_sgb=sum(not r['species_file_found'] for r in audit),
                                      with_named_species=sum(r['named_species_assignment_count'] > 0 for r in audit),
                                      without_named_species=sum(r['named_species_assignment_count'] == 0 for r in audit))
    LOG.info('Actual species-file schema=%s; pending SGB matches=%s', j(schema), j(STATE['speciesfile_counts']))
    return audit


def zhou_candidates(raw):
    candidates = [(raw, 'original_name')]
    formatted = re.sub(r'\s+', ' ', raw.replace('_', ' ')).strip()
    candidates.append((formatted, 'whitespace_and_underscore_formatting'))
    unbracketed = re.sub(r'^\[([^\[\]]+)\](?=\s)', r'\1', formatted)
    candidates.append((unbracketed, 'leading_bracketed_genus_candidate'))
    result, seen = [], set()
    for name, method in candidates:
        if name not in seen:
            result.append((name, method))
            seen.add(name)
    return result


def name_rescues(taxa, top30, species_audit, framework, stage):
    """Bracket removal only generates a query. Acceptance requires unique current
    species-rank evidence in the frozen NCBI taxonomy or the frozen Chen TaxID
    dictionary. Existing confirmed NCBI identities/ties are not reinterpreted.
    """
    old_lookup = {(r['paper_id'], r['raw_taxon']): r for r in taxa}
    queries, query_meta = [], {}
    for review in top30:
        old = old_lookup[('Zhou2022', review['raw_taxon'])]
        if old['ncbi_stage_status'] == 'mapped_species':
            continue
        for name, operation in zhou_candidates(old['raw_taxon']):
            query_id = 'ZC{:03d}'.format(len(queries) + 1)
            queries.append(dict(old, raw_taxon_id=query_id, normalized_candidate=name))
            query_meta[query_id] = dict(route='zhou', original_raw_taxon=old['raw_taxon'], candidate_operation=operation)
    for audit in species_audit:
        old = old_lookup[('Lacruz2026', audit['raw_taxon'])]
        for label in json.loads(audit['species_file_named_species_labels']):
            query_id = 'SC{:03d}'.format(len(queries) + 1)
            queries.append(dict(old, raw_taxon_id=query_id, normalized_candidate=label.replace('_', ' '), lacruz_taxon_class='named_species',
                                source_taxonomy_type='MetaPhlAn4_exact_version_species_file'))
            query_meta[query_id] = dict(route='speciesfile', original_raw_taxon=old['raw_taxon'], sgb_id=audit['sgb_id'], candidate_operation='species_file_named_label')
    taxonomy = framework.Taxonomy(queries)
    resolved, evidence = [], []
    for query in queries:
        row, hits = taxonomy.resolve(query)
        row.update(query_meta[query['raw_taxon_id']])
        resolved.append(row)
        evidence.extend(hits)
    evidence_fields = ['raw_taxon_id', 'paper_id', 'raw_taxon', 'match_tier', 'matched_ncbi_alias', 'name_class',
                       'alias_taxid', 'current_taxid', 'current_rank', 'current_domain']
    write(stage / 'ncbi_candidate_name_evidence.tsv', evidence, evidence_fields)
    by_raw = defaultdict(list)
    for row in resolved:
        if row['route'] == 'zhou':
            by_raw[row['original_raw_taxon']].append(row)
    proposals, decisions = {}, []
    for review in top30:
        raw = review['raw_taxon']
        old = old_lookup[('Zhou2022', raw)]
        candidates = by_raw[raw]
        good = [r for r in candidates if r['ncbi_stage_status'] == 'mapped_species']
        ids = {r['ncbi_taxid'] for r in good}
        decision = dict(raw_taxon=raw, degree_rank=review['degree_rank'], total_unique_affected_edges=review['total_unique_affected_edges'],
                        candidate_name=old['normalized_candidate'], candidate_taxid=old['ncbi_taxid'], candidate_rank=old['ncbi_rank'],
                        ncbi_evidence=old['ncbi_evidence'], old_status=old['mapping_status'], new_status=old['mapping_status'],
                        decision='remain_unresolved', decision_reason=NA, final_gtdb_r226_species=NA,
                        all_candidates=j([dict(candidate_name=r['normalized_candidate'], candidate_taxid=r['ncbi_taxid'],
                                               candidate_rank=r['ncbi_rank'], ncbi_status=r['ncbi_stage_status'],
                                               operation=r['candidate_operation'], evidence=r['ncbi_evidence']) for r in candidates]))
        if old['ncbi_stage_status'] == 'mapped_species':
            decision['decision_reason'] = 'Existing validated NCBI identity retained; prior official/genome ambiguity or absence is not bypassed: ' + old['mapping_reason']
        elif len(ids) > 1:
            decision.update(decision='candidate_not_accepted', candidate_taxid=j(sorted(ids)),
                            decision_reason='Deterministic candidates resolve to multiple species TaxIDs; no choice or expansion allowed')
        elif len(ids) == 1:
            selected = good[0]
            proposals[raw] = selected
            decision.update(candidate_name=selected['normalized_candidate'], candidate_taxid=selected['ncbi_taxid'], candidate_rank='species',
                            ncbi_evidence=selected['ncbi_evidence'], decision='accepted_strict_rescue',
                            decision_reason='New NCBI identity uniquely confirmed by accepted taxonomy/Chen stable TaxID evidence; formatting alone is not evidence')
        else:
            if any(r['ncbi_candidate_taxids'] not in ('[]', NA) for r in candidates):
                decision['decision'] = 'candidate_not_accepted'
            decision['decision_reason'] = 'No deterministic candidate uniquely resolves to an accepted prokaryotic species-rank TaxID'
        decisions.append(decision)
    # Only newly evidenced identities and named SGB candidates enter this step.
    to_gtdb = list(proposals.values()) + [r for r in resolved if r['route'] == 'speciesfile']
    official, forward = framework.load_official()
    framework.apply_official(to_gtdb, official, forward)
    if any(r.get('fallback_eligible') for r in to_gtdb):
        framework.genome_fallback(to_gtdb, taxonomy, stage)
    else:
        LOG.info('No new named candidate requires genome fallback')
    for decision in decisions:
        if decision['raw_taxon'] in proposals:
            row = proposals[decision['raw_taxon']]
            decision.update(new_status=row['mapping_status'], final_gtdb_r226_species=row['gtdb_r226_species'])
    write(stage / 'zhou_rescue_audit.tsv', decisions)
    candidate_fields = ['raw_taxon_id', 'route', 'original_raw_taxon', 'normalized_candidate', 'candidate_operation',
                        'ncbi_name', 'ncbi_taxid', 'ncbi_rank', 'ncbi_mapping_method', 'ncbi_stage_status', 'ncbi_evidence',
                        'gtdb_r226_species', 'gtdb_mapping_method', 'genome_support_top', 'genome_support_total',
                        'mapping_status', 'mapping_reason', 'official_status', 'fallback_eligible']
    write(stage / 'named_candidate_mapping_evidence.tsv', resolved, candidate_fields)
    species_candidates = defaultdict(list)
    for row in to_gtdb:
        if row['route'] == 'speciesfile':
            species_candidates[row['sgb_id']].append(row)
    species_decisions = {}
    for audit in species_audit:
        sgb = audit['sgb_id']
        acceptable = [r for r in species_candidates[sgb] if r['mapping_status'] in OLD_ACCEPTED]
        targets = {r['gtdb_r226_species'] for r in acceptable}
        if len(targets) == 1:
            status, target = 'mapped_speciesfile', next(iter(targets))
            reason = 'All acceptable named-species candidates converge on one R226 species'
        elif len(targets) > 1:
            status, target = 'unresolved_multiple_species_assignments', NA
            reason = 'Acceptable named species point to multiple R226 species; no selection or bridge bypass'
        else:
            status, target = 'needs_sgb2gtdb_route', NA
            reason = 'No usable named-species mapping from the exact-version species file'
        species_decisions[sgb] = dict(speciesfile_status=status, speciesfile_final_r226=target,
                                      speciesfile_reason=reason, speciesfile_candidates=species_candidates[sgb])
    STATE['zhou_rescue_decisions'] = decisions
    STATE['zhou_strict_identity_rescues'] = len(proposals)
    STATE['zhou_new_gtdb_mapped'] = sum(r['mapping_status'] in OLD_ACCEPTED for r in proposals.values())
    STATE['speciesfile_counts'].update(rescued_sgb=sum(r['speciesfile_status'] == 'mapped_speciesfile' for r in species_decisions.values()),
                                       unresolved_sgb=sum(r['speciesfile_status'] != 'mapped_speciesfile' for r in species_decisions.values()))
    LOG.info('Zhou strict NCBI rescues=%d; newly GTDB mapped=%d; species-file rescued SGB=%d',
             len(proposals), STATE['zhou_new_gtdb_mapped'], STATE['speciesfile_counts']['rescued_sgb'])
    return proposals, species_decisions


def sgb_r220_assignments(pending, species_decisions, stage):
    source = defaultdict(list)
    unclassified_markers = Counter()
    with SGB_TABLE.open(encoding='utf-8') as handle:
        for line_no, line in enumerate(handle, 1):
            parts = line.rstrip('\r\n').split('\t')
            check(len(parts) == 2 and re.fullmatch(r'SGB\d+(?:_group)?', parts[0]), 'SGB2GTDB is not the expected observed headerless two-column table')
            if parts[1] in {'Unclassified', 'Unclassified Bacteria'}:
                # These literal values are present in four rows of this pinned
                # file. Preserve the source record; neither supplies a species.
                unclassified_markers[parts[1]] += 1
                lineage = dict.fromkeys('dpcofgs', '')
                if parts[1] == 'Unclassified Bacteria':
                    lineage['d'] = 'Bacteria'
            else:
                lineage = dict((token[0], token[3:]) for token in parts[1].split(';'))
                check(set(lineage) == set('dpcofgs'), 'Unrecognized SGB2GTDB lineage schema')
            source[parts[0]].append(dict(line_number=line_no, raw_line=line.rstrip('\r\n'), lineage=lineage))
    result = []
    for row in pending:
        if species_decisions[row['sgb_id']]['speciesfile_status'] != 'needs_sgb2gtdb_route':
            continue
        hits = source.get(row['sgb_id'], [])
        assignments = {(h['lineage']['d'], h['lineage']['s']) for h in hits if h['lineage']['s']}
        record = dict(sgb_id=row['sgb_id'], raw_taxon=row['raw_taxon'], sgb2gtdb_found=bool(hits), assignment_count=len(hits),
                      domain=NA, phylum=NA, **{'class': NA}, order=NA, family=NA, genus=NA, species=NA, gtdb_r220_species=NA,
                      mapping_status=NA, mapping_reason=NA, source_raw_records=j(hits))
        if len(assignments) > 1:
            record.update(mapping_status='unresolved_multiple_species_assignments', mapping_reason='SGB2GTDB contains distinct R220 species assignments; no one-to-many expansion')
        elif not assignments:
            record.update(mapping_status='unresolved_no_r220_species', mapping_reason='No nonempty R220 species field; genus-level assignment is not a species')
            if len(hits) == 1:
                record.update({name: hits[0]['lineage'][letter] or NA for name, letter in
                               [('domain', 'd'), ('phylum', 'p'), ('class', 'c'), ('order', 'o'), ('family', 'f'), ('genus', 'g')]})
        else:
            domain, species = next(iter(assignments))
            check(domain == row['domain'], 'MetaPhlAn/SGB2GTDB domain mismatch: ' + row['sgb_id'])
            hit = next(h for h in hits if h['lineage']['s'] == species)
            record.update({name: hit['lineage'][letter] or NA for name, letter in
                           [('domain', 'd'), ('phylum', 'p'), ('class', 'c'), ('order', 'o'), ('family', 'f'), ('genus', 'g'), ('species', 's')]})
            record.update(gtdb_r220_species=species, mapping_status='needs_r220_r226_bridge', mapping_reason='R220 assignment is evidence only; requires accession-supported release bridge')
        result.append(record)
    write(stage / 'lacruz_sgb_to_gtdb_r220.tsv', result,
          ['sgb_id', 'raw_taxon', 'sgb2gtdb_found', 'assignment_count', 'domain', 'phylum', 'class', 'order', 'family', 'genus',
           'species', 'gtdb_r220_species', 'mapping_status', 'mapping_reason', 'source_raw_records'])
    STATE['sgb2gtdb_counts'] = dict(examined_sgb=len(result), with_r220_species=sum(r['gtdb_r220_species'] != NA for r in result),
                                   no_r220_species=sum(r['mapping_status'] == 'unresolved_no_r220_species' for r in result),
                                   full_file_unclassified_marker_rows=dict(unclassified_markers))
    LOG.info('SGB2GTDB actual schema: 2 columns; full-file literal unclassified markers=%s; requested SGB=%s', j(unclassified_markers), j(STATE['sgb2gtdb_counts']))
    return result


def accession_key(accession):
    """Only remove one leading RS_/GB_; preserve GCA/GCF and accession version."""
    return accession[3:] if accession.startswith(('RS_', 'GB_')) else accession


def metadata_rows(path):
    with gzip.open(path, 'rt', encoding='utf-8', newline='') as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader)
        check('accession' in header and 'gtdb_taxonomy' in header, 'Missing metadata accession/taxonomy columns')
        a, t = header.index('accession'), header.index('gtdb_taxonomy')
        for index, values in enumerate(reader, 1):
            lineage = values[t]
            species = [token[3:] for token in lineage.split(';') if token.startswith('s__')]
            check(len(species) == 1, 'No unique s__ field in metadata')
            yield values[a], species[0], lineage
            if index % 200000 == 0:
                LOG.info('%s: scanned %d genome records', path.name, index)


def build_bridge(assignments, stage):
    targets = {(r['domain'], r['gtdb_r220_species']) for r in assignments if r['mapping_status'] == 'needs_r220_r226_bridge'}
    LOG.info('R220->R226 bridge: %d unique R220 species; exact normalized genome accession joins only', len(targets))
    r220_genomes, shared, counts = defaultdict(dict), {}, defaultdict(Counter)
    missing_domains = set()
    for domain in sorted({d for d, _ in targets}):
        check(domain in R220, 'Unsupported R220 domain')
        if not R226[domain].is_file():
            check(domain == 'Archaea', 'Missing required bacterial R226 metadata')
            missing_domains.add(domain)
            STATE['warnings'].append('MISSING_R226_ARCHAEA_METADATA')
            continue
        wanted = {s for d, s in targets if d == domain}
        lookup = {}
        for accession, species, lineage in metadata_rows(R220[domain]):
            if species in wanted:  # Equality only within the R220 release.
                normalized = accession_key(accession)
                check(normalized not in lookup, 'Repeated normalized R220 genome accession')
                lookup[normalized] = (domain, species)
                r220_genomes[(domain, species)][normalized] = (accession, lineage)
        LOG.info('%s R220 selected genomes=%d across %d requested species', domain, len(lookup), len(wanted))
        for accession, species, lineage in metadata_rows(R226[domain]):
            normalized = accession_key(accession)
            if normalized in lookup:
                identity = lookup[normalized]
                check((domain, normalized) not in shared, 'Repeated normalized R226 genome accession')
                shared[(domain, normalized)] = (accession, species, lineage)
                if species:
                    counts[identity][species] += 1
    bridge, candidate_rows, genome_rows = [], [], []
    for identity in sorted(targets):
        domain, r220_species = identity
        genomes = r220_genomes[identity]
        shared_ids = [a for a in genomes if (domain, a) in shared]
        support = counts[identity]
        ordered = sorted(support.items(), key=lambda x: (-x[1], x[0]))
        top = ordered[0][1] if ordered else 0
        tops = [s for s, count in ordered if count == top]
        second = ordered[1][1] if len(ordered) > 1 else 0
        row = dict(domain=domain, r220_species=r220_species, r220_total_genomes=len(genomes),
                   shared_genomes_r226=len(shared_ids), classified_shared_genomes=sum(support.values()),
                   r226_candidate_count=len(support), r226_top_species=tops[0] if len(tops) == 1 else NA,
                   r226_top_species_candidates=j(tops), r226_top_count=top, r226_second_count=second,
                   r226_top_proportion=top / len(shared_ids) if shared_ids else NA, r226_top_tie_count=len(tops),
                   mapping_status=NA, mapping_reason=NA, r220_metadata=str(R220[domain]), r226_metadata=str(R226[domain]))
        if domain in missing_domains:
            row.update(mapping_status='pending_missing_r226_archaea_metadata', mapping_reason='MISSING_R226_ARCHAEA_METADATA; no download or bacterial substitution')
        elif not genomes:
            row.update(mapping_status='unresolved_no_r220_genomes', mapping_reason='Exact R220 species has no genome rows in R220 metadata')
        elif not shared_ids:
            row.update(mapping_status='unresolved_no_shared_genomes', mapping_reason='No exact accession-and-version intersection with R226')
        elif not support:
            row.update(mapping_status='unresolved_other', mapping_reason='Shared genomes lack nonempty R226 species assignments')
        elif len(tops) > 1:
            row.update(mapping_status='unresolved_r220_r226_top_tie', mapping_reason='Highest shared-genome count is tied; no choice or expansion')
        else:
            row.update(mapping_status='mapped_r220_r226_bridge', mapping_reason='Unique highest R226 support among exact shared genome accessions; no proportion threshold')
        bridge.append(row)
        for species, count in ordered:
            candidate_rows.append(dict(domain=domain, r220_species=r220_species, r226_species=species,
                                       shared_genome_count=count, is_top=count == top))
        for normalized, (old_accession, old_lineage) in sorted(genomes.items()):
            hit = shared.get((domain, normalized))
            genome_rows.append(dict(domain=domain, r220_species=r220_species, normalized_accession=normalized,
                                    r220_accession=old_accession, r226_accession=hit[0] if hit else NA,
                                    shared=bool(hit), r226_species=hit[1] if hit and hit[1] else NA,
                                    r220_taxonomy=old_lineage, r226_taxonomy=hit[2] if hit else NA))
        check(sum(support.values()) <= len(shared_ids) <= len(genomes), 'Genome bridge support accounting failed')
    write(stage / 'gtdb_r220_to_r226_species_bridge.tsv', bridge,
          ['domain', 'r220_species', 'r220_total_genomes', 'shared_genomes_r226', 'classified_shared_genomes',
           'r226_candidate_count', 'r226_top_species', 'r226_top_species_candidates', 'r226_top_count', 'r226_second_count',
           'r226_top_proportion', 'r226_top_tie_count', 'mapping_status', 'mapping_reason', 'r220_metadata', 'r226_metadata'])
    write(stage / 'gtdb_r220_r226_candidate_genome_support.tsv', candidate_rows,
          ['domain', 'r220_species', 'r226_species', 'shared_genome_count', 'is_top'])
    write(stage / 'gtdb_r220_r226_genome_join_audit.tsv', genome_rows,
          ['domain', 'r220_species', 'normalized_accession', 'r220_accession', 'r226_accession', 'shared', 'r226_species', 'r220_taxonomy', 'r226_taxonomy'])
    STATE['bridge_status_counts'] = dict(Counter(r['mapping_status'] for r in bridge))
    STATE['bridge_species_count'] = len(bridge)
    LOG.info('Genome bridge decisions: %s', j(STATE['bridge_status_counts']))
    return {(r['domain'], r['r220_species']): r for r in bridge}


def final_sgb_mapping(pending, species_audit, species_decisions, assignments, bridges, stage):
    sources = {r['sgb_id']: r for r in species_audit}
    assignment_lookup = {r['sgb_id']: r for r in assignments}
    rows = []
    for manifest in pending:
        sgb = manifest['sgb_id']
        source, decision = sources[sgb], species_decisions[sgb]
        row = dict(manifest, species_file_found=source['species_file_found'], species_file_assignments=source['species_file_species_labels'],
                   speciesfile_final_r226=decision['speciesfile_final_r226'], gtdb_r220_species=NA,
                   r220_total_genomes=NA, shared_genomes_r226=NA, r226_candidate_count=NA, r226_top_species=NA,
                   r226_top_count=NA, r226_second_count=NA, r226_top_proportion=NA, r226_top_tie_count=NA,
                   final_gtdb_r226_species=NA, mapping_method=NA, mapping_status=NA, mapping_reason=NA)
        if decision['speciesfile_status'] == 'mapped_speciesfile':
            row.update(final_gtdb_r226_species=decision['speciesfile_final_r226'], mapping_status='mapped_speciesfile',
                       mapping_method='metaphlan_speciesfile_consensus_to_r226', mapping_reason=decision['speciesfile_reason'])
        elif decision['speciesfile_status'] == 'unresolved_multiple_species_assignments':
            row.update(mapping_status=decision['speciesfile_status'], mapping_reason=decision['speciesfile_reason'])
        else:
            assigned = assignment_lookup[sgb]
            row['gtdb_r220_species'] = assigned['gtdb_r220_species']
            if assigned['mapping_status'] != 'needs_r220_r226_bridge':
                row.update(mapping_status=assigned['mapping_status'], mapping_reason=assigned['mapping_reason'])
            else:
                bridge = bridges[(assigned['domain'], assigned['gtdb_r220_species'])]
                for field in ['r220_total_genomes', 'shared_genomes_r226', 'r226_candidate_count', 'r226_top_species',
                              'r226_top_count', 'r226_second_count', 'r226_top_proportion', 'r226_top_tie_count']:
                    row[field] = bridge[field]
                row.update(mapping_status=bridge['mapping_status'], mapping_reason=bridge['mapping_reason'])
                if bridge['mapping_status'] == 'mapped_r220_r226_bridge':
                    row.update(final_gtdb_r226_species=bridge['r226_top_species'], mapping_method='metaphlan_sgb2gtdb_r220_exact_accession_bridge_to_r226')
        check(row['mapping_status'] != NA, 'Unassigned SGB decision')
        check((row['final_gtdb_r226_species'] != NA) == (row['mapping_status'] in ACCEPTED), 'SGB target/status inconsistency')
        rows.append(row)
    write(stage / 'lacruz_sgb_to_gtdb_r226_final.tsv', rows)
    STATE['sgb_final_counts'] = dict(Counter(r['mapping_status'] for r in rows))
    STATE['sgb_rescued'] = sum(r['mapping_status'] in ACCEPTED for r in rows)
    STATE['sgb_unresolved'] = 138 - STATE['sgb_rescued']
    STATE['sgb_rescue_rate'] = STATE['sgb_rescued'] / 138
    return rows


def overlay_mapping(taxa, proposals, sgb_rows, species_decisions, stage):
    sgb_lookup = {r['raw_taxon']: r for r in sgb_rows}
    taxonomy_fields = ['normalized_candidate', 'ncbi_name', 'ncbi_taxid', 'ncbi_rank', 'ncbi_domain', 'ncbi_mapping_method',
                       'ncbi_stage_status', 'ncbi_name_classes', 'ncbi_evidence', 'ncbi_candidate_taxids',
                       'gtdb_r226_species', 'gtdb_mapping_method', 'genome_support_top', 'genome_support_total', 'genome_support_basis',
                       'mapping_status', 'mapping_reason', 'official_status', 'official_source_name', 'official_target_name',
                       'official_candidate_count', 'official_top_percent', 'official_top_tie_count', 'fallback_eligible',
                       'gtdb_evidence_file', 'genome_metadata_rows', 'genome_unclassified_rows', 'genome_candidate_count',
                       'genome_top_tie_count', 'genome_top_percent']
    result = []
    for old in taxa:
        row = dict(old, stage2b_changed=False, stage2b_method='unchanged', stage2b_evidence=NA,
                   old_mapping_status=old['mapping_status'], new_mapping_status=old['mapping_status'], stage2b_mapping_reason=old['mapping_reason'])
        if old['paper_id'] == 'Zhou2022' and old['raw_taxon'] in proposals:
            proposal = proposals[old['raw_taxon']]
            check(old['mapping_status'] not in OLD_ACCEPTED, 'Attempt to overwrite accepted Zhou mapping')
            for field in taxonomy_fields:
                row[field] = proposal.get(field, NA)
            row.update(stage2b_changed=True, stage2b_method='zhou_strict_ncbi_identity_rescue_then_chen_r226',
                       stage2b_evidence=j(dict(candidate_name=proposal['normalized_candidate'], candidate_operation=proposal['candidate_operation'],
                                               ncbi_evidence=proposal['ncbi_evidence'], mapping_evidence_file=proposal.get('gtdb_evidence_file', NA))),
                       new_mapping_status=proposal['mapping_status'], stage2b_mapping_reason=proposal['mapping_reason'])
        elif old['paper_id'] == 'Lacruz2026' and old['raw_taxon'] in sgb_lookup:
            sgb = sgb_lookup[old['raw_taxon']]
            check(old['mapping_status'] == 'pending_SGB_mapping', 'Only frozen pending SGB labels may be rescued')
            row.update(stage2b_method=sgb['mapping_method'] if sgb['mapping_method'] != NA else 'sgb_audit_no_accepted_rescue',
                       stage2b_evidence=j(dict(sgb_id=sgb['sgb_id'], final_audit='lacruz_sgb_to_gtdb_r226_final.tsv',
                                               species_file_sha256=HASHES[SPECIES], sgb2gtdb_sha256=HASHES[SGB_TABLE],
                                               r220_species=sgb['gtdb_r220_species'], bridge='gtdb_r220_to_r226_species_bridge.tsv')),
                       new_mapping_status=sgb['mapping_status'], stage2b_mapping_reason=sgb['mapping_reason'])
            if sgb['mapping_status'] in ACCEPTED:
                row.update(stage2b_changed=True, gtdb_r226_species=sgb['final_gtdb_r226_species'], gtdb_mapping_method=sgb['mapping_method'],
                           mapping_status=sgb['mapping_status'], mapping_reason=sgb['mapping_reason'],
                           official_status='not_applicable_direct_gtdb_bridge', fallback_eligible=False)
                if sgb['mapping_status'] == 'mapped_r220_r226_bridge':
                    row.update(ncbi_stage_status='not_applicable_genome_accession_bridge', ncbi_mapping_method=NA,
                               genome_support_top=sgb['r226_top_count'], genome_support_total=sgb['shared_genomes_r226'],
                               genome_support_basis='exact shared R220/R226 genome accession-and-version counts',
                               gtdb_evidence_file=str(SGB_TABLE) + ';' + str(R220[sgb['domain']]) + ';' + str(R226[sgb['domain']]))
                else:
                    accepted_candidates = [r for r in species_decisions[sgb['sgb_id']]['speciesfile_candidates'] if r['mapping_status'] in OLD_ACCEPTED]
                    ids = {r['ncbi_taxid'] for r in accepted_candidates}
                    row.update(ncbi_stage_status='verified_speciesfile_named_candidate_consensus',
                               ncbi_evidence=j([dict(ncbi_taxid=r['ncbi_taxid'], ncbi_name=r['ncbi_name'], gtdb_species=r['gtdb_r226_species']) for r in accepted_candidates]),
                               ncbi_mapping_method='exact_version_speciesfile_named_candidate_consensus', gtdb_evidence_file=str(SPECIES))
                    if len(ids) == 1:
                        for field in ['ncbi_taxid', 'ncbi_name', 'ncbi_rank', 'ncbi_domain']:
                            row[field] = accepted_candidates[0][field]
        if not row['stage2b_changed']:
            check(all(row[f] == old[f] for f in STATE['map_fields']), 'Nonrescued core mapping was altered')
        if old['mapping_status'] in OLD_ACCEPTED or old['paper_id'] == 'Guo2023':
            check(not row['stage2b_changed'] and all(row[f] == old[f] for f in STATE['map_fields']), 'Frozen accepted/Guo mapping changed')
        result.append(row)
    check(len(result) == len(taxa) == 1289, 'Taxon mapping row count changed')
    fields = STATE['map_fields'] + ['stage2b_changed', 'stage2b_method', 'stage2b_evidence', 'old_mapping_status', 'new_mapping_status', 'stage2b_mapping_reason']
    write(stage / '03_external_taxa_to_gtdb_r226_stage2b.tsv', result, fields)
    return result


def is_accepted(row):
    return row['new_mapping_status'] in ACCEPTED and row['gtdb_r226_species'] != NA


def rebuild_edges(edges, taxa, before_pairs, before_audit, stage):
    """Rebuild from all 15534 frozen raw edge records, not the old pair table.
    That table is read only for independent before/after invariant checks.
    """
    lookup = {(r['paper_id'], r['raw_taxon']): r for r in taxa}
    groups, audit = defaultdict(list), []
    for raw, before in zip(edges, before_audit):
        endpoints = [lookup[(raw['paper_id'], raw[f])] for f in ['taxon1_raw', 'taxon2_raw']]
        failures = [r for r in endpoints if not is_accepted(r)]
        missing = any(r['new_mapping_status'].startswith('pending_missing_') for r in failures)
        sgb_failure = any(r['paper_id'] == 'Lacruz2026' and r['lacruz_taxon_class'] == 'unnamed_SGB' for r in failures)
        ncbi_failure = any(r['ncbi_stage_status'] == 'unresolved' and r['lacruz_taxon_class'] != 'unnamed_SGB' for r in failures)
        row = dict(raw, stage2_primary_fate=before['primary_fate'], endpoint1_taxon_id=endpoints[0]['raw_taxon_id'], endpoint2_taxon_id=endpoints[1]['raw_taxon_id'],
                   endpoint1_mapping_status=endpoints[0]['new_mapping_status'], endpoint2_mapping_status=endpoints[1]['new_mapping_status'],
                   endpoint1_mapping_reason=endpoints[0]['stage2b_mapping_reason'], endpoint2_mapping_reason=endpoints[1]['stage2b_mapping_reason'],
                   endpoint1_gtdb_species=endpoints[0]['gtdb_r226_species'], endpoint2_gtdb_species=endpoints[1]['gtdb_r226_species'],
                   endpoint1_stage2b_changed=endpoints[0]['stage2b_changed'], endpoint2_stage2b_changed=endpoints[1]['stage2b_changed'],
                   both_endpoints_accepted=not failures, newly_fully_mapped_record=not failures and before['both_endpoints_accepted'] != 'True',
                   any_missing_reference=missing, any_unresolved_sgb=sgb_failure, any_unresolved_ncbi=ncbi_failure,
                   species_a=NA, species_b=NA, pair_key=NA, primary_fate=NA, representative_record_id=NA,
                   collapsed_pair_raw_record_count=NA, collapsed_pair_sign_status=NA, sign_concordance_eligible=False)
        if missing:
            row['primary_fate'] = 'pending_missing_reference'
        elif sgb_failure:
            row['primary_fate'] = 'unresolved_sgb'
        elif ncbi_failure:
            row['primary_fate'] = 'unresolved_ncbi'
        elif failures:
            row['primary_fate'] = 'unresolved_gtdb'
        else:
            a, b = sorted(r['gtdb_r226_species'] for r in endpoints)
            check(' || ' not in a and ' || ' not in b, 'Pair delimiter appears in species name')
            row.update(species_a=a, species_b=b, pair_key=a + ' || ' + b)
            if a == b:
                row['primary_fate'] = 'self_loop_after_mapping'
            else:
                groups[(raw['network_id'], a, b)].append(row)
        audit.append(row)
    final = []
    for (network, a, b), members in sorted(groups.items()):
        first = members[0]
        signs = {r['sign'] for r in members}
        conflict = len(signs) > 1
        check(len({r['primary_status'] for r in members}) == 1, 'Collapse mixed primary/exploratory records')
        final.append(dict(paper_id=first['paper_id'], network_id=network, condition=first['condition'], primary_status=first['primary_status'],
                          source_file=first['source_file'], species_a=a, species_b=b, pair_key=first['pair_key'], raw_record_count=len(members),
                          sign=NA if conflict else next(iter(signs)), sign_status='conflict' if conflict else 'consistent', sign_values=';'.join(sorted(signs)),
                          presence_comparison_eligible=True, sign_concordance_eligible=not conflict, primary_paper_support_eligible=first['primary_status'] == 'primary',
                          representative_record_id=first['raw_record_id'], raw_record_ids=j([r['raw_record_id'] for r in members])))
        for index, row in enumerate(members):
            row.update(primary_fate='duplicate_after_mapping' if index else 'kept_unique_gtdb_pair',
                       representative_record_id=first['raw_record_id'], collapsed_pair_raw_record_count=len(members),
                       collapsed_pair_sign_status='conflict' if conflict else 'consistent', sign_concordance_eligible=not conflict)
    check(len(audit) == 15534 and len({r['raw_record_id'] for r in audit}) == 15534, 'Raw-record audit incomplete')
    check(all(r['primary_fate'] in FATES for r in audit), 'Unassigned/invalid Stage 2B primary fate')
    check(Counter(r['primary_status'] for r in audit) == {'primary': 14836, 'exploratory': 698}, 'Primary/exploratory state changed')
    old_lookup = {(r['network_id'], r['pair_key']): r for r in before_pairs}
    new_lookup = {(r['network_id'], r['pair_key']): r for r in final}
    check(len(new_lookup) == len(final), 'Duplicate final pair within a network')
    check(set(old_lookup).issubset(new_lookup), 'An accepted Stage 2 pair was lost')
    for identity, old in old_lookup.items():
        if old['paper_id'] == 'Guo2023':
            check(all(str(new_lookup[identity][f]) == value for f, value in old.items()), 'Guo final pair mapping changed')
    for network, before_count in BASELINE_PAIRS.items():
        records = [r for r in audit if r['network_id'] == network]
        pairs = [r for r in final if r['network_id'] == network]
        counts = Counter(r['primary_fate'] for r in records)
        n_raw = len(records)
        unresolved = sum(counts[f] for f in FATES[:4])
        fully = sum(r['both_endpoints_accepted'] for r in records)
        row = dict(network_id=network, paper_id=records[0]['paper_id'], primary_status=records[0]['primary_status'],
                   raw_edges=n_raw, raw_unique_taxa=len({r[f] for r in records for f in ['taxon1_raw', 'taxon2_raw']}),
                   fully_mapped_edge_records=fully, **{f: counts[f] for f in FATES},
                   unresolved_mapping_records=unresolved, self_loops_after_gtdb=counts['self_loop_after_mapping'],
                   duplicate_records_collapsed=counts['duplicate_after_mapping'], final_unique_gtdb_pairs=len(pairs),
                   pair_mapping_coverage=len(pairs)/n_raw, sign_conflict_pairs=sum(r['sign_status'] == 'conflict' for r in pairs),
                   newly_fully_mapped_raw_records=sum(r['newly_fully_mapped_record'] for r in records),
                   edge_fate_balance_ok=sum(counts.values()) == n_raw == unresolved + counts['self_loop_after_mapping'] + counts['duplicate_after_mapping'] + len(pairs))
        check(row['edge_fate_balance_ok'] and fully == n_raw - unresolved, 'Network edge-fate balance failed: ' + network)
        check(counts['kept_unique_gtdb_pair'] == len(pairs), 'Representative/final pair count mismatch')
        check(sum(r['raw_record_count'] for r in pairs) == fully - counts['self_loop_after_mapping'], 'Collapse raw-record conservation failed')
        STATE['qc'].append(row)
        STATE['summary'].append(dict(network=network, raw_edges=n_raw, primary_status=row['primary_status'],
                                     Stage2_final_pairs=before_count, Stage2_coverage=before_count/n_raw,
                                     Stage2B_final_pairs=len(pairs), Stage2B_coverage=len(pairs)/n_raw,
                                     absolute_pair_gain=len(pairs)-before_count, percentage_point_gain=100*(len(pairs)-before_count)/n_raw))
        LOG.info('%s: pairs %d -> %d; gains=%d; fates=%s; accounting=PASS', network, before_count, len(pairs), len(pairs)-before_count, j(counts))
    check(len(STATE['qc']) == 8 and sum(r['raw_edges'] for r in STATE['qc']) == 15534, 'Full-network reconciliation failed')
    write(stage / '04_all_external_edges_gtdb_r226_stage2b.tsv', final)
    write(stage / '04_edge_mapping_audit_stage2b.tsv', audit)
    write(stage / 'stage2_vs_stage2b_network_summary.tsv', STATE['summary'])
    write(REPORTS / '02b_taxonomy_rescue_qc.tsv', STATE['qc'])
    STATE['fates'] = dict(Counter(r['primary_fate'] for r in audit))
    STATE['final_pairs'] = len(final)
    STATE['zhou_rescued_raw_records'] = sum(r['newly_fully_mapped_record'] for r in audit if r['paper_id'] == 'Zhou2022')
    return final, audit


def loss_audits(taxa, audit, stage):
    occurrence, network_counts = Counter(), defaultdict(Counter)
    for edge in audit:
        for raw in {edge['taxon1_raw'], edge['taxon2_raw']}:
            identity = (edge['paper_id'], raw)
            occurrence[identity] += 1
            network_counts[identity][edge['network_id']] += 1
    unresolved = []
    for row in taxa:
        if is_accepted(row):
            continue
        identity = (row['paper_id'], row['raw_taxon'])
        unresolved.append(dict(paper_id=row['paper_id'], raw_taxon=row['raw_taxon'], affected_edges=occurrence[identity],
                               per_network_affected_edges=j(network_counts[identity]), ncbi_taxid=row['ncbi_taxid'], ncbi_rank=row['ncbi_rank'],
                               old_mapping_status=row['old_mapping_status'], mapping_status=row['new_mapping_status'],
                               mapping_reason=row['stage2b_mapping_reason'], stage2b_changed=row['stage2b_changed']))
    paper_order = {'Zhou2022': 0, 'Guo2023': 1, 'Lacruz2026': 2}
    unresolved.sort(key=lambda r: (paper_order[r['paper_id']], -r['affected_edges'], r['raw_taxon']))
    write(stage / 'stage2b_unresolved_taxa.tsv', unresolved)
    STATE['remaining_unresolved_by_paper'] = dict(Counter(r['paper_id'] for r in unresolved))
    STATE['remaining_zhou_top'] = [r for r in unresolved if r['paper_id'] == 'Zhou2022'][:15]
    STATE['remaining_zhou_reasons'] = dict(Counter(r['mapping_status'] for r in unresolved if r['paper_id'] == 'Zhou2022'))
    pending = [r for r in audit if r['raw_record_id'] in STATE['pending_record_ids']]
    check(len(pending) == 1241, 'Frozen pending-record cohort changed')
    rows = []
    for network in ['Lacruz2026_' + g for g in GROUPS] + ['ALL', 'ALL_PRIMARY', 'ALL_EXPLORATORY']:
        subset = [r for r in pending if r['network_id'] == network] if network.startswith('Lacruz') else \
                 pending if network == 'ALL' else [r for r in pending if r['primary_status'] == ('primary' if network == 'ALL_PRIMARY' else 'exploratory')]
        c = Counter(r['primary_fate'] for r in subset)
        unresolved_n = sum(c[f] for f in FATES[:4])
        rows.append(dict(network=network, stage2_pending_records=len(subset),
                         kept_unique_gtdb_pair=c['kept_unique_gtdb_pair'], duplicate_after_mapping=c['duplicate_after_mapping'],
                         self_loop_after_mapping=c['self_loop_after_mapping'], still_unresolved=unresolved_n,
                         fate_balance_ok=len(subset) == c['kept_unique_gtdb_pair'] + c['duplicate_after_mapping'] + c['self_loop_after_mapping'] + unresolved_n))
    check(all(r['fate_balance_ok'] for r in rows), 'Pending-SGB raw-record cohort accounting failed')
    check(next(r for r in rows if r['network'] == 'ALL_PRIMARY')['stage2_pending_records'] == 874 and
          next(r for r in rows if r['network'] == 'ALL_EXPLORATORY')['stage2_pending_records'] == 367, 'Pending cohort status changed')
    write(stage / 'lacruz_pending_sgb_edge_fate_summary.tsv', rows)
    STATE['pending_cohort_fates'] = rows
    STATE['accepted_by_paper'] = dict(Counter(r['paper_id'] for r in taxa if is_accepted(r)))


def table(rows, columns):
    def cell(value):
        return str(value).replace('|', '\\|').replace('\n', ' ')
    return '\n'.join(['| ' + ' | '.join(columns) + ' |', '| ' + ' | '.join(['---'] * len(columns)) + ' |'] +
                     ['| ' + ' | '.join(cell(r.get(c, NA)) for c in columns) + ' |' for r in rows])


def report():
    decisions = STATE.get('zhou_rescue_decisions', [])
    rescued = [r for r in decisions if r['decision'] == 'accepted_strict_rescue']
    lines = ['# External Validation Stage 2B taxonomy rescue QC', '',
             'overall_status = **{}**；UTC：{}。'.format(STATE['overall_status'], datetime.now(timezone.utc).isoformat()), '',
             '本阶段仅审查 Zhou 前 30 个 high-degree unresolved taxa、处理 Lacruz 138 个 unnamed SGB，'
             '并从全部 15534 条冻结 standardized raw edge records 重建网络和损失对账。'
             '未联网、下载、修改 Stage 1/2、查询 MACK-DB、做 overlap/source exclusion/sign concordance with MACK-DB/external recurrence。', '',
             '复现：`python3 stage2b/02b_taxonomy_rescue.py`。所有固定文件及冻结产物的 SHA256 见 `02b_input_inventory.tsv`；'
             '前后哈希核对见 `02b_input_integrity_qc.tsv`。', '',
             '## 固定 reference 与 Chen 规则', '',
             table([dict(file=str(p.relative_to(BASE)), expected_sha256=d) for p, d in HASHES.items()], ['file', 'expected_sha256']), '',
             'SGB2GTDB 同时核对了 wheel 整体 SHA256 及 wheel 内部成员：`{}`。'.format(STATE.get('wheel_member', NA)), '',
             '沿用既有 Chen scientific/synonym/equivalent name、merged TaxID、稳定 historical TaxID、domain-aware official 优先、'
             'official ambiguity gate 和 unique-highest genome fallback 规则。'
             '完整 official forward 表沿用冻结 Stage 2 使用的 `official_mapping_v1/01_official_correspondence_all_parsed.tsv`，'
             '用于区分缺失与 top tie；参考脚本只读取，不执行其中的 overlap 代码。没有 fuzzy matching 或 one-to-many expansion。', '',
             '## A. Zhou high-degree unresolved audit', '',
             '重新从 raw edges 计算全部 148 个 unresolved taxa 的 MS/HHC degree 与 affected records。'
             'degree 是该网络内独立 raw 邻居数；affected_edges 是独立 raw_record_id 数；'
             'total_unique_affected_edges 是 MS/HHC record IDs 的并集数。同一个 pair 在不同 network 的 records 保持独立。', '',
             '148 taxa 的 total affected-record 分布：`{}`。'.format(j(STATE.get('zhou_degree_distribution', {}))), '',
             '完整分布及排名见 `zhou_unresolved_degree_audit.tsv`；审查范围固定为下列 top 30：', '',
             table(STATE.get('zhou_top30', []), ['degree_rank', 'raw_taxon', 'MS_degree', 'HHC_degree', 'MS_affected_edges',
                                               'HHC_affected_edges', 'total_unique_affected_edges', 'old_gtdb_status']), '',
             '方括号 genus、空白及 underscore 格式变化只用于生成 candidates。接受必须有 NCBI accepted-name '
             '或既有 Chen stable species TaxID 的唯一 identity 证据；candidate 自身的字符串相似性不是证据。'
             '已确认 NCBI identity 的 official/genome top ties 保持 unresolved，不通过别名改走 fallback。', '',
             '新增 strict NCBI identity rescue taxa：**{}**；其中新增 accepted GTDB mapping：**{}**。'
             '原来未 fully mapped、现在两端都 accepted 的 Zhou raw edge records：**{}**（包含随后 self-loop/collapse 的 records）。'.format(
                 STATE.get('zhou_strict_identity_rescues', 0), STATE.get('zhou_new_gtdb_mapped', 0), STATE.get('zhou_rescued_raw_records', 0)), '',
             table(rescued, ['raw_taxon', 'candidate_name', 'candidate_taxid', 'candidate_rank', 'new_status', 'final_gtdb_r226_species']), '']
    if not rescued and 'zhou_rescue_decisions' in STATE:
        lines += ['No additional strict Zhou taxonomy rescue was supported.', '']
    lines += ['Zhou remaining unresolved taxa：**{}**；原因分布：`{}`。'.format(
                  STATE.get('remaining_unresolved_by_paper', {}).get('Zhou2022', NA), j(STATE.get('remaining_zhou_reasons', {}))), '',
              table(STATE.get('remaining_zhou_top', []), ['raw_taxon', 'affected_edges', 'mapping_status', 'mapping_reason']), '',
              '## B. Lacruz SGB schema、rescue 与 genome bridge', '',
              'species.txt 的前 20 行已写入日志。实际 schema（完整内容另存 `metaphlan_species_file_schema.txt`）：', '',
              '`' + j(STATE.get('species_schema', {})) + '`', '',
              'species-file 匹配及一级 route：`{}`。'.format(j(STATE.get('speciesfile_counts', {}))), '',
              '每条完整匹配记录保存在 `lacruz_sgb_speciesfile_audit.tsv`。若存在多个可接受 named species，'
              '只在其最终 GTDB R226 targets 全部相同时接受；指向多个 species 时 unresolved，不任选其中一个。', '',
              'SGB2GTDB route：`{}`。'.format(j(STATE.get('sgb2gtdb_counts', {}))), '',
              'headerless SGB2GTDB 文件按实际两列解析。空 s__ 字段为 unresolved_no_r220_species，'
              '不从 genus 猜 species。所有 R220 名称只用于选择同一 R220 metadata 中的 genomes，'
              '不与 R226 名称做字符串匹配。', '',
              'R220→R226 bridge 唯一 R220 species 数：{}；bridge 决策分布：`{}`。'.format(
                  STATE.get('bridge_species_count', 0), j(STATE.get('bridge_status_counts', {}))), '',
              'Bridge 分类计数：mapped unique-highest={}；top-tie={}；no-shared-genomes={}；'
              'no-R220-genomes={}；missing-R226-Archaea={}。'.format(*[
                  STATE.get('bridge_status_counts', {}).get(status, 0) for status in
                  ['mapped_r220_r226_bridge', 'unresolved_r220_r226_top_tie', 'unresolved_no_shared_genomes',
                   'unresolved_no_r220_genomes', 'pending_missing_r226_archaea_metadata']]), '',
              'accession normalization 仅删除开头 RS_ 或 GB_，保留 GCF_/GCA_ 和 version。'
              '对每个相关 R220 species 使用其全部 genomes，通过 normalized exact accession inner join R226，'
              '按 shared genomes 统计 R226 species support。unique-highest 才接受，不设置 80%（或 50%）硬阈值。'
              'top proportion 的分母为所有 shared genomes；另记录 classified_shared_genomes。'
              '完整逐 genome 连接证据见 `gtdb_r220_r226_genome_join_audit.tsv`，全部候选支持数另有 audit。', '',
              'archaeal_pending_sgb = **{}**；本地 R226 archaeal metadata 存在：**{}**。'.format(
                  STATE.get('archaeal_pending_sgb', NA), STATE.get('r226_archaea_metadata_exists', False)), '']
    if STATE.get('archaeal_pending_sgb') == 0:
        lines += ['本次 pending SGB 不含 Archaea，无需 archaeal bridge，未标记缺失 archaeal reference。', '']
    lines += ['最终 SGB rescue：**{} / 138**；unresolved：**{}**；rescue rate：**{:.2%}**。'
              '状态分布：`{}`。'.format(STATE.get('sgb_rescued', 0), STATE.get('sgb_unresolved', 138),
                                       STATE.get('sgb_rescue_rate', 0), j(STATE.get('sgb_final_counts', {}))), '',
              '## C. Stage 2 before / Stage 2B after', '',
              table(STATE['summary'], ['network', 'primary_status', 'raw_edges', 'Stage2_final_pairs', 'Stage2_coverage',
                                       'Stage2B_final_pairs', 'Stage2B_coverage', 'absolute_pair_gain', 'percentage_point_gain']), '',
              'Coverage 的分母始终是 Stage 1 固定 raw_edges；percentage_point_gain 为百分点变化。'
              'Stage 2 baseline=10016 network-level unique pairs；Stage 2B={}。不跨 network 去重。'.format(STATE.get('final_pairs', NA)), '',
              'Lacruz MHO 始终 exploratory，primary_paper_support_eligible=False。'
              'sign 全同才保留；冲突设 sign_status=conflict、sign=NA、sign_concordance_eligible=False，presence 保留。', '',
              '## D. 全量互斥 raw-edge fate 与损失对账', '',
              '每条原始记录必须且只能有一个 primary fate。失败优先级为 pending_missing_reference → unresolved_sgb → '
              'unresolved_ncbi → unresolved_gtdb；两个 endpoints 都成功才进入 self-loop 或 network 内 pair collapse。'
              '每个 pair 以输入次序首条 record 作代表，完整保留 member record IDs 和 raw_record_count。', '',
              table(STATE['qc'], ['network_id', 'raw_edges'] + FATES + ['edge_fate_balance_ok']), '',
              '全量 fate 总数：`{}`；总和={}。'.format(j(STATE.get('fates', {})), sum(STATE.get('fates', {}).values())), '',
              'self-loop 和 duplicate collapse 单独列出，不算 taxonomy failure。原 accepted mapping 与既有 pairs 均不能减少；'
              'Guo taxon mapping 与最终 pair 行逐字段保持一致。', '',
              'Stage 2 的 1241 条 pending-SGB raw records 在 Stage 2B 的去向：', '',
              table(STATE.get('pending_cohort_fates', []), ['network', 'stage2_pending_records', 'kept_unique_gtdb_pair',
                                                         'duplicate_after_mapping', 'self_loop_after_mapping', 'still_unresolved', 'fate_balance_ok']), '',
              '上述表为固定 raw-record cohort 的 fate，不能等同于新增 unique pairs；'
              '新增 pairs 以 before/after network summary 的差值为准。', '',
              '## 不变性与自动 QC', '',
              '所有冻结产物与 reference 前后 SHA256 一致：**{}**。'.format(STATE.get('inputs_unchanged', False)), '',
              'Stage 2B taxon table 中，只有 accepted Zhou identity rescue 或 accepted SGB rescue 可改变原有字段。'
              'stage2b_changed 指原有 core mapping 是否改变；未 rescue 的 SGB 保留原始 pending core 字段，'
              '新增 new_mapping_status / stage2b_mapping_reason 记录本次明确的 unresolved 原因，供重建 edge fate 使用。'
              'Guo 与原有 accepted mappings 的原字段逐一检查不变。', '',
              '自动断言：15534 records、Guo 不变、初始 pending SGB=138、accepted mapping 不减少、'
              '每条 edge 一个 fate、8/8 network 对账、primary/exploratory 不变、network 内独立 collapse。'
              '所有输出 TSV 写出后回读逐字段检查。失败时非零退出，不发布 staging 中的 Stage 2B 结果。', '']
    lines += ['- ' + w for w in STATE['warnings']] or ['未发现需要警告的外部资源缺失。']
    if 'error' in STATE:
        lines += ['', '**停止原因：** ' + STATE['error']]
    lines += ['', 'Stage 2B 到此停止；不开展任何下游 overlap、source exclusion 或 external recurrence。', '']
    (REPORTS / '02b_taxonomy_rescue_qc.md').write_text('\n'.join(lines), encoding='utf-8')


def main():
    if sys.flags.optimize:
        raise RuntimeError('Required assertions must not be disabled with -O/-OO')
    OUT.mkdir(exist_ok=True)
    REPORTS.mkdir(exist_ok=True)
    LOG.setLevel(logging.INFO)
    for handler in [logging.FileHandler(REPORTS / '02b_taxonomy_rescue.log', mode='w', encoding='utf-8'), logging.StreamHandler(sys.stdout)]:
        handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))
        LOG.addHandler(handler)
    LOG.info('START offline Stage 2B; frozen Stage 1/2 remain read-only')
    try:
        framework = load_framework()
        snapshot(framework)
        edges, taxa, before_pairs, before_audit = fixed_inputs(framework)
        with tempfile.TemporaryDirectory(prefix='.stage2b-', dir=str(OUT)) as temp:
            stage = Path(temp)
            top30 = zhou_degrees(edges, taxa, stage)
            pending = pending_manifest(edges, taxa, stage)
            species_audit = species_schema_and_matches(pending, stage)
            proposals, species_decisions = name_rescues(taxa, top30, species_audit, framework, stage)
            assignments = sgb_r220_assignments(pending, species_decisions, stage)
            bridges = build_bridge(assignments, stage)
            sgb_final = final_sgb_mapping(pending, species_audit, species_decisions, assignments, bridges, stage)
            new_taxa = overlay_mapping(taxa, proposals, sgb_final, species_decisions, stage)
            final, audit = rebuild_edges(edges, new_taxa, before_pairs, before_audit, stage)
            loss_audits(new_taxa, audit, stage)
            verify_inputs()
            for path in sorted(stage.iterdir()):
                os.replace(str(path), str(OUT / path.name))
        STATE['overall_status'] = 'PASS_WITH_WARNING' if STATE['warnings'] else 'PASS'
        report()
        LOG.info('DONE overall_status=%s; final_pairs=%d; fates=%s; Stage 2B stopped', STATE['overall_status'], STATE['final_pairs'], j(STATE['fates']))
        return 0
    except Exception as error:
        STATE.update(overall_status='FAIL', error='{}: {}'.format(type(error).__name__, error))
        LOG.exception('STOP: %s', STATE['error'])
        if STATE.get('inventory'):
            try:
                verify_inputs()
            except Exception as integrity_error:
                STATE['warnings'].append('Integrity verification failed: ' + str(integrity_error))
        report()
        return 1


if __name__ == '__main__':
    sys.exit(main())
