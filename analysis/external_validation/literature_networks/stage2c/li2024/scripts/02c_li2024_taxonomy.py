#!/usr/bin/env python3
"""Offline Li 2024 Stage 2C. Frozen Stage 2/2B definitions are imported read-only.

Run with python3 -s. Never calls either frozen main(), parses earlier raw
Excel/CSV, reads human_gut.info, or performs any downstream overlap analysis.
NCBI/official/genome mapping uses the frozen implementation. Pair conversion
executes the unchanged computation prefix of Stage 2B rebuild_edges(), omitting
only its old-network-specific assertions, summaries and file writes.
"""

# v0.1: repository-local IO configuration; scientific functions below are unchanged.
import sys as _path_sys
from pathlib import Path as _RepoPath
_repo_root = next(p for p in _RepoPath(__file__).resolve().parents if (p / "lib/mackdb_paths.py").is_file())
_path_sys.path.insert(0, str(_repo_root / "lib"))
from mackdb_paths import resolve_path
import ast
import csv
import hashlib
import importlib.util
import json
import logging
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.dont_write_bytecode = True
BASE = Path(str(resolve_path('literature')))
OUT = BASE / 'stage2c/li2024'
RAW = OUT / 'raw/sparcc_all_species.txt'
REPORT = BASE / 'reports/02c_li2024_taxonomy_qc.md'
PRIOR = BASE / 'stage2b/03_external_taxa_to_gtdb_r226_stage2b.tsv'
NA = 'NA'
QC, INVENTORY, STATS = [], [], {}
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')
LOG = logging.getLogger('stage2c')


def j(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'))


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    with path.open(encoding='utf-8', newline='') as f:
        r = csv.DictReader(f, delimiter='\t')
        return r.fieldnames, list(r)


def write(path, rows, fields=None):
    path = path.resolve()
    assert OUT in path.parents and OUT / 'raw' not in path.parents, path
    fields = fields or list(dict.fromkeys(k for r in rows for k in r))
    assert fields, path
    path.parent.mkdir(parents=True, exist_ok=True)
    encoded = [{k: NA if r.get(k) is None else str(r[k]) for k in fields} for r in rows]
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(encoded)
    assert read(path) == (fields, encoded), 'TSV round-trip: ' + str(path)


def check(name, value, detail=''):
    QC.append(dict(check=name, status='PASS' if value else 'FAIL', detail=detail))
    if not value:
        raise AssertionError(name + ': ' + str(detail))


def load_frameworks():
    path = BASE / 'stage2b/02b_taxonomy_rescue.py'
    spec = importlib.util.spec_from_file_location('frozen_stage2b_li_adapter', path)
    b = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(b)
    b.LOG = LOG
    f = b.load_framework()
    # Keep the actual, unmodified Stage 2B computational statements. The first
    # old-dataset assertion is the boundary before side effects and summaries.
    tree = ast.parse(path.read_text(encoding='utf-8'))
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'rebuild_edges')
    boundary = next(i for i, n in enumerate(node.body)
                    if isinstance(n, ast.Expr) and isinstance(n.value, ast.Call)
                    and isinstance(n.value.func, ast.Name) and n.value.func.id == 'check'
                    and len(n.value.args) > 1 and isinstance(n.value.args[1], ast.Constant)
                    and n.value.args[1].value == 'Raw-record audit incomplete')
    prefix = node.body[:boundary]
    check('frozen_pair_core_boundary_verified', boundary > 5 and
          prefix[-1].end_lineno == node.body[boundary].lineno - 1,
          '{}:{}-{}'.format(path, prefix[0].lineno, prefix[-1].end_lineno))
    STATS['pair_core_lines'] = '{}-{}'.format(prefix[0].lineno, prefix[-1].end_lineno)
    STATS['pair_core_ast_sha256'] = hashlib.sha256(ast.dump(ast.Module(body=prefix, type_ignores=[])).encode()).hexdigest()
    node.name = 'li_pair_core'
    node.body = prefix + [ast.Return(value=ast.Tuple(elts=[ast.Name(id='final', ctx=ast.Load()),
                                                         ast.Name(id='audit', ctx=ast.Load())], ctx=ast.Load()))]
    module = ast.fix_missing_locations(ast.Module(body=[node], type_ignores=[]))
    namespace = dict(b.__dict__)
    exec(compile(module, str(path) + ':readonly_Li_computation_adapter', 'exec'), namespace)
    return f, b, namespace['li_pair_core']


def snapshot():
    expected = {}
    for name in ['01_input_inventory.tsv', '02_input_inventory.tsv', '02b_input_inventory.tsv']:
        for r in read(BASE / 'reports' / name)[1]:
            if r.get('sha256') not in (None, '', NA):
                p = Path(r['absolute_path'])
                check('prior_manifest_hash_consistency:' + p.name,
                      p not in expected or expected[p] == r['sha256'])
                expected[p] = r['sha256']
    paths = set(expected) | {RAW}
    for name in ['parsed', 'mapping', 'scripts', 'stage2b', 'raw', 'reference', 'reports']:
        paths.update(p for p in (BASE / name).rglob('*') if p.is_file() and p != REPORT)
    for p in sorted(paths):
        digest = sha(p)
        scope = 'li2024_raw' if p == RAW else 'frozen_stage1_2_2b_or_reference'
        INVENTORY.append(dict(absolute_path=str(p), scope=scope, size_bytes=p.stat().st_size,
                              sha256_before=digest, prior_manifest_sha256=expected.get(p, NA),
                              matches_prior_manifest=p not in expected or digest == expected[p],
                              sha256_after=NA, unchanged='PENDING'))
    write(OUT / 'qc/li2024_input_integrity.tsv', INVENTORY)
    check('frozen_baseline_matches_prior_manifests', all(r['matches_prior_manifest'] for r in INVENTORY))
    STATS['raw_sha256'] = next(r['sha256_before'] for r in INVENTORY if r['scope'] == 'li2024_raw')
    LOG.info('Read-only snapshot: %d files', len(INVENTORY))


def input_audit():
    flags, v1, v2, missing, seen = Counter(), set(), set(), Counter(), Counter()
    published, self_ids = [], []
    missing_tokens = {'', 'na', 'nan', 'null', 'none', 'n/a'}
    with RAW.open(encoding='utf-8', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fields = reader.fieldnames
        check('raw_schema', fields == ['Var1', 'Var2', 'cor', 'p', 'fdr', 'signature_species_in_figure'], j(fields))
        for row_number, r in enumerate(reader, 2):
            assert None not in r and all(v is not None for v in r.values()), 'Malformed source row'
            flags[r['signature_species_in_figure']] += 1
            v1.add(r['Var1']); v2.add(r['Var2'])
            missing.update(k for k in fields if r[k].strip().lower() in missing_tokens)
            seen[tuple(r[k] for k in fields)] += 1
            if r['Var1'] == r['Var2']:
                self_ids.append(row_number)
            if r['signature_species_in_figure'] == 'Yes':
                published.append(dict(source_row_id=row_number, **r))
    with RAW.open('rb') as f:
        line_count = sum(1 for _ in f)
    taxa = sorted({r[k] for r in published for k in ['Var1', 'Var2']})
    stats = dict(sha256=STATS['raw_sha256'], total_lines_including_header=line_count,
                 total_data_records=sum(flags.values()), field_names=j(fields), unique_Var1=len(v1),
                 unique_Var2=len(v2), all_unique_taxa=len(v1 | v2),
                 signature_species_in_figure_counts=j(flags), missing_value_definition=j(sorted(missing_tokens)),
                 missing_cells=sum(missing.values()), missing_cells_by_field=j({k: missing[k] for k in fields}),
                 raw_self_pair_records=len(self_ids), raw_self_pair_source_rows=j(self_ids),
                 exact_duplicate_groups=sum(n > 1 for n in seen.values()),
                 exact_duplicate_records_excess=sum(n - 1 for n in seen.values()),
                 source_row_id_definition='Physical source file line number, 1-based; header is line 1')
    write(OUT / 'qc/li2024_raw_input_audit.tsv', [dict(metric=k, value=v) for k, v in stats.items()])
    STATS.update(raw_audit=stats, published_edges=len(published), published_taxa=len(taxa))
    write(OUT / 'processed/li2024_published_edges_raw.tsv', published, ['source_row_id'] + fields)
    check('published_edges_exact_742', len(published) == 742, len(published))
    check('published_unique_taxa_exact_60', len(taxa) == 60, len(taxa))
    check('source_row_ids_physical_lines', line_count == sum(flags.values()) + 1, line_count)
    check('published_required_fields_not_missing', all(r[k].strip().lower() not in missing_tokens for r in published for k in fields))
    for r in published:
        for k in ['cor', 'p', 'fdr']:
            assert math.isfinite(float(r[k])), (r['source_row_id'], k)
    signs = Counter(sign(r['cor']) for r in published)
    STATS['published_signs'] = {k: signs[k] for k in ['positive', 'negative', 'zero']}
    check('published_signs_partition', sum(signs.values()) == len(published), j(STATS['published_signs']))
    occurrences = Counter(r[k] for r in published for k in ['Var1', 'Var2'])
    write(OUT / 'processed/li2024_published_taxa_unique.tsv',
          [dict(original_taxon=t, endpoint_occurrences=occurrences[t],
                incident_raw_edges=sum(t in (r['Var1'], r['Var2']) for r in published)) for t in taxa])
    LOG.info('Published extraction: %d records / %d taxa / signs %s', len(published), len(taxa), signs)
    return published, taxa


def sign(value):
    x = float(value)
    return 'positive' if x > 0 else 'negative' if x < 0 else 'zero'


def harmonize(taxa, f, b):
    old_fields, all_prior = read(PRIOR)
    prior = [r for r in all_prior if b.is_accepted(r)]
    by_name, by_id = defaultdict(list), defaultdict(list)
    for r in prior:
        # No bracket removal, punctuation guessing, GTDB-label lookup, or SGB
        # label-derived species identity is used for direct name reuse.
        names = {r['raw_taxon'], r['normalized_candidate'], r['ncbi_name']}
        for n in names - {NA, ''}:
            by_name[f.key(n)].append(r)
        if r['ncbi_stage_status'] == 'mapped_species' and r['ncbi_rank'] == 'species':
            by_id[r['ncbi_taxid']].append(r)
    chosen, routes, reuse_refs, novel = {}, {}, {}, []
    for i, taxon in enumerate(taxa, 1):
        hits = list({r['raw_taxon_id']: r for r in by_name[f.key(taxon)]}.values())
        if hits:
            states = {(r['ncbi_taxid'], r['ncbi_rank'], r['gtdb_r226_species']) for r in hits}
            if len(states) == 1:
                chosen[taxon] = dict(hits[0])
                routes[taxon], reuse_refs[taxon] = 'reused_accepted_name', hits
                continue
            chosen[taxon] = dict(ncbi_stage_status='unresolved', ncbi_taxid=NA, ncbi_rank=NA,
                                 gtdb_r226_species=NA, mapping_status='unresolved_conflicting_evidence',
                                 mapping_reason='Conflicting accepted mappings for deterministic normalized label',
                                 ncbi_evidence=j([r['raw_taxon_id'] for r in hits]))
            routes[taxon] = 'conflicting_prior_evidence'
            continue
        novel.append(dict(raw_taxon_id='LI_T{:03d}'.format(i), paper_id='Li2024', raw_taxon=taxon,
                          normalized_candidate=re.sub(r'\s+', ' ', taxon.replace('_', ' ')).strip(),
                          lacruz_taxon_class='not_applicable'))
    LOG.info('Direct accepted-name reuse: %d; taxa requiring identity resolution: %d', len(reuse_refs), len(novel))
    evidence, to_gtdb = [], []
    if novel:
        # Li labels contain no bracketed genus. The existing deterministic
        # Stage 2B candidate generator yields only original/underscore forms;
        # Stage 2's formatted candidate already covers those identical keys.
        check('novel_candidates_use_existing_format_rules', all(
            {f.key(n) for n, _ in b.zhou_candidates(r['raw_taxon'])} == {f.key(r['normalized_candidate'])}
            for r in novel), 'No new normalization or fuzzy rule')
        taxonomy = f.Taxonomy(novel)
        for q in novel:
            row, hits = taxonomy.resolve(q)
            evidence.extend(hits)
            taxon = q['raw_taxon']
            previous = by_id.get(str(row['ncbi_taxid']), []) if row['ncbi_stage_status'] == 'mapped_species' else []
            targets = {r['gtdb_r226_species'] for r in previous}
            if previous and len(targets) == 1:
                # Keep freshly verified Li identity evidence; reuse every GTDB
                # decision/support field without generating a different target.
                old = previous[0]
                for k in old_fields[old_fields.index('gtdb_r226_species'):old_fields.index('stage2b_changed')]:
                    row[k] = old[k]
                reuse_refs[taxon] = previous
                routes[taxon] = 'reused_accepted_current_ncbi_taxid'
            elif len(targets) > 1:
                row.update(mapping_status='unresolved_conflicting_evidence',
                           mapping_reason='Previously accepted mappings conflict at confirmed current NCBI TaxID')
                routes[taxon] = 'conflicting_prior_evidence'
            else:
                routes[taxon] = 'frozen_stage2_ncbi_official_genome'
                to_gtdb.append(row)
            chosen[taxon] = row
        official, forward = f.load_official()
        f.apply_official(to_gtdb, official, forward)
        if any(r['fallback_eligible'] for r in to_gtdb):
            # Keep all fallback evidence local, with the unchanged evidence
            # writer protected by this Stage 2C-only output guard.
            f.write_tsv = write
            f.genome_fallback(to_gtdb, taxonomy, OUT / 'taxonomy')
    evidence_fields = ['raw_taxon_id', 'paper_id', 'raw_taxon', 'match_tier', 'matched_ncbi_alias',
                       'name_class', 'alias_taxid', 'current_taxid', 'current_rank', 'current_domain']
    write(OUT / 'taxonomy/li2024_ncbi_candidate_evidence.tsv', evidence, evidence_fields)
    mapping, reuse = [], []
    core_fields = old_fields[old_fields.index('ncbi_name'):old_fields.index('stage2b_changed')]
    for i, taxon in enumerate(taxa, 1):
        r = chosen[taxon]
        references = reuse_refs.get(taxon, [])
        row = dict(original_taxon=taxon, normalized_taxon=re.sub(r'\s+', ' ', taxon.replace('_', ' ')).strip(),
                   ncbi_taxid=r.get('ncbi_taxid', NA), ncbi_rank=r.get('ncbi_rank', NA),
                   mapping_method=routes[taxon] + ';' + r.get('ncbi_mapping_method', NA) + ';' + r.get('gtdb_mapping_method', NA),
                   gtdb_r226_species=r['gtdb_r226_species'], mapping_status=r['mapping_status'],
                   evidence_source=';'.join(dict.fromkeys(([str(PRIOR)] if references else []) +
                        [str(f.TAXDUMP / 'names.dmp'), str(f.TAXDUMP / 'nodes.dmp'), str(f.TAXDUMP / 'merged.dmp'),
                         str(f.REF_NAME), str(f.OFFICIAL), str(f.FORWARD)] +
                        ([r['gtdb_evidence_file']] if r.get('gtdb_evidence_file', NA) != NA else []))),
                   evidence_detail=j(dict(reuse_route=routes[taxon], prior_taxon_ids=[x['raw_taxon_id'] for x in references],
                                          ncbi_evidence=r.get('ncbi_evidence', NA), mapping_reason=r['mapping_reason'],
                                          stage2b_evidence=r.get('stage2b_evidence', NA))),
                   reuse_route=routes[taxon], prior_accepted_taxon_ids=j([x['raw_taxon_id'] for x in references]))
        row.update({k: r.get(k, NA) for k in core_fields})
        row.update(raw_taxon_id='LI_T{:03d}'.format(i), raw_taxon=taxon, paper_id='Li2024',
                   new_mapping_status=r['mapping_status'], stage2b_mapping_reason=r['mapping_reason'],
                   stage2b_changed=False, lacruz_taxon_class='not_applicable')
        for old in references:
            identical = all(str(row[k]) == old[k] for k in ['ncbi_taxid', 'ncbi_rank', 'gtdb_r226_species'])
            reuse.append(dict(original_taxon=taxon, prior_taxon_id=old['raw_taxon_id'], prior_paper_id=old['paper_id'],
                              prior_raw_taxon=old['raw_taxon'], reuse_route=routes[taxon],
                              ncbi_taxid=old['ncbi_taxid'], gtdb_r226_species=old['gtdb_r226_species'], unchanged=identical))
            check('reused_accepted_identity_unchanged:' + taxon + ':' + old['raw_taxon_id'], identical)
        mapping.append(row)
    fields = list(mapping[0])
    write(OUT / 'taxonomy/li2024_taxonomy_mapping.tsv', mapping, fields)
    unresolved = [r for r in mapping if not b.is_accepted(r)]
    write(OUT / 'taxonomy/li2024_unresolved_taxa.tsv', unresolved, fields)
    write(OUT / 'qc/li2024_mapping_reuse_audit.tsv', reuse,
          ['original_taxon', 'prior_taxon_id', 'prior_paper_id', 'prior_raw_taxon', 'reuse_route', 'ncbi_taxid', 'gtdb_r226_species', 'unchanged'])
    check('one_final_mapping_per_published_taxon', len(mapping) == len(taxa) == len({r['original_taxon'] for r in mapping}))
    check('all_mappings_terminal', all(r['mapping_status'].startswith(('mapped_', 'unresolved_')) for r in mapping))
    check('accepted_taxa_have_species_identity', all(r['ncbi_rank'] == 'species' and r['gtdb_r226_species'] not in ('', NA)
                                                  for r in mapping if b.is_accepted(r)))
    STATS.update(mapped_taxa=len(mapping) - len(unresolved), unresolved_taxa=len(unresolved),
                 taxonomy_statuses=dict(Counter(r['mapping_status'] for r in mapping)),
                 unresolved_reasons=dict(Counter(r['mapping_status'] for r in unresolved)),
                 reuse_routes=dict(Counter(r['reuse_route'] for r in mapping)),
                 unresolved_details=[{k: r[k] for k in ['original_taxon', 'mapping_status', 'mapping_reason']} for r in unresolved])
    LOG.info('Mapping: %s; routes %s', STATS['taxonomy_statuses'], STATS['reuse_routes'])
    return mapping


def convert(published, mapping, core, b):
    edges = [dict(r, paper_id='Li2024', network_id='Li2024_published', condition='published_figure_network',
                  primary_status='primary', source_file=str(RAW), raw_record_id=str(r['source_row_id']),
                  taxon1_raw=r['Var1'], taxon2_raw=r['Var2'], sign=sign(r['cor'])) for r in published]
    before = [dict(primary_fate='not_applicable_new_network', both_endpoints_accepted='False') for _ in edges]
    pairs, audit = core(edges, mapping, [], before, None)
    for r in audit:
        r.update(gtdb_species_1=r['endpoint1_gtdb_species'], gtdb_species_2=r['endpoint2_gtdb_species'],
                 canonical_species_a=r['species_a'], canonical_species_b=r['species_b'], external_sign=r['sign'],
                 fate_category='unresolved_endpoint' if r['primary_fate'] in b.FATES[:4] else r['primary_fate'])
    by_id = {r['raw_record_id']: r for r in audit}
    raw_fields = ['source_row_id', 'Var1', 'Var2', 'cor', 'p', 'fdr', 'signature_species_in_figure']
    for r in pairs:
        first = by_id[r['representative_record_id']]
        r.update({k: first[k] for k in raw_fields + ['gtdb_species_1', 'gtdb_species_2']})
        r.update(canonical_species_a=r['species_a'], canonical_species_b=r['species_b'], external_sign=r['sign'],
                 raw_value_policy='first input record; no weight aggregation', representative_external_sign=first['external_sign'])
    mapped = [r for r in audit if r['both_endpoints_accepted']]
    write(OUT / 'taxonomy/li2024_mapped_edges_all.tsv', mapped, list(audit[0]))
    write(OUT / 'taxonomy/li2024_gtdb_r226_pairs.tsv', pairs)
    write(OUT / 'qc/li2024_edge_fate.tsv', audit)
    counts = Counter(r['primary_fate'] for r in audit)
    summary = [dict(primary_fate=k, fate_category='unresolved_endpoint' if k in b.FATES[:4] else k,
                    count=counts[k], fraction=counts[k] / len(published)) for k in b.FATES]
    write(OUT / 'qc/li2024_edge_fate_summary.tsv', summary)
    check('exactly_one_fate_per_raw_edge', len(audit) == len(published) == len(by_id) and
          all(r['primary_fate'] in b.FATES for r in audit))
    check('edge_fate_balance', sum(counts.values()) == len(published), j(counts))
    check('mapped_edge_balance', len(mapped) == counts['self_loop_after_mapping'] + counts['duplicate_after_mapping'] + len(pairs))
    check('kept_pairs_valid_unique_canonical', len({(r['species_a'], r['species_b']) for r in pairs}) == len(pairs) and
          all(r['species_a'] not in ('', NA) and r['species_b'] not in ('', NA) and r['species_a'] < r['species_b'] for r in pairs))
    # Independently reconstruct expected pairs and representative order from
    # the published file and mapping table, not from the computed pair rows.
    lookup = {r['original_taxon']: r for r in mapping}
    expected = defaultdict(list)
    expected_fates = {}
    for r in published:
        endpoints = [lookup[r[k]] for k in ['Var1', 'Var2']]
        failures = [x for x in endpoints if not b.is_accepted(x)]
        if failures:
            fate = 'unresolved_ncbi' if any(x['ncbi_stage_status'] == 'unresolved' for x in failures) else 'unresolved_gtdb'
        else:
            pair = tuple(sorted(x['gtdb_r226_species'] for x in endpoints))
            if pair[0] == pair[1]:
                fate = 'self_loop_after_mapping'
            else:
                fate = 'duplicate_after_mapping' if expected[pair] else 'kept_unique_gtdb_pair'
                expected[pair].append(r)
        expected_fates[str(r['source_row_id'])] = fate
    check('independent_all_record_fates', all(expected_fates[r['raw_record_id']] == r['primary_fate'] for r in audit))
    check('independent_canonical_pair_membership', set(expected) == {(r['species_a'], r['species_b']) for r in pairs})
    for pair in pairs:
        members = expected[(pair['species_a'], pair['species_b'])]
        signs = {sign(x['cor']) for x in members}
        assert json.loads(pair['raw_record_ids']) == [str(x['source_row_id']) for x in members]
        assert pair['representative_record_id'] == str(members[0]['source_row_id'])
        assert all(pair[k] == members[0][k] for k in raw_fields)
        assert pair['external_sign'] == (next(iter(signs)) if len(signs) == 1 else NA)
    check('independent_representative_provenance_and_signs', True)
    check('raw_fields_lossless_in_edge_fate', all(all(r[k] == p[k] for k in raw_fields) for r, p in zip(audit, published)))
    STATS.update(fates={k: counts[k] for k in b.FATES}, final_pairs=len(pairs), fully_mapped_edges=len(mapped),
                 final_pair_signs=dict(Counter(r['external_sign'] for r in pairs)),
                 conflict_pairs=sum(r['sign_status'] == 'conflict' for r in pairs))
    LOG.info('Final pairs=%d; fates=%s', len(pairs), counts)


def verify_integrity():
    for r in INVENTORY:
        p = Path(r['absolute_path'])
        r['sha256_after'] = sha(p) if p.is_file() else 'MISSING'
        r['unchanged'] = r['sha256_before'] == r['sha256_after']
    if INVENTORY:
        write(OUT / 'qc/li2024_input_integrity.tsv', INVENTORY)
        check('frozen_and_raw_hashes_unchanged', all(r['unchanged'] for r in INVENTORY), len(INVENTORY))
        paths_before = {r['absolute_path'] for r in INVENTORY if r['scope'] != 'li2024_raw'}
        check('no_new_files_in_frozen_directories', all(str(p) in paths_before for name in
              ['parsed', 'mapping', 'scripts', 'stage2b', 'raw', 'reference', 'reports']
              for p in (BASE / name).rglob('*') if p.is_file() and p != REPORT))
        STATS['frozen_file_count'] = sum(r['scope'] != 'li2024_raw' for r in INVENTORY)


def report(status, error):
    QC.append(dict(check='overall_status', status=status, detail=error or 'All required Stage 2C checks passed'))
    write(OUT / 'qc/li2024_qc_summary.tsv', QC)
    lines = ['# Stage 2C — Li et al. 2024 intake and taxonomy harmonization', '',
             'overall_status = ' + status, '',
             '## 原始文件与提取', '',
             '原始文件：`' + str(RAW) + '`。source_row_id 为原文件物理行号（表头第 1 行）。',
             '唯一提取条件为原值 `signature_species_in_figure == Yes`；没有 correlation/FDR 二次筛选。', '',
             '```json', json.dumps(STATS.get('raw_audit', {}), ensure_ascii=False, indent=2), '```', '',
             'Published raw edges = {}; taxa = {}; positive/negative/zero = {}。'.format(
                 STATS.get('published_edges', NA), STATS.get('published_taxa', NA), j(STATS.get('published_signs', {}))), '',
             '## Taxonomy', '',
             'Mapped taxa = {}/{} ({:.2%}); unresolved = {}。'.format(STATS.get('mapped_taxa', 0),
                 STATS.get('published_taxa', 0), STATS.get('mapped_taxa', 0) / max(1, STATS.get('published_taxa', 0)),
                 STATS.get('unresolved_taxa', NA)),
             '复用路径：`' + j(STATS.get('reuse_routes', {})) + '`。',
             '最终状态：`' + j(STATS.get('taxonomy_statuses', {})) + '`。', '',
             '首先按冻结 Stage 2 格式 key（下划线/空白/大小写）查找 Stage 2B accepted 原名、已验证候选名和 NCBI scientific name；'
             '证据一致时直接复用。新名称调用冻结 Taxonomy.resolve；已确认 current species TaxID 如已有 accepted mapping，优先复用其 GTDB 决策。'
             '其余直接调用冻结 load_official/apply_official/genome_fallback。没有 fuzzy mapping、手工指定、tie bypass 或新 threshold。',
             'NCBI 名称类别、merged TaxID、7 项 historical TaxID、species rank/domain gate 与冻结实现一致。'
             '官方 top tie 禁止 genome fallback；仅无官方 forward 行时，允许原有 metadata genome support 路径。', '',
             'Unresolved 明细：', '', '```json', json.dumps(STATS.get('unresolved_details', []), ensure_ascii=False, indent=2), '```', '',
             '## Pair、符号与 edge fate', '',
             '执行冻结 Stage 2B `rebuild_edges` 的原始计算语句（行 ' + STATS.get('pair_core_lines', NA) + '），'
             '仅去除针对旧 15534 条边/8 个网络的断言及输出尾部；计算代码未经重写。',
             '计算 AST SHA256：`' + STATS.get('pair_core_ast_sha256', NA) + '`。',
             '按精确 GTDB species 字符串排序（保留后缀），按 network + pair 分组。先移除 self-loop；首条原始记录为代表，'
             '其余为 duplicate；不合成权重，raw_record_ids 保留完整来源。组内符号一致时保留，冲突则 external_sign=NA、sign_status=conflict。',
             '`mapped_edges_all` 包含两端 accepted 的所有原始记录，包括之后的 self-loop/duplicate；'
             '`edge_fate` 覆盖全部 published 原始记录；final pair 的 cor/p/fdr 是明确标记的首条代表原值。',
             '失败优先级沿用 pending_missing_reference > unresolved_sgb > unresolved_ncbi > unresolved_gtdb；'
             '细分类均归属 unresolved_endpoint，summary 不加入重叠小计。', '',
             'Fully mapped raw edges = {}; final unique pairs = {}; conflict pairs = {}。'.format(
                 STATS.get('fully_mapped_edges', NA), STATS.get('final_pairs', NA), STATS.get('conflict_pairs', NA)),
             'Final pair signs = `' + j(STATS.get('final_pair_signs', {})) + '`。',
             'Fates = `' + j(STATS.get('fates', {})) + '`；合计 {} / {}。'.format(sum(STATS.get('fates', {}).values()), STATS.get('published_edges', NA)), '',
             '## QC 与冻结完整性', '',
             '运行前后逐文件 SHA256 对照，保护 {} 个冻结文件/参考输入及 Li raw；并核对已有各阶段 input inventory 的历史 SHA256。'
             '详见 qc/li2024_input_integrity.tsv。'.format(STATS.get('frozen_file_count', NA)),
             'qc/li2024_mapping_reuse_audit.tsv 逐条列出复用的 accepted taxon identity/target 对照。'
             '此前 accepted mapping 表全文 hash 不变；该表未被修改或重新生成。',
             '另从 published records + final taxon table 独立重建 pair membership、每条 fate、代表顺序、来源记录、符号和原始字段，并与冻结算法结果对照。',
             '所有逐项 QC 见 qc/li2024_qc_summary.tsv。', '',
             '本阶段到 QC 停止。未读取 human_gut.info，未执行 overlap、source-exclusion、sign concordance、null、enrichment 或 recurrence。']
    if error:
        lines += ['', 'FAIL reason: ' + error]
    lines += ['', '## 新生成正式文件', '']
    lines += ['- `' + str(p) + '`' for p in sorted(OUT.rglob('*')) if p.is_file() and OUT / 'raw' not in p.parents]
    lines += ['- `' + str(REPORT) + '`', '']
    REPORT.write_text('\n'.join(lines), encoding='utf-8')


def main():
    error = ''
    try:
        snapshot()
        f, b, core = load_frameworks()
        published, taxa = input_audit()
        mapping = harmonize(taxa, f, b)
        convert(published, mapping, core, b)
    except Exception as exc:
        error = '{}: {}'.format(type(exc).__name__, exc)
        LOG.exception('Stage 2C failed; dependent processing stopped')
        QC.append(dict(check='execution', status='FAIL', detail=error))
    try:
        verify_integrity()
    except Exception as exc:
        error += '; integrity: ' + str(exc)
    status = 'FAIL' if error or any(r['status'] == 'FAIL' for r in QC) else 'PASS'
    report(status, error)
    print(j(dict(overall_status=status, **STATS)), flush=True)
    return 0 if status == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
