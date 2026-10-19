#!/usr/bin/env python3
"""Phase 2, offline and reproducible: frozen parsed edges -> GTDB R226 pairs.

Use only standard-library Python. Read the two Chen reference scripts as text/AST;
never import or run them, since their top-level code performs prohibited overlap.
Retain their accepted NCBI classes, stable historical TaxIDs, domain-aware forward
official priority, and unique-highest genome fallback. Deliberately retain square
brackets in ALL name matching (the user's stricter Zhou rule). No fuzzy matching,
taxon expansion, raw Excel/CSV parsing, source exclusion, or overlap is performed.
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
import json
import logging
import os
import re
import sys
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(str(resolve_path('external_root')))
BASE = resolve_path('literature')
OUT = BASE / 'mapping'
REPORTS = BASE / 'reports'
EDGES = BASE / 'parsed/all_new_external_edges_raw.tsv'
REF_NAME = ROOT / 'rescue_external_coabundance_taxonomy_v2.py'
REF_GENOME = ROOT / 'gtdb_r226_genome_level_mapping_rescue_v1.py'
TAXDUMP = ROOT / 'taxdump_20240914'
OFFICIAL = ROOT / 'official_mapping_v2_domainaware/03_ncbi_to_gtdb_unique_majority.tsv'
FORWARD = ROOT / 'official_mapping_v1/01_official_correspondence_all_parsed.tsv'
GENOMES = Path(str(resolve_path('gtdb_r226', 'bac120_metadata_r226.tsv.gz')))
VERSION = 'mpa_vJan25_CHOCOPhlAnSGB_202503'
PAPERS = ['Zhou2022', 'Guo2023', 'Lacruz2026']
EXPECTED = {'Zhou2022_MS': 6742, 'Zhou2022_HHC': 5688, 'Guo2023_AggControl': 388,
            'Guo2023_AggMECFS': 358, 'Lacruz2026_MHNO': 690, 'Lacruz2026_MHO': 698,
            'Lacruz2026_MUNO': 514, 'Lacruz2026_MUO': 456}
ACCEPTED = {'mapped_official_unique', 'mapped_official_unique_majority', 'mapped_genome_fallback'}
NAME_CLASSES = {'scientific name', 'synonym', 'equivalent name'}
NA = 'NA'
STATE = {'status': 'RUNNING', 'warnings': [], 'network_qc': [], 'taxon_qc': []}
LOG = logging.getLogger('phase2')


def check(value, message):
    assert value, message


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path):
    with path.open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        return reader.fieldnames, list(reader)


def scalar(value):
    return NA if value is None else str(value)


def write_tsv(path, rows, fields=None):
    check(BASE in path.resolve().parents and BASE / 'raw' not in path.resolve().parents and
          BASE / 'parsed' not in path.resolve().parents, 'Output outside authorized phase-2 area')
    if fields is None:
        fields = list(dict.fromkeys(k for r in rows for k in r))
    check(fields, 'Empty output schema: ' + str(path))
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows({f: scalar(r.get(f)) for f in fields} for r in rows)
    columns, observed = read_tsv(path)
    check(columns == fields and observed == [{f: scalar(r.get(f)) for f in fields} for r in rows],
          'TSV round-trip mismatch: ' + str(path))


def key(name):
    """Chen formatting normalization; brackets and GTDB suffix letters retained."""
    return re.sub(r'\s+', ' ', re.sub(r'^[A-Za-z]__', '', str(name)).replace('_', ' ')).strip().lower()


def relaxed(name):
    """Reference format rescue except bracket removal, which is prohibited here.
    This is deterministic punctuation normalization, NOT approximate matching.
    A unique species TaxID under accepted NCBI name classes remains mandatory.
    """
    value = key(name)
    value = re.sub(r'\bsp\.\s*', 'sp ', value, flags=re.I)
    return re.sub(r'\s+', ' ', re.sub(r'[./-]+', ' ', value)).strip()


def json_text(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'))


def protected_inputs():
    paths = list((BASE / 'parsed').glob('*')) + [BASE / 'scripts/01_parse_external_networks.py']
    paths += [p for p in REPORTS.glob('*') if p.name.startswith(('01_', 'zhou2022_', 'guo2023_', 'lacruz2026_'))]
    paths += [REF_NAME, REF_GENOME, TAXDUMP / 'names.dmp', TAXDUMP / 'nodes.dmp', TAXDUMP / 'merged.dmp', OFFICIAL, FORWARD, GENOMES]
    result = []
    for path in sorted(set(paths)):
        check(path.is_file(), 'Missing protected/input file: ' + str(path))
        result.append(dict(absolute_path=str(path), size_bytes=path.stat().st_size, sha256=sha(path)))
    STATE['input_inventory'] = result
    write_tsv(REPORTS / '02_input_inventory.tsv', result)
    LOG.info('Protected/input snapshot: %d files; phase-1 edge SHA256=%s', len(result), next(r['sha256'] for r in result if r['absolute_path'] == str(EDGES)))


def verify_protected():
    rows = []
    for before in STATE['input_inventory']:
        path = Path(before['absolute_path'])
        after = sha(path)
        rows.append(dict(before, sha256_after=after, unchanged=before['sha256'] == after and before['size_bytes'] == path.stat().st_size))
    write_tsv(REPORTS / '02_input_integrity_qc.tsv', rows)
    check(all(r['unchanged'] for r in rows), 'Protected input/phase-1 result changed')
    STATE['inputs_unchanged'] = True


def load_edges_manifest(stage):
    LOG.info('STEP 1: validate frozen input and build unique (paper, raw_taxon) manifest')
    fields, edges = read_tsv(EDGES)
    check(len(edges) == 15534, 'Frozen input must contain exactly 15534 records')
    check(Counter(r['network_id'] for r in edges) == EXPECTED, 'Frozen network denominators changed')
    check(Counter(r['primary_status'] for r in edges) == {'primary': 14836, 'exploratory': 698}, 'Frozen primary/exploratory counts changed')
    occurrences, networks, statuses, neighbors = Counter(), defaultdict(Counter), defaultdict(Counter), defaultdict(set)
    for index, row in enumerate(edges, 1):
        check(row['primary_status'] == ('exploratory' if row['network_id'] == 'Lacruz2026_MHO' else 'primary'), 'Status mismatch')
        row['raw_record_id'] = 'E{:05d}'.format(index)
        row['input_tsv_row'] = index + 1
        for a, b in [('taxon1_raw', 'taxon2_raw'), ('taxon2_raw', 'taxon1_raw')]:
            identity = (row['paper_id'], row[a])
            occurrences[identity] += 1
            networks[identity][row['network_id']] += 1
            statuses[identity][row['primary_status']] += 1
            neighbors[identity].add(row[b])
    manifest = []
    for i, identity in enumerate(sorted(occurrences, key=lambda k: (PAPERS.index(k[0]), k[1])), 1):
        paper, raw = identity
        record = dict(raw_taxon_id='T{:04d}'.format(i), raw_taxon=raw, paper_id=paper,
                      source_taxonomy_type='NCBI_like_species_name' if paper == 'Zhou2022' else
                      'NCBI_like_underscore_species' if paper == 'Guo2023' else 'MetaPhlAn4_lineage',
                      normalized_candidate=raw if paper == 'Zhou2022' else raw.replace('_', ' '),
                      taxon_class='species_name_label', occurrence_count=occurrences[identity],
                      network_count=len(networks[identity]), network_ids=';'.join(sorted(networks[identity])),
                      per_network_occurrence_count=json_text(networks[identity]),
                      primary_occurrence_count=statuses[identity]['primary'], exploratory_occurrence_count=statuses[identity]['exploratory'],
                      unique_raw_neighbor_degree=len(neighbors[identity]), metaphlan_species_label=NA,
                      metaphlan_sgb=NA, metaphlan_sgb_has_group_suffix=NA, lacruz_taxon_class=NA, metaphlan_database=NA)
        if paper == 'Lacruz2026':
            parts = raw.split('|')
            species = [s[3:] for s in parts if s.startswith('s__')]
            sgb = [s[3:] for s in parts if s.startswith('t__')]
            record.update(metaphlan_database=VERSION, metaphlan_species_label=species[0] if len(species) == 1 else NA,
                          metaphlan_sgb=sgb[0] if len(sgb) == 1 else NA,
                          metaphlan_sgb_has_group_suffix=sgb[0].endswith('_group') if len(sgb) == 1 else NA)
            # Named species are resolved from s__, independently of whether t__
            # contains an SGB group suffix. Keep the entire t__ value as evidence;
            # never expand an SGB group or reinterpret it as a GTDB species.
            if len(species) != 1 or not species[0] or len(sgb) != 1 or not sgb[0]:
                kind, candidate = 'other_unresolved', NA
            elif re.search(r'(?:^|_)SGB\d+(?:_|$)', species[0]):
                kind, candidate = 'unnamed_SGB', NA
                check(any(token == sgb[0] or token + '_group' == sgb[0] for token in species[0].split('_')),
                      'Species/SGB identifiers disagree: ' + raw)
            else:
                kind, candidate = 'named_species', species[0].replace('_', ' ')
            record.update(taxon_class=kind, lacruz_taxon_class=kind, normalized_candidate=candidate)
        elif re.search(r'\b(?:bacterium|sp\.?|strain)\b', record['normalized_candidate'], flags=re.I):
            record['taxon_class'] = 'provisional_or_strain_like_label'
        manifest.append(record)
    check(sum(r['occurrence_count'] for r in manifest) == 2 * len(edges), 'Manifest occurrence accounting failed')
    write_tsv(stage / '02_external_taxon_manifest.tsv', manifest)
    STATE['manifest'] = manifest
    STATE['input_fields'] = fields
    LOG.info('Unique raw taxa by paper: %s', dict(Counter(r['paper_id'] for r in manifest)))
    return edges, manifest


def search_sgb_metadata():
    LOG.info('STEP 2: bounded local search for exact-version MetaPhlAn metadata')
    home = Path(str(resolve_path('legacy_external_inputs')))
    roots = [(home / p, 5) for p in ['metaphlan_databases', 'metaphlan', 'database', 'databases', 'db',
              'new_co_abundance_consistency', 'project/species_co-abund/ncbi_metabolic_work/data', 'project/MetaKSSD/data']]
    roots += [(Path(p), 3) for p in ['/opt/metaphlan', '/data/metaphlan', '/data/databases/metaphlan']]
    for envroot in [home / 'anaconda3', home / '.conda']:
        for lib in [envroot / 'lib'] + list((envroot / 'envs').glob('*/lib')):
            for py in lib.glob('python*'):
                for name in ['site-packages/metaphlan', 'site-packages/metaphlan/metaphlan_databases', 'site-packages/metaphlan/databases']:
                    if (py / name).is_dir():
                        roots.append((py / name, 3))
    audit, matches = [], []
    for root, limit in roots:
        visited, found = 0, []
        if root.is_dir():
            for current, dirs, files in os.walk(root, followlinks=False):
                dirs[:] = [d for d in dirs if not d.startswith('.') and d not in
                           ['raw', 'parsed', 'mapping', 'overlap', '__pycache__', 'node_modules', 'taxdump_20240914']]
                if len(Path(current).relative_to(root).parts) >= limit:
                    dirs[:] = []
                visited += 1
                found.extend(str(Path(current) / name) for name in files if VERSION in name)
        matches.extend(found)
        audit.append(dict(search_root=str(root), exists=root.is_dir(), max_depth=limit,
                          directories_checked=visited, matching_files=json_text(found)))
    write_tsv(REPORTS / '02_metaphlan_metadata_search.tsv', audit)
    STATE['sgb_metadata_matches'] = sorted(set(matches))
    # Never unpickle an unknown file or silently use a different database build.
    # If local evidence changes, inspect its schema/provenance before accepting it.
    check(not matches, 'Exact-version candidate files found; inspect provenance/schema before continuing SGB mapping: ' + json_text(matches))
    STATE['warnings'].append('No exact-version MetaPhlAn metadata found in the explicitly audited local paths; unnamed SGB taxa remain pending_SGB_mapping.')
    LOG.info('Exact-version metadata matches=0; unnamed SGB mapping remains pending')


class Taxonomy:
    def __init__(self, manifest):
        LOG.info('STEP 3: load frozen NCBI nodes/merged and targeted accepted-name aliases')
        self.parents, self.ranks, self.merged, self.domains = {}, {}, {}, {}
        with (TAXDUMP / 'nodes.dmp').open(encoding='utf-8') as handle:
            for line in handle:
                values = [s.strip() for s in line.split('|')]
                tid = int(values[0])
                self.parents[tid], self.ranks[tid] = int(values[1]), sys.intern(values[2])
        with (TAXDUMP / 'merged.dmp').open(encoding='utf-8') as handle:
            for line in handle:
                values = [s.strip() for s in line.split('|')]
                self.merged[int(values[0])] = int(values[1])
        module = ast.parse(REF_NAME.read_text(encoding='utf-8'))
        self.historical = None
        for node in module.body:
            if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'HISTORICAL_TAXID' for t in node.targets):
                self.historical = ast.literal_eval(node.value)
        check(isinstance(self.historical, dict) and len(self.historical) == 7, 'Historical stable TaxID reference changed')
        # Reading the reference does not execute its overlap or file-write code.
        ast.parse(REF_GENOME.read_text(encoding='utf-8'))
        candidates = {r['normalized_candidate'] for r in manifest if r['normalized_candidate'] != NA}
        wanted = [candidates, {key(n) for n in candidates}, {relaxed(n) for n in candidates}]
        self.aliases = [defaultdict(set), defaultdict(set), defaultdict(set)]
        for line in (TAXDUMP / 'names.dmp').open(encoding='utf-8'):
            values = [s.strip() for s in line.split('|')]
            tid, name, category = int(values[0]), values[1], values[3]
            if category not in NAME_CLASSES:
                continue
            keys = (name, key(name), relaxed(name))
            for level, value in enumerate(keys):
                if value in wanted[level]:
                    self.aliases[level][value].add((tid, self.current(tid), category, name))
        needed_ids = {self.current(tid) for tid in self.historical.values()}
        for tier in self.aliases:
            needed_ids.update(hit[1] for hits in tier.values() for hit in hits)
        self.scientific = {}
        with (TAXDUMP / 'names.dmp').open(encoding='utf-8') as handle:
            for line in handle:
                values = [s.strip() for s in line.split('|')]
                tid = int(values[0])
                if tid in needed_ids and values[3] == 'scientific name':
                    self.scientific[tid] = values[1]
        LOG.info('NCBI nodes=%d; merged=%d; candidate-related current TaxIDs=%d', len(self.parents), len(self.merged), len(needed_ids))

    def current(self, tid):
        seen = set()
        while tid in self.merged:
            check(tid not in seen, 'Cycle in merged TaxIDs')
            seen.add(tid)
            tid = self.merged[tid]
        return tid

    def domain(self, tid):
        tid = self.current(tid)
        if tid in self.domains:
            return self.domains[tid]
        original, seen, result = tid, set(), 'other'
        while tid in self.parents and tid not in seen:
            if tid in {2, 2157, 2759}:
                result = {2: 'bacteria', 2157: 'archaea', 2759: 'eukaryota'}[tid]
                break
            seen.add(tid)
            tid = self.parents[tid]
        self.domains[original] = result
        return result

    def valid_species(self, tid):
        return tid in self.scientific and self.ranks.get(tid) == 'species' and self.domain(tid) in {'bacteria', 'archaea'}

    def resolve(self, record):
        row = dict(record, ncbi_name=NA, ncbi_taxid=NA, ncbi_rank=NA, ncbi_domain=NA,
                   ncbi_mapping_method=NA, ncbi_stage_status='unresolved', ncbi_name_classes=NA,
                   ncbi_evidence=NA, ncbi_candidate_taxids=NA, gtdb_r226_species=NA,
                   gtdb_mapping_method=NA, genome_support_top=NA, genome_support_total=NA,
                   genome_support_basis=NA, mapping_status='unresolved_ncbi', mapping_reason=NA)
        if record['lacruz_taxon_class'] == 'unnamed_SGB':
            row.update(mapping_status='pending_SGB_mapping', ncbi_stage_status='pending_sgb',
                       mapping_reason='Exact-version MetaPhlAn SGB metadata not found locally; no string-derived species assignment')
            return row, []
        if record['normalized_candidate'] == NA:
            row.update(mapping_status='unresolved_SGB', mapping_reason='Malformed/absent MetaPhlAn species or SGB field')
            return row, []
        candidate = record['normalized_candidate']
        evidence, selected, selected_hits, selected_method = [], None, [], None
        best_hits = []
        for level, (label, query) in enumerate([('exact_primary_name', candidate), ('format_normalized_primary_name', key(candidate)),
                                                 ('chen_punctuation_rescue_brackets_preserved', relaxed(candidate))]):
            hits = sorted(self.aliases[level].get(query, set()))
            for old, current, name_class, name in hits:
                evidence.append(dict(raw_taxon_id=row['raw_taxon_id'], paper_id=row['paper_id'], raw_taxon=row['raw_taxon'],
                                     match_tier=label, matched_ncbi_alias=name, name_class=name_class,
                                     alias_taxid=old, current_taxid=current, current_rank=self.ranks.get(current, NA),
                                     current_domain=self.domain(current)))
            if hits:
                best_hits = hits
                species_ids = {h[1] for h in hits if self.valid_species(h[1])}
                if len(species_ids) == 1:
                    selected = next(iter(species_ids))
                    selected_hits = [h for h in hits if h[1] == selected]
                    selected_method = label
                # Never use progressively relaxed matching to bypass ambiguity.
                break
        if selected is None and key(candidate) in self.historical:
            historical_tid = self.historical[key(candidate)]
            current = self.current(historical_tid)
            if self.valid_species(current):
                selected, selected_method = current, 'chen_historical_stable_taxid_rescue'
                row['ncbi_evidence'] = 'Reference HISTORICAL_TAXID: {} -> {} -> current {}'.format(candidate, historical_tid, current)
        if selected is not None:
            row.update(ncbi_name=self.scientific[selected], ncbi_taxid=selected, ncbi_rank='species',
                       ncbi_domain=self.domain(selected), ncbi_mapping_method=selected_method, ncbi_stage_status='mapped_species',
                       ncbi_name_classes=';'.join(sorted({h[2] for h in selected_hits})) or NA,
                       ncbi_candidate_taxids=json_text(sorted({h[1] for h in best_hits})) if best_hits else str(selected),
                       mapping_status='awaiting_gtdb', mapping_reason='Unique accepted prokaryotic species-rank NCBI TaxID')
            if selected_hits:
                row['ncbi_evidence'] = json_text([dict(alias=h[3], name_class=h[2], original_taxid=h[0], current_taxid=h[1]) for h in selected_hits])
                if any(h[0] != h[1] for h in selected_hits):
                    row['ncbi_mapping_method'] += '+merged_taxid'
        else:
            ids = {h[1] for h in best_hits}
            row['ncbi_candidate_taxids'] = json_text(sorted(ids))
            row['ncbi_evidence'] = json_text([dict(alias=h[3], name_class=h[2], taxid=h[1]) for h in best_hits])
            if ids and all(self.ranks.get(tid) != 'species' for tid in ids):
                row.update(mapping_status='unresolved_not_species_rank', mapping_reason='Exact/format NCBI hit does not resolve to species rank; no promotion from strain/genus')
                if len(ids) == 1:
                    tid = next(iter(ids))
                    row.update(ncbi_taxid=tid, ncbi_name=self.scientific.get(tid, NA), ncbi_rank=self.ranks.get(tid, NA), ncbi_domain=self.domain(tid))
            elif len({h[1] for h in best_hits if self.valid_species(h[1])}) > 1:
                row['mapping_reason'] = 'Multiple accepted prokaryotic species TaxIDs; ambiguity not expanded'
            elif ids:
                row['mapping_reason'] = 'No uniquely valid prokaryotic species TaxID (outside scope or missing current scientific name)'
            else:
                row['mapping_reason'] = 'No accepted NCBI scientific/synonym/equivalent name or validated Chen stable-TaxID rescue'
        return row, evidence


def load_official():
    LOG.info('STEP 4: domain-aware official mapping and complete forward ambiguity gate')
    official, forward = {}, {}
    for path, target in [(OFFICIAL, official), (FORWARD, forward)]:
        _, rows = read_tsv(path)
        for row in rows:
            if row['direction'] != 'ncbi_to_gtdb' or row['rank'] != 'species' or row['domain'] not in {'bacteria', 'archaea'}:
                continue
            identity = (row['domain'], key(row['source_key']))
            check(identity not in target, 'Duplicate official source key: ' + str(identity))
            target[identity] = row
    for identity, row in official.items():
        check(identity in forward and forward[identity]['mapping_status'] == 'unique_top_valid', 'Accepted official entry contradicts full forward evidence')
        check(key(row['target_name']) == key(forward[identity]['target_name']), 'Official target disagreement')
        check(float(row['top_percent']) > float(row['second_percent']), 'Official top not uniquely highest')
    return official, forward


def apply_official(rows, official, forward):
    fallback_taxids = set()
    for row in rows:
        row.update(official_status=NA, official_source_name=NA, official_target_name=NA,
                   official_candidate_count=NA, official_top_percent=NA, official_top_tie_count=NA,
                   fallback_eligible=False, gtdb_evidence_file=NA)
        if row['ncbi_stage_status'] != 'mapped_species':
            continue
        identity = (row['ncbi_domain'], key(row['ncbi_name']))
        full = forward.get(identity)
        if full:
            row.update(official_status=full['mapping_status'], official_source_name=full['source_name'],
                       official_target_name=full['target_name'] or NA,
                       official_candidate_count=full['candidate_count'], official_top_percent=full['top_percent'],
                       official_top_tie_count=full['top_tie_count'], genome_support_total=int(full['source_genomes']),
                       genome_support_basis='official_source_genomes; exact top genome count not supplied',
                       gtdb_evidence_file=str(OFFICIAL if identity in official else FORWARD))
            if full['mapping_status'] == 'top_tie':
                row.update(mapping_status='unresolved_official_top_tie', mapping_reason='Official forward top tie; genome fallback is prohibited')
            elif full['mapping_status'] == 'unique_top_valid' and full['target_name']:
                selected = official.get(identity, full)
                row.update(gtdb_r226_species=selected['target_name'],
                           gtdb_mapping_method='official_domainaware_forward' if identity in official else 'official_forward_unique_top_rescue',
                           mapping_status='mapped_official_unique' if int(selected['candidate_count']) == 1 else 'mapped_official_unique_majority',
                           mapping_reason='Official unique highest support accepted without a 50% or 80% threshold')
            else:
                row.update(mapping_status='unresolved_official_nonaccepted', mapping_reason='Official nonaccepted forward row; fallback prohibited')
        else:
            row.update(official_status='absent', fallback_eligible=True,
                       mapping_status='awaiting_genome_fallback', mapping_reason='No official forward row for resolved current NCBI scientific name')
            fallback_taxids.add(row['ncbi_taxid'])
    return fallback_taxids


def genome_fallback(rows, taxonomy, stage):
    """Count ALL bacterial metadata genomes for each fallback-eligible stable
    NCBI species TaxID. Never consult metadata for a taxon with an official tie.
    Keep exact GTDB spelling/suffixes, count candidates, require unique top only.
    """
    wanted = {r['ncbi_taxid'] for r in rows if r['fallback_eligible'] and r['ncbi_domain'] == 'bacteria'}
    LOG.info('STEP 5: genome fallback for %d bacterial stable species TaxIDs; full bac120 metadata scan', len(wanted))
    counts, total, unclassified, evidence = defaultdict(Counter), Counter(), Counter(), []
    with gzip.open(GENOMES, 'rt', encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        check({'accession', 'ncbi_species_taxid', 'gtdb_taxonomy'}.issubset(reader.fieldnames), 'Genome metadata columns missing')
        for index, record in enumerate(reader, 1):
            value = record['ncbi_species_taxid']
            if value.isdigit():
                original_tid = int(value)
                tid = taxonomy.current(original_tid)
                if tid in wanted:
                    total[tid] += 1
                    species = [s[3:] for s in record['gtdb_taxonomy'].split(';') if s.startswith('s__')]
                    check(len(species) == 1, 'Malformed GTDB species metadata field')
                    target = species[0]
                    if target:
                        counts[tid][target] += 1
                    else:
                        unclassified[tid] += 1
                    evidence.append(dict(ncbi_species_taxid=tid, metadata_ncbi_species_taxid=original_tid,
                                         genome_accession=record['accession'], gtdb_r226_species=target or NA,
                                         metadata_file=str(GENOMES)))
            if index % 100000 == 0:
                LOG.info('Genome metadata rows scanned=%d; matching genomes=%d', index, sum(total.values()))
    STATE['genome_rows_scanned'] = index
    check(len({r['genome_accession'] for r in evidence}) == len(evidence), 'Repeated genome accession in fallback metadata evidence')
    candidates = []
    for tid in sorted(wanted):
        ordered = sorted(counts[tid].items(), key=lambda x: (-x[1], x[0]))
        for target, support in ordered:
            candidates.append(dict(ncbi_species_taxid=tid, gtdb_candidate=target, genome_support=support,
                                   classified_genomes=sum(counts[tid].values()), metadata_rows=total[tid],
                                   unclassified_genomes=unclassified[tid], is_top=support == ordered[0][1]))
    write_tsv(stage / '03_genome_fallback_candidate_support.tsv', candidates,
              ['ncbi_species_taxid', 'gtdb_candidate', 'genome_support', 'classified_genomes', 'metadata_rows', 'unclassified_genomes', 'is_top'])
    write_tsv(stage / '03_genome_fallback_genome_evidence.tsv', evidence,
              ['ncbi_species_taxid', 'metadata_ncbi_species_taxid', 'genome_accession', 'gtdb_r226_species', 'metadata_file'])
    for row in rows:
        if not row['fallback_eligible']:
            continue
        tid = row['ncbi_taxid']
        c = counts[tid]
        maximum = max(c.values()) if c else 0
        tops = sorted(t for t, n in c.items() if n == maximum)
        row.update(genome_support_top=maximum, genome_support_total=sum(c.values()),
                   genome_support_basis='exact counts of classified bac120 genomes at current stable NCBI species TaxID',
                   genome_metadata_rows=total[tid], genome_unclassified_rows=unclassified[tid],
                   genome_candidate_count=len(c), genome_top_tie_count=len(tops),
                   genome_top_percent=100.0 * maximum / sum(c.values()) if c else NA,
                   gtdb_evidence_file=str(GENOMES))
        if row['ncbi_domain'] != 'bacteria':
            row.update(mapping_status='unresolved_genome_metadata_domain_unavailable',
                       mapping_reason='Only bac120 metadata is supplied; no archaeal genome fallback')
        elif not c:
            row.update(mapping_status='unresolved_genome_no_support', mapping_reason='No classified GTDB species support for this stable TaxID in supplied R226 metadata')
        elif len(tops) > 1:
            row.update(mapping_status='unresolved_genome_top_tie', mapping_reason='Equal highest genome support across {} GTDB species; no expansion'.format(len(tops)))
        else:
            row.update(gtdb_r226_species=tops[0], gtdb_mapping_method='genome_unique_highest_species_taxid',
                       mapping_status='mapped_genome_fallback', mapping_reason='Unique highest classified genome count; no 50%/80% threshold')
    LOG.info('Taxon mapping statuses: %s', dict(Counter(r['mapping_status'] for r in rows)))


def convert_edges(edges, taxa, stage):
    LOG.info('STEP 6: endpoint conversion, mutually exclusive fates, per-network collapse')
    lookup = {(r['paper_id'], r['raw_taxon']): r for r in taxa}
    check(len(lookup) == len(taxa), 'Taxon mapping is not unique per paper/raw label')
    target_spellings = defaultdict(set)
    for row in taxa:
        if row['mapping_status'] in ACCEPTED:
            check(row['ncbi_stage_status'] == 'mapped_species', 'GTDB accepted without species-rank NCBI identity')
            check(row['gtdb_r226_species'] not in ('', NA), 'Accepted mapping missing target')
            target_spellings[key(row['gtdb_r226_species'])].add(row['gtdb_r226_species'])
            if row['mapping_status'] == 'mapped_genome_fallback':
                check(row['official_status'] == 'absent' and row['genome_top_tie_count'] == 1, 'Fallback bypassed official ambiguity')
        else:
            check(row['gtdb_r226_species'] == NA, 'Unresolved mapping has a GTDB target')
    check(all(len(v) == 1 for v in target_spellings.values()), 'Inconsistent exact GTDB spellings for a single reference key')
    audit, groups = [], defaultdict(list)
    for raw in edges:
        endpoints = [lookup[(raw['paper_id'], raw[f])] for f in ['taxon1_raw', 'taxon2_raw']]
        accepted = [r['mapping_status'] in ACCEPTED for r in endpoints]
        pending = any(r['mapping_status'] == 'pending_SGB_mapping' for r in endpoints)
        ncbi_failure = any(r['ncbi_stage_status'] == 'unresolved' for r in endpoints)
        gtdb_failure = any(r['ncbi_stage_status'] == 'mapped_species' and r['mapping_status'] not in ACCEPTED for r in endpoints)
        record = dict(raw, endpoint1_taxon_id=endpoints[0]['raw_taxon_id'], endpoint2_taxon_id=endpoints[1]['raw_taxon_id'],
                      endpoint1_ncbi_stage=endpoints[0]['ncbi_stage_status'], endpoint2_ncbi_stage=endpoints[1]['ncbi_stage_status'],
                      endpoint1_mapping_status=endpoints[0]['mapping_status'], endpoint2_mapping_status=endpoints[1]['mapping_status'],
                      endpoint1_mapping_reason=endpoints[0]['mapping_reason'], endpoint2_mapping_reason=endpoints[1]['mapping_reason'],
                      endpoint1_gtdb_species=endpoints[0]['gtdb_r226_species'], endpoint2_gtdb_species=endpoints[1]['gtdb_r226_species'],
                      both_endpoints_accepted=all(accepted), any_pending_sgb=pending,
                      any_ncbi_failure=ncbi_failure, any_gtdb_failure=gtdb_failure,
                      species_a=NA, species_b=NA, pair_key=NA, primary_fate=NA,
                      representative_record_id=NA, collapsed_pair_raw_record_count=NA,
                      collapsed_pair_sign_status=NA, sign_concordance_eligible=False)
        # Mutually exclusive precedence handles mixed failures without counting an
        # edge twice. Boolean endpoint failure flags above remain available.
        if pending:
            record['primary_fate'] = 'pending_sgb'
        elif ncbi_failure:
            record['primary_fate'] = 'unresolved_ncbi'
        elif gtdb_failure:
            record['primary_fate'] = 'unresolved_gtdb'
        else:
            check(all(accepted), 'An edge endpoint has no defined failure/acceptance state')
            a, b = sorted(r['gtdb_r226_species'] for r in endpoints)
            check(' || ' not in a and ' || ' not in b, 'Pair-key delimiter present in GTDB name')
            record.update(species_a=a, species_b=b, pair_key=a + ' || ' + b)
            if a == b:
                record['primary_fate'] = 'self_loop_after_mapping'
            else:
                groups[(raw['network_id'], a, b)].append(record)
        audit.append(record)
    final = []
    for (network, a, b), members in sorted(groups.items()):
        representative = members[0]  # Input order only; never choose by weight.
        signs = {r['sign'] for r in members}
        status = 'consistent' if len(signs) == 1 else 'conflict'
        check(len({r['paper_id'] for r in members}) == len({r['primary_status'] for r in members}) == 1,
              'Network collapse mixed papers or primary/exploratory statuses')
        row = dict(paper_id=representative['paper_id'], network_id=network, condition=representative['condition'],
                   primary_status=representative['primary_status'], source_file=representative['source_file'],
                   species_a=a, species_b=b, pair_key=a + ' || ' + b, raw_record_count=len(members),
                   sign=next(iter(signs)) if len(signs) == 1 else NA, sign_status=status,
                   sign_values=';'.join(sorted(signs)), presence_comparison_eligible=True,
                   sign_concordance_eligible=status != 'conflict',
                   primary_paper_support_eligible=representative['primary_status'] == 'primary',
                   representative_record_id=representative['raw_record_id'],
                   raw_record_ids=json_text([r['raw_record_id'] for r in members]))
        final.append(row)
        for i, record in enumerate(members):
            record.update(primary_fate='kept_unique_gtdb_pair' if i == 0 else 'duplicate_after_mapping',
                          representative_record_id=row['representative_record_id'], collapsed_pair_raw_record_count=len(members),
                          collapsed_pair_sign_status=status, sign_concordance_eligible=row['sign_concordance_eligible'])
    allowed_fates = {'kept_unique_gtdb_pair', 'unresolved_ncbi', 'unresolved_gtdb', 'pending_sgb',
                     'self_loop_after_mapping', 'duplicate_after_mapping'}
    check(len(audit) == 15534 and len({r['raw_record_id'] for r in audit}) == 15534, 'Edge audit must cover each input record exactly once')
    check(all(r['primary_fate'] in allowed_fates for r in audit), 'Unassigned edge fate')
    check(len({(r['network_id'], r['pair_key']) for r in final}) == len(final), 'Duplicate final network pair')
    for network in EXPECTED:
        original = [r for r in audit if r['network_id'] == network]
        collapsed = [r for r in final if r['network_id'] == network]
        fates = Counter(r['primary_fate'] for r in original)
        n_unresolved = sum(fates[f] for f in ['unresolved_ncbi', 'unresolved_gtdb', 'pending_sgb'])
        fully_mapped = sum(r['both_endpoints_accepted'] for r in original)
        qc = dict(paper_id=original[0]['paper_id'], network_id=network, primary_status=original[0]['primary_status'],
                  raw_edges=len(original), raw_unique_taxa=len({r[f] for r in original for f in ['taxon1_raw', 'taxon2_raw']}),
                  fully_mapped_edge_records=fully_mapped, edge_records_lost_ncbi=fates['unresolved_ncbi'],
                  edge_records_lost_gtdb=fates['unresolved_gtdb'], edge_records_pending_sgb=fates['pending_sgb'],
                  unresolved_mapping_records=n_unresolved, self_loops_after_gtdb=fates['self_loop_after_mapping'],
                  duplicate_records_collapsed=fates['duplicate_after_mapping'], final_unique_gtdb_pairs=len(collapsed),
                  pair_mapping_coverage=len(collapsed)/len(original), sign_conflict_pairs=sum(r['sign_status'] == 'conflict' for r in collapsed),
                  any_ncbi_failure_records=sum(r['any_ncbi_failure'] for r in original),
                  any_gtdb_failure_records=sum(r['any_gtdb_failure'] for r in original),
                  any_pending_sgb_records=sum(r['any_pending_sgb'] for r in original),
                  edge_fate_balance_ok=len(original) == n_unresolved + fates['self_loop_after_mapping'] + fates['duplicate_after_mapping'] + len(collapsed))
        check(qc['edge_fate_balance_ok'], 'Edge-fate balance failed: ' + network)
        check(fully_mapped == qc['self_loops_after_gtdb'] + qc['duplicate_records_collapsed'] + len(collapsed), 'Fully-mapped reconciliation failed')
        check(sum(r['raw_record_count'] for r in collapsed) == fully_mapped - qc['self_loops_after_gtdb'], 'Collapsed raw-count reconciliation failed')
        check(fates['kept_unique_gtdb_pair'] == len(collapsed), 'Representative fate count differs from final pairs')
        STATE['network_qc'].append(qc)
        LOG.info('%s: raw=%d; pending=%d; NCBI-loss=%d; GTDB-loss=%d; self-loop=%d; collapse=%d; final=%d; coverage=%.4f; balance=PASS',
                 network, len(original), qc['edge_records_pending_sgb'], qc['edge_records_lost_ncbi'], qc['edge_records_lost_gtdb'],
                 qc['self_loops_after_gtdb'], qc['duplicate_records_collapsed'], len(collapsed), qc['pair_mapping_coverage'])
    write_tsv(stage / '04_edge_mapping_audit.tsv', audit)
    write_tsv(stage / '04_all_external_edges_gtdb_r226.tsv', final)
    write_tsv(REPORTS / '02_taxonomy_mapping_qc.tsv', STATE['network_qc'])
    STATE['edge_fates'] = dict(Counter(r['primary_fate'] for r in audit))
    STATE['pending_affected'] = dict(total=sum(r['any_pending_sgb'] for r in audit),
                                    primary=sum(r['any_pending_sgb'] and r['primary_status'] == 'primary' for r in audit),
                                    MHO_exploratory=sum(r['any_pending_sgb'] and r['primary_status'] == 'exploratory' for r in audit))
    return final, audit


def taxon_audit(taxa):
    for paper in PAPERS:
        part = [r for r in taxa if r['paper_id'] == paper]
        ncbi = Counter(r['ncbi_stage_status'] for r in part)
        mapped = sum(r['mapping_status'] in ACCEPTED for r in part)
        STATE['taxon_qc'].append(dict(paper_id=paper, unique_raw_taxa=len(part),
                                     ncbi_mapped_taxa=ncbi['mapped_species'], ncbi_unresolved_taxa=ncbi['unresolved'],
                                     pending_sgb_taxa=ncbi['pending_sgb'], gtdb_mapped_taxa=mapped,
                                     gtdb_unresolved_after_ncbi_taxa=ncbi['mapped_species']-mapped,
                                     total_unmapped_taxa_including_pending=len(part)-mapped,
                                     mapped_official_unique=sum(r['mapping_status'] == 'mapped_official_unique' for r in part),
                                     mapped_official_unique_majority=sum(r['mapping_status'] == 'mapped_official_unique_majority' for r in part),
                                     mapped_genome_fallback=sum(r['mapping_status'] == 'mapped_genome_fallback' for r in part)))
        check(ncbi['mapped_species'] + ncbi['unresolved'] + ncbi['pending_sgb'] == len(part), 'Taxon-stage balance failed')
    write_tsv(REPORTS / '02_taxon_stage_summary.tsv', STATE['taxon_qc'])
    unresolved = []
    for r in taxa:
        if r['mapping_status'] not in ACCEPTED:
            unresolved.append(dict(paper_id=r['paper_id'], taxon=r['raw_taxon'], raw_taxon_id=r['raw_taxon_id'],
                                   affected_edges=r['occurrence_count'], unique_raw_neighbor_degree=r['unique_raw_neighbor_degree'],
                                   primary_affected_edges=r['primary_occurrence_count'], exploratory_affected_edges=r['exploratory_occurrence_count'],
                                   network_count=r['network_count'], network_ids=r['network_ids'],
                                   per_network_affected_edges=r['per_network_occurrence_count'],
                                   normalized_candidate=r['normalized_candidate'], ncbi_name=r['ncbi_name'],
                                   ncbi_taxid=r['ncbi_taxid'], ncbi_rank=r['ncbi_rank'],
                                   mapping_status=r['mapping_status'], mapping_reason=r['mapping_reason']))
    unresolved.sort(key=lambda r: (PAPERS.index(r['paper_id']), -r['affected_edges'], r['taxon'], r['mapping_reason']))
    write_tsv(REPORTS / '02_unresolved_taxa.tsv', unresolved)
    STATE['high_degree_unresolved'] = [{k: r[k] for k in ['paper_id', 'taxon', 'affected_edges', 'mapping_status']} for r in sorted(unresolved, key=lambda r: -r['affected_edges'])[:20]]
    STATE['ncbi_methods'] = dict(Counter(r['ncbi_mapping_method'] for r in taxa if r['ncbi_stage_status'] == 'mapped_species'))
    STATE['mapping_status_counts'] = dict(Counter(r['mapping_status'] for r in taxa))
    STATE['lacruz_classes'] = dict(Counter(r['lacruz_taxon_class'] for r in taxa if r['paper_id'] == 'Lacruz2026'))


def md_table(rows, columns):
    def escaped(value):
        return scalar(value).replace('|', '\\|').replace('\n', ' ')
    return '\n'.join(['| ' + ' | '.join(columns) + ' |', '| ' + ' | '.join(['---'] * len(columns)) + ' |'] +
                     ['| ' + ' | '.join(escaped(r.get(c)) for c in columns) + ' |' for r in rows])


def write_report():
    lines = ['# 第二阶段：external taxon harmonization / GTDB R226 QC', '',
             '状态：**{}**；UTC：{}。'.format(STATE['status'], datetime.now(timezone.utc).isoformat()), '',
             '唯一网络输入为 `parsed/all_new_external_edges_raw.tsv`，固定 15534 records '
             '(primary=14836，exploratory=698)。未重读原始 Excel/CSV，未筛选原始 edges，未修改第一阶段结果或 Chen/Wallen 文件。'
             '未联网、下载、使用 MACK-DB、做 overlap/source exclusion/external recurrence。', '',
             '复现：`python3 scripts/02_harmonize_external_taxonomy.py`。仅依赖 Python 标准库。', '',
             '## 固定 framework 与证据', '',
             '- 只读取两份指定参考脚本；从 AST 提取原有 7 条 HISTORICAL_TAXID，不执行参考脚本中的 overlap 代码。',
             '- 使用 taxdump_20240914 的 names/nodes/merged，接受 scientific name、synonym、equivalent name。'
             '原名称精确匹配优先，其次确定性大小写/空白/下划线格式归一化，最后 Chen 的 sp. 与点/斜线/连字符格式 rescue。'
             '所有层级保留方括号；不删除 [] 强行匹配，不使用 fuzzy matching。',
             '- 必须唯一解析到 bacterial/archaeal species-rank TaxID。merged TaxID 按冻结 taxdump 跟踪；'
             '不把 strain、genus 或其他非 species-rank 命中提升为 species，不猜测未命名标签。',
             '- official 优先使用指定 domain-aware unique-majority 表。为区分“缺失”与“被拒绝”，'
             '另只读参考脚本原本使用的 `official_mapping_v1/01_official_correspondence_all_parsed.tsv` 完整 forward 表。'
             '该补充本地证据文件记录在 input inventory；official top tie/nonaccepted 永不进入 genome fallback。',
             '- 仅当前 NCBI scientific name 缺少 official forward source row 时，按稳定 species TaxID 统计全部匹配 bac120 genomes。'
             'genome-level unique-highest 才接受，top tie unresolved；不设 50%/80% 硬阈值、不做 one-to-many expansion。',
             '- GTDB 输出保留官方物种名拼写及 `_A`、`_B` 等后缀。无向 pair 用 species_a <= species_b 排序，'
             'pair_key=`species_a || species_b`，仅在 network 内 collapse。',
             '- official 表只提供 source_genomes 与四舍五入的 top_percent，未提供精确 top genome count。'
             '因此 official rows 的 genome_support_top=NA、genome_support_total=source_genomes；'
             '不从百分比反推整数。fallback rows 的两项支持数则为 metadata 中的精确计数，来源由 genome_support_basis 标明。', '',
             '## Taxon manifest 与阶段计数', '',
             'manifest 每个 (paper_id, raw_taxon) 一行；raw_taxon 完整保留。occurrence_count 为 endpoint 出现次数，'
             'network_ids / per_network_occurrence_count 记录各网络来源。Guo 仅先做 `_`→空格；'
             'Lacruz 保留完整 lineage，另提取 s__、t__ 字段。', '',
             md_table(STATE['taxon_qc'], ['paper_id', 'unique_raw_taxa', 'ncbi_mapped_taxa', 'ncbi_unresolved_taxa',
                                       'pending_sgb_taxa', 'gtdb_mapped_taxa', 'gtdb_unresolved_after_ncbi_taxa',
                                       'total_unmapped_taxa_including_pending']), '',
             '`gtdb_unresolved_after_ncbi_taxa` 只统计已通过 NCBI stage 的 taxa；最后一列另含 NCBI unresolved 和 pending SGB。'
             '各 study 的 unique raw taxa 不等于去重后的 NCBI/GTDB biological species 数。', '',
             md_table(STATE['taxon_qc'], ['paper_id', 'mapped_official_unique', 'mapped_official_unique_majority', 'mapped_genome_fallback']), '',
             'NCBI mapping methods：`{}`。'.format(json_text(STATE.get('ncbi_methods', {}))), '',
             'Mapping statuses：`{}`。'.format(json_text(STATE.get('mapping_status_counts', {}))), '',
             '## Lacruz species labels 与 pending SGB', '',
             '指定 MetaPhlAn build：`' + VERSION + '`。只在 `02_metaphlan_metadata_search.tsv` 列出的有限本地目录与深度检查文件名，'
             '没有全服务器搜索、换版本或下载 metadata。', '',
             'Taxon classes：`{}`。'.format(json_text(STATE.get('lacruz_classes', {}))), '',
             '`unnamed_SGB` 包括 GGB…_SGB…、命名 genus_SGB… 和 unclassified_SGB… 等所有明确 SGB placeholder。'
             '`named_species` 表示 species field 中非 SGB placeholder 的名称标签，包括 provisional/strain-like 名称；'
             '并不预先保证 species rank，仍必须通过 NCBI species-rank 检查。'
             't__ 的 `_group` 后缀完整保留并单独标记；有明确 s__ 名称的记录正常进入 NCBI 核验，'
             '不因 t__ 后缀额外丢弃，也不展开 SGB group。结构缺失/不一致归 other_unresolved 或终止核查。', '',
             'Pending 影响（每条 edge 仅计一次，即使两个 endpoints 都 pending）：`{}`。'.format(json_text(STATE.get('pending_affected', {}))), '',
             '没有可靠同版本 metadata 时保留 pending_SGB_mapping，未创建推测性的 SGB→species mapping。此部分到此停止。', '',
             '## 逐 network taxonomy loss、self-loop、duplicate collapse', '',
             md_table(STATE['network_qc'], ['network_id', 'primary_status', 'raw_edges', 'raw_unique_taxa',
                                          'fully_mapped_edge_records', 'edge_records_lost_ncbi', 'edge_records_lost_gtdb',
                                          'edge_records_pending_sgb', 'self_loops_after_gtdb', 'duplicate_records_collapsed',
                                          'final_unique_gtdb_pairs', 'pair_mapping_coverage', 'sign_conflict_pairs']), '',
             'pair_mapping_coverage = final_unique_gtdb_pairs / raw_edges。fully_mapped_edge_records 包含随后成为 self-loop '
             '或被 collapse 的 records；这些不计作 taxonomy failure。', '',
             '互斥 fate 优先级固定为 **pending_sgb → unresolved_ncbi → unresolved_gtdb → self_loop_after_mapping → '
             'kept_unique_gtdb_pair / duplicate_after_mapping**。同一 edge 有多个失败端点时，只按此顺序计一次；'
             '审计表仍保存两个 endpoint status/reason 及 any_*_failure 标记，额外 any_* 汇总不是互斥数，不能相加。', '',
             '每个网络必须满足：raw_edges = unresolved_mapping_records + self_loops_after_gtdb + '
             'duplicate_records_collapsed + final_unique_gtdb_pairs。', '',
             md_table(STATE['network_qc'], ['network_id', 'raw_edges', 'unresolved_mapping_records', 'self_loops_after_gtdb',
                                          'duplicate_records_collapsed', 'final_unique_gtdb_pairs', 'edge_fate_balance_ok']), '',
             '全量互斥 edge fates：`{}`。'.format(json_text(STATE.get('edge_fates', {}))), '',
             'collapse 选择输入次序中的首条 record 作为审计代表；保留 raw_record_count 与完整 raw_record_ids。'
             '不合成 association weight。sign 全同保留原 sign；任何不一致保留 pair presence，设置 sign_status=conflict、'
             'sign=NA、sign_concordance_eligible=False。', '',
             'Lacruz MHO 仍为 exploratory，primary_paper_support_eligible=False；其独有 pairs 不能增加后续 primary paper support。'
             '所有其他网络保持 primary。本阶段不计算跨网络 recurrence。', '',
             '## 高 degree unresolved taxa', '',
             '完整清单：`02_unresolved_taxa.tsv`，按 paper、affected_edges 降序、taxon、mapping_reason 排列。'
             '每个 taxon 的 affected_edges 是涉及该 taxon 的记录数；不同 taxon 的该数值不可直接相加为损失 edge 数。', '',
             md_table(STATE.get('high_degree_unresolved', []), ['paper_id', 'taxon', 'affected_edges', 'mapping_status']), '',
             '## 完整性与停止范围', '',
             '输入及第一阶段产物 SHA256 前后不变：**{}**。'.format(STATE.get('inputs_unchanged', False)), '',
             '所有 TSV 写出后逐字段回读核对；逐 edge/逐 network 的记录数、collapse membership、互斥 fate 均有 assert。'
             '不允许 optimized Python 绕过断言；失败时非零退出，不发布 staging mapping 文件。', '']
    lines += ['- ' + w for w in STATE['warnings']]
    if 'error' in STATE:
        lines += ['', '**停止原因：** ' + STATE['error']]
    lines += ['', '本阶段完成后停止：未开始 MACK-DB overlap、source exclusion 或 external recurrence。', '']
    (REPORTS / '02_taxonomy_mapping_qc.md').write_text('\n'.join(lines), encoding='utf-8')


def main():
    if sys.flags.optimize:
        raise RuntimeError('Optimized Python is prohibited: required assertions must remain active')
    OUT.mkdir(exist_ok=True)
    REPORTS.mkdir(exist_ok=True)
    LOG.setLevel(logging.INFO)
    for handler in [logging.FileHandler(REPORTS / '02_harmonize_external_taxonomy.log', mode='w', encoding='utf-8'),
                    logging.StreamHandler(sys.stdout)]:
        handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))
        LOG.addHandler(handler)
    LOG.info('START phase 2; fixed BASE=%s; no phase-1/raw edits; offline', BASE)
    try:
        protected_inputs()
        with tempfile.TemporaryDirectory(prefix='.phase2-', dir=str(OUT)) as temp:
            stage = Path(temp)
            edges, manifest = load_edges_manifest(stage)
            search_sgb_metadata()
            taxonomy = Taxonomy(manifest)
            taxa, evidence = [], []
            for row in manifest:
                resolved, hits = taxonomy.resolve(row)
                taxa.append(resolved)
                evidence.extend(hits)
            LOG.info('NCBI stages: %s', dict(Counter(r['ncbi_stage_status'] for r in taxa)))
            write_tsv(stage / '03_ncbi_name_evidence.tsv', evidence)
            official, forward = load_official()
            apply_official(taxa, official, forward)
            genome_fallback(taxa, taxonomy, stage)
            write_tsv(stage / '03_external_taxa_to_gtdb_r226.tsv', taxa)
            taxon_audit(taxa)
            convert_edges(edges, taxa, stage)
            LOG.info('STEP 7: verify protected hashes and publish validated phase-2 outputs')
            verify_protected()
            for source in sorted(stage.iterdir()):
                os.replace(str(source), str(OUT / source.name))
        STATE['status'] = 'PASS_WITH_PENDING_SGB' if STATE['pending_affected']['total'] else 'PASS'
        write_report()
        LOG.info('DONE %s; raw=15534; final within-network unique pairs=%d; every network fate balance passed',
                 STATE['status'], sum(r['final_unique_gtdb_pairs'] for r in STATE['network_qc']))
        return 0
    except Exception as error:
        STATE.update(status='FAIL', error='{}: {}'.format(type(error).__name__, error))
        LOG.exception('STOP phase 2: %s', STATE['error'])
        if STATE.get('input_inventory'):
            try:
                verify_protected()
            except Exception as integrity_error:
                STATE['warnings'].append('Protected-input verification failed: ' + str(integrity_error))
        write_report()
        return 1


if __name__ == '__main__':
    sys.exit(main())
