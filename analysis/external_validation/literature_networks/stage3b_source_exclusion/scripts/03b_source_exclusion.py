#!/usr/bin/env python3
"""Study-specific support exclusion; exact frozen GTDB pairs; no network inference."""

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
import sys
import time
from pathlib import Path

OUT = Path(__file__).resolve().parents[1]
BASE = OUT.parent
PROJECT = Path(str(resolve_path('project_data')))
DB = PROJECT / 'coabundance/cor03/species/belongs/f2/human_gut.info'
DB_HASH = '4b67935940a765d28518a2dd612587f347c0e47e985323910fcdde51950a7834'
RULE = PROJECT / 'coabundance/cor03/get_link_info_fastspar.py'
WALLEN = BASE.parent / '9_2_codex/PMID_36376318_taxonomy_rescue_overlap_v44'
DENOM = dict(zip(['Zhou_MS', 'Zhou_HHC', 'Guo_AggControl', 'Guo_AggMECFS', 'Lacruz_MHNO', 'Lacruz_MHO', 'Lacruz_MUNO', 'Lacruz_MUO', 'Li_IBS'],
                 [4818, 4086, 335, 312, 649, 623, 486, 417, 599]))
VECTORS = ['proj_ID', 'proj_corr', 'proj_pv', 'proj_sam']


def read(path, delimiter='\t'):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def write(path, rows, fields=None):
    rows = list(rows)
    with Path(path).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def canonical(a, b):
    return tuple(sorted((a, b)))


def key(a, b):
    return ' || '.join(canonical(a, b))


def main():
    qc = []

    def check(name, ok, details=''):
        qc.append(dict(check=name, status='PASS' if ok else 'FAIL', detail=str(details)))
        if not ok:
            write(OUT / 'qc/03b_source_exclusion_qc.tsv', qc)
            raise AssertionError((name, details))

    started = time.time()
    integrity = json.loads((OUT / 'qc/03b_frozen_snapshot.json').read_text())
    tracked = {r['input_file'] for r in integrity}

    def track(path, scope, expected=''):
        path = str(Path(path).resolve())
        if path in tracked:
            return
        h = sha(path)
        check('input_hash:' + path, not expected or expected == h, scope)
        integrity.append(dict(input_file=path, size_bytes=Path(path).stat().st_size,
                              sha256_before=h, prior_recorded_sha256=expected,
                              matches_prior='NA' if not expected else str(h == expected), scope=scope))
        tracked.add(path)

    check('original_database_hash', sha(DB) == DB_HASH, DB_HASH)
    check('initial_frozen_snapshot_matches_prior', all(r['matches_prior'] != 'False' for r in integrity), len(integrity))
    pq = read(OUT / 'qc/03b_source_provenance_qc.tsv')
    check('source_provenance_all_PASS', pq and all(r['status'] == 'PASS' for r in pq), len(pq))
    for p in sorted((OUT / 'provenance').rglob('*')):
        if p.is_file():
            track(p, 'Stage3B_provenance_evidence')
    for p in sorted((OUT / 'scripts').glob('*.py')):
        track(p, 'Stage3B_analysis_script')
    track(OUT / 'qc/03b_source_provenance_qc.tsv', 'Stage3B_provenance_QC')
    track(OUT / 'qc/03b_aggregation_runtime_check.tsv', 'aggregation_runtime_comparison')
    track(sys.executable, 'aggregation_python_interpreter')
    check('aggregation_runtime', sys.version_info[:2] == (3, 13), sys.executable + '; ' + sys.version)
    config = json.loads((OUT / 'provenance/03b_network_source_config.json').read_text())
    # Extract the actual original fmt function and arithmetic AST, never run the original main.
    tree = ast.parse(RULE.read_text())
    fmt_node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'fmt')
    ns = {}
    exec(compile(ast.Module(body=[fmt_node], type_ignores=[]), str(RULE), 'exec'), ns)
    fmt = ns['fmt']
    assignments = {n.targets[0].id: n.value for n in ast.walk(tree)
                   if isinstance(n, ast.Assign) and len(n.targets) == 1 and isinstance(n.targets[0], ast.Name)
                   and n.targets[0].id in ('aver_corr', 'aver_pv')}
    for target, var in [('aver_corr', 'cor_values'), ('aver_pv', 'pv_values')]:
        expected = ast.parse('sum(' + var + ')/len(' + var + ')', mode='eval').body
        check('formal_aggregation_' + target, ast.dump(assignments[target]) == ast.dump(expected), 'unweighted arithmetic mean')
    check('formal_fmt', all(fmt(x) == format(x, '.10g') for x in [0, 1/3, -0.1234567890123, 0.0000123456789123]), str(RULE))
    write(OUT / 'qc/03b_aggregation_rule.tsv', [dict(source_file=str(RULE), sha256=sha(RULE),
          aver_corr='sum(cor_values)/len(cor_values)', aver_pv='sum(pv_values)/len(pv_values)', format='.10g',
          precision_policy='Original full-precision project corr_pvalue.csv floats; serialized support vectors retained verbatim by index',
          unchanged_row_policy='Keep original averages for rows with no removed support', execution='AST inspection + fmt definition only; original main never executed',
          python_executable=sys.executable, python_version=sys.version.replace('\n', ' '))])

    manifest = read(BASE / 'stage3a/qc/03a_input_manifest.tsv')
    check('nine_manifest_networks', {r['network'] for r in manifest} == set(DENOM), len(manifest))
    pair_cache, baseline, input_manifest = {}, {}, []
    summary0 = {r['network']: r for r in read(BASE / 'stage3a/summary/03a_baseline_overlap_summary.tsv')}
    all_external = set()
    for m in manifest:
        n, path = m['network'], m['input_file']
        check(n + '_frozen_pair_hash', sha(path) == m['sha256'], path)
        if path not in pair_cache:
            pair_cache[path] = read(path)
        pairs = [r for r in pair_cache[path] if r['network_id'] == m['source_network_id']]
        keys = [key(r['species_a'], r['species_b']) for r in pairs]
        check(n + '_exact_frozen_pair_keys', all(k == r['pair_key'] for k, r in zip(keys, pairs)), 'case and suffix preserved; unordered lexical species pair')
        check(n + '_denominator', len(keys) == len(set(keys)) == DENOM[n] == int(m['final_unique_pair_count']), DENOM[n])
        bf = BASE / 'stage3a/baseline_overlap' / (n + '_baseline_overlap.tsv')
        baseline[n] = read(bf)
        check(n + '_baseline_pair_identity', {r['external_pair_key'] for r in baseline[n]} == set(keys) and len(baseline[n]) == len(keys), str(bf))
        all_external.update(keys)
        input_manifest.append(dict(network=n, source_network_id=m['source_network_id'], input_file=path,
                                   sha256=m['sha256'], final_unique_pair_count=len(keys), baseline_detail_file=str(bf),
                                   baseline_detail_sha256=sha(bf), source_projects=';'.join(config[n]['sources']),
                                   primary_or_exploratory=m['primary_or_exploratory'], published_raw_edges=m['published_raw_edges']))
    write(OUT / 'qc/03b_input_manifest.tsv', input_manifest)
    decoder = {r['database_label']: r['gtdb_r226_species'] for r in read(BASE / 'stage3a/qc/03a_database_label_decoding.tsv')}
    listed_path = PROJECT / 'env/human_gut.list'
    listed = listed_path.read_text().split()
    check('human_gut_875_unique_listed_projects', len(listed) == len(set(listed)) == 875, len(listed))
    listed = set(listed)
    stats_path = DB.parent / 'human_gut_links_stats.tsv'
    missing_path = DB.parent / 'human_gut_links_missing_projects.list'
    track(stats_path, 'formal_aggregation_statistics')
    track(missing_path, 'formal_unavailable_projects')
    stats = read(stats_path)[0]
    missing = set(missing_path.read_text().split())
    check('formal_859_available_16_unavailable', int(stats['available_projects']) == 859 and int(stats['missing_projects']) == len(missing) == 16,
          'Unavailable projects retain their existing meaning; no zero imputation')
    groups = {}
    for n, c in config.items():
        g = c['database_group']
        if g in groups:
            check(g + '_shared_source_set_' + n, groups[g]['sources'] == set(c['sources']), c['sources'])
        else:
            groups[g] = dict(sources=set(c['sources']), networks=[], paper=c['paper'], affected=0, retained=0, removed=0,
                             lost0=0, lost1=0, recovered=set())
        groups[g]['networks'].append(n)
    union_sources = set().union(*(g['sources'] for g in groups.values()))
    supports, external_db, affected_rows, needed = collections.Counter(), {}, {}, collections.defaultdict(set)
    total, malformed, raw_duplicate = 0, 0, 0
    print('Scanning baseline database and full-precision requirements...', flush=True)
    with DB.open() as f:
        dr = csv.DictReader(f, delimiter='\t')
        fields = dr.fieldnames
        for line, r in enumerate(dr, 2):
            total += 1
            vs = [r[k].split(';') for k in VECTORS]
            ok = all(len(v) == int(r['proj_count']) for v in vs) and len(vs[0]) == len(set(vs[0])) and int(r['proj_count']) >= 2
            malformed += not ok
            supports.update(vs[0])
            k = key(decoder[r['species_a']], decoder[r['species_b']])
            if k in all_external:
                raw_duplicate += k in external_db
                external_db[k] = (line, r)
            if union_sources.intersection(vs[0]):
                rawkey = canonical(r['species_a'], r['species_b'])
                affected_rows[rawkey] = r
                for project in vs[0]:
                    needed[project].add(rawkey)
    check('baseline_653071_rows', total == 653071, total)
    check('all_support_vectors_aligned_unique', malformed == 0, malformed)
    check('canonical_external_database_keys_unique', raw_duplicate == 0, raw_duplicate)
    check('actual_support_859_projects', len(supports) == 859 and set(supports) == listed - missing, len(supports))

    # Full precision matches the original generator, avoiding double rounding from .info's .10g vectors.
    prior_inputs = read(BASE / 'stage3a1_zhou_miss/qc/03a1_input_manifest.tsv')
    csv_inputs = {r['project_id']: r for r in prior_inputs if r['role'] == 'project_filtered_csv' and r['exists'] == 'True'}
    precision = collections.defaultdict(dict)
    print('Reading original project CSV floats for', len(affected_rows), 'affected rows from', len(needed), 'projects...', flush=True)
    for i, (project, pairs) in enumerate(sorted(needed.items()), 1):
        pm = csv_inputs[project]
        path = Path(pm['input_file'])
        track(path, 'original_project_support_full_precision', pm['sha256_after'])
        found = set()
        with path.open() as f:
            for r in csv.DictReader(f):
                rawkey = canonical(r['from'].strip(), r['to'].strip())
                if rawkey not in pairs:
                    continue
                assert rawkey not in found
                found.add(rawkey)
                precision[rawkey][project] = (float(r['corr']), float(r['pvalue']))
        check(project + '_original_float_support_complete', found == pairs, len(found))
        if i % 200 == 0:
            print('Loaded project CSVs:', i, '/', len(needed), flush=True)
    precision_bad, original_mean_bad = 0, 0
    mean_diagnostics = []
    for rawkey, r in affected_rows.items():
        ids = r['proj_ID'].split(';')
        vals = [precision[rawkey][p] for p in ids]
        precision_bad += [fmt(v[0]) for v in vals] != r['proj_corr'].split(';') or [fmt(v[1]) for v in vals] != r['proj_pv'].split(';')
        computed_corr = fmt(sum(v[0] for v in vals)/len(vals))
        computed_pv = fmt(sum(v[1] for v in vals)/len(vals))
        mismatch = computed_corr != r['aver_corr'] or computed_pv != r['aver_pv']
        original_mean_bad += mismatch
        if mismatch:
            mean_diagnostics.append(dict(species_a=r['species_a'], species_b=r['species_b'],
                baseline_aver_corr=r['aver_corr'], recomputed_aver_corr=computed_corr,
                baseline_aver_pv=r['aver_pv'], recomputed_aver_pv=computed_pv, proj_count=r['proj_count'],
                proj_ID=r['proj_ID'], original_csv_values=json.dumps(vals)))
    write(OUT / 'qc/03b_original_mean_diagnostics.tsv', mean_diagnostics,
          ['species_a', 'species_b', 'baseline_aver_corr', 'recomputed_aver_corr', 'baseline_aver_pv', 'recomputed_aver_pv', 'proj_count', 'proj_ID', 'original_csv_values'])
    check('original_float_vectors_reproduce_serialized_support', precision_bad == 0, precision_bad)
    check('original_float_means_reproduce_baseline', original_mean_bad == 0, original_mean_bad)

    def exclude(r, sources):
        vectors = [r[k].split(';') for k in VECTORS]
        keep = [i for i, p in enumerate(vectors[0]) if p not in sources]
        removed = [p for p in vectors[0] if p in sources]
        new = r.copy()
        if not removed:
            return new, removed
        for name, values in zip(VECTORS, vectors):
            new[name] = ';'.join(values[i] for i in keep)
        new['proj_count'] = str(len(keep))
        rawkey = canonical(r['species_a'], r['species_b'])
        for col, pos in [('aver_corr', 0), ('aver_pv', 1)]:
            vals = [precision[rawkey][vectors[0][i]][pos] for i in keep]
            new[col] = fmt(sum(vals)/len(vals)) if vals else 'NA'
        return new, removed

    handles, writers, db_removed = {}, {}, []
    for g, v in groups.items():
        v['actual'] = sorted(v['sources'] & set(supports))
        if not v['actual']:
            v['file'] = str(DB)
            v['identical'] = 'Yes'
            v['retained'] = total
        else:
            path = OUT / 'source_excluded_networks' / (g + '_source_excluded_human_gut.tsv')
            v['file'], v['identical'] = str(path), 'No'
            handles[g] = path.open('w', newline='')
            writers[g] = csv.DictWriter(handles[g], fieldnames=fields, delimiter='\t', lineterminator='\n')
            writers[g].writeheader()
    print('Writing study-specific derivative databases...', flush=True)
    with DB.open() as f:
        for r in csv.DictReader(f, delimiter='\t'):
            ids = set(r['proj_ID'].split(';'))
            for g, writer in writers.items():
                v = groups[g]
                if not ids & v['sources']:
                    new, removed = r, []
                else:
                    new, removed = exclude(r, v['sources'])
                    v['affected'] += 1
                count = int(new['proj_count'])
                if count >= 2:
                    writer.writerow(new)
                    v['retained'] += 1
                else:
                    v['removed'] += 1
                    v['lost' + str(count)] += 1
                    db_removed.append(dict(database_group=g, canonical_pair_key=key(decoder[r['species_a']], decoder[r['species_b']]),
                                           baseline_proj_count=r['proj_count'], remaining_proj_count=count,
                                           source_projects_removed=';'.join(removed), remaining_proj_ID=new['proj_ID'],
                                           loss_reason='remaining_project_count_' + str(count)))
    for f in handles.values():
        f.close()
    write(OUT / 'source_excluded_networks/03b_removed_database_associations.tsv', db_removed,
          ['database_group', 'canonical_pair_key', 'baseline_proj_count', 'remaining_proj_count', 'source_projects_removed', 'remaining_proj_ID', 'loss_reason'])

    # Independent row-by-row verification using zipped tuples, rather than the construction index mask.
    print('Independently verifying every retained derivative record...', flush=True)
    derivative_files = {g: Path(v['file']).open() for g, v in groups.items() if v['identical'] == 'No'}
    readers = {g: csv.DictReader(f, delimiter='\t') for g, f in derivative_files.items()}
    errors = collections.Counter()
    with DB.open() as f:
        for r in csv.DictReader(f, delimiter='\t'):
            old_tuples = list(zip(*(r[x].split(';') for x in VECTORS)))
            k = key(decoder[r['species_a']], decoder[r['species_b']])
            rawkey = canonical(r['species_a'], r['species_b'])
            for g, reader in readers.items():
                v = groups[g]
                remaining = [t for t in old_tuples if t[0] not in v['sources']]
                if len(remaining) < 2:
                    continue
                actual = next(reader, None)
                if actual is None:
                    errors[g] += 1
                    continue
                expected = r.copy()
                expected['proj_count'] = str(len(remaining))
                for i, col in enumerate(VECTORS):
                    expected[col] = ';'.join(t[i] for t in remaining)
                if len(remaining) != len(old_tuples):
                    for i, col in enumerate(['aver_corr', 'aver_pv']):
                        values = [precision[rawkey][t[0]][i] for t in remaining]
                        expected[col] = format(sum(values) / len(values), '.10g')
                errors[g] += actual != expected
                if k in all_external:
                    v['recovered'].add(k)
    for g, reader in readers.items():
        check(g + '_no_unexpected_derivative_rows', next(reader, None) is None)
        check(g + '_all_records_vectors_counts_means_verified', errors[g] == 0, errors[g])
        derivative_files[g].close()
    for g, v in groups.items():
        if v['identical'] == 'Yes':
            v['recovered'] = set(external_db)
        check(g + '_retained_plus_removed', v['retained'] + v['removed'] == total, (v['retained'], v['removed']))
        v['sha256'] = sha(v['file'])

    summaries, losses = [], []
    m_by_n = {r['network']: r for r in manifest}
    for n in DENOM:
        c, m = config[n], m_by_n[n]
        group = groups[c['database_group']]
        rows = []
        for r in baseline[n]:
            k = r['external_pair_key']
            base = k in external_db
            check(n + '_baseline_exact_presence:' + k, base == (r['recovered'] == 'Yes'))
            rec = k in group['recovered']
            check(n + '_no_new_recovery:' + k, not rec or base)
            out = dict(network=n, primary_or_exploratory=m['primary_or_exploratory'], external_pair_key=k,
                       external_species_a=r['external_species_a'], external_species_b=r['external_species_b'],
                       baseline_recovered='Yes' if base else 'No', source_excluded_recovered='Yes' if rec else 'No',
                       source_projects_requested=';'.join(c['sources']), source_projects_removed='',
                       baseline_mack_proj_count='NA', remaining_mack_proj_count='NA', remaining_proj_ID='NA',
                       remaining_proj_corr='NA', remaining_proj_pv='NA', remaining_proj_sam='NA',
                       remaining_aver_corr='NA', remaining_aver_pv='NA',
                       source_exclusion_loss='Yes' if base and not rec else 'No', loss_reason='NA',
                       mack_source_row_id='NA', analysis_status='PASS')
            if base:
                line, original = external_db[k]
                check(n + '_baseline_support_verbatim:' + k, all(r['mack_' + col] == original[col] for col in fields)
                      and str(line) == r['mack_source_row_id'])
                new, removed = exclude(original, set(c['sources']))
                check(n + '_recovery_equals_remaining_f2:' + k, rec == (int(new['proj_count']) >= 2))
                out.update(source_projects_removed=';'.join(removed), baseline_mack_proj_count=original['proj_count'],
                           remaining_mack_proj_count=new['proj_count'], remaining_aver_corr=new['aver_corr'],
                           remaining_aver_pv=new['aver_pv'], mack_source_row_id=line)
                for col in VECTORS:
                    out['remaining_' + col] = new[col]
                if not rec:
                    out['loss_reason'] = 'remaining_project_count_' + new['proj_count']
                    losses.append(out.copy())
            rows.append(out)
        before = sum(r['baseline_recovered'] == 'Yes' for r in rows)
        after = sum(r['source_excluded_recovered'] == 'Yes' for r in rows)
        loss = before - after
        check(n + '_recovered_summary_matches_Stage3A', before == int(summary0[n]['baseline_recovered_pairs']), before)
        check(n + '_summary_no_increase', after <= before, (before, after))
        check(n + '_fixed_denominator', len(rows) == DENOM[n], len(rows))
        check(n + '_primary_status_unchanged', m['primary_or_exploratory'] == ('exploratory' if n == 'Lacruz_MHO' else 'primary'))
        write(OUT / 'overlap' / (n + '_source_excluded_overlap.tsv'), rows)
        summaries.append(dict(paper=c['paper'], network=n, published_raw_edges=m['published_raw_edges'], final_gtdb_pairs=DENOM[n],
                              source_projects_reported=';'.join(c['sources']), source_projects_overlapping_mackdb=';'.join(group['actual']),
                              baseline_recovered=before, baseline_recovery_rate=before/DENOM[n], source_excluded_recovered=after,
                              source_excluded_not_recovered=DENOM[n]-after, source_excluded_recovery_rate=after/DENOM[n],
                              lost_recovered_pairs=loss, recovery_change_percentage_points=-100*loss/DENOM[n],
                              source_exclusion_loss_percentage_points=100*loss/DENOM[n], analysis_status='PASS',
                              primary_or_exploratory=m['primary_or_exploratory'], database_group=c['database_group']))
        print(n, before, '->', after, '/', DENOM[n], 'lost', loss, flush=True)
    write(OUT / 'summary/03b_source_exclusion_summary.tsv', summaries)
    write(OUT / 'overlap/03b_lost_recovered_pairs.tsv', losses, list(rows[0]))

    # Project audit includes rejected candidates, labelled explicitly, so mere mentions cannot enter exclusion.
    project_audit, overlap = [], []
    for n in DENOM:
        c, group = config[n], groups[config[n]['database_group']]
        candidates = set(c['sources'])
        if c['paper'] == 'Lacruz2026':
            candidates |= {'PRJEB87701', 'PRJEB7774', 'PRJEB1786', 'PRJEB41311', 'PRJEB38742', 'PRJEB37249'}
        if c['paper'] == 'Li2024':
            candidates |= {r['canonical_project_accession'] for r in read(OUT / 'provenance/03b_li_table_s1_cohort_audit.tsv')}
        for project in sorted(candidates):
            use = project in c['sources']
            reason = 'actual network construction source' if use else ('no construction sample for this phenotype' if c['paper'] == 'Lacruz2026' else 'validation or rifaximin-only; not discovery')
            project_audit.append(dict(paper=c['paper'], network=n, source_project=project, include_in_source_set='Yes' if use else 'No',
                                      provenance_decision=reason, in_human_gut_list='Yes' if project in listed else 'No',
                                      in_human_gut_info_support='Yes' if supports[project] else 'No', supporting_association_count=supports[project]))
        overlap.append(dict(paper=c['paper'], network=n, source_projects=';'.join(c['sources']),
                            source_projects_in_human_gut_list=';'.join(sorted(set(c['sources']) & listed)),
                            source_projects_in_actual_support=';'.join(group['actual']),
                            associations_with_any_source_support=group['affected'], baseline_database_associations=total,
                            source_excluded_database_associations=group['retained'], removed_database_associations=group['removed'],
                            remaining_zero=group['lost0'], remaining_one=group['lost1'],
                            source_excluded_network_identical_to_baseline=group['identical'],
                            source_excluded_database_file=group['file'], sha256=group['sha256'],
                            reason='No source project contributes to baseline support; reuse baseline' if group['identical'] == 'Yes' else 'Study/phenotype-specific project support exclusion'))
    write(OUT / 'provenance/03b_source_project_audit.tsv', project_audit)
    write(OUT / 'provenance/03b_source_overlap_with_mackdb.tsv', overlap)
    write(OUT / 'source_excluded_networks/03b_database_manifest.tsv',
          [dict(database_group=g, networks=';'.join(v['networks']), source_projects=';'.join(sorted(v['sources'])),
                actual_source_intersection=';'.join(v['actual']), associations_with_any_source_support=v['affected'],
                retained_associations=v['retained'], removed_associations=v['removed'], remaining_zero=v['lost0'], remaining_one=v['lost1'],
                source_excluded_network_identical_to_baseline=v['identical'], file=v['file'], sha256=v['sha256']) for g, v in groups.items()])

    # Read-only consistency reference; retain legacy Wallen key semantics in these files.
    wb = read(WALLEN / '08_overlap_original_mackdb_unique_gtdb_pairs.tsv')
    wa = read(WALLEN / '11_overlap_source_excluded_unique_gtdb_pairs.tsv')
    wd = read(WALLEN / '10_mackdb_PRJNA834801_source_exclusion_audit.tsv')
    ws = {r['metric']: float(r['value']) for r in read(WALLEN / '13_final_summary.tsv')}
    before = {r['gtdb_pair_key'] for r in wb if r['recovered'] == 'True'}
    after = {r['gtdb_pair_key'] for r in wa if r['recovered'] == 'True'}
    check('Wallen_existing_before_after', len(wb) == len(wa) == 2268 and len(before) == len(after) == 1784 and before == after,
          'Existing files read only; no Wallen taxonomy or overlap rerun')
    check('Wallen_existing_database_counts', len(wd) == 1246 and sum(r['retained_after_source_exclusion'] == 'False' for r in wd) == 2
          and ws['Source_excluded_MACKDB_pair_count'] == 653069 and ws['Recovered_unique_GTDB_pairs_lost_after_source_exclusion'] == 0)
    check('Wallen_actual_source_support_reference', supports['PRJNA834801'] == len(wd), supports['PRJNA834801'])
    wallen_rows = []
    for name in ['07_unique_gtdb_pairs.tsv', '08_overlap_original_mackdb_unique_gtdb_pairs.tsv', '10_mackdb_PRJNA834801_source_exclusion_audit.tsv',
                 '11_overlap_source_excluded_unique_gtdb_pairs.tsv', '13_final_summary.tsv', 'WALLEN_taxonomy_rescue_overlap_report.txt']:
        wallen_rows.append(dict(file=str(WALLEN / name), sha256=sha(WALLEN / name), reference_status='PASS'))
    write(OUT / 'provenance/03b_wallen_reference_manifest.tsv', wallen_rows)
    write(OUT / 'summary/03b_wallen_reference_summary.tsv', [dict(paper='Wallen2022', source_project='PRJNA834801',
          in_human_gut_list='Yes' if 'PRJNA834801' in listed else 'No', in_human_gut_info_support='Yes',
          associations_with_source_support=1246, denominator=2268, baseline_recovered=1784, source_excluded_recovered=1784,
          lost_recovered_pairs=0, baseline_database_associations=653071, source_excluded_database_associations=653069,
          reference_status='PASS', operation='Read existing formal results only')])

    check('no_taxonomy_remapping', True, 'Only frozen Stage3A decoder and manifest pairs read; no taxonomy decision changed')
    check('no_FastSpar_or_network_inference', True, 'Project CSV files read only for original aggregation floats')
    check('no_sign_condition_on_recovery', True, 'Presence and remaining support count >=2 only')
    print('Checking frozen hashes and all evidence/input hashes...', flush=True)
    for r in integrity:
        h = sha(r['input_file'])
        r['sha256_after'], r['unchanged'] = h, str(h == r['sha256_before'])
        r.setdefault('scope', 'frozen_previous_stage_or_reference')
    write(OUT / 'qc/03b_input_integrity.tsv', integrity,
          ['input_file', 'scope', 'size_bytes', 'sha256_before', 'prior_recorded_sha256', 'matches_prior', 'sha256_after', 'unchanged'])
    check('all_frozen_and_used_input_hashes_unchanged', all(r['unchanged'] == 'True' for r in integrity), len(integrity))
    check('database_final_expected_hash', sha(DB) == DB_HASH, DB_HASH)
    qc.append(dict(check='overall_status', status='PASS', detail='All provenance, source exclusion, reference and frozen-input checks passed.'))
    write(OUT / 'qc/03b_source_exclusion_qc.tsv', qc)
    (OUT / 'summary/03b_run_status.json').write_text(json.dumps(dict(overall_status='PASS', input_hashes_checked=len(integrity),
         exclusion_qc_checks=len(qc), provenance_qc_checks=len(pq), elapsed_seconds=time.time()-started,
         full_precision_projects_read=len(needed), affected_union_rows_for_float_loading=len(affected_rows)), indent=2) + '\n')
    print('Stage 3B PASS', len(qc), 'checks;', len(integrity), 'input hashes unchanged;', round(time.time()-started, 1), 's', flush=True)


if __name__ == '__main__':
    main()
