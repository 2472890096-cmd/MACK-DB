#!/usr/bin/env python3
"""Recompute concordance from an already unified-R six-table candidate snapshot.

Delegates all scientific set calculations to the unchanged accepted release.
This entry never applies the historical four-project replacement a second time.
It does not infer project eligibility or scientific run success from candidates.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

CORE = Path(__file__).with_name('concordance_release.py')
CORE_SHA256 = '4866756e5494048188d69409ef95c2e20f2194fcb98c2756d4ed9540186e95c1'
if hashlib.sha256(CORE.read_bytes()).hexdigest() != CORE_SHA256:
    raise ValueError('Accepted concordance core SHA256 mismatch')
import concordance_release as core


def read(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write(path, rows, fields):
    with Path(path).open('w') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidates-dir', type=Path, required=True,
                        help='FINAL, already R-corrected six key_*.tsv files')
    parser.add_argument('--r-candidates', type=Path, required=True,
                        help='Accepted61 rows; assertion only, never replacement')
    parser.add_argument('--scopes-dir', type=Path, required=True)
    parser.add_argument('--expected-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), 'Use a new output directory')
    sources, environments, input_files = {}, {}, []
    for level in ['species', 'genus']:
        for method, suffix in zip(core.methods, ['', '_epi', '_dki']):
            path = args.candidates_dir / f'key_{level}{suffix}.tsv'
            rows = read(path)
            values = {(r['bioproject'], r['name'], r['environment']) for r in rows}
            require(len(values) == len(rows), f'Duplicate candidate records: {path.name}')
            for project, taxon, environment in values:
                require(all((project, taxon, environment)) and environment.lower() != 'unknown',
                        'Empty name/project/environment or unknown environment')
                require(environments.setdefault(project, environment) == environment,
                        'Inconsistent project environment: ' + project)
            sources[level, method] = values
            input_files.append(path)
    replacement = read(args.r_candidates)
    accepted = {(r['bioproject'], r['name'], r['environment']) for r in replacement}
    require(len(replacement) == len(accepted) == 61, 'Accepted R reference must contain61 unique rows')
    require({p for p, _, _ in accepted} == core.projects, 'Unexpected R reference projects')
    actual = {r for r in sources['species', 'DKI'] if r[0] in core.projects}
    require(actual == accepted,
            'Species DKI is not already the accepted unified-R snapshot; use the historical entry for old data')
    reference_files = [args.expected_dir / 'FIGURE3_UNIFIED_R_22_VALUES.tsv',
                       args.expected_dir / 'FIGURE3_UNIFIED_R_SENSITIVITY.tsv']
    expected_rows = [r for p in reference_files for r in read(p)]
    expected = {(int(r['scope']), r['panel'], r['category'], r['taxonomic_level']):
                int(r['unified_R_value']) for r in expected_rows}
    require(len(expected_rows) == len(expected) == 66, 'Expected66 unique main/sensitivity values')
    output_rows, support, coverage, details, universes = [], [], [], [], {}
    for count, filename in core.scopes.items():
        path = args.scopes_dir / filename
        items = path.read_text().split()
        universe = set(items)
        require(len(items) == len(universe) == count, 'Scope count/duplicates')
        require('PRJEB51353' not in universe and core.projects <= universe, 'Unexpected scope membership')
        universes[count] = universe
        selected = {key: {r for r in values if r[0] in universe}
                    for key, values in sources.items()}
        values, contexts = core.calc(selected)
        require(len(values) == 22, 'Expected22 values per scope')
        for key, value in values.items():
            require(expected[(count, *key)] == value, f'Reference mismatch: {count}, {key}, {value}')
            panel, category, level = key
            output_rows.append(dict(scope=count, panel=panel, category=category,
                                    taxonomic_level=level, value=value))
        for level in ['species', 'genus']:
            for project in sorted(universe):
                counts = {m: sum(p == project for p, _, _ in selected[level, m]) for m in core.methods}
                coverage.append(dict(scope=count, taxonomic_level=level, bioproject=project,
                                     environment=environments.get(project, ''), **counts,
                                     candidate_state='HAS_CANDIDATES' if any(counts.values()) else 'NO_CANDIDATES',
                                     eligibility='ACCEPTED_SCOPE_LIST_NOT_INFERRED_FROM_CANDIDATES'))
            for mode in ['Broad', 'Strict']:
                for (context, taxon), row in sorted(contexts[level, mode].items()):
                    item = dict(scope=count, taxonomic_level=level, mode=mode, name=taxon,
                                environment=row['environment'], bioproject=context if mode == 'Strict' else '',
                                combination=core.cats[row['mask']], method_count=row['mask'].bit_count(),
                                **{m + '_projects': ';'.join(sorted(row['sources'][m])) for m in core.methods})
                    details.append(item)
                    if mode == 'Strict' and row['mask'] == 7:
                        support.append(item)
        input_files.append(path)
    require(universes[1568] < universes[1570] < universes[1571], 'Scope nesting')
    args.output.mkdir(parents=True)
    write(args.output / 'FIGURE3_ALL_SCOPES_22_VALUES.tsv', output_rows, list(output_rows[0]))
    write(args.output / 'ALL_COMBINATIONS.tsv', details, list(details[0]))
    write(args.output / 'STRICT_ALL3_NAMES_AND_SUPPORT.tsv', support, list(support[0]))
    write(args.output / 'PROJECT_SCOPE_AND_CANDIDATES.tsv', coverage, list(coverage[0]))
    inputs = [dict(file=p.name, bytes=p.stat().st_size, sha256=hashlib.sha256(p.read_bytes()).hexdigest())
              for p in input_files + [args.r_candidates] + reference_files]
    (args.output / 'VALIDATION.json').write_text(json.dumps(dict(
        status='PASS', input_state='FINAL_UNIFIED_R', replacement_applied=False,
        exact_values=66, independent_set_checks='inside unchanged core.calc',
        scope_from_candidates=False, scientific_inference_run=False, inputs=inputs), indent=2) + '\n')
    print('PASS:66 values from FINAL unified-R candidates; no replacement applied')


if __name__ == '__main__':
    main()
