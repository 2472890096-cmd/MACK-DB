#!/usr/bin/env python3
"""Read-only, versioned source/data-boundary guard for the phase-one candidate.
The original84 manifest remains immutable historical evidence; the current
manifest protects every released script, including helpers and this guard.
"""
from pathlib import Path
import collections
import csv
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
HISTORICAL = {
    'docs/provenance/ACCEPTED_84_SCRIPTS_ORIGINAL.tsv': '9c7bd6967fada1e2e64cbabd16914cabc1ecbdb243442de936e9c99fc9e1d800',
    'docs/provenance/BASELINE_93_SCRIPTS.tsv': 'e88dc91d2683c30cbaf68342bed38ca37bdbf33fca5eeda5470660cba4074ff5',
    'docs/provenance/BASELINE_151_TRACKED_FILES.tsv': 'cb16d3f603741ed4894179a3417a32f29194ae6727d418766fa1e890b89460f3',
}


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def included_file(path):
    rel = path.relative_to(ROOT)
    return path.is_file() and not any(x in {'.git', '__pycache__'} for x in rel.parts) and not (
        rel.parts[0] in {'data', 'work', 'external', 'audit', 'models'} or
        rel.parts[:2] == ('tests', 'results'))


def main():
    for path, digest in HISTORICAL.items():
        require(sha(ROOT/path) == digest, 'Historical manifest changed: '+path)
    old = read(ROOT/'docs/provenance/ACCEPTED_84_SCRIPTS_ORIGINAL.tsv')
    require(len(old) == len({r['path'] for r in old}) == 84, 'Historical84 count')
    require(collections.Counter(r['purpose'] for r in old) ==
            {'KEEP_MAIN':35, 'KEEP_SUPPLEMENTARY':45, 'KEEP_TEST':4}, 'Historical roles')
    release = json.loads((ROOT/'config/release_integrity.json').read_text())
    scripts = read(ROOT/'docs/RELEASE_SCRIPT_MANIFEST.tsv')
    require(len(scripts) == len({r['path'] for r in scripts}) == release['script_count'],
            'Current script count/duplicates')
    for row in scripts:
        require(sha(ROOT/row['path']) == row['sha256'], 'Current script changed: '+row['path'])
    discovered = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*')
                  if included_file(p) and p.suffix.lower() in {'.py', '.r', '.sh'}}
    require(discovered == {r['path'] for r in scripts}, 'Unregistered or missing script')
    disposition = read(ROOT/'docs/provenance/ORIGINAL_93_DISPOSITION.tsv')
    base = {r['path']:r for r in read(ROOT/'docs/provenance/BASELINE_93_SCRIPTS.tsv')}
    require(len(disposition) == 93 and {r['original_path'] for r in disposition} == set(base),
            'Original93 coverage')
    require(sum(r['action'] == 'REMOVED_FROM_DEFAULT' for r in disposition) == 17,
            'Reviewed removal count')
    for row in disposition:
        require(row['original_sha256'] == base[row['original_path']]['sha256'],
                'Original source hash mismatch: '+row['original_path'])
        if row['action'] == 'REMOVED_FROM_DEFAULT':
            require(not (ROOT/row['original_path']).exists(),
                    'Removed default script unexpectedly present: '+row['original_path'])
        else:
            require((ROOT/row['candidate_path']).is_file(),
                    'Kept original script missing: '+row['candidate_path'])
            # Phase1 authorizes a test import literal and the version guard;
            # every other retained original script must keep its original bytes.
            if row['original_path'] not in {
                'tests/run_small_tests.py', 'tools/verify_private_checkout.py'
            }:
                require(sha(ROOT/row['candidate_path']) == base[row['original_path']]['sha256'],
                        'Retained original source changed: '+row['original_path'])
    files = read(ROOT/'docs/RELEASE_FILE_MANIFEST.tsv')
    require(len(files) == len({r['path'] for r in files}) == release['file_count'],
            'Current file manifest count/duplicates')
    for row in files:
        require(sha(ROOT/row['path']) == row['sha256'], 'Release file changed: '+row['path'])
    aliases = [
        ('pipeline/epi_v3/core/get_keystone_EPI_single_safe_qecdf_v3_fast.R',
         'pipeline/epi_v3/genus/scripts/get_keystone_EPI_single_safe_qecdf_v3_fast.R'),
        ('pipeline/epi_v3/species/scripts/workflow.py',
         'pipeline/epi_v3/genus/scripts/species_production_workflow.py'),
    ]
    for left, right in aliases:
        require(sha(ROOT/left) == sha(ROOT/right), 'EPI compatibility copies disagree')
    availability = read(ROOT/'docs/FIGURE3_INPUT_AVAILABILITY.tsv')
    require(len(availability) == 70, 'Frozen input availability count')
    supplied = [r for r in availability if r['decision'] == 'INCLUDE_PRIVATE']
    require(len(supplied) == 6, 'Approved frozen subset count')
    for row in supplied:
        require(sha(ROOT/row['path']) == row['sha256'], 'Frozen reference changed')
    for row in availability:
        if row['decision'] != 'INCLUDE_PRIVATE':
            require(not (ROOT/row['path']).exists(), 'Held input unexpectedly bundled: '+row['path'])
    scopes = {}
    for n, name in {
        1571:'COMMON_PROJECTS_SCIENTIFICALLY_USABLE.txt',
        1570:'COMMON_PROJECTS_SCIENTIFICALLY_USABLE_EXCLUDING_UNVERIFIED_DKI.txt',
        1568:'COMMON_PROJECTS_STRICT_SUCCESS.txt',
    }.items():
        values = (ROOT/'frozen_inputs/figure3/scopes'/name).read_text().split()
        scopes[n] = set(values)
        require(len(values) == len(scopes[n]) == n, 'Scope count/duplicates')
        require('PRJEB51353' not in scopes[n], 'Unexpected scope project')
    require(scopes[1568] < scopes[1570] < scopes[1571], 'Scope nesting')
    deps = [r for r in read(ROOT/'DEPENDENCY_RESOLUTION.tsv')
            if r['dependency_type'] == 'RETAINED_REPOSITORY_CODE']
    for row in deps:
        require((ROOT/row['dependency']).is_file(), 'Local dependency missing: '+row['dependency'])
    for name in ['lib/mackdb_paths.py', 'tools/epi_project.py', 'tools/network_matrix_io.py',
                 'tests/run_small_tests.py', 'tests/unit/epi.R', 'tests/unit/dki.R',
                 'tests/unit/zipi.R', 'pipeline/epi_v3/species/audit/config.json',
                 'pipeline/epi_v3/genus/audit/config.json',
                 'pipeline/epi_v3/genus/audit/species_reference_config.json']:
        require((ROOT/name).is_file(), 'Runtime/test helper missing: '+name)
    for row in read(ROOT/'tests/fixtures/FIXTURE_MANIFEST.tsv'):
        require(sha(ROOT/row['path']) == row['sha256'], 'Synthetic fixture changed')
    print(json.dumps(dict(status='PASS', release=release['release_id'],
                          historical84_verified=True, original93_accounted=True,
                          scripts=len(scripts), protected_files=len(files),
                          removed_default_scripts=17, approved_frozen_files=6,
                          held_frozen_files=64, scientific_computation_executed=False), indent=2))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError) as exc:
        print('FAIL: '+str(exc), file=sys.stderr)
        sys.exit(1)
