#!/usr/bin/env python3
"""Read-only integrity and local-file dependency checks for the private snapshot."""
from pathlib import Path
import collections
import csv
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    scripts = read(ROOT/'docs/ACCEPTED_84_SCRIPTS.tsv')
    require(len(scripts) == len({r['path'] for r in scripts}) == 84, 'Fixed script count')
    counts = collections.Counter(r['purpose'] for r in scripts)
    require(counts == {'KEEP_MAIN':35, 'KEEP_SUPPLEMENTARY':45, 'KEEP_TEST':4}, 'Purpose counts')
    for row in scripts:
        require((ROOT/row['path']).is_file() and sha(ROOT/row['path']) == row['sha256'],
                'Accepted script changed/missing: '+row['path'])
    aliases = [
        ('pipeline/epi_v3/core/get_keystone_EPI_single_safe_qecdf_v3_fast.R',
         'pipeline/epi_v3/genus/scripts/get_keystone_EPI_single_safe_qecdf_v3_fast.R'),
        ('pipeline/epi_v3/species/scripts/workflow.py',
         'pipeline/epi_v3/genus/scripts/species_production_workflow.py'),
    ]
    for left, right in aliases:
        require(sha(ROOT/left) == sha(ROOT/right), 'Compatibility copies disagree')
    availability = read(ROOT/'docs/FIGURE3_INPUT_AVAILABILITY.tsv')
    require(len(availability) == 70, 'Figure3 availability count')
    included = [r for r in availability if r['decision']=='INCLUDE_PRIVATE']
    require(len(included) == 6, 'Approved Figure3 subset count')
    for row in included:
        require(sha(ROOT/row['path']) == row['sha256'], 'Frozen approved reference mismatch')
    for row in availability:
        if row['decision'] != 'INCLUDE_PRIVATE':
            require(not (ROOT/row['path']).exists(), 'Authorization-held data unexpectedly bundled: '+row['path'])
    scopefiles = {
        1571:'COMMON_PROJECTS_SCIENTIFICALLY_USABLE.txt',
        1570:'COMMON_PROJECTS_SCIENTIFICALLY_USABLE_EXCLUDING_UNVERIFIED_DKI.txt',
        1568:'COMMON_PROJECTS_STRICT_SUCCESS.txt',
    }
    scopes = {}
    for n, name in scopefiles.items():
        values = (ROOT/'frozen_inputs/figure3/scopes'/name).read_text().split()
        scopes[n] = set(values)
        require(len(values) == len(scopes[n]) == n, 'Scope count/duplicates')
        require('PRJEB51353' not in scopes[n], 'Unexpected Figure3 project')
    require(scopes[1568] < scopes[1570] < scopes[1571], 'Scope nesting')
    deps = [r for r in read(ROOT/'DEPENDENCY_RESOLUTION.tsv')
            if r['dependency_type']=='RETAINED_REPOSITORY_CODE']
    for row in deps:
        require((ROOT/row['dependency']).is_file(), 'Declared local dependency missing: '+row['dependency'])
    for name in ['lib/mackdb_paths.py','tools/epi_project.py','tools/network_matrix_io.py',
                 'tests/run_small_tests.py','tests/unit/epi.R','tests/unit/dki.R','tests/unit/zipi.R',
                 'pipeline/epi_v3/species/audit/config.json','pipeline/epi_v3/genus/audit/config.json',
                 'pipeline/epi_v3/genus/audit/species_reference_config.json']:
        require((ROOT/name).is_file(), 'Runtime/test helper missing: '+name)
    fixtures = read(ROOT/'tests/fixtures/FIXTURE_MANIFEST.tsv')
    for row in fixtures:
        require(sha(ROOT/row['path']) == row['sha256'], 'Synthetic fixture mismatch: '+row['path'])
    print(json.dumps(dict(status='PASS',accepted_scripts=84,main=35,local_dependency_edges=len(deps),
                          synthetic_fixture_files=len(fixtures),figure3_approved_files=6,
                          figure3_held_files=64,scientific_computation_executed=False), indent=2))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError) as exc:
        print('FAIL: '+str(exc), file=sys.stderr)
        sys.exit(1)
