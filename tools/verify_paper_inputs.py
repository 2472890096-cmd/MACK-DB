#!/usr/bin/env python3
"""Read-only SHA256/size verification of an extracted paper-input package."""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--root', type=Path, required=True,
                        help='Directory against which archive_path is resolved')
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    seen, total = set(), 0
    with args.manifest.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if not rows:
        raise ValueError('Empty input manifest')
    for row in rows:
        relative = Path(row['archive_path'] if 'archive_path' in row else row['path'])
        if relative.is_absolute() or '..' in relative.parts or relative.as_posix() in seen:
            raise ValueError('Unsafe or duplicate archive path')
        seen.add(relative.as_posix())
        path = (root / relative).resolve(strict=True)
        if root not in path.parents or not path.is_file():
            raise ValueError('Manifest path escapes package root or is not a file')
        expected_bytes = int(row['bytes'] if 'bytes' in row else row['size_bytes'])
        if path.stat().st_size != expected_bytes:
            raise ValueError('Size mismatch: ' + relative.as_posix())
        digest = hashlib.sha256()
        with path.open('rb') as handle:
            for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
                digest.update(block)
        if digest.hexdigest() != row['sha256']:
            raise ValueError('SHA256 mismatch: ' + relative.as_posix())
        total += expected_bytes
    print(json.dumps(dict(status='PASS', files=len(seen), bytes=total,
                          data_modified=False, scientific_computation=False)))


if __name__ == '__main__':
    main()
