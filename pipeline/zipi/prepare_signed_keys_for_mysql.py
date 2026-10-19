#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import csv
from pathlib import Path

SPECIES_DIR = Path("work/zipi/species_all_signed/keyspecies")
GENUS_DIR = Path("work/zipi/genus_all_signed/keygenus")
OUT_DIR = Path("work/zipi/mysql_import")

COLUMNS = [
    "environment", "bioproject", "name", "degree",
    "betweenness_centrality", "closeness_centrality",
    "eigenvector_centrality", "Abundance", "module", "z", "p",
    "taxa_roles", "positive_degree", "negative_degree",
    "within_positive_degree", "within_negative_degree",
    "between_positive_degree", "between_negative_degree",
    "positive_strength", "negative_strength_abs",
    "Domain", "Phylum", "Class", "order_name",
    "Family", "Genus", "Species"
]

ALIASES = {
    "order_name": ["order_name", "Order", "order"],
    "taxa_roles": ["taxa_roles", "taxa_role", "roles", "role"],
}


def get_value(row, target):
    if target in row:
        return (row.get(target) or "").strip()
    for candidate in ALIASES.get(target, []):
        if candidate in row:
            return (row.get(candidate) or "").strip()
    return ""


def merge_level(level, root, output):
    if not root.is_dir():
        raise FileNotFoundError(f"{level}目录不存在：{root}")

    key_files = sorted(root.glob("*/*_key.txt"))
    merged = []
    seen = {}
    conflicts = []
    empty_files = 0
    exact_duplicates = 0

    for key_file in key_files:
        environment = key_file.parent.name
        project = key_file.name.removesuffix("_key.txt")

        with key_file.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames:
                empty_files += 1
                continue

            valid_rows = 0
            for source in reader:
                if not any((v or "").strip() for v in source.values()):
                    continue

                name = get_value(source, "name")
                if not name:
                    continue

                row = {"environment": environment, "bioproject": project}
                for column in COLUMNS[2:]:
                    row[column] = get_value(source, column)

                key = (environment, project, name)
                if key in seen:
                    if seen[key] == row:
                        exact_duplicates += 1
                    else:
                        conflicts.append((environment, project, name, str(key_file)))
                    continue

                seen[key] = row
                merged.append(row)
                valid_rows += 1

            if valid_rows == 0:
                empty_files += 1

    merged.sort(key=lambda r: (r["environment"], r["bioproject"], r["name"]))

    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(merged)

    conflict_path = output.with_name(output.stem + "_conflicts.tsv")
    with conflict_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["environment", "bioproject", "name", "source_file"])
        writer.writerows(conflicts)

    print(
        f"{level}: key文件={len(key_files)}, 空文件={empty_files}, "
        f"输出行={len(merged)}, 完全重复={exact_duplicates}, 冲突={len(conflicts)}"
    )
    print(f"输出：{output}")
    print(f"冲突报告：{conflict_path}")

    if conflicts:
        raise RuntimeError(f"{level}发现冲突记录，请检查：{conflict_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-dir", type=Path, default=SPECIES_DIR)
    parser.add_argument("--genus-dir", type=Path, default=GENUS_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--level", choices=["species", "genus", "both"], default="both")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.level in ("species", "both"):
        merge_level("species", args.species_dir, args.output_dir / "key_species_new.tsv")

    if args.level in ("genus", "both"):
        merge_level("genus", args.genus_dir, args.output_dir / "key_genus_new.tsv")


if __name__ == "__main__":
    main()
