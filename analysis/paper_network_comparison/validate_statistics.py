#!/usr/bin/env python3
"""Compare fresh Figure 4 statistics to the 40 final display values."""
import argparse
import csv
import json
import math
from pathlib import Path


def read(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def unique_index(records, fields, expected_count):
    result = {}
    if len(records) != expected_count:
        raise ValueError(f"Expected {expected_count} actual records; received {len(records)}")
    for row in records:
        key = tuple(row[field] for field in fields)
        if key in result:
            raise ValueError("Duplicate actual statistics key: " + repr(key))
        result[key] = row
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--a", type=Path, required=True)
    p.add_argument("--b", type=Path, required=True)
    p.add_argument("--c", type=Path, required=True)
    p.add_argument("--d", type=Path, required=True)
    p.add_argument("--reference", type=Path,
                   default=Path(__file__).parent / "reference/Figure4_expected_values.tsv")
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    if args.output.exists():
        raise ValueError("Refusing to overwrite validation output")
    aa = {key[0]: row for key, row in unique_index(
        read(args.a / "Figure4A_fixed_order_recomputed.tsv"), ["study"], 6).items()}
    bb = {key[0]: row for key, row in unique_index(
        read(args.b / "07_statistics/stage5B2_primary_summary.tsv"), ["network"], 5).items()}
    cc = {key[0]: row for key, row in unique_index(
        read(args.c / "data/stage6B_all_comparison_summary.tsv"), ["comparison_id"], 5).items()}
    dd = unique_index(read(args.d / "species_genus_overlap_comparison.tsv"),
                      ["comparison_id", "scope"], 20)
    expected = read(args.reference)
    keys = [(r["panel"], r["comparison_id"], r["series"]) for r in expected]
    if len(keys) != len(set(keys)) or len(keys) != 40:
        raise ValueError("Reference must contain exactly 40 unique display entries")
    checks = []
    for row in expected:
        panel, cid, series = row["panel"], row["comparison_id"], row["series"]
        numerator = denominator = ""
        if panel == "A":
            value = aa[cid]["cumulative_newly_recovered_pairs"]
            numerator = value
            denominator = aa[cid]["cumulative_reported_union_pairs"]
        elif panel == "B":
            if series == "Species":
                record = bb[cid]
                value = record["published_edge_recovery"]
                numerator, denominator = record["shared_edges"], record["published_eligible_edges"]
            else:
                record = dd[cid, "COMMON_SPECIES"]
                value = record["genus_published_pair_overlap"]
                numerator, denominator = record["genus_shared_pairs"], record["genus_published_pairs"]
        elif panel == "C":
            field = {"Published": "eligible_published_edges",
                     "MACK-DB complete project": "project_edges_all_input_taxa",
                     "MACK-DB common species": "project_edges_within_common_taxa",
                     "Shared": "shared_edges"}[series]
            value = cc[cid][field]
        elif panel == "D":
            scope = "ALL_MAPPED_SPECIES" if series == "MACK-DB complete project" else "COMMON_SPECIES"
            field = {"Published": "genus_published_pairs",
                     "MACK-DB complete project": "genus_mackdb_pairs",
                     "MACK-DB comparable set": "genus_mackdb_pairs",
                     "Shared": "genus_shared_pairs"}[series]
            value = dd[cid, scope][field]
        else:
            raise ValueError("Unknown panel")
        difference = abs(float(value) - float(row["expected_value"]))
        ok = math.isfinite(float(value)) and difference <= 1e-15
        for field, actual in [("numerator", numerator), ("denominator", denominator)]:
            if row[field] and int(row[field]) != int(actual):
                ok = False
        checks.append(dict(panel=panel, comparison_id=cid, series=series,
                           expected_value=row["expected_value"], actual_value=value,
                           absolute_difference=difference, status="PASS" if ok else "FAIL"))
    args.output.mkdir(parents=True)
    with (args.output / "Figure4_40_values_comparison.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(checks[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(checks)
    state = dict(status="PASS" if all(r["status"] == "PASS" for r in checks) else "FAIL",
                 checked_values=len(checks), maximum_absolute_difference=max(r["absolute_difference"] for r in checks),
                 renderer_executed=False)
    (args.output / "VALIDATION.json").write_text(json.dumps(state, indent=2) + "\n")
    print(json.dumps(state))
    if state["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
