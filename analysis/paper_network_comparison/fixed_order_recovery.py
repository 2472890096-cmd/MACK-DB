#!/usr/bin/env python3
"""Recheck Figure 4A from frozen publication pairs and project-support records.

This transparent fixed-order set check does not rerun the historical 720-order
producer. Input is its pre-existing processed lookup, not sequencing reads.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

ORDER = ["Chen", "Wallen", "Li", "Guo", "Lacruz", "Zhou"]


def rows(path):
    with path.open(newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def sha(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def pair(row):
    a, b = row["canonical_species_a"], row["canonical_species_b"]
    if not a < b:
        raise ValueError("Noncanonical or self pair")
    return a, b


def archived_steps(path):
    saved_rows = [row for row in rows(path)
                  if row["analysis"] == "primary_all_known_source_excluded"
                  and row["order_id"] == "22"]
    saved = {int(row["position"]): row for row in saved_rows}
    if len(saved_rows) != len(saved) or set(saved) != set(range(1, 7)):
        raise ValueError("Missing/duplicate archived fixed-order steps")
    return saved


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairsets-dir", type=Path, required=True)
    parser.add_argument("--support-lookup", type=Path, required=True)
    parser.add_argument("--source-exclusion", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True,
                        help="Saved 720-order table; only order 22 is read for comparison")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("Output already exists")
    files = [args.support_lookup, args.source_exclusion, args.archive]
    sets = {}
    for study in ORDER:
        path = args.pairsets_dir / (study + "_frozen_pairs.tsv")
        records = list(rows(path))
        sets[study] = {pair(row) for row in records}
        if len(records) != len(sets[study]):
            raise ValueError("Duplicate publication pair: " + study)
        files.append(path)
    exclusions = {row["source_BioProject"] for row in rows(args.source_exclusion)
                  if row["excluded_from_primary"] == "True"}
    if len(exclusions) != 9:
        raise ValueError("Expected the frozen nine-source exclusion set")
    published_union = set.union(*sets.values())
    present, recovered = set(), set()
    for row in rows(args.support_lookup):
        key = pair(row)
        if key in present or key not in published_union:
            raise ValueError("Duplicate or extraneous original-database lookup pair")
        present.add(key)
        ids = row["original_proj_ID"].split(";")
        if len(ids) != len(set(ids)) or len(ids) != int(row["original_proj_count"]):
            raise ValueError("Invalid original project-support vector")
        kept = [project for project in ids if project not in exclusions]
        if kept != row["remaining_proj_ID"].split(";") and not (
                not kept and row["remaining_proj_ID"] == ""):
            raise ValueError("Frozen exclusion-vector disagreement")
        if len(kept) != int(row["remaining_proj_count"]):
            raise ValueError("Frozen exclusion-count disagreement")
        accepted = len(kept) >= 2
        if accepted != (row["in_primary_D0"] == "True"):
            raise ValueError("Frozen source-filtered membership disagreement")
        if accepted:
            recovered.add(key)
    saved = archived_steps(args.archive)
    union, found = set(), set()
    result = []
    for position, study in enumerate(ORDER, 1):
        newly = (sets[study] & recovered) - found
        union.update(sets[study])
        found.update(newly)
        record = dict(position=position, study=study,
                      newly_recovered_pairs=len(newly),
                      cumulative_newly_recovered_pairs=len(found),
                      cumulative_reported_union_pairs=len(union),
                      cumulative_union_recovery_fraction=len(found) / len(union))
        old = saved[position]
        if old["study"] != study or any(int(old[field]) != record[field] for field in [
                "newly_recovered_pairs", "cumulative_newly_recovered_pairs",
                "cumulative_reported_union_pairs"]):
            raise ValueError("First differing fixed-order step: " + str(position))
        result.append(record)
    args.output.mkdir(parents=True)
    with (args.output / "Figure4A_fixed_order_recomputed.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(result)
    state = dict(status="PASS", validation="FIXED_ORDER_SET_RECOMPUTATION_FROM_PROCESSED_SUPPORT_INPUT",
                 historical_720_orders_rerun=False, original_inference_rerun=False,
                 publication_union=len(published_union), recovered_union=len(recovered),
                 source_exclusions=sorted(exclusions), all_six_archived_steps_match=True,
                 unlisted_pairs="Absent from the audited original f>=2 lookup; sub-f2 support unknown",
                 inputs=[dict(path=str(p.resolve()), bytes=p.stat().st_size, sha256=sha(p)) for p in files])
    (args.output / "VALIDATION.json").write_text(json.dumps(state, indent=2) + "\n")
    print(json.dumps({k: state[k] for k in ["status", "publication_union", "recovered_union"]}))


if __name__ == "__main__":
    main()
