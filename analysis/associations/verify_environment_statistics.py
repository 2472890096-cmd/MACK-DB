#!/usr/bin/env python3
"""Compare recomputed Figure 2 statistics with the fixed summary/plot values.

This is numerical validation, not a plotting program. The original final
renderer is included unchanged in original/ and hash checked solely to reuse its actual
selected_environments() function (the renderer main() is never called).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path

RENDERER_SHA256 = "50dda8e6c2169aa7fd864a8d68f230ac77535529d4b79201a09ff1929956c14a"
DEFAULT_SELECTION_SOURCE = Path(__file__).resolve().parent / "original" / "plot_coabundance_panel_ABCDE_ARIAL.py"
COUNTS = {"fastspar_edges", "cclasso_edges", "common_edges", "fastspar_only", "cclasso_only", "union_edges"}
FIELDS = COUNTS | {"jaccard", "common_fraction_fastspar", "common_fraction_cclasso", "sign_agreement", "common_corr_spearman", "median_abs_corr_diff"}


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def compare(actual, expected, tolerance):
    if actual == expected:
        return True, 0.0
    if actual == "" or expected == "":
        return False, None
    a, b = float(actual), float(expected)
    if math.isnan(a) and math.isnan(b):
        return True, 0.0
    diff = abs(a - b)
    return math.isfinite(diff) and diff <= tolerance, diff


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--reference-summary", type=Path, required=True)
    parser.add_argument("--plot-values", type=Path, required=True)
    parser.add_argument("--selection-source", type=Path, default=DEFAULT_SELECTION_SOURCE,
                        help="Optional renderer override; the frozen SHA256 is still enforced")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--absolute-tolerance", type=float, default=1e-15)
    args = parser.parse_args()
    if hashlib.sha256(args.selection_source.read_bytes()).hexdigest() != RENDERER_SHA256:
        parser.error("Final renderer SHA256 differs from the frozen source")
    if args.output.exists():
        parser.error("Refusing an existing output directory")
    if args.absolute_tolerance < 0:
        parser.error("Tolerance cannot be negative")
    actual = read(args.summary)
    expected = read(args.reference_summary)
    keyed = {(r["rank"], r["env"]): r for r in actual}
    fixed = {(r["rank"], r["env"]): r for r in expected}
    if len(actual) != 90 or len(keyed) != 90 or set(keyed) != set(fixed):
        raise ValueError("Expected identical 90 rank-environment keys")
    comparisons = []
    for key, row in keyed.items():
        for field in sorted(FIELDS):
            ok, difference = compare(row[field], fixed[key][field], args.absolute_tolerance)
            comparisons.append(dict(panel="ALL_ENVIRONMENTS", rank=key[0], environment=key[1], field=field,
                                    actual=row[field], expected=fixed[key][field], absolute_difference=difference, match=ok))
    # Reuse the archived selection rule rather than substitute a new top-N rule.
    spec = importlib.util.spec_from_file_location("mackdb_final_figure2_renderer", args.selection_source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    frame = module.pd.read_csv(args.summary, sep="\t", float_precision="round_trip")
    ab, c, de = module.selected_environments(frame)
    envs = {"A": ab, "B": ab, "C": c, "D": de, "E": de}
    plot = read(args.plot_values)
    plot_keys = {(r["panel"], r["rank"], r["env"]) for r in plot}
    expected_keys = {(panel, rank, env) for panel, names in envs.items()
                     for rank in (["species"] if panel == "A" else ["genus"] if panel == "B" else ["species", "genus"])
                     for env in names}
    if len(plot) != 120 or len(plot_keys) != 120 or plot_keys != expected_keys:
        raise ValueError("Final plot rank/environment selection differs from recomputed selection")
    for row in plot:
        panel = row["panel"]
        fields = ["fastspar_only", "common_edges", "cclasso_only"] if panel in "AB" else {
            "C": ["jaccard"], "D": ["sign_agreement"], "E": ["common_corr_spearman"]}[panel]
        for field in fields:
            value = keyed[row["rank"], row["env"]][field]
            ok, difference = compare(value, row[field], args.absolute_tolerance)
            comparisons.append(dict(panel=panel, rank=row["rank"], environment=row["env"], field=field,
                                    actual=value, expected=row[field], absolute_difference=difference, match=ok))
    args.output.mkdir(parents=True, exist_ok=False)
    with (args.output / "VALUE_COMPARISONS.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(comparisons[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(comparisons)
    result = dict(status="PASS" if all(r["match"] for r in comparisons) else "FAIL",
                  summary_rows=len(actual), plot_rows=len(plot),
                  comparison_cells=len(comparisons), panel_numeric_cells=sum(r["panel"] != "ALL_ENVIRONMENTS" for r in comparisons),
                  mismatched_cells=sum(not r["match"] for r in comparisons),
                  maximum_absolute_difference=max((r["absolute_difference"] for r in comparisons if r["absolute_difference"] is not None), default=0),
                  missing_spearman_rows=sum(r["common_corr_spearman"] == "" for r in actual),
                  plot_rendered=False, source_selection_function="selected_environments",
                  input_sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [args.summary, args.reference_summary, args.plot_values, args.selection_source]})
    (args.output / "VALIDATION.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    if result["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
