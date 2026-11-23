#!/usr/bin/env python3
"""Run the archived Figure 2 statistics without its historical plotting paths.

The original source is included unchanged in original/. Its SHA256 is checked
before importing its actual read_info(),
summarize_env(), build_summaries() and global_summary() implementations. This
adapter does not replace their scientific calculations or render final artwork.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

SOURCE_SHA256 = "48920f30e40f7212fb73c3394b6718a5634acd6d5d5f4554fa5d6f97f29dfb9c"
DEFAULT_SOURCE = Path(__file__).resolve().parent / "original" / "plot_species_genus_coabundance_complete.py"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--statistics-source", type=Path, default=DEFAULT_SOURCE,
                        help="Optional source override; the frozen SHA256 is still enforced")
    for name in ("species-fastspar", "genus-fastspar", "species-cclasso", "genus-cclasso"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True,
                        help="New output directory; existing directories are refused")
    args = parser.parse_args()
    source = args.statistics_source.resolve(strict=True)
    if sha256(source) != SOURCE_SHA256:
        parser.error("Statistics source SHA256 does not match the archived producer")
    dirs = {
        ("species", "FastSpar"): args.species_fastspar.resolve(strict=True),
        ("genus", "FastSpar"): args.genus_fastspar.resolve(strict=True),
        ("species", "CCLasso"): args.species_cclasso.resolve(strict=True),
        ("genus", "CCLasso"): args.genus_cclasso.resolve(strict=True),
    }
    if any(not path.is_dir() for path in dirs.values()):
        parser.error("All four association inputs must be directories")
    environments = [{p.stem for p in path.glob("*.info")} for path in dirs.values()]
    if any(len(s) != 45 for s in environments) or any(s != environments[0] for s in environments):
        parser.error("The fixed paper input requires the same 45 environments in all four directories")
    output = args.output.resolve()
    if output.exists():
        parser.error("Refusing an existing output directory")
    # Check all inputs before importing the original module; no import-time main().
    input_manifest = []
    for (rank, method), directory in dirs.items():
        for path in sorted(directory.glob("*.info")):
            input_manifest.append(dict(rank=rank, method=method, environment=path.stem,
                                       path=str(path), bytes=path.stat().st_size,
                                       sha256=sha256(path)))
    spec = importlib.util.spec_from_file_location("mackdb_archived_environment_statistics", source)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load the archived statistics source")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    output.mkdir(parents=True, exist_ok=False)
    start = time.perf_counter()
    summary, support = module.build_summaries(dirs)
    global_summary = module.global_summary(summary)
    summary["env_display"] = summary["env"].map(module.disp_env)
    summary.to_csv(output / "environment_method_summary.tsv", sep="\t", index=False)
    support.to_csv(output / "support_count_distribution.tsv", sep="\t", index=False)
    global_summary.to_csv(output / "global_rank_summary.tsv", sep="\t", index=False)
    with (output / "INPUT_SHA256.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(input_manifest[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(input_manifest)
    result = dict(status="STATISTICS_COMPLETED", statistics_source=str(source),
                  statistics_source_sha256=SOURCE_SHA256, input_files=len(input_manifest),
                  input_bytes=sum(r["bytes"] for r in input_manifest),
                  environment_rank_rows=len(summary), statistics_seconds=time.perf_counter() - start,
                  python=sys.version, plot_rendered=False,
                  validation_boundary="Real environment association tables to statistics; means are read from aver_corr")
    (output / "EXECUTION.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
