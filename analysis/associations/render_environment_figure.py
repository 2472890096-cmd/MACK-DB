#!/usr/bin/env python3
"""Render the final Figure 2 from its fixed summary using the original source.

Only output/font paths and the original --summary argument are configured.
No statistics, thresholds, labels or selection functions are substituted.
The held renderer and separately licensed Arial fonts are not bundled.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

RENDERER_SHA256 = "50dda8e6c2169aa7fd864a8d68f230ac77535529d4b79201a09ff1929956c14a"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--renderer-source", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--font-dir", type=Path, required=True,
                        help="Separately obtained licensed Arial TTF directory")
    parser.add_argument("--output", type=Path, required=True,
                        help="New isolated output directory")
    args = parser.parse_args()
    source = args.renderer_source.resolve(strict=True)
    summary = args.summary.resolve(strict=True)
    fonts = args.font_dir.resolve(strict=True)
    if hashlib.sha256(source.read_bytes()).hexdigest() != RENDERER_SHA256:
        parser.error("Final renderer SHA256 differs from the frozen source")
    required_fonts = [fonts / n for n in ("Arial.TTF", "Arialbd.TTF", "Ariali.TTF", "Arialbi.TTF")]
    if any(not p.is_file() for p in required_fonts):
        parser.error("The original renderer requires all four named Arial TTF files")
    output = args.output.resolve()
    if output.exists():
        parser.error("Refusing an existing output directory")
    spec = importlib.util.spec_from_file_location("mackdb_final_figure2_renderer", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    output.mkdir(parents=True, exist_ok=False)
    module.OUTDIR = output
    module.FONT_DIR = fonts
    original_argv = sys.argv
    start = time.perf_counter()
    try:
        sys.argv = [str(source), "--summary", str(summary)]
        module.main()
    finally:
        sys.argv = original_argv
    record = dict(status="RENDER_COMPLETED", seconds=time.perf_counter() - start,
                  renderer_source_sha256=RENDERER_SHA256,
                  summary_sha256=hashlib.sha256(summary.read_bytes()).hexdigest(),
                  font_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in required_fonts},
                  changed_configuration=["OUTDIR", "FONT_DIR", "--summary"],
                  new_association_statistics_computed=False)
    (output / "RENDER_EXECUTION.json").write_text(json.dumps(record, indent=2) + "\n")


if __name__ == "__main__":
    main()
