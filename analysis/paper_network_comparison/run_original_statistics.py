#!/usr/bin/env python3
"""Execute frozen Figure 4 producers with explicit read-only input/output roots.

Scientific functions are compiled unchanged from SHA-pinned source files. Only
module-level path/state initialization is supplied here. No inference, training,
or null simulation is invoked. The accepted original sources are bundled under
original/ and are checked against their approved SHA256 identities.
The B source location/data-root initializer was adapted without changing its functions.
"""
import argparse
import ast
import csv
import hashlib
import json
import shutil
import sys
import time
import traceback
from collections import defaultdict
from pathlib import Path

PINS = {
    "b": "1ccc1687dfaeb7bea2faa8c3c3c8f61899c29f9d35b1c99212ed4d026d58b30b",
    "c": "2b704d1fc35d16411995f904e03a8312a6cef100d4c889c9b779ea6eb02e8c86",
    "d": "c4e898425a991c429100002e83e6d0c41f128dbc9a8abc9522f02bdb6f718d25",
    "c_common": "aa33e14b18bb89e23626f352c0d0de7bf8aee520321d3a851237d9669366b10c",
    "d_common": "1a698efe6b84e5309a09bc251c1dd2d986b73a8c5fdddf611a182b74a43df2ce",
}


def default_source(mode):
    """Resolve the bundled producer without using author-server paths."""
    here = Path(__file__).resolve().parent
    return {
        "b": here / "original/matched_subgroups/stage5b2_analysis.py",
        "c": here / "original/complete_project/compare_networks.py",
        "d": here / "original/genus_groups/genus_overlap.py",
    }[mode]


def digest(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def read(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_functions(path, expected, namespace):
    """Keep all function ASTs intact; omit original top-level path/IO actions."""
    if digest(path) != expected:
        raise ValueError(f"Unapproved source version: {path}")
    tree = ast.parse(path.read_text(), filename=str(path))
    nodes = [node for node in tree.body if isinstance(node, (ast.Import, ast.ImportFrom,
                                                            ast.FunctionDef))
             and not (isinstance(node, ast.ImportFrom)
                      and node.module in {"common", "mackdb_paths"})]
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), "exec"), namespace)
    return tree


def literal_assignment(tree, name, namespace):
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == name
                                                for t in node.targets):
            return eval(compile(ast.Expression(node.value), "<frozen initializer>", "eval"), namespace)
    raise ValueError("Missing frozen initializer: " + name)


def resolve_config(path):
    cfg = json.loads(path.read_text())
    result = {}
    for key, value in cfg.items():
        if key == "path_rebase":
            result[key] = {old.rstrip("/"): str((Path(new) if Path(new).is_absolute()
                                                 else path.parent / new).resolve())
                           for old, new in value.items()}
        else:
            result[key] = (Path(value) if Path(value).is_absolute() else path.parent / value).resolve()
    return result


def install_io_rebase(namespace, mappings):
    """Adapt path-bearing IO records, preserving the original input file bytes.

    Archived manifest keys and Run-map path fields must move along with data.
    Expected digest strings and scientific fields are not changed. This wrapper
    does not execute any top-level original module IO.
    """
    ordered = sorted(mappings.items(), key=lambda item: len(item[0]), reverse=True)
    already_local = {new for _, new in ordered}
    already_local.update(str(namespace[key]) for key in ["ROOT", "OUT"] if key in namespace)
    # Source provenance is already a real local file. A broad archived data
    # prefix (for example the author's project directory) may also contain the
    # checkout; do not redirect this one explicitly known code file into data.
    # This is an exact-file exception, never an exemption for the whole repo.
    code_files = {str(Path(namespace["__file__"]).resolve())} if "__file__" in namespace else set()

    def translate(value):
        value = str(value)
        if value in code_files:
            return value
        if any(value == root or value.startswith(root + "/") for root in already_local):
            return value
        for old, new in ordered:
            if value == old or value.startswith(old + "/"):
                return new + value[len(old):]
        return value

    if "read" in namespace:
        original_read = namespace["read"]

        def rebase_read(path):
            records = original_read(Path(translate(path)))
            return [{key: translate(value) if isinstance(value, str) else value
                     for key, value in row.items()} for row in records]

        namespace["read"] = rebase_read
    for name in ["reg", "register"]:
        if name in namespace:
            original = namespace[name]

            def wrapped(path, *args, _original=original, **kwargs):
                return _original(Path(translate(path)), *args, **kwargs)

            namespace[name] = wrapped
    return translate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["b", "c", "d"], required=True)
    parser.add_argument("--config", type=Path, required=True,
                        help="JSON: base, project_inputs; optional genus_root")
    parser.add_argument("--source", type=Path,
                        help="Optional override of the bundled SHA-pinned producer")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cfg = resolve_config(args.config.resolve())
    mappings = cfg.get("path_rebase", {})
    out = args.output.resolve()
    if out.exists():
        raise ValueError("Output already exists; select a new directory")
    base = cfg["base"]
    project = cfg["project_inputs"]
    b5 = base / "9_22_stage5B_same_data_cross_method"
    n5 = b5 / "04_network_comparison"
    follow = base / "10_13_group_meeting_followup"
    ext = base / "9_16_external_validation"
    s6 = follow / "02_stage6B_project_network_comparison"
    genus = cfg.get("genus_root", follow / "figure7_group_genus_check_20261029")
    source = (args.source or default_source(args.mode)).resolve()
    # Check the producer before any output creation.
    if digest(source) != PINS[args.mode]:
        raise ValueError("Source SHA256 does not match the accepted producer")
    ns = dict(__name__="frozen_statistics", __file__=str(source), ROOT=out, OUT=out,
              BASE=base, FOLLOW=follow, EXT=ext, B5=b5, N5=n5, S6=s6,
              PROJECT=project, NA="NA", HIST=defaultdict(list), QC=[], Q=[],
              MAN={}, INPUTS={}, DATA={})
    out.mkdir(parents=True)
    start = time.monotonic()
    state = dict(mode=args.mode, source_sha256=PINS[args.mode],
                 source=str(source), output=str(out), status="RUNNING",
                 inference_executed=False, null_simulations_executed=False,
                 function_ast_modified=False)
    try:
        for name in ["data", "qc", "logs", "00_inputs", "01_published_networks",
                     "02_fastspar_networks", "03_common_universe", "04_overlap",
                     "05_sign_concordance", "07_statistics", "genus"]:
            (out / name).mkdir()
        if args.mode == "b":
            tree = load_functions(source, PINS["b"], ns)
            install_io_rebase(ns, mappings)
            ns.update(STAGE=b5, PREP=b5 / "02_prepare_fastspar_input",
                      FSROOT=b5 / "03_fastspar_result", SPECROOT=project / "GDTB226_species_abundance_result",
                      GUOWB=ext / "raw/guo2023/1-s2.0-S193131282300029X-mmc8.xlsx",
                      WALLENWB=b5 / "00_metadata_cache/Wallen_supplementary_datasets.xlsx",
                      WMAP=base / "figure56/PMID_36376318_taxonomy_rescue_overlap_v44/02_taxon_mapping_audit.tsv",
                      GMAP=ext / "stage2b/03_external_taxa_to_gtdb_r226_stage2b.tsv")
            ns["NETWORKS"] = literal_assignment(tree, "NETWORKS", ns)
            for name in ["gtdb_r226_species_vocabulary.tsv", "reference_manifest.json"]:
                shutil.copyfile(n5 / "00_inputs" / name, out / "00_inputs" / name)
        else:
            helper = source.with_name("common.py")
            tree = load_functions(helper, PINS[args.mode + "_common"], ns)
            translate = install_io_rebase(ns, mappings)
            state["common_sha256"] = PINS[args.mode + "_common"]
            if args.mode == "c":
                ns["RES"] = b5 / "01_sample_rescue"
                ns["HFILES"] = literal_assignment(tree, "HFILES", ns)
                for path, key, hashkey in ns["HFILES"]:
                    for row in ns["read"](path):
                        ns["HIST"][row[key]].append((str(path), row[hashkey]))
                snap = ns["RES"] / "qc/input_snapshot.json"
                ns["sp"] = snap
                for path, row in json.loads(snap.read_text()).items():
                    if row.get("exists") and row.get("sha256"):
                        ns["HIST"][translate(path)].append((str(snap), row["sha256"]))
                pids = ["PRJNA751448", "PRJNA834801"]
                ns["PIDS"] = pids
                ns["ABUNDANCE"] = {p + "_abundance.matrix": project / "GDTB226_species_abundance_result/abundance" / (p + "_abundance.matrix") for p in pids}
                ns["FASTSPAR"] = {p + "_" + kind + ".tsv": project / "r226_fastspar_species_result" / (p + "_" + kind + ".tsv") for p in pids for kind in ["correlation", "pvalues", "covariance"]}
            else:
                for path, key, hashkey in [(n5 / "00_inputs/input_manifest.tsv", "path", "sha256_before"),
                                           (n5 / "qc/stage5B2_delivery_manifest.tsv", "path", "sha256"),
                                           (s6 / "stage6B_input_manifest.tsv", "input_path", "sha256_before")]:
                    for row in ns["read"](path):
                        ns["HIST"][row[key]].append(row[hashkey])
                dataset_source = cfg.get("comparison_datasets", genus / "data/comparison_datasets.json")
                datasets = json.loads(dataset_source.read_text())
                for dataset in datasets:
                    for key, value in dataset.items():
                        if key.endswith("_path"):
                            value = Path(translate(value))
                            dataset[key] = str(value if value.is_absolute() else (dataset_source.parent / value).resolve())
                (out / "data/comparison_datasets.json").write_text(json.dumps(datasets, indent=2) + "\n")
            load_functions(source, PINS[args.mode], ns)
            # Producer imports cannot replace the helper IO wrappers: its only
            # common import was deliberately excluded by load_functions.
        state["io_path_rebase_entries"] = len(mappings)
        ns["main"]()
        state["status"] = "PASS"
    except Exception as exc:
        state["status"] = "FAIL"
        state["error"] = repr(exc)
        (out / "ADAPTER_FAILURE.txt").write_text(traceback.format_exc())
        raise
    finally:
        state["wall_seconds"] = time.monotonic() - start
        (out / "ADAPTER_EXECUTION.json").write_text(json.dumps(state, indent=2) + "\n")


if __name__ == "__main__":
    main()
