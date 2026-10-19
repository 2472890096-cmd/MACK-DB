#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Final Chen + Wallen Null A / Null B analysis
--------------------------------------------
One script runs:
  Chen   Null A x1000
  Chen   Null B x1000
  Wallen Null A x1000
  Wallen Null B x1000

Final observed networks:
  Chen   2733 unique GTDB pairs, 2045 recovered after source exclusion
  Wallen 2268 unique GTDB pairs, 1784 recovered after source exclusion

Null A:
  Same external species set and same edge count; sample edges uniformly
  without replacement from all possible undirected pairs.

Null B:
  Exact-degree-preserving double-edge swaps. Each replicate starts from the
  observed external network and performs 50E SUCCESSFUL rewiring operations,
  rejecting self-loops and duplicate edges.

Empirical one-sided P:
  (1 + number(null recovery >= observed recovery)) / (1000 + 1)

The script is restartable. Existing replicate TSVs are reused automatically.
"""

from __future__ import annotations

# v0.1: repository-local IO configuration; scientific functions below are unchanged.
import sys as _path_sys
from pathlib import Path as _RepoPath
_repo_root = next(p for p in _RepoPath(__file__).resolve().parents if (p / "lib/mackdb_paths.py").is_file())
_path_sys.path.insert(0, str(_repo_root / "lib"))
from mackdb_paths import resolve_path

import argparse
import math
import os
import random
import re
from dataclasses import dataclass
from itertools import combinations
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(str(resolve_path('external_root')))

MACKDB = Path(
    str(resolve_path('environment_network'))
)

CHEN_EXTERNAL = (
    ROOT / "chen_postrescue_source_exclusion_overlap_v45"
    / "01_final_external_unique_gtdb_pairs.tsv"
)

WALLEN_EXTERNAL = (
    ROOT / "9_2_codex"
    / "PMID_36376318_taxonomy_rescue_overlap_v44"
    / "07_unique_gtdb_pairs.tsv"
)

DEFAULT_OUTDIR = ROOT / "chen_wallen_null_models_final_v1"

EXPECTED = {
    "Chen": {
        "edges": 2733,
        "recovered": 2045,
        "source_ref_pairs": 653068,
        "exclude": {"PRJNA319574", "PRJNA942468"},
    },
    "Wallen": {
        "edges": 2268,
        "recovered": 1784,
        "source_ref_pairs": 653069,
        "exclude": {"PRJNA834801"},
    },
}

DEFAULT_REPS = 1000
DEFAULT_SWAPS_PER_EDGE = 50
DEFAULT_SEED = 20260910


def die(msg):
    raise RuntimeError(msg)


def clean(x):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return ""
    s = str(x).strip()
    return "" if s.casefold() in {"", "nan", "none", "null", "na"} else s


def norm_gtdb(x):
    s = clean(x)
    if not s:
        return ""
    if ";" in s and "s__" in s:
        m = re.search(r"(?:^|;)s__([^;]*)", s)
        if m:
            s = m.group(1)
    s = s.strip("\"'")
    s = re.sub(r"^s__", "", s, flags=re.I)
    s = re.sub(r"^\d+_", "", s)  # remove MACK-DB internal numeric prefix
    s = s.replace("_", " ")
    s = re.sub(r"\s+", " ", s).strip()
    return s.casefold()


def canonical_name_pair(a, b):
    a = norm_gtdb(a)
    b = norm_gtdb(b)
    if not a or not b or a == b:
        return None
    return (a, b) if a <= b else (b, a)


def project_ids(x):
    s = clean(x)
    if not s:
        return []
    hits = re.findall(r"\bPRJ(?:NA|EB|DB)\d+\b", s, flags=re.I)
    if hits:
        return list(dict.fromkeys(h.upper() for h in hits))
    vals = [v.strip().upper() for v in re.split(r"[;,|]+", s) if v.strip()]
    return list(dict.fromkeys(vals))


def require_file(path, label):
    if not Path(path).is_file():
        die(f"{label} not found: {path}")


def encode_edge(u, v, n):
    if u > v:
        u, v = v, u
    return u * n + v


def load_external(path, study):
    require_file(path, f"{study} external network")
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    if not {"gtdb_a", "gtdb_b"}.issubset(df.columns):
        die(f"{study} external file must contain gtdb_a and gtdb_b.")
    pairs = []
    for _, r in df.iterrows():
        p = canonical_name_pair(r["gtdb_a"], r["gtdb_b"])
        if p:
            pairs.append(p)
    return sorted(set(pairs))


def load_mackdb():
    require_file(MACKDB, "MACK-DB human_gut.info")
    df = pd.read_csv(
        MACKDB, sep="\t", dtype=str, low_memory=False,
        usecols=["species_a", "species_b", "proj_ID"]
    )
    df["_pair"] = [
        canonical_name_pair(a, b)
        for a, b in zip(df["species_a"], df["species_b"])
    ]
    return df


def build_source_excluded_reference(db, excluded_projects):
    retained = set()
    affected = 0
    removed = 0
    per_project = {p: 0 for p in sorted(excluded_projects)}

    for _, r in db.iterrows():
        p = r["_pair"]
        if p is None:
            continue

        projects = project_ids(r["proj_ID"])
        present = [x for x in projects if x in excluded_projects]

        if present:
            affected += 1
            for x in present:
                per_project[x] += 1

        remaining = [x for x in projects if x not in excluded_projects]
        if len(set(remaining)) >= 2:
            retained.add(p)
        else:
            removed += 1

    return retained, {
        "affected_rows": affected,
        "removed_rows": removed,
        "retained_pairs": len(retained),
        "per_project": per_project,
    }


@dataclass
class Study:
    name: str
    nodes: list
    edges: list
    observed_keys: set
    ref_keys: set
    possible_keys: list
    n: int
    E: int
    observed_recovered: int
    observed_rate: float


def prepare_study(name, name_pairs, source_ref):
    nodes = sorted({x for p in name_pairs for x in p})
    idx = {x: i for i, x in enumerate(nodes)}
    n = len(nodes)

    edges = []
    obs = set()
    for a, b in name_pairs:
        u, v = idx[a], idx[b]
        if u > v:
            u, v = v, u
        edges.append((u, v))
        obs.add(encode_edge(u, v, n))

    if len(edges) != len(obs):
        die(f"{name}: duplicate observed edges after encoding.")

    ref_keys = set()
    for a, b in source_ref:
        if a in idx and b in idx:
            ref_keys.add(encode_edge(idx[a], idx[b], n))

    observed_recovered = sum(k in ref_keys for k in obs)
    possible = [encode_edge(u, v, n) for u, v in combinations(range(n), 2)]

    return Study(
        name=name,
        nodes=nodes,
        edges=edges,
        observed_keys=obs,
        ref_keys=ref_keys,
        possible_keys=possible,
        n=n,
        E=len(edges),
        observed_recovered=observed_recovered,
        observed_rate=observed_recovered / len(edges),
    )


def load_done(path, n_reps):
    if not path.is_file():
        return pd.DataFrame()
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        return pd.DataFrame()
    if "replicate" not in df.columns:
        return pd.DataFrame()
    df = df.drop_duplicates("replicate", keep="last")
    df = df[(df["replicate"] >= 1) & (df["replicate"] <= n_reps)]
    return df.sort_values("replicate")


def append_row(path, row):
    exists = path.is_file() and path.stat().st_size > 0
    pd.DataFrame([row]).to_csv(
        path, sep="\t", index=False, mode="a", header=not exists
    )


def run_null_a(study, n_reps, seed_base, outfile):
    old = load_done(outfile, n_reps)
    done = set(old["replicate"].astype(int)) if len(old) else set()
    print(f"[{study.name} Null A] resume {len(done)}/{n_reps}", flush=True)

    for rep in range(1, n_reps + 1):
        if rep in done:
            continue
        rng = random.Random(seed_base + rep)
        sampled = rng.sample(study.possible_keys, study.E)
        recovered = sum(k in study.ref_keys for k in sampled)
        append_row(outfile, {
            "replicate": rep,
            "seed": seed_base + rep,
            "recovered_count": recovered,
            "recovery_rate": recovered / study.E,
        })
        if rep % 100 == 0 or rep == n_reps:
            print(f"[{study.name} Null A] {rep}/{n_reps}", flush=True)

    final = load_done(outfile, n_reps)
    if len(final) != n_reps:
        die(f"{study.name} Null A incomplete: {len(final)}/{n_reps}")
    return final


# multiprocessing globals for Null B
_B_U = _B_V = _B_BASE_KEYS = _B_REF = None
_B_N = _B_E = _B_TARGET = _B_SEED = None


def init_b_worker(edges, ref_keys, n, target_swaps, seed_base):
    global _B_U, _B_V, _B_BASE_KEYS, _B_REF
    global _B_N, _B_E, _B_TARGET, _B_SEED

    _B_U = [u for u, _ in edges]
    _B_V = [v for _, v in edges]
    _B_N = n
    _B_E = len(edges)
    _B_BASE_KEYS = {encode_edge(u, v, n) for u, v in edges}
    _B_REF = set(ref_keys)
    _B_TARGET = target_swaps
    _B_SEED = seed_base


def null_b_worker(rep):
    rng = random.Random(_B_SEED + rep)

    U = _B_U.copy()
    V = _B_V.copy()
    edge_set = _B_BASE_KEYS.copy()

    accepted = 0
    attempts = 0
    max_attempts = max(_B_TARGET * 1000, _B_TARGET + 1_000_000)

    while accepted < _B_TARGET:
        attempts += 1
        if attempts > max_attempts:
            raise RuntimeError(
                f"replicate {rep}: exceeded max attempts after "
                f"{accepted}/{_B_TARGET} accepted swaps"
            )

        i = rng.randrange(_B_E)
        j = rng.randrange(_B_E - 1)
        if j >= i:
            j += 1

        a, b = U[i], V[i]
        c, d = U[j], V[j]

        # require four distinct nodes
        if a == c or a == d or b == c or b == d:
            continue

        if rng.getrandbits(1):
            x1, y1 = a, d
            x2, y2 = c, b
        else:
            x1, y1 = a, c
            x2, y2 = b, d

        if x1 > y1:
            x1, y1 = y1, x1
        if x2 > y2:
            x2, y2 = y2, x2

        if x1 == y1 or x2 == y2:
            continue

        new1 = x1 * _B_N + y1
        new2 = x2 * _B_N + y2

        if new1 == new2 or new1 in edge_set or new2 in edge_set:
            continue

        old1 = a * _B_N + b
        old2 = c * _B_N + d

        edge_set.remove(old1)
        edge_set.remove(old2)
        edge_set.add(new1)
        edge_set.add(new2)

        U[i], V[i] = x1, y1
        U[j], V[j] = x2, y2

        accepted += 1

    recovered = sum(k in _B_REF for k in edge_set)

    return {
        "replicate": rep,
        "seed": _B_SEED + rep,
        "recovered_count": recovered,
        "recovery_rate": recovered / _B_E,
        "accepted_swaps": accepted,
        "attempts": attempts,
        "swap_acceptance_rate": accepted / attempts,
    }


def run_null_b(study, n_reps, swaps_per_edge, seed_base, workers, outfile):
    old = load_done(outfile, n_reps)
    done = set(old["replicate"].astype(int)) if len(old) else set()
    missing = [r for r in range(1, n_reps + 1) if r not in done]

    target = swaps_per_edge * study.E
    print(
        f"[{study.name} Null B] resume {len(done)}/{n_reps}; "
        f"E={study.E}; accepted swaps/rep={target:,}; workers={workers}",
        flush=True
    )

    if missing:
        with Pool(
            processes=workers,
            initializer=init_b_worker,
            initargs=(study.edges, study.ref_keys, study.n, target, seed_base)
        ) as pool:
            completed = len(done)
            for row in pool.imap_unordered(null_b_worker, missing, chunksize=1):
                append_row(outfile, row)
                completed += 1
                if completed % 25 == 0 or completed == n_reps:
                    print(
                        f"[{study.name} Null B] {completed}/{n_reps}",
                        flush=True
                    )

    final = load_done(outfile, n_reps)
    if len(final) != n_reps:
        die(f"{study.name} Null B incomplete: {len(final)}/{n_reps}")
    return final


def summarize(study, null_name, reps, n_reps):
    counts = pd.to_numeric(reps["recovered_count"], errors="raise").to_numpy(float)
    rates = counts / study.E

    mean_count = float(np.mean(counts))
    sd_count = float(np.std(counts, ddof=1))
    mean_rate = float(np.mean(rates))
    sd_rate = float(np.std(rates, ddof=1))

    q025_count = float(np.quantile(counts, 0.025))
    q975_count = float(np.quantile(counts, 0.975))
    q025_rate = float(np.quantile(rates, 0.025))
    q975_rate = float(np.quantile(rates, 0.975))

    enrichment = study.observed_rate / mean_rate if mean_rate > 0 else float("inf")
    z = (
        (study.observed_recovered - mean_count) / sd_count
        if sd_count > 0 else float("inf")
    )

    n_ge = int(np.sum(counts >= study.observed_recovered))
    p_emp = (1 + n_ge) / (n_reps + 1)

    return {
        "study": study.name,
        "null_model": null_name,
        "external_species_count": len(study.nodes),
        "observed_edge_count": study.E,
        "observed_recovered_count": study.observed_recovered,
        "observed_recovery_rate": study.observed_rate,
        "null_replicates": n_reps,
        "null_mean_recovered_count": mean_count,
        "null_sd_recovered_count": sd_count,
        "null_2.5pct_recovered_count": q025_count,
        "null_97.5pct_recovered_count": q975_count,
        "null_mean_recovery_rate": mean_rate,
        "null_sd_recovery_rate": sd_rate,
        "null_2.5pct_recovery_rate": q025_rate,
        "null_97.5pct_recovery_rate": q975_rate,
        "observed_to_null_mean_enrichment": enrichment,
        "z_score": z,
        "null_replicates_ge_observed": n_ge,
        "empirical_one_sided_P": p_emp,
    }


def pct(x):
    return f"{100*x:.2f}%"


def write_report(outdir, chen, wallen, source_stats, summary, reps, swaps_per_edge):
    lines = []
    lines.append("FINAL CHEN + WALLEN NULL MODEL REPORT")
    lines.append("=" * 90)
    lines.append("")
    lines.append("Observed source-excluded networks")
    lines.append(
        f"Chen: {chen.observed_recovered}/{chen.E} = {pct(chen.observed_rate)}"
    )
    lines.append(
        f"Wallen: {wallen.observed_recovered}/{wallen.E} = {pct(wallen.observed_rate)}"
    )
    lines.append("")
    lines.append("Source-excluded MACK-DB references")
    for study in ["Chen", "Wallen"]:
        s = source_stats[study]
        lines.append(
            f"{study}: retained_pairs={s['retained_pairs']}; "
            f"affected_rows={s['affected_rows']}; "
            f"removed_rows={s['removed_rows']}"
        )
        for p, n in s["per_project"].items():
            lines.append(f"  {p}: affected_rows={n}")

    lines.append("")
    lines.append("Randomization settings")
    lines.append(f"Replicates per null model: {reps}")
    lines.append(
        "Null A: same external species set and edge count; uniform sampling "
        "without replacement from all possible undirected pairs."
    )
    lines.append(
        f"Null B: exact-degree-preserving rewiring; {swaps_per_edge}E successful "
        "double-edge swaps per replicate; every replicate restarts from observed."
    )
    lines.append(
        f"Chen Null B swaps/replicate: {swaps_per_edge * chen.E:,}"
    )
    lines.append(
        f"Wallen Null B swaps/replicate: {swaps_per_edge * wallen.E:,}"
    )

    lines.append("")
    lines.append("Results")
    for _, r in summary.iterrows():
        lines.append("")
        lines.append(f"{r['study']} {r['null_model']}")
        lines.append(
            f"Observed: {int(r['observed_recovered_count'])}/"
            f"{int(r['observed_edge_count'])} = "
            f"{pct(r['observed_recovery_rate'])}"
        )
        lines.append(
            f"Null mean: {r['null_mean_recovered_count']:.3f} pairs = "
            f"{pct(r['null_mean_recovery_rate'])}"
        )
        lines.append(
            f"Null SD: {r['null_sd_recovered_count']:.3f} pairs = "
            f"{pct(r['null_sd_recovery_rate'])}"
        )
        lines.append(
            f"95% empirical interval: "
            f"{r['null_2.5pct_recovered_count']:.3f}-"
            f"{r['null_97.5pct_recovered_count']:.3f} pairs; "
            f"{pct(r['null_2.5pct_recovery_rate'])}-"
            f"{pct(r['null_97.5pct_recovery_rate'])}"
        )
        lines.append(
            f"Observed/null enrichment: "
            f"{r['observed_to_null_mean_enrichment']:.6f}"
        )
        lines.append(f"z-score: {r['z_score']:.6f}")
        lines.append(
            f"Null replicates >= observed: "
            f"{int(r['null_replicates_ge_observed'])}/{reps}"
        )
        lines.append(
            f"Empirical one-sided P: {r['empirical_one_sided_P']:.6f}"
        )

    lines.append("")
    lines.append("Interpretation")
    lines.append(
        "Null A controls for the external species set and total edge count. "
        "Null B additionally preserves the exact degree of every external species. "
        "The new results replace all old null-model values based on the former "
        "Chen 1,808-edge and Wallen 1,568-edge mapped networks."
    )

    (outdir / "NULL_MODEL_FINAL_REPORT.txt").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main():
    default_workers = min(16, max(1, (os.cpu_count() or 2) - 1))

    ap = argparse.ArgumentParser()
    ap.add_argument("--replicates", type=int, default=DEFAULT_REPS)
    ap.add_argument("--swaps-per-edge", type=int, default=DEFAULT_SWAPS_PER_EDGE)
    ap.add_argument("--workers", type=int, default=default_workers)
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    ap.add_argument("--outdir", default=str(DEFAULT_OUTDIR))
    ap.add_argument("--allow-baseline-mismatch", action="store_true")
    args = ap.parse_args()

    if args.replicates <= 0 or args.swaps_per_edge <= 0 or args.workers <= 0:
        die("replicates, swaps-per-edge and workers must all be > 0")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 100, flush=True)
    print("FINAL CHEN + WALLEN NULL A / NULL B ANALYSIS", flush=True)
    print("=" * 100, flush=True)
    print(f"replicates/model={args.replicates}", flush=True)
    print(f"swaps-per-edge={args.swaps_per_edge}", flush=True)
    print(f"workers={args.workers}", flush=True)
    print(f"seed={args.seed}", flush=True)
    print(f"outdir={outdir}", flush=True)

    chen_pairs = load_external(CHEN_EXTERNAL, "Chen")
    wallen_pairs = load_external(WALLEN_EXTERNAL, "Wallen")

    print(f"Chen external unique GTDB pairs={len(chen_pairs)}", flush=True)
    print(f"Wallen external unique GTDB pairs={len(wallen_pairs)}", flush=True)

    db = load_mackdb()
    print(f"MACK-DB rows loaded={len(db):,}", flush=True)

    source_refs = {}
    source_stats = {}
    for study in ["Chen", "Wallen"]:
        ref, stats = build_source_excluded_reference(
            db, EXPECTED[study]["exclude"]
        )
        source_refs[study] = ref
        source_stats[study] = stats
        print(
            f"{study} source-excluded MACK-DB={len(ref):,}; "
            f"affected={stats['affected_rows']:,}; "
            f"removed={stats['removed_rows']:,}",
            flush=True
        )

    chen = prepare_study("Chen", chen_pairs, source_refs["Chen"])
    wallen = prepare_study("Wallen", wallen_pairs, source_refs["Wallen"])

    mismatches = []
    for s in [chen, wallen]:
        exp = EXPECTED[s.name]
        if s.E != exp["edges"]:
            mismatches.append(
                f"{s.name} edges {s.E} != expected {exp['edges']}"
            )
        if s.observed_recovered != exp["recovered"]:
            mismatches.append(
                f"{s.name} observed recovered {s.observed_recovered} "
                f"!= expected {exp['recovered']}"
            )
        if source_stats[s.name]["retained_pairs"] != exp["source_ref_pairs"]:
            mismatches.append(
                f"{s.name} source ref {source_stats[s.name]['retained_pairs']} "
                f"!= expected {exp['source_ref_pairs']}"
            )

    if mismatches and not args.allow_baseline_mismatch:
        die(
            "FINAL BASELINE SANITY CHECK FAILED:\n- "
            + "\n- ".join(mismatches)
        )

    print(
        f"Chen observed: species={len(chen.nodes)}, edges={chen.E}, "
        f"recovered={chen.observed_recovered}, rate={chen.observed_rate:.4%}",
        flush=True
    )
    print(
        f"Wallen observed: species={len(wallen.nodes)}, edges={wallen.E}, "
        f"recovered={wallen.observed_recovered}, rate={wallen.observed_rate:.4%}",
        flush=True
    )

    baseline = []
    for s in [chen, wallen]:
        st = source_stats[s.name]
        baseline.append({
            "study": s.name,
            "external_species_count": len(s.nodes),
            "external_edge_count": s.E,
            "observed_recovered_count": s.observed_recovered,
            "observed_recovery_rate": s.observed_rate,
            "source_excluded_MACKDB_pairs": st["retained_pairs"],
            "MACKDB_rows_affected": st["affected_rows"],
            "MACKDB_rows_removed_after_remaining_BioProjects_ge2":
                st["removed_rows"],
            "excluded_projects": ";".join(sorted(EXPECTED[s.name]["exclude"])),
        })
    pd.DataFrame(baseline).to_csv(
        outdir / "00_baseline_summary.tsv", sep="\t", index=False
    )

    seeds = {
        ("Chen", "A"): args.seed + 100000,
        ("Chen", "B"): args.seed + 200000,
        ("Wallen", "A"): args.seed + 300000,
        ("Wallen", "B"): args.seed + 400000,
    }

    chen_a = run_null_a(
        chen, args.replicates, seeds[("Chen", "A")],
        outdir / "01_Chen_NullA_replicates.tsv"
    )
    chen_b = run_null_b(
        chen, args.replicates, args.swaps_per_edge,
        seeds[("Chen", "B")], args.workers,
        outdir / "02_Chen_NullB_replicates.tsv"
    )
    wallen_a = run_null_a(
        wallen, args.replicates, seeds[("Wallen", "A")],
        outdir / "03_Wallen_NullA_replicates.tsv"
    )
    wallen_b = run_null_b(
        wallen, args.replicates, args.swaps_per_edge,
        seeds[("Wallen", "B")], args.workers,
        outdir / "04_Wallen_NullB_replicates.tsv"
    )

    summary = pd.DataFrame([
        summarize(chen, "Null A", chen_a, args.replicates),
        summarize(chen, "Null B", chen_b, args.replicates),
        summarize(wallen, "Null A", wallen_a, args.replicates),
        summarize(wallen, "Null B", wallen_b, args.replicates),
    ])
    summary.to_csv(
        outdir / "05_null_model_summary.tsv", sep="\t", index=False
    )

    manuscript = []
    for s in [chen, wallen]:
        a = summary[
            (summary["study"] == s.name) &
            (summary["null_model"] == "Null A")
        ].iloc[0]
        b = summary[
            (summary["study"] == s.name) &
            (summary["null_model"] == "Null B")
        ].iloc[0]
        manuscript.append({
            "study": s.name,
            "external_unique_GTDB_pairs": s.E,
            "observed_recovered": s.observed_recovered,
            "observed_recovery_percent": 100 * s.observed_rate,
            "Null_A_mean_recovery_percent": 100 * a["null_mean_recovery_rate"],
            "Null_B_mean_recovery_percent": 100 * b["null_mean_recovery_rate"],
            "empirical_P_Null_A": a["empirical_one_sided_P"],
            "empirical_P_Null_B": b["empirical_one_sided_P"],
            "Null_A_enrichment": a["observed_to_null_mean_enrichment"],
            "Null_B_enrichment": b["observed_to_null_mean_enrichment"],
            "Null_A_z": a["z_score"],
            "Null_B_z": b["z_score"],
        })
    pd.DataFrame(manuscript).to_csv(
        outdir / "06_manuscript_table_values.tsv", sep="\t", index=False
    )

    source_rows = []
    for study in ["Chen", "Wallen"]:
        st = source_stats[study]
        row = {
            "study": study,
            "excluded_projects": ";".join(sorted(EXPECTED[study]["exclude"])),
            "affected_rows": st["affected_rows"],
            "removed_rows_after_remaining_BioProjects_ge2": st["removed_rows"],
            "source_excluded_MACKDB_pairs": st["retained_pairs"],
        }
        for p, n in st["per_project"].items():
            row[f"affected_rows_{p}"] = n
        source_rows.append(row)
    pd.DataFrame(source_rows).to_csv(
        outdir / "07_source_exclusion_summary.tsv", sep="\t", index=False
    )

    write_report(
        outdir, chen, wallen, source_stats, summary,
        args.replicates, args.swaps_per_edge
    )

    print("", flush=True)
    print("=" * 100, flush=True)
    print("FINAL NULL MODEL RESULTS", flush=True)
    print("=" * 100, flush=True)
    for _, r in summary.iterrows():
        print(
            f"{r['study']} {r['null_model']}: "
            f"observed={int(r['observed_recovered_count'])}/"
            f"{int(r['observed_edge_count'])} "
            f"({100*r['observed_recovery_rate']:.2f}%); "
            f"null_mean={100*r['null_mean_recovery_rate']:.2f}%; "
            f"95%={100*r['null_2.5pct_recovery_rate']:.2f}-"
            f"{100*r['null_97.5pct_recovery_rate']:.2f}%; "
            f"enrichment={r['observed_to_null_mean_enrichment']:.4f}; "
            f"z={r['z_score']:.4f}; "
            f"P={r['empirical_one_sided_P']:.6f}",
            flush=True
        )

    print("", flush=True)
    print("Important outputs:", flush=True)
    print(outdir / "05_null_model_summary.tsv", flush=True)
    print(outdir / "06_manuscript_table_values.tsv", flush=True)
    print(outdir / "NULL_MODEL_FINAL_REPORT.txt", flush=True)


if __name__ == "__main__":
    main()
