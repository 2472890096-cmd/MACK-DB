#!/usr/bin/env python3

from pathlib import Path
import csv
import os
import sys

COR_CUTOFF = 0.3
PVALUE_CUTOFF = 0.0001

CONFIGS = [
    {
        "rank": "species",
        "base": None,
        "pattern": "*_fastCCLasso_freq_sp_pairs.tsv",
        "suffix": "_fastCCLasso_freq_sp_pairs.tsv",
        "csv_suffix": "_fastCCLasso_freq_sp_pairs.csv",
    },
    {
        "rank": "genus",
        "base": None,
        "pattern": "*_fastCCLasso_freq_genus_pairs.tsv",
        "suffix": "_fastCCLasso_freq_genus_pairs.tsv",
        "csv_suffix": "_fastCCLasso_freq_genus_pairs.csv",
    },
]

def atomic_replace(path, writer):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        writer(tmp)
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()

def update_status_file(path, n_edges):
    if not path.exists():
        return

    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        fieldnames = reader.fieldnames

    if not fieldnames or "selected_edges" not in fieldnames:
        return

    for row in rows:
        row["selected_edges"] = str(n_edges)

    def write_tmp(tmp):
        with tmp.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=fieldnames,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)

    atomic_replace(path, write_tmp)

def update_summary_file(path, edge_counts):
    if not path.exists():
        return

    try:
        with path.open("r", newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
            fieldnames = reader.fieldnames
    except Exception:
        return

    if (
        not fieldnames
        or "project" not in fieldnames
        or "selected_edges" not in fieldnames
    ):
        return

    changed = False

    for row in rows:
        project = row.get("project", "")
        if project in edge_counts:
            row["selected_edges"] = str(edge_counts[project])
            changed = True

    if not changed:
        return

    def write_tmp(tmp):
        with tmp.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=fieldnames,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)

    atomic_replace(path, write_tmp)

def process_config(cfg):
    rank = cfg["rank"]
    base = cfg["base"]

    tsv_dir = base / "cor-pv_tsvformat"
    csv_dir = base / "cor-pv_csvformat"
    links_dir = base / "proj_links"
    status_dir = base / "project_status"

    files = sorted(tsv_dir.glob(cfg["pattern"]))

    print("=" * 72)
    print(f"{rank}: found {len(files)} pair TSV files")
    print("=" * 72)

    edge_counts = {}
    total_old = 0
    total_new = 0

    report_rows = []

    for i, tsv_path in enumerate(files, 1):
        project = tsv_path.name[:-len(cfg["suffix"])]

        kept = []
        old_count = 0

        with tsv_path.open("r", newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")

            required = {"from", "to", "corr", "pvalue"}

            if not reader.fieldnames or not required.issubset(reader.fieldnames):
                raise RuntimeError(
                    f"Unexpected columns in {tsv_path}: {reader.fieldnames}"
                )

            for row in reader:
                old_count += 1

                try:
                    corr = float(row["corr"])
                    pvalue = float(row["pvalue"])
                except (ValueError, TypeError):
                    continue

                if (
                    abs(corr) >= COR_CUTOFF
                    and pvalue <= PVALUE_CUTOFF
                ):
                    kept.append({
                        "from": row["from"],
                        "to": row["to"],
                        "corr": row["corr"],
                        "pvalue": row["pvalue"],
                    })

        new_count = len(kept)

        total_old += old_count
        total_new += new_count
        edge_counts[project] = new_count

        # TSV
        def write_tsv(tmp):
            with tmp.open("w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(
                    f,
                    fieldnames=["from", "to", "corr", "pvalue"],
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(kept)

        atomic_replace(tsv_path, write_tsv)

        # CSV
        csv_path = csv_dir / f"{project}{cfg['csv_suffix']}"

        def write_csv(tmp):
            with tmp.open("w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(
                    f,
                    fieldnames=["from", "to", "corr", "pvalue"],
                    delimiter=",",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(kept)

        atomic_replace(csv_path, write_csv)

        # links.txt，与原R脚本相同：每条边两个taxon，无表头
        links_path = links_dir / f"{project}_links.txt"

        def write_links(tmp):
            with tmp.open("w", newline="", encoding="utf-8") as f:
                for row in kept:
                    a = row["from"]
                    b = row["to"]
                    x, y = sorted((a, b))
                    f.write(f"{x}\t{y}\n")

        atomic_replace(links_path, write_links)

        # 更新单project状态中的 selected_edges
        update_status_file(
            status_dir / f"{project}.tsv",
            new_count
        )

        report_rows.append({
            "project": project,
            "old_edges_cor0p1": old_count,
            "new_edges_cor0p3": new_count,
        })

        if i % 100 == 0 or i == len(files):
            print(
                f"[{rank}] {i}/{len(files)} "
                f"processed; current={project}; "
                f"{old_count} -> {new_count}"
            )

    # 更新group summary中的selected_edges
    group_dir = base / "group_summaries"

    if group_dir.exists():
        for path in group_dir.glob("project_timing_and_status_*.tsv"):
            update_summary_file(path, edge_counts)

    # 更新可能存在的总summary
    for path in [
        base / "project_timing_and_status_all.tsv",
        base / "project_timing_and_status.tsv",
    ]:
        update_summary_file(path, edge_counts)

    report_file = base / "refilter_cor0p3_summary.tsv"

    with report_file.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "project",
                "old_edges_cor0p1",
                "new_edges_cor0p3",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(report_rows)

    print()
    print(f"{rank} finished")
    print(f"old |cor|>=0.1 edges : {total_old}")
    print(f"new |cor|>=0.3 edges : {total_new}")
    print(f"report               : {report_file}")
    print()

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Apply the accepted |r|>=0.3 and P<=1e-4 filter to an explicit working copy; writes in place.")
    parser.add_argument("--level", choices=["species", "genus"], required=True)
    parser.add_argument("--work-root", type=Path, required=True)
    args = parser.parse_args()
    if not (args.work_root / "cor-pv_tsvformat").is_dir():
        parser.error("work-root must contain cor-pv_tsvformat")
    cfg = dict(next(c for c in CONFIGS if c['rank'] == args.level))
    cfg['base'] = args.work_root.resolve()
    process_config(cfg)

if __name__ == '__main__':
    main()
