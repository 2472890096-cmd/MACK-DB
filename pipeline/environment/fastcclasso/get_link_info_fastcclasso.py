#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
from pathlib import Path
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Generate 11-column environment-level "
            "fastCCLasso link info table."
        )
    )

    p.add_argument(
        "-i", "--input",
        required=True,
        help="Filtered .belonging file."
    )

    p.add_argument(
        "-d", "--pair-dir",
        required=True,
        help=(
            "Directory containing fastCCLasso "
            "PROJECT pairs TSV files."
        )
    )

    p.add_argument(
        "-a", "--abundance-dir",
        required=True,
        help=(
            "Directory containing "
            "PROJECT_abundance.matrix."
        )
    )

    p.add_argument(
        "-o", "--out",
        required=True,
        help="Output directory."
    )

    p.add_argument(
        "-n", "--name",
        required=True,
        help="Output filename."
    )

    p.add_argument(
        "-e", "--environment",
        required=True,
        help=(
            "Environment name, e.g. human_gut."
        )
    )

    p.add_argument(
        "--rank",
        choices=["species", "genus"],
        required=True
    )

    return p.parse_args()


def canonical_pair(a, b):
    return tuple(sorted((a.strip(), b.strip())))


def fmt(x):
    # 避免把1e-8这样的P值round成0.00000
    return f"{x:.10g}"


def get_pair_suffix(rank):
    if rank == "species":
        return "_fastCCLasso_freq_sp_pairs.tsv"

    return "_fastCCLasso_freq_genus_pairs.tsv"


def get_environment_parts(env):
    if "_" in env:
        super_env, sub_env = env.split("_", 1)
    else:
        super_env = env
        sub_env = "-"

    return super_env, sub_env


def main():
    args = parse_args()

    input_file = Path(args.input).resolve()
    pair_dir = Path(args.pair_dir).resolve()
    abundance_dir = Path(
        args.abundance_dir
    ).resolve()
    out_dir = Path(args.out).resolve()

    out_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    suffix = get_pair_suffix(args.rank)

    # -------------------------------------------------
    # 读取 belonging
    # -------------------------------------------------
    records = []

    needed_by_project = defaultdict(set)

    with input_file.open(
        "r",
        encoding="utf-8"
    ) as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        required_cols = {
            "species_a",
            "species_b",
            "project_count",
            "project_ID",
        }

        if (
            reader.fieldnames is None
            or not required_cols.issubset(
                reader.fieldnames
            )
        ):
            raise RuntimeError(
                f"Unexpected belonging columns: "
                f"{reader.fieldnames}"
            )

        for row in reader:
            a = row["species_a"]
            b = row["species_b"]

            if (
                a == "NA"
                and b == "NA"
            ):
                continue

            pair = canonical_pair(a, b)

            projects = [
                x
                for x in row["project_ID"].split(";")
                if x
            ]

            expected_count = int(
                row["project_count"]
            )

            if expected_count != len(projects):
                raise RuntimeError(
                    f"project_count mismatch for "
                    f"{a} -- {b}"
                )

            records.append(
                (pair, projects)
            )

            for project in projects:
                needed_by_project[
                    project
                ].add(pair)

    unique_projects = sorted(
        needed_by_project
    )

    print(
        f"Need {len(records)} links "
        f"from {len(unique_projects)} projects."
    )

    # -------------------------------------------------
    # 每个project只读取一次
    # 而且只缓存belonging真正需要的边
    # -------------------------------------------------
    pair_cache = defaultdict(dict)
    sample_cache = {}

    for idx, project in enumerate(
        unique_projects,
        1
    ):
        abundance_file = (
            abundance_dir /
            f"{project}_abundance.matrix"
        )

        if not abundance_file.exists():
            raise FileNotFoundError(
                f"Missing abundance file: "
                f"{abundance_file}"
            )

        with abundance_file.open(
            "r",
            encoding="utf-8"
        ) as f:
            header = f.readline().rstrip("\n")

        sample_count = (
            len(header.split("\t")) - 1
        )

        sample_cache[project] = str(
            sample_count
        )

        pair_file = (
            pair_dir /
            f"{project}{suffix}"
        )

        if not pair_file.exists():
            raise FileNotFoundError(
                f"Missing pair file: "
                f"{pair_file}"
            )

        needed = needed_by_project[
            project
        ]

        with pair_file.open(
            "r",
            encoding="utf-8"
        ) as f:

            reader = csv.DictReader(
                f,
                delimiter="\t"
            )

            required_pair_cols = {
                "from",
                "to",
                "corr",
                "pvalue",
            }

            if (
                reader.fieldnames is None
                or not required_pair_cols.issubset(
                    reader.fieldnames
                )
            ):
                raise RuntimeError(
                    f"Unexpected columns in "
                    f"{pair_file}: "
                    f"{reader.fieldnames}"
                )

            for row in reader:
                pair = canonical_pair(
                    row["from"],
                    row["to"]
                )

                if pair not in needed:
                    continue

                pair_cache[
                    project
                ][pair] = (
                    float(row["corr"]),
                    float(row["pvalue"]),
                )

        missing_edges = (
            needed -
            set(pair_cache[project])
        )

        if missing_edges:
            examples = list(
                missing_edges
            )[:5]

            raise RuntimeError(
                f"{project}: "
                f"{len(missing_edges)} links exist "
                f"in links.txt/belonging but cannot "
                f"be found in pair TSV. "
                f"Examples: {examples}"
            )

        if (
            idx % 20 == 0
            or idx == len(unique_projects)
        ):
            print(
                f"Cached {idx}/"
                f"{len(unique_projects)} projects"
            )

    # -------------------------------------------------
    # 输出11列info
    # -------------------------------------------------
    super_env, sub_env = (
        get_environment_parts(
            args.environment
        )
    )

    output_file = (
        out_dir /
        args.name
    )

    with output_file.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as out:

        writer = csv.writer(
            out,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writerow([
            "species_a",
            "species_b",
            "aver_corr",
            "aver_pv",
            "proj_count",
            "proj_ID",
            "proj_corr",
            "proj_pv",
            "proj_sam",
            "super_env",
            "sub_env",
        ])

        for pair, projects in records:
            cor_values = []
            pv_values = []
            sample_values = []

            for project in projects:
                corr, pv = (
                    pair_cache[
                        project
                    ][pair]
                )

                cor_values.append(corr)
                pv_values.append(pv)

                sample_values.append(
                    sample_cache[project]
                )

            aver_corr = (
                sum(cor_values) /
                len(cor_values)
            )

            aver_pv = (
                sum(pv_values) /
                len(pv_values)
            )

            writer.writerow([
                pair[0],
                pair[1],
                fmt(aver_corr),
                fmt(aver_pv),
                len(projects),
                ";".join(projects),
                ";".join(
                    fmt(x)
                    for x in cor_values
                ),
                ";".join(
                    fmt(x)
                    for x in pv_values
                ),
                ";".join(
                    sample_values
                ),
                super_env,
                sub_env,
            ])

    print(
        f"INFO written: {output_file}"
    )


if __name__ == "__main__":
    main()
