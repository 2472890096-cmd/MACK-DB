#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import math
from pathlib import Path
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(
        description="Aggregate fastCCLasso links across projects in one environment."
    )

    p.add_argument(
        "-l", "--list",
        required=True,
        help="Environment project list, one project ID per line."
    )

    p.add_argument(
        "-d", "--dir",
        required=True,
        help="Directory containing PROJECT_links.txt files."
    )

    p.add_argument(
        "-o", "--out",
        required=True,
        help="Output directory."
    )

    p.add_argument(
        "-n", "--name",
        required=True,
        help="Output basename."
    )

    p.add_argument(
        "-f", "--filter",
        type=float,
        default=2,
        help=(
            "<=1: minimum fraction of analyzable projects; "
            ">1: minimum absolute project count."
        )
    )

    p.add_argument(
        "-s", "--switch",
        default="T",
        help="T = output filtered belonging file."
    )

    p.add_argument(
        "--denominator",
        choices=["available", "listed"],
        default="available",
        help=(
            "available: successfully analyzed projects with links file; "
            "listed: all projects in environment list."
        )
    )

    return p.parse_args()


def unique_preserve_order(values):
    seen = set()
    out = []

    for x in values:
        if x and x not in seen:
            seen.add(x)
            out.append(x)

    return out


def canonical_pair(a, b):
    return tuple(sorted((a.strip(), b.strip())))


def write_belonging(path, data):
    with path.open("w", encoding="utf-8") as out:
        out.write(
            "species_a\tspecies_b\tproject_count\tproject_ID\n"
        )

        for (a, b), projects in sorted(data.items()):
            out.write(
                f"{a}\t{b}\t{len(projects)}\t"
                + ";".join(projects)
                + "\n"
            )


def main():
    args = parse_args()

    list_file = Path(args.list).resolve()
    links_dir = Path(args.dir).resolve()
    out_dir = Path(args.out).resolve()

    out_dir.mkdir(parents=True, exist_ok=True)

    projects = unique_preserve_order([
        x.strip()
        for x in list_file.read_text(
            encoding="utf-8"
        ).splitlines()
    ])

    link_projects = defaultdict(list)

    available_projects = []
    missing_projects = []

    for project in projects:
        link_file = links_dir / f"{project}_links.txt"

        if not link_file.exists():
            missing_projects.append(project)
            continue

        # 即使links文件为空，也代表该project成功分析，
        # 因此计入环境分母。
        available_projects.append(project)

        seen_in_project = set()

        with link_file.open("r", encoding="utf-8") as f:
            for line_no, line in enumerate(f, 1):
                line = line.rstrip("\n")

                if not line:
                    continue

                parts = line.split("\t")

                if len(parts) < 2:
                    raise RuntimeError(
                        f"Malformed links file: "
                        f"{link_file}:{line_no}"
                    )

                pair = canonical_pair(parts[0], parts[1])

                # 同一个project中同一条边只计1次
                if pair in seen_in_project:
                    continue

                seen_in_project.add(pair)
                link_projects[pair].append(project)

    if args.denominator == "available":
        denominator = len(available_projects)
    else:
        denominator = len(projects)

    if denominator == 0:
        raise RuntimeError(
            f"No analyzable projects for {list_file}"
        )

    raw_file = out_dir / f"{args.name}.belonging"
    write_belonging(raw_file, link_projects)

    if args.filter <= 1:
        required = math.ceil(
            denominator * args.filter - 1e-12
        )
        required = max(required, 1)
        criterion = (
            f">={args.filter:.6g} of projects "
            f"(required={required}/{denominator})"
        )
    else:
        required = math.ceil(args.filter)
        criterion = (
            f">={required} projects"
        )

    filtered = {
        pair: projs
        for pair, projs in link_projects.items()
        if len(projs) >= required
    }

    filtered_file = (
        out_dir /
        f"{args.name}_filtered.belonging"
    )

    if args.switch.upper() == "T":
        write_belonging(
            filtered_file,
            filtered
        )

    missing_file = (
        out_dir /
        f"{args.name}_missing_projects.list"
    )

    missing_file.write_text(
        "\n".join(missing_projects)
        + ("\n" if missing_projects else ""),
        encoding="utf-8"
    )

    stats_file = (
        out_dir /
        f"{args.name}_stats.tsv"
    )

    with stats_file.open(
        "w",
        encoding="utf-8"
    ) as out:
        out.write(
            "listed_projects\t"
            "available_projects\t"
            "missing_projects\t"
            "denominator\t"
            "required_projects\t"
            "unfiltered_links\t"
            "filtered_links\n"
        )

        out.write(
            f"{len(projects)}\t"
            f"{len(available_projects)}\t"
            f"{len(missing_projects)}\t"
            f"{denominator}\t"
            f"{required}\t"
            f"{len(link_projects)}\t"
            f"{len(filtered)}\n"
        )

    print("=" * 70)
    print("Environment list :", list_file)
    print("Listed projects  :", len(projects))
    print("Available        :", len(available_projects))
    print("Missing          :", len(missing_projects))
    print("Denominator      :", denominator)
    print("Filter           :", criterion)
    print("All links        :", len(link_projects))
    print("Filtered links   :", len(filtered))
    print("Output           :", filtered_file)
    print("=" * 70)


if __name__ == "__main__":
    main()
