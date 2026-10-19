#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extract Zi-Pi candidate keystone taxa from species/genus node tables
and organize project_key.txt files by environment.

Accepted taxa_roles values, case-insensitive:
- Module hub / Module hubs
- Connector / Connectors
- Network hub / Network hubs

The output keeps the original header and all original columns.
A project listed in multiple environments is written to every matching
environment directory and recorded in a duplicate mapping report.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


DEFAULT_ENV_DIR = Path("data/env")
DEFAULT_SPECIES_NODE_DIR = Path(
    "work/zipi/"
    "two_pairs_test/species_all_signed/node_tables"
)
DEFAULT_SPECIES_OUT_DIR = Path(
    "work/zipi/"
    "two_pairs_test/species_all_signed/keyspecies"
)
DEFAULT_GENUS_NODE_DIR = Path(
    "work/zipi/"
    "two_pairs_test/genus_all_signed/node_tables"
)
DEFAULT_GENUS_OUT_DIR = Path(
    "work/zipi/"
    "two_pairs_test/genus_all_signed/keygenus"
)

ACCEPTED_ROLES = {
    "module hub",
    "module hubs",
    "connector",
    "connectors",
    "network hub",
    "network hubs",
}


@dataclass
class LevelSummary:
    level: str
    node_files: int = 0
    projects_mapped: int = 0
    projects_unmapped: int = 0
    projects_multiple_envs: int = 0
    key_files_written: int = 0
    projects_with_keys: int = 0
    projects_without_keys: int = 0
    key_rows: int = 0
    missing_role_column: int = 0
    malformed_files: int = 0


def normalize_project_id(value: str) -> str:
    """Strip whitespace and a possible UTF-8 BOM."""
    return value.strip().lstrip("\ufeff")


def normalize_role(value: str) -> str:
    """Normalize spaces, case, underscores, and hyphens in taxa role labels."""
    text = value.strip().lower().replace("_", " ").replace("-", " ")
    return " ".join(text.split())


def read_environment_mapping(
    env_dir: Path,
) -> tuple[list[str], dict[str, list[str]], list[tuple[str, str]]]:
    """
    Return:
      environment names,
      project -> environments mapping,
      duplicate (project, environments joined) records.
    """
    if not env_dir.is_dir():
        raise FileNotFoundError(f"环境目录不存在：{env_dir}")

    list_files = sorted(env_dir.glob("*.list"))
    if not list_files:
        raise FileNotFoundError(f"环境目录中没有找到 .list 文件：{env_dir}")

    project_to_envs: dict[str, list[str]] = defaultdict(list)
    env_names: list[str] = []

    for list_file in list_files:
        env_name = list_file.stem
        env_names.append(env_name)

        with list_file.open("r", encoding="utf-8-sig") as handle:
            for line in handle:
                project = normalize_project_id(line)
                if not project or project.startswith("#"):
                    continue
                if env_name not in project_to_envs[project]:
                    project_to_envs[project].append(env_name)

    duplicate_records = [
        (project, ",".join(envs))
        for project, envs in sorted(project_to_envs.items())
        if len(envs) > 1
    ]

    return env_names, dict(project_to_envs), duplicate_records


def prepare_environment_directories(
    out_dir: Path,
    env_names: Iterable[str],
    clean: bool,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    if clean:
        # Only remove environment subdirectories and prior report TSVs.
        for env_name in env_names:
            env_path = out_dir / env_name
            if env_path.exists():
                shutil.rmtree(env_path)
        for report_file in out_dir.glob("_*.tsv"):
            report_file.unlink()

    for env_name in env_names:
        (out_dir / env_name).mkdir(parents=True, exist_ok=True)


def detect_role_column(fieldnames: list[str] | None) -> str | None:
    if not fieldnames:
        return None

    normalized = {
        "".join(name.strip().lower().split()).replace("_", ""): name
        for name in fieldnames
    }

    # Current node files use taxa_roles. Keep a few compatible variants.
    for candidate in ("taxaroles", "taxarole", "roles", "role"):
        if candidate in normalized:
            return normalized[candidate]
    return None


def extract_key_rows(node_file: Path) -> tuple[list[str], list[dict[str, str]], str]:
    """
    Return fieldnames, filtered rows, and role column name.
    Raise ValueError for malformed files.
    """
    with node_file.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("缺少表头")

        role_column = detect_role_column(reader.fieldnames)
        if role_column is None:
            raise KeyError("未找到 taxa_roles 列")

        rows: list[dict[str, str]] = []
        for row in reader:
            role = normalize_role(row.get(role_column, "") or "")
            if role in ACCEPTED_ROLES:
                rows.append(row)

        return reader.fieldnames, rows, role_column


def write_key_file(
    output_file: Path,
    fieldnames: list[str],
    rows: list[dict[str, str]],
) -> None:
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_tsv(path: Path, header: list[str], rows: Iterable[Iterable[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            writer.writerow(list(row))


def process_level(
    level: str,
    node_dir: Path,
    out_dir: Path,
    env_names: list[str],
    project_to_envs: dict[str, list[str]],
    duplicate_records: list[tuple[str, str]],
    write_empty: bool,
    clean: bool,
) -> LevelSummary:
    summary = LevelSummary(level=level)
    prepare_environment_directories(out_dir, env_names, clean=clean)

    node_files = sorted(node_dir.glob("*_node.txt"))
    if not node_dir.is_dir():
        raise FileNotFoundError(f"{level} node目录不存在：{node_dir}")
    if not node_files:
        raise FileNotFoundError(f"{level} node目录中没有 *_node.txt：{node_dir}")

    unmatched_rows: list[tuple[str, str]] = []
    multiple_env_rows: list[tuple[str, str]] = []
    project_rows: list[tuple[str, str, int, int, str]] = []
    role_count_rows: list[tuple[str, str, str, int]] = []

    for node_file in node_files:
        summary.node_files += 1
        project = node_file.name.removesuffix("_node.txt")
        envs = project_to_envs.get(project, [])

        if not envs:
            summary.projects_unmapped += 1
            unmatched_rows.append((project, str(node_file)))
            continue

        summary.projects_mapped += 1
        if len(envs) > 1:
            summary.projects_multiple_envs += 1
            multiple_env_rows.append((project, ",".join(envs)))

        try:
            fieldnames, key_rows, role_column = extract_key_rows(node_file)
        except KeyError as exc:
            summary.missing_role_column += 1
            project_rows.append(
                (project, ",".join(envs), 0, 0, f"ERROR:{exc}")
            )
            continue
        except Exception as exc:
            summary.malformed_files += 1
            project_rows.append(
                (project, ",".join(envs), 0, 0, f"ERROR:{type(exc).__name__}:{exc}")
            )
            continue

        role_counts: dict[str, int] = defaultdict(int)
        for row in key_rows:
            role_counts[row.get(role_column, "")] += 1

        for original_role, count in sorted(role_counts.items()):
            role_count_rows.append(
                (project, ",".join(envs), original_role, count)
            )

        key_count = len(key_rows)
        summary.key_rows += key_count

        if key_count > 0:
            summary.projects_with_keys += 1
        else:
            summary.projects_without_keys += 1

        written = 0
        if key_count > 0 or write_empty:
            for env_name in envs:
                output_file = out_dir / env_name / f"{project}_key.txt"
                write_key_file(output_file, fieldnames, key_rows)
                written += 1
                summary.key_files_written += 1

        project_rows.append(
            (
                project,
                ",".join(envs),
                key_count,
                written,
                "OK" if key_count > 0 else ("EMPTY_WRITTEN" if write_empty else "NO_KEY_TAXA"),
            )
        )

    # Reports stored directly under keyspecies/keygenus root.
    write_tsv(
        out_dir / "_project_extraction_summary.tsv",
        ["project", "environments", "key_rows", "key_files_written", "status"],
        project_rows,
    )
    write_tsv(
        out_dir / "_role_counts.tsv",
        ["project", "environments", "taxa_roles", "count"],
        role_count_rows,
    )
    write_tsv(
        out_dir / "_unmatched_projects.tsv",
        ["project", "node_file"],
        unmatched_rows,
    )
    write_tsv(
        out_dir / "_multiple_environment_projects.tsv",
        ["project", "environments"],
        multiple_env_rows,
    )
    write_tsv(
        out_dir / "_environment_mapping_duplicates_all.tsv",
        ["project", "environments"],
        duplicate_records,
    )
    write_tsv(
        out_dir / "_run_summary.tsv",
        ["metric", "value"],
        [
            ("level", level),
            ("environment_directories", len(env_names)),
            ("node_files", summary.node_files),
            ("projects_mapped", summary.projects_mapped),
            ("projects_unmapped", summary.projects_unmapped),
            ("projects_multiple_environments", summary.projects_multiple_envs),
            ("projects_with_key_taxa", summary.projects_with_keys),
            ("projects_without_key_taxa", summary.projects_without_keys),
            ("key_rows_total", summary.key_rows),
            ("key_files_written", summary.key_files_written),
            ("missing_taxa_roles_column", summary.missing_role_column),
            ("malformed_files", summary.malformed_files),
            ("write_empty_files", write_empty),
        ],
    )

    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "筛选species和genus node表中的Module hubs、Connectors和"
            "Network hubs，并按环境组织project_key.txt。"
        )
    )
    parser.add_argument("--env-dir", type=Path, default=DEFAULT_ENV_DIR)
    parser.add_argument(
        "--species-node-dir", type=Path, default=DEFAULT_SPECIES_NODE_DIR
    )
    parser.add_argument(
        "--species-out-dir", type=Path, default=DEFAULT_SPECIES_OUT_DIR
    )
    parser.add_argument(
        "--genus-node-dir", type=Path, default=DEFAULT_GENUS_NODE_DIR
    )
    parser.add_argument(
        "--genus-out-dir", type=Path, default=DEFAULT_GENUS_OUT_DIR
    )
    parser.add_argument(
        "--level",
        choices=("species", "genus", "both"),
        default="both",
        help="处理species、genus或两者，默认both。",
    )
    parser.add_argument(
        "--write-empty",
        action="store_true",
        help="没有关键类群时也写出仅含表头的project_key.txt。",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="运行前清空目标环境子目录中的旧结果。",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()

    try:
        env_names, project_to_envs, duplicate_records = read_environment_mapping(
            args.env_dir
        )
    except Exception as exc:
        print(f"错误：读取环境映射失败：{exc}", file=sys.stderr)
        return 1

    print(f"识别环境文件：{len(env_names)} 个")
    if len(env_names) != 45:
        print(
            f"警告：检测到 {len(env_names)} 个 .list 文件，不是预期的45个。",
            file=sys.stderr,
        )
    print(f"环境映射中的唯一project：{len(project_to_envs)} 个")
    print(f"属于多个环境的project：{len(duplicate_records)} 个")

    summaries: list[LevelSummary] = []

    try:
        if args.level in ("species", "both"):
            summaries.append(
                process_level(
                    level="species",
                    node_dir=args.species_node_dir,
                    out_dir=args.species_out_dir,
                    env_names=env_names,
                    project_to_envs=project_to_envs,
                    duplicate_records=duplicate_records,
                    write_empty=args.write_empty,
                    clean=args.clean,
                )
            )

        if args.level in ("genus", "both"):
            summaries.append(
                process_level(
                    level="genus",
                    node_dir=args.genus_node_dir,
                    out_dir=args.genus_out_dir,
                    env_names=env_names,
                    project_to_envs=project_to_envs,
                    duplicate_records=duplicate_records,
                    write_empty=args.write_empty,
                    clean=args.clean,
                )
            )
    except Exception as exc:
        print(f"错误：处理失败：{type(exc).__name__}: {exc}", file=sys.stderr)
        return 1

    print("\n处理完成：")
    for s in summaries:
        print(
            f"{s.level}: node文件={s.node_files}, "
            f"有环境映射={s.projects_mapped}, "
            f"无环境映射={s.projects_unmapped}, "
            f"含关键类群project={s.projects_with_keys}, "
            f"无关键类群project={s.projects_without_keys}, "
            f"关键行总数={s.key_rows}, "
            f"写出key文件={s.key_files_written}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
