#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Prepare species- and genus-level DKI keystone taxa for MySQL import.

For each project:
  1. Read results/Keystoneness_structural_median_final.csv
  2. Sort taxa by Ks_pred in descending order
  3. Let N be the number of valid unique taxa
  4. Candidate pool size:
         K = max(3, ceil(0.01 * N))
     capped at N when N < 3
  5. Compute project-wide mean and SD from all valid Ks_pred
  6. Retain only taxa within the top-K pool satisfying:
         Ks_pred >= mean + 2 * SD
  7. Add environment from data/env/*.list
  8. Write MySQL-ready TSV:
         name, bioproject, environment, ks_score, rank_in_project

Default SD mode is sample SD (N-1), matching R sd(); use --sd-mode population\nfor population SD (N). Species and genus use the same selection rule. The first input column may be\nnamed "species" even for genus-level result files.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


DEFAULT_SPECIES_ROOT = Path(
    "work/dki/species"
)
DEFAULT_GENUS_ROOT = Path(
    "work/dki/genus"
)
DEFAULT_ENV_DIR = Path(
    "data/env"
)
DEFAULT_OUTPUT_DIR = Path(
    "work/candidates/"
    "two_pairs_test/mysql_import_dki"
)

RESULT_RELATIVE_PATH = Path("results/Keystoneness_structural_median_final.csv")


@dataclass
class LevelStats:
    level: str
    project_dirs: int = 0
    result_files_found: int = 0
    missing_result_files: int = 0
    malformed_projects: int = 0
    unmatched_environment_projects: int = 0
    multiple_environment_projects: int = 0
    projects_with_candidates: int = 0
    projects_without_candidates: int = 0
    unique_project_taxa: int = 0
    mysql_rows: int = 0


def normalize_project(value: str) -> str:
    return value.strip().lstrip("\ufeff")


def read_environment_mapping(
    env_dir: Path,
) -> Tuple[List[str], Dict[str, List[str]]]:
    if not env_dir.is_dir():
        raise FileNotFoundError(f"环境目录不存在：{env_dir}")

    list_files = sorted(env_dir.glob("*.list"))
    if not list_files:
        raise FileNotFoundError(f"没有找到环境 .list 文件：{env_dir}")

    env_names: List[str] = []
    project_to_envs: Dict[str, List[str]] = defaultdict(list)

    for list_file in list_files:
        environment = list_file.stem
        env_names.append(environment)

        with list_file.open("r", encoding="utf-8-sig") as handle:
            for line in handle:
                project = normalize_project(line)
                if not project or project.startswith("#"):
                    continue
                if environment not in project_to_envs[project]:
                    project_to_envs[project].append(environment)

    return env_names, dict(project_to_envs)


def resolve_columns(fieldnames: Sequence[str] | None) -> Tuple[str, str]:
    if not fieldnames:
        raise ValueError("CSV缺少表头")

    stripped = {name.strip(): name for name in fieldnames}
    lower_map = {name.strip().lower(): name for name in fieldnames}

    taxon_candidates = (
        "species",
        "genus",
        "name",
        "taxon",
        "taxa",
    )
    score_candidates = (
        "ks_pred",
        "k_s_pred",
        "ks",
        "ks_score",
        "keystoneness",
    )

    taxon_col = next(
        (lower_map[x] for x in taxon_candidates if x in lower_map),
        None,
    )
    score_col = next(
        (lower_map[x] for x in score_candidates if x in lower_map),
        None,
    )

    if taxon_col is None:
        raise ValueError(
            f"未找到分类单元列。实际表头：{list(fieldnames)}"
        )
    if score_col is None:
        raise ValueError(
            f"未找到Ks_pred列。实际表头：{list(fieldnames)}"
        )

    return taxon_col, score_col


def read_project_scores(
    csv_path: Path,
) -> Tuple[List[Tuple[str, float]], int, int]:
    """
    Return:
      sorted unique (taxon, score),
      invalid row count,
      duplicate taxon count.

    For duplicated taxa, keep the highest Ks_pred and record the duplicate count.
    """
    best_score: Dict[str, float] = {}
    invalid_rows = 0
    duplicate_taxa = 0

    with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        taxon_col, score_col = resolve_columns(reader.fieldnames)

        for row in reader:
            taxon = (row.get(taxon_col) or "").strip()
            score_text = (row.get(score_col) or "").strip()

            if not taxon or not score_text:
                invalid_rows += 1
                continue

            try:
                score = float(score_text)
            except ValueError:
                invalid_rows += 1
                continue

            if not math.isfinite(score):
                invalid_rows += 1
                continue

            if taxon in best_score:
                duplicate_taxa += 1
                if score > best_score[taxon]:
                    best_score[taxon] = score
            else:
                best_score[taxon] = score

    ordered = sorted(
        best_score.items(),
        key=lambda x: (-x[1], x[0]),
    )
    return ordered, invalid_rows, duplicate_taxa


def calculate_threshold(
    scores: Sequence[float],
    sd_mode: str,
) -> Tuple[float, float, float]:
    if not scores:
        raise ValueError("没有有效Ks_pred")

    mean_value = statistics.fmean(scores)

    if len(scores) >= 2:
        if sd_mode == "sample":
            sd_value = statistics.stdev(scores)   # denominator N-1
        else:
            sd_value = statistics.pstdev(scores)  # denominator N
    else:
        sd_value = 0.0

    threshold = mean_value + 2.0 * sd_value
    return mean_value, sd_value, threshold


def write_tsv(
    path: Path,
    header: Sequence[str],
    rows: Iterable[Sequence[object]],
) -> int:
    count = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writerow(list(header))
        for row in rows:
            writer.writerow(list(row))
            count += 1
    return count


def process_level(
    level: str,
    root: Path,
    out_dir: Path,
    project_to_envs: Dict[str, List[str]],
    sd_mode: str,
) -> LevelStats:
    if not root.is_dir():
        raise FileNotFoundError(f"{level}根目录不存在：{root}")

    stats = LevelStats(level=level)

    project_dirs = sorted(
        path for path in root.iterdir()
        if path.is_dir() and path.name.startswith("PRJ")
    )
    stats.project_dirs = len(project_dirs)

    mysql_rows: List[Tuple[str, str, str, float, int]] = []
    selected_audit_rows: List[Tuple[object, ...]] = []
    topk_audit_rows: List[Tuple[object, ...]] = []
    project_summary_rows: List[Tuple[object, ...]] = []
    missing_rows: List[Tuple[str, str]] = []
    malformed_rows: List[Tuple[str, str, str]] = []
    unmatched_env_rows: List[Tuple[str, str]] = []
    multi_env_rows: List[Tuple[str, str]] = []
    no_candidate_rows: List[Tuple[object, ...]] = []

    seen_mysql = set()

    for project_dir in project_dirs:
        project = project_dir.name
        result_file = project_dir / RESULT_RELATIVE_PATH

        if not result_file.is_file():
            stats.missing_result_files += 1
            missing_rows.append((project, str(result_file)))
            project_summary_rows.append(
                (
                    project, 0, 0, "", "", "", 0, "",
                    "MISSING_RESULT_FILE",
                )
            )
            continue

        stats.result_files_found += 1

        try:
            ordered, invalid_rows, duplicate_taxa = read_project_scores(
                result_file
            )
            if not ordered:
                raise ValueError("没有有效分类单元及Ks_pred")

            scores = [score for _, score in ordered]
            n_taxa = len(ordered)
            k_raw = max(3, math.ceil(0.01 * n_taxa))
            k_value = min(n_taxa, k_raw)

            mean_value, sd_value, threshold = calculate_threshold(
                scores,
                sd_mode,
            )

            top_k = ordered[:k_value]
            selected = [
                (rank, taxon, score)
                for rank, (taxon, score) in enumerate(top_k, start=1)
                if score >= threshold
            ]
        except Exception as exc:
            stats.malformed_projects += 1
            malformed_rows.append(
                (project, str(result_file), f"{type(exc).__name__}: {exc}")
            )
            project_summary_rows.append(
                (
                    project, 0, 0, "", "", "", 0, "",
                    f"MALFORMED:{type(exc).__name__}",
                )
            )
            continue

        envs = project_to_envs.get(project, [])
        if not envs:
            stats.unmatched_environment_projects += 1
            unmatched_env_rows.append((project, str(result_file)))
            project_summary_rows.append(
                (
                    project,
                    n_taxa,
                    k_value,
                    mean_value,
                    sd_value,
                    threshold,
                    len(selected),
                    "",
                    "NO_ENVIRONMENT_MAPPING",
                )
            )
            continue

        if len(envs) > 1:
            stats.multiple_environment_projects += 1
            multi_env_rows.append((project, ";".join(envs)))

        if selected:
            stats.projects_with_candidates += 1
        else:
            stats.projects_without_candidates += 1
            no_candidate_rows.append(
                (
                    project,
                    ";".join(envs),
                    n_taxa,
                    k_value,
                    mean_value,
                    sd_value,
                    threshold,
                )
            )

        stats.unique_project_taxa += len(selected)

        selected_set = {taxon for _, taxon, _ in selected}

        for rank, (taxon, score) in enumerate(top_k, start=1):
            topk_audit_rows.append(
                (
                    project,
                    taxon,
                    rank,
                    score,
                    n_taxa,
                    k_value,
                    mean_value,
                    sd_value,
                    threshold,
                    int(score >= threshold),
                    ";".join(envs),
                )
            )

        for rank, taxon, score in selected:
            selected_audit_rows.append(
                (
                    project,
                    taxon,
                    score,
                    rank,
                    n_taxa,
                    k_value,
                    mean_value,
                    sd_value,
                    threshold,
                    ";".join(envs),
                )
            )

            for environment in envs:
                row = (
                    taxon,
                    project,
                    environment,
                    score,
                    rank,
                )
                if row not in seen_mysql:
                    seen_mysql.add(row)
                    mysql_rows.append(row)

        project_summary_rows.append(
            (
                project,
                n_taxa,
                k_value,
                mean_value,
                sd_value,
                threshold,
                len(selected),
                ";".join(envs),
                "OK" if selected else "NO_CANDIDATES",
                invalid_rows,
                duplicate_taxa,
            )
        )

    mysql_rows.sort(key=lambda x: (x[2], x[1], x[4], x[0]))
    selected_audit_rows.sort(key=lambda x: (x[0], x[3], x[1]))
    topk_audit_rows.sort(key=lambda x: (x[0], x[2], x[1]))

    mysql_path = out_dir / f"key_{level}_dki_new.tsv"
    selected_audit_path = out_dir / f"key_{level}_dki_selected_audit.tsv"
    topk_audit_path = out_dir / f"key_{level}_dki_topk_audit.tsv"
    project_summary_path = out_dir / f"key_{level}_dki_project_summary.tsv"
    missing_path = out_dir / f"key_{level}_dki_missing_results.tsv"
    malformed_path = out_dir / f"key_{level}_dki_malformed_projects.tsv"
    unmatched_path = out_dir / f"key_{level}_dki_unmatched_environment.tsv"
    multi_env_path = out_dir / f"key_{level}_dki_multiple_environments.tsv"
    no_candidate_path = out_dir / f"key_{level}_dki_no_candidates.tsv"
    run_summary_path = out_dir / f"key_{level}_dki_run_summary.tsv"

    stats.mysql_rows = write_tsv(
        mysql_path,
        [
            "name",
            "bioproject",
            "environment",
            "ks_score",
            "rank_in_project",
        ],
        mysql_rows,
    )

    write_tsv(
        selected_audit_path,
        [
            "bioproject",
            "name",
            "ks_score",
            "rank_in_project",
            "N",
            "K",
            "mean_Ks",
            "sample_SD_Ks",
            "threshold_mean_plus_2SD",
            "environments",
        ],
        selected_audit_rows,
    )

    write_tsv(
        topk_audit_path,
        [
            "bioproject",
            "name",
            "rank_in_project",
            "ks_score",
            "N",
            "K",
            "mean_Ks",
            "sample_SD_Ks",
            "threshold_mean_plus_2SD",
            "passes_threshold",
            "environments",
        ],
        topk_audit_rows,
    )

    write_tsv(
        project_summary_path,
        [
            "bioproject",
            "N",
            "K",
            "mean_Ks",
            "sample_SD_Ks",
            "threshold_mean_plus_2SD",
            "selected_count",
            "environments",
            "status",
            "invalid_input_rows",
            "duplicate_taxa_rows",
        ],
        project_summary_rows,
    )

    write_tsv(
        missing_path,
        ["bioproject", "expected_file"],
        missing_rows,
    )

    write_tsv(
        malformed_path,
        ["bioproject", "file", "error"],
        malformed_rows,
    )

    write_tsv(
        unmatched_path,
        ["bioproject", "result_file"],
        unmatched_env_rows,
    )

    write_tsv(
        multi_env_path,
        ["bioproject", "environments"],
        multi_env_rows,
    )

    write_tsv(
        no_candidate_path,
        [
            "bioproject",
            "environments",
            "N",
            "K",
            "mean_Ks",
            "sample_SD_Ks",
            "threshold_mean_plus_2SD",
        ],
        no_candidate_rows,
    )

    write_tsv(
        run_summary_path,
        ["metric", "value"],
        [
            ("level", level),
            ("input_root", str(root)),
            ("project_directories", stats.project_dirs),
            ("result_files_found", stats.result_files_found),
            ("missing_result_files", stats.missing_result_files),
            ("malformed_projects", stats.malformed_projects),
            (
                "unmatched_environment_projects",
                stats.unmatched_environment_projects,
            ),
            (
                "multiple_environment_projects",
                stats.multiple_environment_projects,
            ),
            ("projects_with_candidates", stats.projects_with_candidates),
            ("projects_without_candidates", stats.projects_without_candidates),
            ("unique_project_taxa", stats.unique_project_taxa),
            ("mysql_rows", stats.mysql_rows),
            ("K_formula", "max(3, ceil(0.01*N)), capped at N"),
            ("threshold_scope", "all valid project taxa"),
            (
                "SD_definition",
                "sample SD, denominator N-1"
                if sd_mode == "sample"
                else "population SD, denominator N",
            ),
            ("selection_rule", "rank<=K and Ks_pred>=mean+2*SD"),
        ],
    )

    print(
        f"{level}: project_dirs={stats.project_dirs}, "
        f"results={stats.result_files_found}, "
        f"missing={stats.missing_result_files}, "
        f"malformed={stats.malformed_projects}, "
        f"projects_with_candidates={stats.projects_with_candidates}, "
        f"projects_without_candidates={stats.projects_without_candidates}, "
        f"unique_project_taxa={stats.unique_project_taxa}, "
        f"MySQL rows={stats.mysql_rows}"
    )
    print(f"{level} MySQL导入文件：{mysql_path}")
    print(f"{level}项目审计文件：{project_summary_path}")

    return stats


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "按 top-K 与 mean+2SD 规则筛选species/genus DKI关键类群，"
            "并生成MySQL导入TSV。"
        )
    )
    parser.add_argument(
        "--species-root",
        type=Path,
        default=DEFAULT_SPECIES_ROOT,
    )
    parser.add_argument(
        "--genus-root",
        type=Path,
        default=DEFAULT_GENUS_ROOT,
    )
    parser.add_argument(
        "--env-dir",
        type=Path,
        default=DEFAULT_ENV_DIR,
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
    )
    parser.add_argument(
        "--level",
        choices=("species", "genus", "both"),
        default="both",
    )
    parser.add_argument(
        "--sd-mode",
        choices=("sample", "population"),
        default="sample",
        help=(
            "标准差定义。sample使用N-1分母，与R的sd()一致；"
            "population使用N分母。默认sample。"
        ),
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    env_names, project_to_envs = read_environment_mapping(args.env_dir)
    print(f"环境列表：{len(env_names)} 个")
    if len(env_names) != 45:
        print(
            f"警告：检测到 {len(env_names)} 个环境 .list 文件，不是预期的45个。",
            file=sys.stderr,
        )
    print(f"环境映射唯一project：{len(project_to_envs)} 个")

    if args.level in ("species", "both"):
        process_level(
            "species",
            args.species_root,
            args.output_dir,
            project_to_envs,
            args.sd_mode,
        )

    if args.level in ("genus", "both"):
        process_level(
            "genus",
            args.genus_root,
            args.output_dir,
            project_to_envs,
            args.sd_mode,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
