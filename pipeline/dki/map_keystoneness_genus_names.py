#!/usr/bin/env python3

import argparse
import os
import re
from pathlib import Path

import pandas as pd


def clean_name(value):
    if pd.isna(value):
        return ""

    value = str(value).strip()

    # 例如：2075_Bacteroides -> Bacteroides
    value = re.sub(r"^\d+_", "", value)

    return value


def load_genus_mapping(genus_file):
    df = pd.read_csv(
        genus_file,
        sep=r"\s+",
        header=None,
        engine="python",
        dtype=str,
    )

    if df.shape[1] < 2:
        raise ValueError(f"Invalid genus matrix: {genus_file}")

    mapping = {}

    for _, row in df.iterrows():
        idx_raw = row.iloc[0]
        name_raw = row.iloc[1]

        if pd.isna(idx_raw):
            continue

        # genus.matrix 通常是 0-based；
        # Ks 文件中的 species ID 通常是 1-based。
        idx = int(float(str(idx_raw)))
        mapping[idx + 1] = clean_name(name_raw)

    return mapping


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--project-list", required=True)
    parser.add_argument("--pred-base", required=True)
    parser.add_argument("--genus-dir", required=True)

    parser.add_argument(
        "--input-name",
        default="Keystoneness_structural_median.csv",
    )

    parser.add_argument(
        "--output-name",
        default="Keystoneness_structural_median_final.csv",
    )

    parser.add_argument("--status-file", required=True)

    args = parser.parse_args()

    pred_base = Path(args.pred_base)
    genus_dir = Path(args.genus_dir)
    status_file = Path(args.status_file)

    status_file.parent.mkdir(parents=True, exist_ok=True)

    with open(args.project_list, "r") as handle:
        projects = [line.strip() for line in handle if line.strip()]

    with open(status_file, "w") as status_handle:
        status_handle.write(
            "project\tstatus\tmapped_n\tunmapped_n\t"
            "input_file\tgenus_file\toutput_file\tmessage\n"
        )

        for index, project in enumerate(projects, start=1):
            result_dir = pred_base / project / "results"

            input_file = result_dir / args.input_name
            output_file = result_dir / args.output_name

            # 优先查找 <PROJECT>_genus.matrix
            genus_file = genus_dir / f"{project}_genus.matrix"

            # 兼容个别仍使用 species 后缀的文件
            if not genus_file.exists():
                alternative = genus_dir / f"{project}_species.matrix"

                if alternative.exists():
                    genus_file = alternative

            print(
                f"[{index}/{len(projects)}] {project}",
                flush=True,
            )

            if not input_file.exists() or input_file.stat().st_size == 0:
                print("  MISSING_KS", flush=True)

                status_handle.write(
                    f"{project}\tMISSING_KS\t0\t0\t"
                    f"{input_file}\t{genus_file}\t{output_file}\t"
                    f"missing median Ks file\n"
                )
                status_handle.flush()
                continue

            if not genus_file.exists() or genus_file.stat().st_size == 0:
                print("  MISSING_GENUS_MAPPING", flush=True)

                status_handle.write(
                    f"{project}\tMISSING_GENUS_MAPPING\t0\t0\t"
                    f"{input_file}\t{genus_file}\t{output_file}\t"
                    f"missing genus matrix\n"
                )
                status_handle.flush()
                continue

            try:
                mapping = load_genus_mapping(genus_file)

                ks = pd.read_csv(input_file)

                if "species" not in ks.columns:
                    raise ValueError(
                        f"missing species ID column: {input_file}"
                    )

                if "Ks_pred" not in ks.columns:
                    if "median_Ks_pred" in ks.columns:
                        ks = ks.rename(
                            columns={
                                "median_Ks_pred": "Ks_pred"
                            }
                        )
                    else:
                        raise ValueError(
                            f"missing Ks_pred column: {input_file}"
                        )

                species_id = pd.to_numeric(
                    ks["species"],
                    errors="coerce",
                )

                ks_pred = pd.to_numeric(
                    ks["Ks_pred"],
                    errors="coerce",
                )

                result = pd.DataFrame(
                    {
                        "species_id": species_id,
                        "Ks_pred": ks_pred,
                    }
                )

                result = result.dropna(
                    subset=["species_id", "Ks_pred"]
                ).copy()

                result["species_id"] = (
                    result["species_id"].astype(int)
                )

                # 映射数字 ID 到 genus 名称
                result["species"] = (
                    result["species_id"].map(mapping)
                )

                mapped_n = int(result["species"].notna().sum())
                unmapped_n = int(result["species"].isna().sum())

                # 未映射编号保留，避免丢失 Ks
                result["species"] = result["species"].fillna(
                    "UNMAPPED_"
                    + result["species_id"].astype(str)
                )

                # 按 Ks_pred 从大到小排列
                result = result.sort_values(
                    by=["Ks_pred", "species_id"],
                    ascending=[False, True],
                    kind="stable",
                ).reset_index(drop=True)

                # 最终严格只保留两列
                result = result[
                    ["species", "Ks_pred"]
                ]

                # 原子写入，避免出现残缺文件
                tmp_file = output_file.with_name(
                    output_file.name + ".tmp"
                )

                result.to_csv(tmp_file, index=False)
                os.replace(tmp_file, output_file)

                print(
                    f"  DONE: rows={len(result)}, "
                    f"mapped={mapped_n}, "
                    f"unmapped={unmapped_n}",
                    flush=True,
                )

                status_handle.write(
                    f"{project}\tDONE\t{mapped_n}\t"
                    f"{unmapped_n}\t{input_file}\t"
                    f"{genus_file}\t{output_file}\tok\n"
                )
                status_handle.flush()

            except Exception as exc:
                print(f"  FAILED: {exc}", flush=True)

                status_handle.write(
                    f"{project}\tFAILED\t0\t0\t"
                    f"{input_file}\t{genus_file}\t"
                    f"{output_file}\t{repr(exc)}\n"
                )
                status_handle.flush()

    print(
        f"Status written to: {status_file}",
        flush=True,
    )


if __name__ == "__main__":
    main()
