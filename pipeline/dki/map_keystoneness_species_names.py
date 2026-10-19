#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import re


def clean_species_name(x):
    """
    清洗 species.matrix 里的物种名。
    例如:
    2075_Amylibacter_sp900197625 -> Amylibacter_sp900197625
    """
    if pd.isna(x):
        return ""

    x = str(x).strip()
    x = re.sub(r"^\d+_", "", x)
    return x


def load_mapping_from_data_mapping(mapping_file):
    """
    优先读取 generate_dki_input_from_abundance.py 生成的:
    DKI/<PROJECT>/data/species_mapping.tsv

    期望列:
    species_id_1based    otu_id    species_name

    species_id_1based 是 DKI 结果里的 species 编号。
    """
    df = pd.read_csv(mapping_file, sep="\t", dtype=str)

    if "species_id_1based" not in df.columns:
        raise ValueError(f"species_mapping.tsv missing species_id_1based: {mapping_file}")

    mapping = {}

    for _, row in df.iterrows():
        sid_raw = row.get("species_id_1based", "")

        if pd.isna(sid_raw) or str(sid_raw).strip() == "":
            continue

        sid = int(float(str(sid_raw)))

        species_name = row.get("species_name", "")
        otu_id = row.get("otu_id", "")

        name = clean_species_name(species_name)

        if name == "" or name.lower() == "nan":
            name = clean_species_name(otu_id)

        mapping[sid] = name

    return mapping


def load_mapping_from_species_matrix(species_file):
    """
    回退读取:
    abundance/species/<PROJECT>_species.matrix

    常见格式:
    0    0_Bacteroides_xxx
    1    1_Faecalibacterium_xxx

    species.matrix 里通常是 0-based；
    Keystoneness_structural_median.csv 里的 species 通常是 1-based；
    所以这里用 idx + 1 作为映射键。
    """
    df = pd.read_csv(
        species_file,
        sep=r"\s+",
        header=None,
        engine="python",
        dtype=str
    )

    if df.shape[1] < 2:
        raise ValueError(f"Invalid species matrix: {species_file}")

    mapping = {}

    for _, row in df.iterrows():
        idx_raw = row.iloc[0]
        raw_name = row.iloc[1]

        if pd.isna(idx_raw):
            continue

        idx = int(float(str(idx_raw)))
        mapping[idx + 1] = clean_species_name(raw_name)

    return mapping


def main():
    parser = argparse.ArgumentParser(
        description="Map numeric species IDs in Keystoneness_structural_median.csv to species names."
    )

    parser.add_argument("--project-list", required=True)
    parser.add_argument("--dki-base", required=True)
    parser.add_argument("--pred-base", required=True)
    parser.add_argument("--species-dir", required=True)
    parser.add_argument("--input-name", default="Keystoneness_structural_median.csv")
    parser.add_argument("--output-name", default="Keystoneness_structural_median_final.csv")
    parser.add_argument("--status-file", required=True)

    args = parser.parse_args()

    project_list = Path(args.project_list)
    dki_base = Path(args.dki_base)
    pred_base = Path(args.pred_base)
    species_dir = Path(args.species_dir)
    status_file = Path(args.status_file)

    status_file.parent.mkdir(parents=True, exist_ok=True)

    with open(project_list, "r") as f:
        projects = [x.strip() for x in f if x.strip()]

    with open(status_file, "w") as sf:
        sf.write(
            "project\tstatus\tinput_file\toutput_file\tmapped_n\tunmapped_n\tmapping_source\tmessage\n"
        )

        for i, project in enumerate(projects, start=1):
            result_dir = pred_base / project / "results"
            input_file = result_dir / args.input_name
            output_file = result_dir / args.output_name

            data_mapping_file = dki_base / project / "data" / "species_mapping.tsv"
            species_file = species_dir / f"{project}_species.matrix"

            print(f"[{i}/{len(projects)}] {project}", flush=True)

            if not input_file.exists() or input_file.stat().st_size == 0:
                print(f"  ⚠️ missing input: {input_file}", flush=True)
                sf.write(
                    f"{project}\tMISSING_INPUT\t{input_file}\t{output_file}\t0\t0\tNA\tmissing Keystoneness file\n"
                )
                continue

            try:
                if data_mapping_file.exists() and data_mapping_file.stat().st_size > 0:
                    mapping = load_mapping_from_data_mapping(data_mapping_file)
                    mapping_source = str(data_mapping_file)
                elif species_file.exists() and species_file.stat().st_size > 0:
                    mapping = load_mapping_from_species_matrix(species_file)
                    mapping_source = str(species_file)
                else:
                    print(f"  ⚠️ missing mapping file", flush=True)
                    sf.write(
                        f"{project}\tMISSING_MAPPING\t{input_file}\t{output_file}\t0\t0\tNA\tmissing species_mapping.tsv and species.matrix\n"
                    )
                    continue

                df = pd.read_csv(input_file)

                if "species" not in df.columns:
                    print(f"  ❌ no species column", flush=True)
                    sf.write(
                        f"{project}\tNO_SPECIES_COLUMN\t{input_file}\t{output_file}\t0\t0\t{mapping_source}\tmissing species column\n"
                    )
                    continue

                # 原始 species 编号
                species_id = pd.to_numeric(df["species"], errors="coerce").astype("Int64")

                # 映射成物种名
                species_name = species_id.map(mapping)

                mapped_n = int(species_name.notna().sum())
                unmapped_n = int(species_name.isna().sum())

                df = df.copy()

                # 保留原始 species 编号
                df.insert(0, "species_id", species_id)

                # species 列替换成真实物种名
                fallback_names = "UNMAPPED_" + species_id.astype(str)
                df["species"] = species_name
                df["species"] = df["species"].fillna(fallback_names)

                # 和原 Keystone_calculate.py 保持一致：删除 Ks_true
                if "Ks_true" in df.columns:
                    df = df.drop(columns=["Ks_true"])

                # 和原 Keystone_calculate.py 保持一致：按 Ks_pred 从大到小排序
                if "Ks_pred" in df.columns:
                    df["Ks_pred"] = pd.to_numeric(df["Ks_pred"], errors="coerce")
                    df = df.sort_values(by="Ks_pred", ascending=False)
                elif "median_Ks_pred" in df.columns:
                    df["median_Ks_pred"] = pd.to_numeric(df["median_Ks_pred"], errors="coerce")
                    df = df.sort_values(by="median_Ks_pred", ascending=False)

                # 调整列顺序
                cols = ["species_id", "species"] + [
                    c for c in df.columns if c not in ["species_id", "species"]
                ]
                df = df[cols]

                df.to_csv(output_file, index=False)

                sf.write(
                    f"{project}\tDONE\t{input_file}\t{output_file}\t{mapped_n}\t{unmapped_n}\t{mapping_source}\tok\n"
                )

                print(f"  ✅ mapped={mapped_n}, unmapped={unmapped_n}", flush=True)

            except Exception as e:
                print(f"  ❌ failed: {e}", flush=True)
                sf.write(
                    f"{project}\tFAILED\t{input_file}\t{output_file}\t0\t0\tNA\t{repr(e)}\n"
                )

    print(f"Status written to: {status_file}", flush=True)


if __name__ == "__main__":
    main()
