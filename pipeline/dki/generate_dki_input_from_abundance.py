#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def read_species_map(species_file):
    species_map = {}
    if not species_file.exists():
        return species_map

    df = pd.read_csv(species_file, sep=r"\s+", header=None, engine="python")
    if df.shape[1] >= 2:
        for _, row in df.iterrows():
            species_map[str(row.iloc[0])] = str(row.iloc[1])
    return species_map


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True)
    parser.add_argument("--abundance-dir", required=True)
    parser.add_argument("--species-dir", required=True)
    parser.add_argument("--out-base", required=True)
    parser.add_argument("--min-abundance", type=float, default=0.0)
    parser.add_argument("--chunk-size", type=int, default=100)
    args = parser.parse_args()

    project = args.project
    abundance_file = Path(args.abundance_dir) / f"{project}_abundance.matrix"
    species_file = Path(args.species_dir) / f"{project}_species.matrix"
    out_dir = Path(args.out_base) / project / "data"
    out_dir.mkdir(parents=True, exist_ok=True)

    if not abundance_file.exists():
        raise FileNotFoundError(f"Missing abundance file: {abundance_file}")

    print(f"[{project}] reading abundance: {abundance_file}", flush=True)

    # abundance.matrix: rows = species, columns = samples
    df = pd.read_csv(abundance_file, sep=r"\s+", engine="python")
    if df.shape[1] < 2:
        raise ValueError(f"Invalid abundance matrix: {abundance_file}")

    otu_col = df.columns[0]
    otu_ids = df[otu_col].astype(str).tolist()
    sample_names = [str(x) for x in df.columns[1:]]

    mat = df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=np.float64)
    mat[mat < 0] = 0.0

    # 删除全零 species / 全零 samples
    species_sum = mat.sum(axis=1)
    sample_sum = mat.sum(axis=0)

    keep_species = species_sum > 0
    keep_samples = sample_sum > 0

    mat = mat[keep_species, :][:, keep_samples]
    kept_otu_ids = [otu_ids[i] for i, keep in enumerate(keep_species) if keep]
    kept_sample_names = [sample_names[i] for i, keep in enumerate(keep_samples) if keep]

    if mat.shape[0] == 0 or mat.shape[1] == 0:
        raise ValueError(f"No nonzero species/sample left for {project}")

    # Ptrain 正确方向：species × samples
    col_sum = mat.sum(axis=0, keepdims=True)
    P = np.divide(mat, col_sum, out=np.zeros_like(mat), where=col_sum > 0)

    n_species, n_samples = P.shape
    print(f"[{project}] Ptrain shape: species={n_species}, samples={n_samples}", flush=True)

    # 保存 Ptrain.csv：species × samples
    ptrain_file = out_dir / "Ptrain.csv"
    np.savetxt(ptrain_file, P, delimiter=",", fmt="%.10g")

    # 保存 species/sample 映射
    species_map = read_species_map(species_file)

    with open(out_dir / "species_mapping.tsv", "w") as f:
        f.write("species_id_1based\totu_id\tspecies_name\n")
        for j, otu in enumerate(kept_otu_ids, start=1):
            f.write(f"{j}\t{otu}\t{species_map.get(str(otu), '')}\n")

    with open(out_dir / "sample_mapping.tsv", "w") as f:
        f.write("sample_id_1based\tsample_name\n")
        for i, s in enumerate(kept_sample_names, start=1):
            f.write(f"{i}\t{s}\n")

    # 构建 remove-one case 列表
    case_sample = []
    case_species = []

    for sample_idx in range(n_samples):
        present_species = np.flatnonzero(P[:, sample_idx] > args.min_abundance)
        for sp_idx in present_species:
            case_sample.append(sample_idx)
            case_species.append(sp_idx)

    case_sample = np.asarray(case_sample, dtype=np.int64)
    case_species = np.asarray(case_species, dtype=np.int64)
    n_cases = len(case_sample)

    print(f"[{project}] Ptest shape will be: species={n_species}, removal_cases={n_cases}", flush=True)

    # 保存 species_id/sample_id：每一行对应 Ptest 的一列
    with open(out_dir / "species_id.csv", "w") as f_sid:
        for sp_idx in case_species:
            f_sid.write(f"{sp_idx + 1}\n")

    with open(out_dir / "sample_id.csv", "w") as f_samid:
        for sample_idx in case_sample:
            f_samid.write(f"{sample_idx + 1}\n")

    # 生成 Ptest.csv：species × removal_cases
    ptest_file = out_dir / "Ptest.csv"

    removed_abund = P[case_species, case_sample]
    denom = 1.0 - removed_abund
    valid = denom > 1e-15

    with open(ptest_file, "w") as f_ptest:
        for sp in range(n_species):
            vals = np.zeros(n_cases, dtype=np.float64)

            if np.any(valid):
                vals[valid] = P[sp, case_sample[valid]] / denom[valid]

            # 被移除的物种在该 removal case 中为 0
            vals[case_species == sp] = 0.0

            # 对于 denom <= 0 的情况，整列保持 0
            vals[~valid] = 0.0

            np.savetxt(f_ptest, vals.reshape(1, -1), delimiter=",", fmt="%.10g")

            if (sp + 1) % 500 == 0:
                print(f"[{project}] wrote Ptest rows: {sp + 1}/{n_species}", flush=True)

    print(f"[{project}] done: {out_dir}", flush=True)


if __name__ == "__main__":
    main()
