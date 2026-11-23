#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Create five publication-ready figures from complete species/genus
FastSpar and FastCCLasso environment-level *.info files.

Figures
-------
1. Environment-level association composition
   (FastSpar-only / common / FastCCLasso-only), species and genus.
2. FastSpar vs FastCCLasso association counts across environments,
   species and genus.
3. Species-vs-genus comparison of FastSpar–FastCCLasso Jaccard overlap.
4. Agreement among common edges:
   A) sign agreement; B) Spearman correlation of average r.
5. BioProject-support long-tail distribution (CCDF), species and genus.

Important
---------
A-B and B-A are treated as the same undirected edge.
Project-level network density is NOT reconstructed from belongs/f2 files,
because f2 already filters edges at the environment level.
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ENV_DISPLAY = {
    "human_gut": "human gut",
    "human_oral": "human oral",
    "human_skin": "human skin",
    "human_airway": "human airway",
    "human_lung": "human lung",
    "human_vagina": "human vagina",
    "mouse_feces": "mouse feces",
    "mouse_cecum": "mouse cecum",
    "pig_feces": "pig feces",
    "bovine_feces": "bovine feces",
    "bovine_rumen": "bovine rumen",
    "chicken_feces": "chicken feces",
    "chicken_cecum": "chicken cecum",
    "horse_feces": "horse feces",
    "canine_feces": "canine feces",
    "cat_fecal": "cat fecal",
    "rat_feces": "rat feces",
    "vole_feces": "vole feces",
    "panda_feces": "panda feces",
    "oceanic_seawater": "oceanic seawater",
    "oceanic_sediment": "oceanic sediment",
    "oceanic_coral-reef": "oceanic coral reef",
    "Fresh-Water_river-or-lake": "freshwater river/lake",
    "Fresh-Water_drinking-water": "drinking water",
    "Fresh-Water_sediment": "freshwater sediment",
    "soil_general": "general soil",
    "soil_rhizosphere": "rhizosphere soil",
    "soil_forest": "forest soil",
    "soil_grassland": "grassland soil",
    "soil_rainforest": "rainforest soil",
    "soil_desert": "desert soil",
    "soil_tundra": "tundra soil",
    "soil_wetland": "wetland soil",
    "sewage_general-wastewater": "general wastewater",
    "sewage_activated-sludge": "activated sludge",
    "sewage_anaeroic-digestor": "anaerobic digestor",
    "sewage_bioreactor": "sewage bioreactor",
    "sewage_hospital-wastewater": "hospital wastewater",
    "Natural-Environment_air": "natural air",
    "Natural-Environment_hot-spring": "hot spring",
    "food_cheese": "cheese",
    "food_milk": "milk",
    "food_pickle": "pickle",
    "honeybee_gut": "honeybee gut",
    "Ancient-Human_Dental-Calculus": "ancient dental calculus",
}


def disp_env(env: str) -> str:
    return ENV_DISPLAY.get(env, env.replace("_", " "))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument(
        "--species-fastspar",
        default="/home/wangmengyang/project/coabundance/cor03/species/belongs/f2",
    )
    p.add_argument(
        "--genus-fastspar",
        default="/home/wangmengyang/project/coabundance/cor03/genus/belongs/f2",
    )
    p.add_argument(
        "--species-cclasso",
        default="/mnt/new_data/all_fastccl_species_count_GTDB226/belongs/f2",
    )
    p.add_argument(
        "--genus-cclasso",
        default="/mnt/new_data/all_fastccl_genus_count_GTDB226/belongs/f2",
    )
    p.add_argument(
        "--outdir",
        default="/home/wangmengyang/project/coabundance/cor03/species_genus_method_figures",
    )
    p.add_argument("--top-n", type=int, default=15)
    p.add_argument("--dpi", type=int, default=600)
    return p.parse_args()


def setup_style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.labelsize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


C_FASTSPAR = "#2C7FB8"
C_COMMON = "#F28E2B"
C_CCLASSO = "#36A148"
C_SPECIES = "#5086C4"
C_GENUS = "#F0A19A"


def save_figure(fig: plt.Figure, base: Path, dpi: int) -> None:
    fig.savefig(str(base) + ".png", dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(str(base) + ".pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(str(base) + ".svg", bbox_inches="tight", facecolor="white")


REQUIRED_COLS = {"species_a", "species_b", "aver_corr", "aver_pv", "proj_count"}


def read_info(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, engine="python")
    missing = REQUIRED_COLS - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing columns: {sorted(missing)}")

    df = df.copy()
    df["species_a"] = df["species_a"].astype(str).str.strip()
    df["species_b"] = df["species_b"].astype(str).str.strip()
    df["aver_corr"] = pd.to_numeric(df["aver_corr"], errors="coerce")
    df["aver_pv"] = pd.to_numeric(df["aver_pv"], errors="coerce")
    df["proj_count"] = pd.to_numeric(df["proj_count"], errors="coerce").fillna(0).astype(int)

    a = df["species_a"].to_numpy(dtype=str)
    b = df["species_b"].to_numpy(dtype=str)
    edge_a = np.where(a <= b, a, b)
    edge_b = np.where(a <= b, b, a)
    df["edge_key"] = pd.Series(edge_a, index=df.index) + "\t" + pd.Series(edge_b, index=df.index)

    if df["edge_key"].duplicated().any():
        warnings.warn(f"{path.name}: duplicated undirected edges found; keeping first.")
        df = df.drop_duplicates("edge_key", keep="first")

    return df


def list_info(directory: Path) -> Dict[str, Path]:
    if not directory.is_dir():
        raise FileNotFoundError(directory)
    return {p.stem: p for p in sorted(directory.glob("*.info"))}


def safe_spearman(x: pd.Series, y: pd.Series) -> float:
    z = pd.DataFrame({"x": x, "y": y}).dropna()
    if len(z) < 3 or z["x"].nunique() < 2 or z["y"].nunique() < 2:
        return np.nan
    return float(z["x"].corr(z["y"], method="spearman"))


def summarize_env(env: str, rank: str, fs: pd.DataFrame, cc: pd.DataFrame) -> dict:
    fs_keys = set(fs["edge_key"])
    cc_keys = set(cc["edge_key"])
    common_keys = fs_keys & cc_keys
    union_keys = fs_keys | cc_keys

    merged = fs[["edge_key", "aver_corr"]].merge(
        cc[["edge_key", "aver_corr"]],
        on="edge_key",
        how="inner",
        suffixes=("_fs", "_cc"),
    )

    if len(merged):
        valid = (
            merged["aver_corr_fs"].notna()
            & merged["aver_corr_cc"].notna()
            & merged["aver_corr_fs"].ne(0)
            & merged["aver_corr_cc"].ne(0)
        )
        same = (
            np.sign(merged.loc[valid, "aver_corr_fs"])
            == np.sign(merged.loc[valid, "aver_corr_cc"])
        )
        sign_agreement = float(same.mean()) if len(same) else np.nan
        spearman = safe_spearman(merged["aver_corr_fs"], merged["aver_corr_cc"])
        median_abs_diff = float((merged["aver_corr_fs"] - merged["aver_corr_cc"]).abs().median())
    else:
        sign_agreement = np.nan
        spearman = np.nan
        median_abs_diff = np.nan

    common = len(common_keys)
    union = len(union_keys)

    return {
        "rank": rank,
        "env": env,
        "fastspar_edges": len(fs_keys),
        "cclasso_edges": len(cc_keys),
        "common_edges": common,
        "fastspar_only": len(fs_keys - cc_keys),
        "cclasso_only": len(cc_keys - fs_keys),
        "union_edges": union,
        "jaccard": common / union if union else np.nan,
        "common_fraction_fastspar": common / len(fs_keys) if fs_keys else np.nan,
        "common_fraction_cclasso": common / len(cc_keys) if cc_keys else np.nan,
        "sign_agreement": sign_agreement,
        "common_corr_spearman": spearman,
        "median_abs_corr_diff": median_abs_diff,
    }


def build_summaries(dirs: Dict[Tuple[str, str], Path]):
    maps = {k: list_info(v) for k, v in dirs.items()}
    all_envs = sorted(set().union(*[set(m.keys()) for m in maps.values()]))

    summary_rows: List[dict] = []
    support_rows: List[dict] = []

    for rank in ("species", "genus"):
        fs_map = maps[(rank, "FastSpar")]
        cc_map = maps[(rank, "CCLasso")]

        for env in all_envs:
            if env not in fs_map or env not in cc_map:
                warnings.warn(f"{rank}/{env}: missing one method; skipped.")
                continue

            fs = read_info(fs_map[env])
            cc = read_info(cc_map[env])
            summary_rows.append(summarize_env(env, rank, fs, cc))

            for method, frame in (("FastSpar", fs), ("CCLasso", cc)):
                vc = frame["proj_count"].value_counts().sort_index()
                for k, n in vc.items():
                    support_rows.append({
                        "rank": rank,
                        "method": method,
                        "env": env,
                        "proj_count": int(k),
                        "edge_count": int(n),
                    })

    return pd.DataFrame(summary_rows), pd.DataFrame(support_rows)


def global_summary(summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for rank, g in summary.groupby("rank"):
        fs = int(g["fastspar_edges"].sum())
        cc = int(g["cclasso_edges"].sum())
        common = int(g["common_edges"].sum())
        union = int(g["union_edges"].sum())
        rows.append({
            "rank": rank,
            "environment_count": g["env"].nunique(),
            "fastspar_environment_edges": fs,
            "cclasso_environment_edges": cc,
            "common_environment_edges": common,
            "union_environment_edges": union,
            "environment_edge_jaccard": common / union if union else np.nan,
        })
    return pd.DataFrame(rows)


# ------------------------------------------------------------------
# Figure 1
# ------------------------------------------------------------------

def figure1(summary: pd.DataFrame, outdir: Path, top_n: int, dpi: int):
    totals = summary.groupby("env")["union_edges"].sum().sort_values(ascending=False)
    envs = list(reversed(totals.head(top_n).index.tolist()))

    fig, axes = plt.subplots(
        1, 2,
        figsize=(13.4, max(6.0, 0.40 * top_n + 2.0)),
        sharey=True,
    )

    for ax, rank, panel in zip(axes, ("species", "genus"), ("A", "B")):
        g = summary[(summary["rank"] == rank) & summary["env"].isin(envs)].set_index("env").reindex(envs)
        y = np.arange(len(envs))
        fs = g["fastspar_only"].to_numpy()
        co = g["common_edges"].to_numpy()
        cc = g["cclasso_only"].to_numpy()

        ax.barh(y, fs, color=C_FASTSPAR, label="FastSpar only", height=0.68)
        ax.barh(y, co, left=fs, color=C_COMMON, label="Common", height=0.68)
        ax.barh(y, cc, left=fs + co, color=C_CCLASSO, label="FastCCLasso only", height=0.68)

        ax.set_yticks(y)
        ax.set_yticklabels([disp_env(e) for e in envs])
        ax.set_xlabel("Number of environment-level associations")
        ax.set_title(rank.capitalize())
        ax.grid(axis="x", linestyle="--", alpha=0.20)
        ax.set_axisbelow(True)
        ax.text(-0.08, 1.03, panel, transform=ax.transAxes, fontsize=16, fontweight="bold")

    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.0))
    fig.subplots_adjust(left=0.20, right=0.98, bottom=0.10, top=0.90, wspace=0.18)
    save_figure(fig, outdir / "Figure_1_environment_association_composition", dpi)
    plt.close(fig)


# ------------------------------------------------------------------
# Figure 2
# ------------------------------------------------------------------

def annotate_top(ax, g: pd.DataFrame, n: int = 9):
    top = g.nlargest(n, "union_edges")
    offsets = [(7, 6), (7, -9), (-7, 8), (-7, -11), (9, 0), (-9, 0), (7, 13), (7, -15), (-7, 15)]
    texts = []

    for i, (_, row) in enumerate(top.iterrows()):
        dx, dy = offsets[i % len(offsets)]
        if row["env"] == "human_gut":
            dx, dy = 9, 13
        elif row["env"] == "oceanic_seawater":
            dx, dy = 9, -15

        texts.append(
            ax.annotate(
                disp_env(row["env"]),
                xy=(row["fastspar_edges"] + 1, row["cclasso_edges"] + 1),
                xytext=(dx, dy),
                textcoords="offset points",
                fontsize=8,
                ha="left" if dx >= 0 else "right",
                va="center",
            )
        )

    try:
        from adjustText import adjust_text
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", lw=0.5, alpha=0.5))
    except Exception:
        pass


def figure2(summary: pd.DataFrame, outdir: Path, dpi: int):
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.4))

    for ax, rank, panel in zip(axes, ("species", "genus"), ("A", "B")):
        g = summary[summary["rank"] == rank].copy()
        x = g["fastspar_edges"].to_numpy() + 1
        y = g["cclasso_edges"].to_numpy() + 1
        max_common = max(1, int(g["common_edges"].max()))
        sizes = 24 + 100 * np.sqrt(g["common_edges"].to_numpy() / max_common)

        ax.scatter(x, y, s=sizes, alpha=0.72, color=C_SPECIES if rank == "species" else C_GENUS, edgecolor="none")
        lo = min(x.min(), y.min())
        hi = max(x.max(), y.max())
        ax.plot([lo, hi], [lo, hi], "--", color="0.4", lw=1.1)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("FastSpar associations + 1")
        ax.set_ylabel("FastCCLasso associations + 1")
        ax.set_title(rank.capitalize())
        ax.grid(which="major", linestyle="--", alpha=0.18)
        annotate_top(ax, g, 9)
        ax.text(-0.10, 1.03, panel, transform=ax.transAxes, fontsize=16, fontweight="bold")

    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.13, top=0.90, wspace=0.22)
    save_figure(fig, outdir / "Figure_2_fastspar_vs_cclasso_counts", dpi)
    plt.close(fig)


# ------------------------------------------------------------------
# Paired rank helper
# ------------------------------------------------------------------

def paired_rank(ax, data: pd.DataFrame, envs: Sequence[str], metric: str, xlabel: str, xlim: Optional[Tuple[float, float]] = None):
    sp = data[data["rank"] == "species"].set_index("env")[metric]
    ge = data[data["rank"] == "genus"].set_index("env")[metric]
    y = np.arange(len(envs))

    for i, env in enumerate(envs):
        a, b = sp.get(env, np.nan), ge.get(env, np.nan)
        if pd.notna(a) and pd.notna(b):
            ax.plot([a, b], [i, i], color="0.75", lw=1.2, zorder=1)

    ax.scatter([sp.get(e, np.nan) for e in envs], y, s=48, color=C_SPECIES, label="Species", zorder=3)
    ax.scatter([ge.get(e, np.nan) for e in envs], y, s=48, color=C_GENUS, label="Genus", zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([disp_env(e) for e in envs])
    ax.set_xlabel(xlabel)
    ax.grid(axis="x", linestyle="--", alpha=0.20)
    ax.set_axisbelow(True)
    if xlim:
        ax.set_xlim(*xlim)


# ------------------------------------------------------------------
# Figure 3
# ------------------------------------------------------------------

def figure3(summary: pd.DataFrame, outdir: Path, top_n: int, dpi: int):
    selected = summary.groupby("env")["union_edges"].sum().nlargest(top_n).index
    envs = (
        summary[summary["env"].isin(selected)]
        .groupby("env")["jaccard"]
        .mean()
        .sort_values()
        .index.tolist()
    )

    fig, ax = plt.subplots(figsize=(8.4, max(5.4, 0.39 * len(envs) + 1.7)))
    xmax = summary["jaccard"].max()
    paired_rank(ax, summary, envs, "jaccard", "FastSpar–FastCCLasso Jaccard overlap", (0, min(1.0, max(0.05, xmax * 1.12))))
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=2, frameon=False)
    fig.subplots_adjust(left=0.27, right=0.97, bottom=0.11, top=0.88)
    save_figure(fig, outdir / "Figure_3_species_genus_jaccard", dpi)
    plt.close(fig)


# ------------------------------------------------------------------
# Figure 4
# ------------------------------------------------------------------

def figure4(summary: pd.DataFrame, outdir: Path, top_n: int, dpi: int):
    selected = summary.groupby("env")["common_edges"].sum().nlargest(top_n).index
    envs = (
        summary[summary["env"].isin(selected)]
        .groupby("env")["sign_agreement"]
        .mean()
        .sort_values()
        .index.tolist()
    )

    fig, axes = plt.subplots(1, 2, figsize=(12.5, max(5.4, 0.39 * len(envs) + 1.7)), sharey=True)

    paired_rank(axes[0], summary, envs, "sign_agreement", "Sign agreement among common edges", (0, 1.02))
    axes[0].text(-0.12, 1.03, "A", transform=axes[0].transAxes, fontsize=16, fontweight="bold")

    paired_rank(axes[1], summary, envs, "common_corr_spearman", "Spearman correlation of average r", (-1, 1))
    axes[1].text(-0.10, 1.03, "B", transform=axes[1].transAxes, fontsize=16, fontweight="bold")

    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.0))
    fig.subplots_adjust(left=0.23, right=0.98, bottom=0.10, top=0.90, wspace=0.18)
    save_figure(fig, outdir / "Figure_4_common_edge_direction_strength_agreement", dpi)
    plt.close(fig)


# ------------------------------------------------------------------
# Figure 5
# ------------------------------------------------------------------

def ccdf(support_dist: pd.DataFrame, rank: str, method: str) -> pd.DataFrame:
    g = support_dist[(support_dist["rank"] == rank) & (support_dist["method"] == method)]
    exact = g.groupby("proj_count")["edge_count"].sum().sort_index()
    if exact.empty:
        return pd.DataFrame(columns=["proj_count", "ccdf_edges"])

    max_k = int(exact.index.max())
    s = pd.Series(0, index=np.arange(1, max_k + 1), dtype=np.int64)
    s.loc[exact.index.astype(int)] = exact.values.astype(np.int64)
    c = s.iloc[::-1].cumsum().iloc[::-1]
    return pd.DataFrame({"proj_count": c.index.astype(int), "ccdf_edges": c.values.astype(np.int64)})


def figure5(support_dist: pd.DataFrame, outdir: Path, dpi: int):
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.2))

    for ax, rank, panel in zip(axes, ("species", "genus"), ("A", "B")):
        for method, color, display_label in (("FastSpar", C_FASTSPAR, "FastSpar"), ("CCLasso", C_CCLASSO, "FastCCLasso")):
            d = ccdf(support_dist, rank, method)
            ax.plot(d["proj_count"], d["ccdf_edges"], lw=1.7, color=color, label=display_label)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Minimum number of supporting BioProjects")
        ax.set_ylabel("Associations with support ≥ x")
        ax.set_title(rank.capitalize())
        ax.grid(which="major", linestyle="--", alpha=0.20)
        ax.text(-0.10, 1.03, panel, transform=ax.transAxes, fontsize=16, fontweight="bold")

    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.0))
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.13, top=0.88, wspace=0.22)
    save_figure(fig, outdir / "Figure_5_bioproject_support_ccdf", dpi)
    plt.close(fig)


def main():
    args = parse_args()
    setup_style()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    dirs = {
        ("species", "FastSpar"): Path(args.species_fastspar),
        ("genus", "FastSpar"): Path(args.genus_fastspar),
        ("species", "CCLasso"): Path(args.species_cclasso),
        ("genus", "CCLasso"): Path(args.genus_cclasso),
    }

    print("Reading complete species/genus FastSpar and FastCCLasso data ...")
    summary, support_dist = build_summaries(dirs)
    gsum = global_summary(summary)

    summary["env_display"] = summary["env"].map(disp_env)
    summary.to_csv(outdir / "environment_method_summary.tsv", sep="\t", index=False)
    support_dist.to_csv(outdir / "support_count_distribution.tsv", sep="\t", index=False)
    gsum.to_csv(outdir / "global_rank_summary.tsv", sep="\t", index=False)

    print("\nGlobal environment-level summary:")
    print(gsum.to_string(index=False))

    print("\nDrawing Figure 1 ...")
    figure1(summary, outdir, args.top_n, args.dpi)
    print("Drawing Figure 2 ...")
    figure2(summary, outdir, args.dpi)
    print("Drawing Figure 3 ...")
    figure3(summary, outdir, args.top_n, args.dpi)
    print("Drawing Figure 4 ...")
    figure4(summary, outdir, args.top_n, args.dpi)
    print("Drawing Figure 5 ...")
    figure5(support_dist, outdir, args.dpi)

    print("\nDone.")
    print("Output directory:", outdir)
    print("NOTE: project-level network density was not reconstructed from f2 files.")


if __name__ == "__main__":
    main()
