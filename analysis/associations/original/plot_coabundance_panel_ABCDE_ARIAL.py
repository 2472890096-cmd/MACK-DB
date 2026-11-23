#!/usr/bin/env python3
"""Draw Figure 2 A–E directly from frozen environment_method_summary.tsv.

Source plotting/statistics script:
/home/wangmengyang/project/coabundance/cor03/species_genus_method_figures/plot_species_genus_coabundance_complete.py

This script reads saved values only. It does not read raw .info networks,
infer associations, refit a model, or recalculate agreement statistics.
Selection and ordering match the source figure1(), figure3(), and figure4().
All panels are drawn on one canvas; PDF/SVG contain no embedded bitmaps.
"""

from pathlib import Path
import argparse
import hashlib
import json
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator, ScalarFormatter, StrMethodFormatter
import numpy as np
import pandas as pd
from PIL import Image

SOURCE = Path('/home/wangmengyang/project/coabundance/cor03/species_genus_method_figures')
DEFAULT_SUMMARY = SOURCE / 'environment_method_summary.tsv'
OUTDIR = Path(__file__).resolve().parent
STEM = 'Figure_2_coabundance_panel_ABCDE_ARIAL_600dpi'
FONT_DIR = Path('/home/wangmengyang/figure_revision_resources/fonts/Arial')
PREVIOUS_SCRIPT = Path('/home/wangmengyang/project/coabundance_figures/panel_600dpi_20261111/plot_coabundance_panel_ABCDE.py')
C_FASTSPAR = '#2C7FB8'
C_COMMON = '#F28E2B'
C_CCLASSO = '#36A148'
C_SPECIES = '#5086C4'
C_GENUS = '#F0A19A'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def selected_environments(summary, top_n=15):
    totals = summary.groupby('env')['union_edges'].sum().sort_values(ascending=False)
    ab = list(reversed(totals.head(top_n).index.tolist()))
    selected_c = summary.groupby('env')['union_edges'].sum().nlargest(top_n).index
    c = (summary[summary['env'].isin(selected_c)].groupby('env')['jaccard']
         .mean().sort_values().index.tolist())
    selected_de = summary.groupby('env')['common_edges'].sum().nlargest(top_n).index
    de = (summary[summary['env'].isin(selected_de)].groupby('env')['sign_agreement']
          .mean().sort_values().index.tolist())
    return ab, c, de


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary', type=Path, default=DEFAULT_SUMMARY)
    args = parser.parse_args()
    before = sha(args.summary)
    summary = pd.read_csv(args.summary, sep='\t', float_precision='round_trip')
    assert not summary.duplicated(['rank', 'env']).any()
    assert summary.groupby('rank').size().to_dict() == {'genus': 45, 'species': 45}
    ab_envs, c_envs, de_envs = selected_environments(summary)
    env_names = summary.drop_duplicates('env').set_index('env')['env_display'].to_dict()
    plotted_rows, checks = [], []

    for filename in ['Arial.TTF', 'Arialbd.TTF', 'Ariali.TTF', 'Arialbi.TTF']:
        font_manager.fontManager.addfont(FONT_DIR / filename)
    actual_fonts = {}
    for style, weight in [('normal', 'normal'), ('normal', 'bold'), ('italic', 'normal')]:
        properties = font_manager.FontProperties(family='Arial', style=style, weight=weight)
        resolved = Path(font_manager.findfont(properties, fallback_to_default=False))
        assert resolved.parent == FONT_DIR and font_manager.FontProperties(fname=resolved).get_name() == 'Arial'
        actual_fonts[f'{style}_{weight}'] = str(resolved)
    plt.rcParams.update({
        'font.family': 'Arial',
        'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial',
        'mathtext.it': 'Arial:italic', 'mathtext.bf': 'Arial:bold',
        'mathtext.bfit': 'Arial:italic:bold', 'mathtext.sf': 'Arial',
        'mathtext.tt': 'Arial', 'mathtext.cal': 'Arial',
        'mathtext.fallback': None, 'mathtext.default': 'regular', 'font.size': 7.5,
        'axes.labelsize': 8, 'axes.titlesize': 9,
        'xtick.labelsize': 6.5, 'ytick.labelsize': 6.8,
        'legend.fontsize': 7.5, 'axes.linewidth': 0.6,
        'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
        'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
        'axes.spines.top': False, 'axes.spines.right': False,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
        'savefig.facecolor': 'white',
    })
    # 180 x 190 mm at the final manuscript size; all text uses this canvas.
    fig = plt.figure(figsize=(180/25.4, 190/25.4))
    positions = {
        'A': [0.170, 0.490, 0.355, 0.425],
        'B': [0.620, 0.490, 0.355, 0.425],
        'C': [0.170, 0.090, 0.265, 0.265],
        'D': [0.635, 0.090, 0.150, 0.265],
        'E': [0.835, 0.090, 0.140, 0.265],
    }
    axes = {letter: fig.add_axes(pos) for letter, pos in positions.items()}

    for letter, rank in [('A', 'species'), ('B', 'genus')]:
        ax = axes[letter]
        rows = summary[summary['rank'] == rank].set_index('env').loc[ab_envs]
        y = np.arange(len(rows))
        fs = rows['fastspar_only'].to_numpy()
        common = rows['common_edges'].to_numpy()
        cc = rows['cclasso_only'].to_numpy()
        for widths, left, color in [(fs, np.zeros(len(fs)), C_FASTSPAR),
                                    (common, fs, C_COMMON),
                                    (cc, fs+common, C_CCLASSO)]:
            bars = ax.barh(y, widths, left=left, height=0.68, color=color)
            np.testing.assert_array_equal([b.get_width() for b in bars], widths)
        ax.set_yticks(y, [env_names[e] for e in ab_envs])
        ax.set_ylim(-0.7, len(rows)-0.3)
        if letter == 'B':
            ax.tick_params(axis='y', labelleft=False)
        ax.set_title(rank.capitalize(), pad=5)
        ax.set_xlabel('Number of associations\nin each environment', labelpad=4)
        ax.set_xlim(left=0)
        locator = MultipleLocator(200_000 if rank == 'species' else 100_000)
        ax.xaxis.set_major_locator(locator)
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((0, 0))
        ax.xaxis.set_major_formatter(formatter)
        ax.xaxis.get_offset_text().set_fontsize(6.5)
        plotted_rows.append(rows.reset_index().assign(panel=letter))

    def paired_plot(letter, envs, field, xlabel, xlim, ticks, show_y=True):
        ax = axes[letter]
        sp = summary[summary['rank'] == 'species'].set_index('env').loc[envs]
        ge = summary[summary['rank'] == 'genus'].set_index('env').loc[envs]
        y = np.arange(len(envs))
        a, b = sp[field].to_numpy(), ge[field].to_numpy()
        assert np.isfinite(a).all() and np.isfinite(b).all()
        for yy, av, bv in zip(y, a, b):
            ax.plot([av, bv], [yy, yy], color='0.75', lw=0.65, zorder=1)
        # Separate identical-sized markers, using the source colors and values.
        for vals, color in [(a, C_SPECIES), (b, C_GENUS)]:
            points = ax.scatter(vals, y, s=10, color=color, linewidths=0, zorder=3)
            np.testing.assert_array_equal(points.get_offsets()[:, 0], vals)
        ax.set_yticks(y, [env_names[e] for e in envs])
        ax.set_ylim(-0.7, len(envs)-0.3)
        ax.tick_params(axis='y', labelleft=show_y)
        ax.set_xlabel(xlabel, labelpad=4)
        ax.set_xlim(*xlim)
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(StrMethodFormatter('{x:g}'))
        plotted_rows.extend([sp.reset_index().assign(panel=letter),
                             ge.reset_index().assign(panel=letter)])
        checks.append({'panel': letter, 'field': field, 'all_coordinates_match_saved_values': True})

    xmax = min(1.0, max(0.05, summary['jaccard'].max()*1.12))
    paired_plot('C', c_envs, 'jaccard', 'FastSpar–FastCCLasso\nJaccard overlap',
                (0, xmax), [0, .2, .4, .6])
    paired_plot('D', de_envs, 'sign_agreement', 'Sign agreement\namong common edges',
                (0, 1.02), [0, .5, 1])
    paired_plot('E', de_envs, 'common_corr_spearman', 'Spearman correlation\nof mean $\\mathit{r}$',
                (-1, 1), [-1, 0, 1], show_y=False)

    for ax in axes.values():
        ax.grid(axis='x', linestyle='--', linewidth=0.45, alpha=0.20)
        ax.set_axisbelow(True)
    for letter, x, y in [('A', .045, .940), ('B', .580, .940),
                         ('C', .045, .378), ('D', .510, .378), ('E', .810, .378)]:
        fig.text(x, y, letter, fontsize=11, fontweight='bold', va='bottom')

    fig.legend(handles=[Patch(color=C_FASTSPAR, label='FastSpar only'),
                        Patch(color=C_COMMON, label='Common'),
                        Patch(color=C_CCLASSO, label='FastCCLasso only')],
               loc='upper center', bbox_to_anchor=(.56, .992), ncol=3,
               frameon=False, handlelength=1.4, columnspacing=1.6)
    fig.legend(handles=[Line2D([], [], marker='o', color=C_SPECIES, lw=0, markersize=3.3, label='Species'),
                        Line2D([], [], marker='o', color=C_GENUS, lw=0, markersize=3.3, label='Genus')],
               loc='center', bbox_to_anchor=(.565, .415), ncol=2,
               frameon=False, handlelength=1, columnspacing=1.5)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    clipped = []
    for text in fig.findobj(matplotlib.text.Text):
        if not text.get_visible() or not text.get_text():
            continue
        bb = text.get_window_extent(renderer)
        # Off-axis automatic tick labels beyond x limits are not rendered.
        if text in [t for ax in axes.values() for t in ax.get_xticklabels()]:
            owner = next(ax for ax in axes.values() if text in ax.get_xticklabels())
            if not owner.get_xlim()[0] <= text.get_position()[0] <= owner.get_xlim()[1]:
                continue
        if bb.x0 < 0 or bb.y0 < 0 or bb.x1 > fig.bbox.width or bb.y1 > fig.bbox.height:
            clipped.append(text.get_text())
    assert not clipped, f'Text outside canvas: {clipped}'
    for letter, ax in axes.items():
        labels = [t for t in ax.get_yticklabels() if t.get_visible()]
        boxes = [t.get_window_extent(renderer) for t in labels]
        assert all(not first.overlaps(second) for first, second in zip(boxes, boxes[1:])), letter

    # Cross-panel labels and legends must stay separate at the final print size.
    layout_artists = list(fig.texts) + list(fig.legends)
    for ax in axes.values():
        layout_artists.extend([ax.xaxis.label, ax.title, ax.xaxis.get_offset_text()])
    layout_artists = [a for a in layout_artists
                      if not isinstance(a, matplotlib.text.Text) or a.get_text()]
    overlap_pairs = []
    for i, first in enumerate(layout_artists):
        for second in layout_artists[i+1:]:
            if first.get_window_extent(renderer).overlaps(second.get_window_extent(renderer)):
                overlap_pairs.append([str(first), str(second)])
    assert not overlap_pairs, f'Labels or legends overlap: {overlap_pairs}'

    for ext in ['png', 'pdf', 'svg', 'tiff']:
        kwargs = {'dpi': 600, 'facecolor': 'white'}
        if ext == 'tiff':
            kwargs['pil_kwargs'] = {'compression': 'tiff_lzw'}
        fig.savefig(OUTDIR / f'{STEM}.{ext}', **kwargs)
    plt.close(fig)
    pd.concat(plotted_rows, ignore_index=True).to_csv(
        OUTDIR/'panel_plot_values.tsv', sep='\t', index=False)
    assert sha(args.summary) == before
    svg = ET.parse(OUTDIR/f'{STEM}.svg').getroot()
    assert not any(n.tag.endswith('image') for n in svg.iter()), 'SVG must be entirely vector.'
    with Image.open(OUTDIR/f'{STEM}.png') as im:
        pixel_size, dpi = list(im.size), list(im.info['dpi'])
    audit = {
        'source_script': str(SOURCE/'plot_species_genus_coabundance_complete.py'),
        'previous_plotting_script': str(PREVIOUS_SCRIPT),
        'font_family': 'Arial', 'resolved_fonts': actual_fonts,
        'E_x_axis': 'Spearman correlation of mean r', 'E_r_is_italic': True,
        'input_summary': str(args.summary), 'input_sha256_unchanged': before,
        'input_rows': len(summary), 'original_summary_statistics_recalculated': False,
        'top_n': 15, 'panel_order_top_to_bottom': {
            'A_B': list(reversed(ab_envs)), 'C': list(reversed(c_envs)),
            'D_E': list(reversed(de_envs))},
        'A_B_x_axis': 'Number of associations in each environment',
        'figure_size_mm': [180, 190], 'png_pixels': pixel_size, 'png_dpi': dpi,
        'svg_embedded_images': 0, 'text_inside_canvas': True,
        'axis_labels_panel_letters_and_legends_do_not_overlap': True,
        'environment_labels_do_not_overlap': True, 'saved_value_checks': checks,
        'outputs': [str(OUTDIR/f'{STEM}.{e}') for e in ['png','pdf','svg','tiff']],
    }
    (OUTDIR/'render_audit.json').write_text(json.dumps(audit, indent=2)+'\n')
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
