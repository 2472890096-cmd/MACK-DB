#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Update MACK-DB Figure 3 to the 1571-project, EPI v3, unified DKI-R release.

Only the heights of the 22 EXISTING bars and their text labels are eligible for
change. All existing SVG layout, axes, limits, colors, typography, legend,
UpSet dots, and UpSet connectors are retained character-for-character.

The user-specified previous result has 22 counts, 10 changes, 12 unchanged.
Any older vector Figure 3 with the same bar/text IDs is also supported;
its actual input counts are recorded in the audit report rather than silently
misrepresented as the immediately preceding dataset.

Inputs, in order of precedence:
    --input-svg /path/to/Figure3_latest.svg
    --docx /path/to/manuscript.docx
    bundled Figure3_EPIv3_SOURCE.svg alongside this script

Outputs: .svg (pure original vector markup with only edited bar/label groups),
         .pdf, .png (4200 px wide, tagged as 600dpi), audit TSV/JSON.

Requires: Python>=3.8; cairoSVG and Pillow for PDF/PNG (SVG itself only stdlib).
Use --svg-only if image-rendering dependencies are not installed.
"""

import argparse
import hashlib
import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from zipfile import ZipFile

# (panel, logical name, bar SVG patch ID, label SVG text ID,
#  previous user-spec value, final official value)
SPECS = [
    ("A", "Zi-Pi only",       3,  8, 4576, 4575),
    ("A", "EPI only",         4,  9, 3605, 3604),
    ("A", "DKI only",         5, 10, 4780, 4786),
    ("A", "Zi-Pi + EPI",      6, 11,  878,  879),
    ("A", "Zi-Pi + DKI",      7, 12, 1140, 1141),
    ("A", "EPI + DKI",        8, 13,  355,  356),
    ("A", "All three",        9, 14,  435,  434),
    ("B", "Zi-Pi only",      15, 26, 2108, 2108),
    ("B", "EPI only",        16, 27, 1260, 1260),
    ("B", "DKI only",        17, 28, 1826, 1826),
    ("B", "Zi-Pi + EPI",     18, 29,  330,  330),
    ("B", "Zi-Pi + DKI",     19, 30,  621,  621),
    ("B", "EPI + DKI",       20, 31,  161,  161),
    ("B", "All three",       21, 32,  273,  273),
    ("C", "Species >=2",    27, 47, 2490, 2491),
    ("C", "Species all3",   28, 48,  417,  416),
    ("C", "Genus >=2",      29, 49,  982,  982),
    ("C", "Genus all3",     30, 50,  243,  243),
    ("D", "Species >=2",    36, 64,  972,  970),
    ("D", "Species all3",   37, 65,   25,   25),
    ("D", "Genus >=2",      38, 66,  399,  399),
    ("D", "Genus all3",     39, 67,   19,   19),
]

RECT = re.compile(
    r"^M(?P<x1>\S+) (?P<bottom1>\S+) "
    r"(?P<x2>\S+) (?P<bottom2>\S+) "
    r"(?P<x3>\S+) (?P<top1>\S+) "
    r"(?P<x4>\S+) (?P<top2>\S+)Z$"
)


def digest(data):
    if isinstance(data, str):
        data = data.encode("utf-8")
    return hashlib.sha256(data).hexdigest()


def identify_source(args):
    if args.input_svg:
        p = Path(args.input_svg).expanduser().resolve()
        return str(p), p.read_bytes()
    if args.docx:
        p = Path(args.docx).expanduser().resolve()
        hits = []
        with ZipFile(p, "r") as archive:
            for member in archive.namelist():
                if member.startswith("word/media/") and member.lower().endswith(".svg"):
                    data = archive.read(member)
                    if (b"Species (same environment)" in data and
                        b"Strict consensus within the same BioProject" in data):
                        hits.append((member, data))
        if len(hits) != 1:
            raise ValueError(f"Expected one Figure 3 SVG inside {p}; found {len(hits)}. "
                             "Pass --input-svg to select the correct source.")
        return f"{p}::{hits[0][0]}", hits[0][1]
    bundled = Path(__file__).with_name("Figure3_EPIv3_SOURCE.svg")
    if not bundled.exists():
        # For installation next to a user's existing old Figure3.
        bundled = Path(__file__).with_name("Figure3_EPIv3.svg")
    if not bundled.exists():
        raise FileNotFoundError(
            "No source SVG was found beside this script. Supply --input-svg "
            "<your latest Figure3.svg> or --docx <your manuscript.docx>, "
            "or unzip the bundle with Figure3_EPIv3_SOURCE.svg."
        )
    return str(bundled), bundled.read_bytes()


def exactly_one(pattern, text, what):
    matches = list(re.finditer(pattern, text, re.DOTALL))
    if len(matches) != 1:
        raise ValueError(f"{what}: expected exactly one match, got {len(matches)}")
    return matches[0]


def fnum(n):
    return ("%.6f" % n).rstrip("0").rstrip(".")


def bar_group_regex(pid):
    return (rf'(?P<before><g id="patch_{pid}"><path d=")'
            rf'(?P<geometry>[^"]+)'
            rf'(?P<after>"[^>]*?/></g>)')


def label_group_regex(tid):
    return (rf'(?P<before><g id="text_{tid}"><text\b[^>]*?\by=")'
            rf'(?P<y>[-+\d.eE]+)'
            rf'(?P<middle>"[^>]*>)'
            rf'(?P<number>\d+)'
            rf'(?P<after></text></g>)')


def read_labels(svg):
    return [int(exactly_one(label_group_regex(tid), svg, f"text_{tid}").group("number"))
            for _, _, _, tid, _, _ in SPECS]


def update_one(svg, item):
    panel, role, pid, tid, previous, target = item
    bar = exactly_one(bar_group_regex(pid), svg, f"panel {panel}, patch_{pid}")
    lab = exactly_one(label_group_regex(tid), svg, f"panel {panel}, text_{tid}")
    geom = RECT.fullmatch(bar.group("geometry"))
    if geom is None:
        raise ValueError(f"patch_{pid} is not an expected rectangular bar")
    d = geom.groupdict()
    x1, x2, x3, x4 = (float(d[k]) for k in ("x1", "x2", "x3", "x4"))
    b1, b2, t1, t2 = (float(d[k]) for k in ("bottom1", "bottom2", "top1", "top2"))
    if (abs(b1-b2) > 1e-4 or abs(t1-t2) > 1e-4
        or abs(x1-x4) > 1e-4 or abs(x2-x3) > 1e-4 or not b1 > t1):
        raise ValueError(f"Bar geometry failed validation: patch_{pid}")
    old = int(lab.group("number"))
    if old <= 0 or target <= 0:
        raise ValueError(f"Invalid count in {panel}/{role}: {old}->{target}")
    old_height = b1-t1
    if old == target:
        return svg, {
            "panel": panel, "role": role, "bar_patch_id": pid,
            "label_text_id": tid, "previous_spec": previous,
            "source_count": old, "final_count": target,
            "source_to_final_changed": False,
            "previous_to_final_changed": previous != target,
            "source_bar_height": old_height, "final_bar_height": old_height,
        }
    new_height = old_height * target / old
    new_top = b1-new_height
    label_y = float(lab.group("y")) + (new_top-t1)
    new_geom = (f'M{d["x1"]} {d["bottom1"]} '
                f'{d["x2"]} {d["bottom2"]} '
                f'{d["x3"]} {fnum(new_top)} '
                f'{d["x4"]} {fnum(new_top)}Z')
    new_bar = bar.group("before") + new_geom + bar.group("after")
    svg = svg[:bar.start()] + new_bar + svg[bar.end():]
    lab = exactly_one(label_group_regex(tid), svg, f"panel {panel}, text_{tid} again")
    new_label = (lab.group("before") + fnum(label_y) + lab.group("middle")
                 + str(target) + lab.group("after"))
    svg = svg[:lab.start()] + new_label + svg[lab.end():]
    return svg, {
        "panel": panel, "role": role, "bar_patch_id": pid,
        "label_text_id": tid, "previous_spec": previous,
        "source_count": old, "final_count": target,
        "source_to_final_changed": True,
        "previous_to_final_changed": previous != target,
        "source_bar_height": old_height, "final_bar_height": new_height,
    }


def mask_target_groups(svg):
    """Remove the 22 bar and 22 label groups for byte-exact untouched QA."""
    masked = svg
    for panel, role, pid, tid, old, new in SPECS:
        bar = exactly_one(bar_group_regex(pid), masked, f"mask patch_{pid}")
        masked = masked[:bar.start()] + f"[EXISTING_BAR_{pid}]" + masked[bar.end():]
        label = exactly_one(label_group_regex(tid), masked, f"mask text_{tid}")
        masked = masked[:label.start()] + f"[EXISTING_LABEL_{tid}]" + masked[label.end():]
    return masked


def render_svg(svg_bytes, pdf_path, png_path, png_width, dpi):
    try:
        import cairosvg
        from PIL import Image
    except ImportError as exc:
        raise RuntimeError(
            "Rendering needs cairosvg and pillow. Run: "
            "python3 -m pip install cairosvg pillow ; "
            "or rerun with --svg-only to create the final SVG without conversions."
        ) from exc
    from io import BytesIO
    cairosvg.svg2pdf(bytestring=svg_bytes, write_to=str(pdf_path))
    pngbuf = cairosvg.svg2png(bytestring=svg_bytes, output_width=png_width)
    im = Image.open(BytesIO(pngbuf)).convert("RGB")
    im.save(png_path, format="PNG", dpi=(dpi, dpi), optimize=True)
    if im.width != png_width:
        raise AssertionError(f"Wrong PNG width: {im.width} != {png_width}")
    return [im.width, im.height]


def run(args):
    assert len(SPECS) == 22
    assert sum(p != f for _,_,_,_,p,f in SPECS) == 10
    source_name, source_bytes = identify_source(args)
    original = source_bytes.decode("utf-8")
    ET.fromstring(source_bytes)
    for check in ('id="axes_1"', 'id="axes_2"', 'id="axes_3"',
                  'id="axes_4"', 'id="axes_5"', 'id="axes_6"'):
        if check not in original:
            raise ValueError("Source SVG is not the matching four-panel Figure 3")
    source_counts = read_labels(original)
    previous_counts = [v[4] for v in SPECS]
    if args.require_previous and source_counts != previous_counts:
        raise ValueError("--require-previous was given, but the source counts do "
                         "not equal the user's supplied previous 22 values.")

    new_svg = original
    rows = []
    for item in SPECS:
        new_svg, audit = update_one(new_svg, item)
        rows.append(audit)

    # Strongest preservation assertion: every character outside the existing
    # 22 bar groups and 22 label groups must be identical to the source.
    if mask_target_groups(original) != mask_target_groups(new_svg):
        raise AssertionError("SVG markup outside bar/label groups was changed!")
    ET.fromstring(new_svg.encode("utf-8"))
    if read_labels(new_svg) != [r[5] for r in SPECS]:
        raise AssertionError("Final bar labels do not match the 22 specified values")
    for row in rows:
        pid = row["bar_patch_id"]
        m = exactly_one(bar_group_regex(pid), new_svg, f"final patch_{pid}")
        rect = RECT.fullmatch(m.group("geometry"))
        p = rect.groupdict()
        final_height = float(p["bottom1"]) - float(p["top1"])
        if abs(final_height - row["final_bar_height"]) > 1e-4:
            raise AssertionError(f"Bar height mismatch in {pid}")
    outdir = Path(args.output_dir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    stem = args.stem
    svg_path = outdir / (stem + ".svg")
    if args.input_svg and svg_path.resolve() == Path(args.input_svg).expanduser().resolve():
        raise ValueError("Output must not overwrite input SVG")
    svg_path.write_text(new_svg, encoding="utf-8")

    pdf_path = outdir / (stem + ".pdf")
    png_path = outdir / (stem + ".png")
    dims = None
    if not args.svg_only:
        dims = render_svg(new_svg.encode("utf-8"), pdf_path, png_path,
                          args.png_width, args.dpi)

    fields = list(rows[0])
    tsv = outdir / "Figure3_final_22_value_audit.tsv"
    with tsv.open("w", encoding="utf-8") as f:
        f.write("\t".join(fields) + "\n")
        for r in rows:
            f.write("\t".join(str(r[k]) for k in fields) + "\n")
    report = {
        "figure": "MACK-DB Figure 3",
        "release": "MACKDB_DKI_R_20261012_v1",
        "epi_species_and_genus": "v3",
        "projects_main_analysis": 1571,
        "strict_all3_species": 25,
        "strict_all3_genera": 19,
        "data_unit_A_B": "environment-taxon combinations",
        "data_unit_C_D": "unique taxa",
        "source_svg": source_name,
        "source_sha256": digest(source_bytes),
        "final_svg_sha256": digest(new_svg),
        "input_counts": source_counts,
        "previous_user_spec_counts": previous_counts,
        "final_counts": [r[5] for r in SPECS],
        "source_was_immediately_previous_spec": source_counts == previous_counts,
        "changed_relative_to_previous_spec": 10,
        "unchanged_relative_to_previous_spec": 12,
        "changed_relative_to_svg_source": sum(r["source_to_final_changed"] for r in rows),
        "unchanged_relative_to_svg_source": sum(not r["source_to_final_changed"] for r in rows),
        "preservation_audit": "PASS - outside original bar+label groups identical",
        "png_px": dims,
        "png_dpi": args.dpi if dims else None,
        "exported_files": [str(svg_path)] + ([] if args.svg_only else [str(pdf_path),str(png_path)]),
    }
    (outdir / "Figure3_final_source_QC.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print("SUCCESS - MACK-DB Figure 3 final DKI unified-R update")
    print("Source:", source_name)
    print("Input is the immediately preceding version:", source_counts==previous_counts)
    if source_counts != previous_counts:
        print("NOTE: supplied source is an earlier Figure 3 version. Styling is preserved;")
        print("      targets are final. Exactly 10 values differ from the USER-SPECIFIED")
        print("      previous dataset, but more may change vs this earlier SVG input.")
    print("22 final values; 10 changed and 12 unchanged vs user-specified previous version")
    for panel in "ABCD":
        print(panel + ": " + ", ".join(str(r[5]) for r in SPECS if r[0]==panel))
    print("All other SVG elements unchanged: axes, fonts, color, legend, dots/lines")
    print("SVG:", svg_path)
    if dims:
        print("PDF:", pdf_path)
        print("PNG:", png_path, f"({dims[0]}x{dims[1]}, {args.dpi} dpi)")
    print("QC:", outdir / "Figure3_final_source_QC.json")


def main():
    p=argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    g=p.add_mutually_exclusive_group()
    g.add_argument("--input-svg", type=Path, help="Path to your latest original vector Figure 3 SVG")
    g.add_argument("--docx", type=Path, help="Manuscript DOCX containing original vector Figure 3")
    p.add_argument("--output-dir", type=Path, default=Path("figure3_DKI_R_final_20261012"))
    p.add_argument("--stem", default="Figure3_FINAL_1571_EPIv3_DKI_R")
    p.add_argument("--png-width", type=int, default=4200)
    p.add_argument("--dpi", type=int, default=600)
    p.add_argument("--svg-only", action="store_true")
    p.add_argument("--require-previous", action="store_true",
                   help="Only accept source SVG with exactly the 22 previous counts specified by user")
    args=p.parse_args()
    if args.png_width < 1000 or args.dpi < 100:
        p.error("--png-width must be >=1000 and --dpi >=100")
    try:
        run(args)
    except Exception as e:
        print(f"ERROR: {e}",file=sys.stderr)
        sys.exit(1)

if __name__=='__main__':
    main()
