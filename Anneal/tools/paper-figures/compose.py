#!/usr/bin/env python3
"""
Side-by-side figure composer for the two-panel paper figures.

Three paper figures are two-panel composites whose panels are produced by
different analysis pipelines, so they cannot be emitted by a single plotting
script. This utility assembles them deterministically from the committed source
figures, preserving vector quality for the .pdf and compositing the .png:

  fig3 = campaign_results/km_survival_primary  +  absorption_results/k_abs_vs_N
  fig5 = phase_results/cp2_rose                +  absorption_results/absorption_phase_rose
  fig8 = reduced_results/cp2_portrait_A0.2     +  reduced_results/cp2_portrait_A0.5

The shorter panel is vertically centered. Run via tools/paper-figures/run_all.py,
or standalone:

  python3 tools/paper-figures/compose.py LEFT_STEM RIGHT_STEM DEST_STEM
    (stems are paths relative to the repo root, without extension)

Dependencies: pypdf (vector PDF merge) and Pillow (PNG). See requirements.txt.
"""
from __future__ import annotations

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def compose_pdf(left: str, right: str, dest: str, gap: float = 10.0) -> None:
    from pypdf import PdfReader, PdfWriter, Transformation, PageObject

    a = PdfReader(left).pages[0]
    b = PdfReader(right).pages[0]
    aw, ah = float(a.mediabox.width), float(a.mediabox.height)
    bw, bh = float(b.mediabox.width), float(b.mediabox.height)
    H = max(ah, bh)
    W = aw + bw + gap
    page = PageObject.create_blank_page(width=W, height=H)
    page.merge_transformed_page(a, Transformation().translate(0, (H - ah) / 2))
    page.merge_transformed_page(b, Transformation().translate(aw + gap, (H - bh) / 2))
    w = PdfWriter()
    w.add_page(page)
    with open(dest, "wb") as f:
        w.write(f)


def compose_png(left: str, right: str, dest: str, gap: int = 16) -> None:
    from PIL import Image

    a = Image.open(left).convert("RGBA")
    b = Image.open(right).convert("RGBA")
    H = max(a.height, b.height)
    W = a.width + b.width + gap
    canvas = Image.new("RGBA", (W, H), (255, 255, 255, 255))
    canvas.paste(a, (0, (H - a.height) // 2), a)
    canvas.paste(b, (a.width + gap, (H - b.height) // 2), b)
    canvas.convert("RGB").save(dest, dpi=(a.info.get("dpi", (150, 150))))


def compose(left_stem: str, right_stem: str, dest_stem: str) -> None:
    """left_stem/right_stem/dest_stem are repo-root-relative, no extension."""
    L = os.path.join(ROOT, left_stem)
    R = os.path.join(ROOT, right_stem)
    D = os.path.join(ROOT, dest_stem)
    compose_pdf(L + ".pdf", R + ".pdf", D + ".pdf")
    if os.path.exists(L + ".png") and os.path.exists(R + ".png"):
        compose_png(L + ".png", R + ".png", D + ".png")
    print(f"  {os.path.basename(dest_stem)} <- {left_stem} + {right_stem}")


# Paper two-panel figures: dest <- (left source, right source), repo-relative stems.
PAIR_MAP = {
    "paper_figures/fig3": (
        "campaign_results/km_survival_primary",
        "absorption_results/k_abs_vs_N",
    ),
    "paper_figures/fig5": (
        "phase_results/cp2_rose",
        "absorption_results/absorption_phase_rose",
    ),
    "paper_figures/fig8": (
        "reduced_results/cp2_portrait_A0.2",
        "reduced_results/cp2_portrait_A0.5",
    ),
}


def compose_all() -> None:
    print("\n# composing two-panel paper figures (fig3, fig5, fig8)")
    for dest, (left, right) in PAIR_MAP.items():
        compose(left, right, dest)


if __name__ == "__main__":
    if len(sys.argv) == 4:
        compose(sys.argv[1], sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 1:
        compose_all()
    else:
        sys.exit("usage: compose.py [LEFT_STEM RIGHT_STEM DEST_STEM]")
