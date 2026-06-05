"""
Shared matplotlib style for the Chaos (AIP) manuscript figures.

Single source of truth for fonts, sizes, line weights, and a colorblind-safe
palette (Okabe & Ito / Wong, Nat. Methods 8, 441 (2011)). Every paper figure is
built through `apply_style()` + `PALETTE`, so the whole figure set can be
regenerated with one consistent look. The palette is mirrored in
figures.config.json (the committed source); this module reads it from there so
the two never drift.

Deterministic: no randomness, no environment-dependent fonts (DejaVu Sans ships
with matplotlib), no interactive backend.
"""
from __future__ import annotations

import json
import os

import matplotlib

matplotlib.use("Agg")  # headless, deterministic raster/vector backend
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))


def load_config(path=None):
    path = path or os.path.join(HERE, "figures.config.json")
    with open(path) as f:
        return json.load(f)


_CFG = load_config()
PALETTE = dict(_CFG["style"]["palette"])
DPI = int(_CFG["dpi"])

# Semantic colour roles, drawn from the one palette, so both figures speak the
# same visual language.
ROLES = {
    "black": PALETTE["black"],
    "grey": PALETTE["grey"],
    "primary": PALETTE["blue"],        # main signal min(R1,R2)
    "faint": PALETTE["grey"],          # backdrop R1/R2
    "graze_thresh": PALETTE["vermillion"],   # theta = 0.85
    "recover_thresh": PALETTE["orange"],     # recovery = 0.80
    "graze_event": PALETTE["orange"],
    "absorption": PALETTE["purple"],
    "supervisor": PALETTE["grey"],
    "sn": PALETTE["blue"],             # saddle-node
    "hopf": PALETTE["vermillion"],     # Hopf
    "homoclinic": PALETTE["green"],    # homoclinic
    "tb": PALETTE["black"],            # Takens-Bogdanov
    "corner_stable": PALETTE["green"],
    "corner_transient": PALETTE["vermillion"],
}


def apply_style():
    """Apply the shared rcParams. Idempotent; call at the top of each figure."""
    base = int(_CFG["style"]["base_fontsize"])
    plt.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": DPI,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03,
        "pdf.fonttype": 42,          # embed TrueType (editable text in the PDF)
        "ps.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": [_CFG["style"]["font_family"], "DejaVu Sans"],
        "font.size": base,
        "axes.titlesize": base + 1,
        "axes.labelsize": base,
        "xtick.labelsize": base - 1,
        "ytick.labelsize": base - 1,
        "legend.fontsize": base - 1,
        "axes.linewidth": 0.8,
        "lines.linewidth": 1.4,
        "lines.markersize": 5,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "legend.frameon": False,
        "figure.autolayout": False,
    })


def savefig(fig, stem):
    """Save `fig` to <output_dir>/<stem>.{pdf,png} at the configured dpi.
    Returns the list of written paths (absolute)."""
    out_dir = os.path.join(ROOT, _CFG["output_dir"])
    os.makedirs(out_dir, exist_ok=True)
    paths = []
    for ext in _CFG["formats"]:
        p = os.path.join(out_dir, f"{stem}.{ext}")
        fig.savefig(p)
        paths.append(p)
    return paths
