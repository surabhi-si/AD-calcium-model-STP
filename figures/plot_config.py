"""
Master style & condition registry for the AD calcium paper figures.

Every figure script in this folder should import its condition colors and
publication rcParams from here instead of redefining them locally. That way
a color only ever gets set in one place, and every figure that references a
condition by key stays visually consistent with every other figure.

Typical usage in a notebook/script::

    from plot_config import CONDITIONS, get_condition, apply_style, DPI
    from plot_config import strip_spines, panel_labeler, label_panel
    from plot_config import set_n_ticks, add_zoom_inset, save_figure

    apply_style()
    cond = get_condition("ad_5x_ryr")
    ax.plot(t, y, color=cond.color, linestyle=cond.linestyle,
            alpha=cond.alpha, zorder=cond.zorder, label=cond.label)

To add a new experimental condition, add one entry to CONDITIONS below and
every script that looks it up by key picks up a consistent color and label.

Note on multi-panel figures built across several cells: Jupyter's inline
backend auto-displays *and closes* the active figure at the end of any cell
that touches it, even without an explicit plt.show(). If you build one
`fig` incrementally (e.g. one cell per panel), everything after the first
such cell silently vanishes from the final output. Keep `fig, axes =
plt.subplots(...)` through the final `plt.show()` in one contiguous cell —
cells that only load data or set config beforehand are fine to keep
separate, since they never touch the figure.
"""

from __future__ import annotations

import string
from dataclasses import dataclass
from typing import Iterable, Optional, Sequence

import numpy as np
from matplotlib import rcParams
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

# ============================================================
# 1. Colorblind-safe palette (Okabe & Ito, 2008), plus a couple of named
#    custom colors requested for specific conditions below. Prefer the
#    Okabe-Ito entries when adding a new condition; the custom entries are
#    a deliberate departure (a specific paper-figure color request), kept
#    visually distinct from "bluish_green" (used elsewhere for VDCC=80
#    arrowheads) rather than re-certified colorblind-safe.
# ============================================================
PALETTE = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "bluish_green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "reddish_purple": "#CC79A7",
    "brick_red": "#B22222",     # custom: requested for ad_5x_ryr
    "greenish_blue": "#0089AA",  # custom: requested for ad_half_ryr
    "gray": "#808080",          # custom: control_er_blocked -- a lighter
                                 # tint of "black" (control), so the two
                                 # read as the same background, perturbed
    "salmon": "#D89090",        # custom: ad_er_blocked -- a lighter tint
                                 # of "brick_red" (ad_5x_ryr), same logic
}


# ============================================================
# 2. Experimental conditions
# ============================================================
@dataclass
class Condition:
    label: str  # legend / display text
    color: str  # hex code or matplotlib color name
    linestyle: str = "-"
    marker: str = "o"  # for scatter/errorbar-style figures
    alpha: float = 1.0
    zorder: int = 2


CONDITIONS: dict[str, Condition] = {
    "control": Condition(
        label="Control",
        color=PALETTE["black"],
        marker="o",
        alpha=0.8,
    ),
    "control_er_blocked": Condition(
        label="Control + ER blocked",
        color=PALETTE["gray"],
        marker="v",
    ),
    "ad": Condition(
        label="AD",
        color=PALETTE["vermillion"],
        marker="o",
    ),
    "ad_5x_ryr": Condition(
        label="AD + 5× RyR",
        color=PALETTE["brick_red"],
        marker="s",
        zorder=3,  # draw the headline AD-overexpression trace on top
    ),
    "ad_3x_ryr": Condition(
        label="AD + 3× RyR",
        color=PALETTE["bluish_green"],
        marker="D",
    ),
    "ad_third_ryr": Condition(
        label=r"AD + 1/3× RyR",
        color=PALETTE["sky_blue"],
        marker="v",
        alpha=0.9,
    ),
    "ad_half_ryr": Condition(
        label="AD + 1/2× RyR",
        color=PALETTE["greenish_blue"],
        marker="^",
        alpha=0.8,
    ),
    "ad_er_blocked": Condition(
        label="AD + ER blocked",
        color=PALETTE["salmon"],
        marker="P",
    ),
    "ad_ps_blocked": Condition(
        label="PS-blocked",
        color=PALETTE["reddish_purple"],
        marker="X",
        zorder=3,  # draw on top of control where traces nearly overlap -- matches ad_5x_ryr's precedent
    ),
    "ad_half_ryr_er250": Condition(
        # RyR under-expression WITHOUT ER calcium overload (ER Ca2+ held at
        # Control's ~250uM, unlike ad_half_ryr which is at AD's ~750uM).
        label=r"Experimental AD (1/2$\times$ RyR; ER Ca$^{2+}$=250$\mu$M)",
        color=PALETTE["blue"],
        marker="^",
    ),
    "ad_half_ryr_er250_er_blocked": Condition(
        # Same as ad_half_ryr_er250 but with ER release pharmacologically blocked.
        label=r"Experimental AD (1/2$\times$ RyR; ER Ca$^{2+}$=250$\mu$M; ER blocked)",
        color=PALETTE["sky_blue"],
        marker="P",
    ),
    "dobrunz_stevens_1997": Condition(
        # Published experimental benchmark (Dobrunz & Stevens 1997, Neuron
        # 18:995-1008, Fig 3A), not a model condition.
        label="Dobrunz & Stevens (1997)",
        color=PALETTE["orange"],
        linestyle="-",
    ),
}


def get_condition(key: str) -> Condition:
    """Look up a condition by key, with a helpful error on typos."""
    try:
        return CONDITIONS[key]
    except KeyError as exc:
        raise KeyError(
            f"Unknown condition {key!r}. Available conditions: "
            f"{sorted(CONDITIONS)}. Add new conditions to CONDITIONS in "
            "plot_config.py."
        ) from exc


# ============================================================
# 3. Publication rcParams
# ============================================================
DPI = 600

STYLE = {
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 12,
    "axes.linewidth": 1.2,
    "lines.linewidth": 2.2,
    "xtick.major.width": 1.2,
    "ytick.major.width": 1.2,
    "xtick.major.size": 4,
    "ytick.major.size": 4,
    "pdf.fonttype": 42,  # keep text editable in Illustrator/Inkscape
    "ps.fonttype": 42,
    "svg.fonttype": "none",
}


def apply_style(overrides: Optional[dict] = None) -> None:
    """Apply the shared publication rcParams. Call once per notebook/script."""
    rcParams.update({**STYLE, **(overrides or {})})


# ============================================================
# 4. Small figure-building helpers shared across scripts
# ============================================================
def strip_spines(ax, spines: Sequence[str] = ("top", "right"), linewidth: Optional[float] = None) -> None:
    """Drop the given spines and apply outward-facing ticks."""
    for s in spines:
        ax.spines[s].set_visible(False)
    ax.tick_params(direction="out")
    if linewidth is not None:
        for s in ax.spines.values():
            s.set_linewidth(linewidth)


def panel_labeler(start: str = "a") -> Iterable[str]:
    """Yield 'a', 'b', 'c', ... for panel labels; start=... to resume mid-alphabet."""
    letters = string.ascii_lowercase
    i = letters.index(start)
    while True:
        yield letters[i]
        i += 1


def label_panel(ax, letter: str, x: float = -0.18, y: float = 1.05, **kwargs) -> None:
    """Draw a bold panel letter (e.g. 'a') in the top-left of an axes."""
    kwargs.setdefault("fontsize", 13)
    kwargs.setdefault("fontweight", "bold")
    ax.text(x, y, letter, transform=ax.transAxes, **kwargs)


def _format_tick(value: float) -> str:
    """3 significant figures, but never scientific notation: Python's "g"
    format (the previous approach here) switches to exponential once a
    value's magnitude reaches ~1000, which read fine for calcium/Pr-scale
    axes but produced "1.47e+03" once a panel's values ran into the
    thousands (e.g. a release-event count)."""
    if abs(value - round(value)) < 1e-6:
        return str(int(round(value)))
    if value == 0:
        return "0"
    exponent = int(np.floor(np.log10(abs(value))))
    decimals = max(0, 2 - exponent)
    text = f"{value:.{decimals}f}"
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    return text


def set_n_ticks(ax, vmin: float, vmax: float, n: int = 3, axis: str = "x") -> np.ndarray:
    """Set exactly n evenly spaced, cleanly formatted ticks spanning [vmin, vmax]."""
    ticks = np.linspace(vmin, vmax, n)
    labels = [_format_tick(t) for t in ticks]
    if axis == "x":
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
    else:
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)
    return ticks


def set_nice_ticks(ax, vmin: float, vmax: float, n: int = 3, axis: str = "x") -> np.ndarray:
    """Like set_n_ticks, but snaps to round numbers (matplotlib's own
    "nice number" step algorithm: steps of 1/2/5/10 x a power of ten)
    instead of exact fractions of [vmin, vmax] -- e.g. 0, 0.5, 1, 1.5
    rather than -0.08, 0.8, 1.69. Useful when the data range itself isn't
    already close to round numbers; set_n_ticks is the simpler default
    elsewhere.

    MaxNLocator's `nbins` is a ceiling, not a target -- asking for exactly
    n often returns fewer once ticks outside [vmin, vmax] are dropped, so
    this requests n+1 to reliably land close to n."""
    from matplotlib.ticker import MaxNLocator

    locator = MaxNLocator(nbins=n + 1, steps=[1, 2, 5, 10])
    ticks = np.array([t for t in locator.tick_values(vmin, vmax) if vmin <= t <= vmax])
    labels = [_format_tick(t) for t in ticks]
    if axis == "x":
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
    else:
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)
    return ticks


def add_zoom_inset(
    ax,
    traces: Sequence[tuple[np.ndarray, np.ndarray, Condition]],
    window: tuple[float, float],
    rect: tuple[float, float, float, float],
    *,
    ylim: Optional[tuple[float, float]] = None,
    n_xticks: int = 2,
    show_yticklabels: bool = False,
    connector: bool = True,
    loc1: int = 3,
    loc2: int = 4,
    title: Optional[str] = None,
    title_fontsize: float = 8,
    tick_fontsize: float = 8,
):
    """
    Add a small inset axes to `ax` that zooms into `window` (in x-data units).

    `traces` is a list of (x, y, Condition) tuples, one per condition, drawn
    with that condition's color/linestyle/alpha/zorder. `rect` places the
    inset in axes-fraction coordinates: (x0, y0, width, height).

    The inset shares `ax`'s y-limits by default (pass `ylim` to override),
    which keeps zoomed pulses visually comparable to the full trace and to
    each other. Only the two window boundary values are ticked on the
    x-axis; set `connector=False` to skip the dashed box linking the inset
    back to its source region on `ax`.
    """
    inset = ax.inset_axes(rect)

    for x, y, cond in traces:
        m = (x >= window[0]) & (x <= window[1])
        inset.plot(
            x[m], y[m],
            color=cond.color,
            linestyle=cond.linestyle,
            alpha=cond.alpha,
            zorder=cond.zorder,
        )

    inset.set_xlim(*window)
    inset.set_ylim(ylim if ylim is not None else ax.get_ylim())

    set_n_ticks(inset, window[0], window[1], n=n_xticks, axis="x")
    if not show_yticklabels:
        inset.set_yticklabels([])

    inset.tick_params(labelsize=tick_fontsize, length=3, width=1.0, direction="out")
    strip_spines(inset, linewidth=1.0)

    if title:
        inset.set_title(title, fontsize=title_fontsize, pad=2)

    if connector:
        mark_inset(ax, inset, loc1=loc1, loc2=loc2, fc="none", ec="0.5", lw=0.8, ls="--", zorder=1)

    return inset


def save_figure(fig, path: str, dpi: float = DPI, **kwargs) -> None:
    """Save at publication DPI with tight bounding box by default."""
    kwargs.setdefault("bbox_inches", "tight")
    fig.savefig(path, dpi=dpi, **kwargs)
