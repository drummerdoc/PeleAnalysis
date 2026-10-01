#!/usr/bin/env python3
"""
Plot the bin table written by amrJointFit (its outfile=).

Two panels over the same coordinate plane:

  left   the volume each bin holds -- the joint distribution of the two
         progress variables, weighted by volume rather than by cell count, so
         the refinement pattern does not tilt it. For a flame obeying a single
         progress variable this collapses onto a curve; the width of the fan
         off that curve IS the departure, and it is what the thermodiffusive
         instability produces.

  right  the conditional mean of a response over the same plane -- psi, the
         best possible function of these two coordinates. Its level curves are
         the answer to "do heat release contours follow the progress-variable
         contours": level curves running along a coordinate axis mean the
         response is controlled by that variable alone there, diagonal ones
         mean it is controlled by the pair.

  usage:  plotJointFit.py jt.dat                       first response, auto
          plotJointFit.py jt.dat --resp HeatRelease
          plotJointFit.py jt.dat --out fig.png --log-volume
          plotJointFit.py jt.dat --list

Reads only the header line for names, so it adapts to whatever coords and
responses the run used. Needs a two-coordinate table; a three-coordinate run
has to be marginalised first (or plotted a slice at a time), and the script
says so rather than guessing.
"""
import argparse
import sys

try:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, LogNorm, TwoSlopeNorm
except ImportError as e:          # noqa: BLE001
    sys.exit("needs numpy and matplotlib: %s" % e)

# One hue, light to dark, for magnitude. Two sequential panels side by side get
# different hues so they are never confused for the same scale. Lightness is
# monotone in OKLab across every ramp, which is the property that makes a
# sequential map readable and the property a rainbow map lacks.
BLUE = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
        "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281",
        "#0d366b"]
ORANGE = ["#fde3d4", "#fbc9ae", "#f7ab85", "#f18d5f", "#eb6834", "#d4541f",
          "#b04418", "#8a3513", "#62250d"]
# Diverging, for a response that changes sign: warm and cool poles that read as
# opposite, with a neutral grey midpoint so zero reads as nothing.
RED = ["#f7d6d6", "#f0a9a9", "#e87b7b", "#e34948", "#c23333", "#9c2626",
       "#761b1b"]
GREY_MID = "#f0efec"

INK = "#0b0b0b"
INK_MUTED = "#898781"
SURFACE = "#fcfcfb"
GRID = "#e1e0d9"

CMAP_VOL = LinearSegmentedColormap.from_list("vol", BLUE)
CMAP_RESP = LinearSegmentedColormap.from_list("resp", ORANGE)
CMAP_DIV = LinearSegmentedColormap.from_list(
    "div", list(reversed(BLUE)) + [GREY_MID] + RED)


def read_table(path):
    names, rows = None, []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                if names is None:
                    names = line.lstrip("#").split()
                continue
            if line.strip():
                rows.append([float(x) for x in line.split()])
    if names is None or not rows:
        sys.exit("%s: no header or no data rows" % path)
    return names, np.asarray(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--resp", default=None, help="response to plot (default: the first)")
    ap.add_argument("--out", default=None, help="output image (default: <table>.png)")
    ap.add_argument("--log-volume", action="store_true",
                    help="log colour scale for the volume panel")
    ap.add_argument("--log-resp", action="store_true",
                    help="log colour scale for the response panel")
    ap.add_argument("--panel", choices=("mean", "relsd"), default="mean",
                    help="right panel: the conditional mean (default), or the "
                         "conditional standard deviation relative to it -- the "
                         "local scatter the coordinates fail to explain, which "
                         "is where the parametrisation breaks down")
    ap.add_argument("--list", action="store_true", help="list the columns and exit")
    ap.add_argument("--mean-floor", type=float, default=1e-3,
                    help="relsd only: mask bins whose |mean| is below this "
                         "fraction of the largest, where a ratio to it is "
                         "meaningless (default 1e-3)")
    ap.add_argument("--clip", type=float, default=99.0,
                    help="relsd only: clip the colour scale at this percentile "
                         "so a few extreme bins do not flatten the rest")
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()

    names, data = read_table(args.table)
    try:
        ivol = names.index("volume")
    except ValueError:
        sys.exit("no 'volume' column; is this an amrJointFit table?")
    coords = names[:ivol]
    resps = [n[:-5] for n in names[ivol + 2:] if n.endswith("_mean")]

    if args.list:
        print("coordinates:", ", ".join(coords))
        print("responses  :", ", ".join(resps) if resps else "(none)")
        print("rows       :", len(data))
        return

    if len(coords) != 2:
        sys.exit("this table has %d coordinates (%s); the plot needs exactly 2. "
                 "Re-run amrJointFit with two coords, or marginalise first."
                 % (len(coords), ", ".join(coords)))
    if not resps:
        sys.exit("no response columns in the table")
    resp = args.resp or resps[0]
    if resp not in resps:
        sys.exit("no response %r; have %s" % (resp, ", ".join(resps)))
    imean = names.index(resp + "_mean")
    isd = names.index(resp + "_sd")

    # Occupied bins only; empty ones are absent from the file, so the grid is
    # rebuilt from the bin centres that are present and the rest left blank.
    # Blank is the honest rendering of "no data here" -- zero is not.
    xs, ys = np.unique(data[:, 0]), np.unique(data[:, 1])
    ix = np.searchsorted(xs, data[:, 0])
    iy = np.searchsorted(ys, data[:, 1])
    vol = np.full((ys.size, xs.size), np.nan)
    psi = np.full((ys.size, xs.size), np.nan)
    vol[iy, ix] = data[:, ivol]
    if args.panel == "relsd":
        # Local relative scatter. R^2 over the whole domain is flattered by all
        # the quiescent cells the coordinates trivially place; this is the same
        # question asked bin by bin, so it says WHERE the pair fails rather than
        # how much of the domain it gets right for free.
        m = np.abs(data[:, imean])
        # Relative scatter is meaningless where the mean passes through zero --
        # the burnt edge and the unburnt corner, where the response itself is
        # ~0. Those bins divide by nothing and would otherwise set the colour
        # scale for the whole panel. Mask them and say how many.
        floor = args.mean_floor * m.max() if m.size else 0.0
        keep = m > floor
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.where(keep, data[:, isd]/np.where(keep, m, 1.0), np.nan)
        psi[iy, ix] = rel
        nmask = int((~keep).sum())
        rlabel = "sd / |mean| of %s" % resp
    else:
        psi[iy, ix] = data[:, imean]
        nmask = 0
        rlabel = "mean %s" % resp

    occupied = np.isfinite(vol).sum()
    frac = data[:, ivol].sum()

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), constrained_layout=True)
    fig.patch.set_facecolor(SURFACE)

    panels = [
        (axes[0], vol, CMAP_VOL, "volume per bin", args.log_volume, False),
        (axes[1], psi, CMAP_RESP, rlabel, args.log_resp, args.panel == "mean"),
    ]
    for ax, field, cmap, label, uselog, signed in panels:
        finite = field[np.isfinite(field)]
        norm = None
        if signed and finite.size and finite.min() < 0 < finite.max():
            # Sign change present: diverging about zero, so the midpoint reads
            # as nothing and the two signs read as opposite.
            cmap = CMAP_DIV
            lim = max(abs(finite.min()), abs(finite.max()))
            norm = TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim)
        elif label.startswith("sd /") and finite.size:
            from matplotlib.colors import Normalize
            norm = Normalize(vmin=0.0, vmax=float(np.percentile(finite, args.clip)))
        elif uselog and finite.size:
            pos = finite[finite > 0]
            if pos.size:
                norm = LogNorm(vmin=pos.min(), vmax=pos.max())
        cmap = cmap.copy()
        cmap.set_bad(SURFACE)

        im = ax.pcolormesh(xs, ys, field, cmap=cmap, norm=norm,
                           shading="nearest", rasterized=True)
        cb = fig.colorbar(im, ax=ax, pad=0.02)
        cb.set_label(label, color=INK, fontsize=9)
        cb.ax.tick_params(colors=INK_MUTED, labelsize=8)
        cb.outline.set_visible(False)

        ax.set_facecolor(SURFACE)
        ax.set_xlabel(coords[0], color=INK, fontsize=10)
        ax.set_ylabel(coords[1], color=INK, fontsize=10)
        ax.tick_params(colors=INK_MUTED, labelsize=8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("bottom", "left"):
            ax.spines[s].set_color(GRID)
        ax.grid(True, color=GRID, linewidth=0.5, alpha=0.6)
        ax.set_axisbelow(True)

    extra = ("" if not nmask else
             ",  %d bin%s masked where |mean| < %g of peak"
             % (nmask, "" if nmask == 1 else "s", args.mean_floor))
    fig.suptitle("%s   (%d occupied bins, %.6g total volume%s)"
                 % (args.table, occupied, frac, extra),
                 color=INK, fontsize=11)
    out = args.out or (args.table + ".png")
    fig.savefig(out, dpi=args.dpi, facecolor=SURFACE)
    print("wrote %s  (%d occupied bins)" % (out, occupied))


if __name__ == "__main__":
    main()
