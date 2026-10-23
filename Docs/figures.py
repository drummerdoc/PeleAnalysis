#!/usr/bin/env python3
"""
Figures for Docs/ThermodiffusiveLens.tex, built from the amrJointFit tables.

Everything here reads a .dat written by amrJointFit and nothing else, so each
figure can be traced back to a table on disk exactly as the numbers in the text
can. Output is vector PDF at a size meant for \\includegraphics at
0.8\\textwidth in a one-column article.

    figures.py <dir with the .dat files> <output dir>
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# Categorical slots in fixed order, light surface. Validated: worst adjacent
# CVD dE 9.1, normal-vision dE 22.9. Every series is also directly labelled,
# which is what the contrast warning on slots 3 and 4 obligates.
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
TU, YU = 693.0, 1.325792e-2

plt.rcParams.update({
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9.5,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.6,
    "xtick.color": INK2, "ytick.color": INK2,
    "axes.labelcolor": INK, "text.color": INK,
    "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    "legend.frameon": False,
})


def style(ax, grid="y"):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid:
        ax.grid(True, axis=grid, color="#e6e5e1", linewidth=0.5, zorder=0)
        ax.set_axisbelow(True)


def load(path):
    with open(path) as fh:
        names = fh.readline()[1:].split()
    d = np.loadtxt(path, ndmin=2)
    return {n: d[:, i] for i, n in enumerate(names)}, names


def centreline(T, Y, v):
    Ts, Ys = np.unique(T), np.unique(Y)
    gi, gj = np.searchsorted(Ts, T), np.searchsorted(Ys, Y)
    V = np.zeros((len(Ts), len(Ys)))
    np.add.at(V, (gi, gj), v)
    tot = V.sum(axis=1)
    with np.errstate(invalid="ignore"):
        YL = (V * Ys).sum(axis=1) / tot
    YL[tot <= 0] = np.nan
    return Ts, Ys, gi, gj, V, YL


# ---------------------------------------------------------------- figure 1
def fig_plane(src, out):
    c, _ = load(os.path.join(src, "fit_TY_L1_fine.dat"))
    T, Y, v = c["temp"], c["Y(H2)"], c["volume"]
    q = -c["I_R(H2)_mean"]
    Ts, Ys, gi, gj, V, YL = centreline(T, Y, v)

    Q = np.full((len(Ts), len(Ys)), np.nan)
    Q[gi, gj] = np.maximum(q, 1e-30)

    fig, ax = plt.subplots(figsize=(5.0, 3.5))
    style(ax, grid=None)
    # One hue, light to dark, for a magnitude. Log, because the consumption
    # density spans five decades across the occupied band.
    m = ax.pcolormesh(Ts, Ys * 1e3, np.ma.masked_invalid(Q).T,
                      cmap="Blues", norm=LogNorm(vmin=1e1, vmax=np.nanmax(Q)),
                      shading="nearest", rasterized=True)
    cb = fig.colorbar(m, ax=ax, pad=0.02)
    cb.set_label(r"consumption density $\langle -\dot\omega_{\rm H_2}\rangle$"
                 "\n" r"[kg m$^{-3}$ s$^{-1}$]", fontsize=8)
    cb.outline.set_linewidth(0.5)

    ax.plot(Ts, YL * 1e3, color=INK, lw=1.4, ls="--", zorder=4)
    ax.annotate("band centreline " r"$Y_L(T)$", xy=(905, 7.3), fontsize=8,
                color=INK, ha="left")

    Tst, Yst = 1445.6, 1.34672e-3
    ax.plot([Tst, Ts.max()], [Yst * 1e3, Yst * 1e3], color=C[1], lw=1.8, zorder=5)
    ax.plot([Tst, Tst], [Yst * 1e3, Ys.max() * 1e3], color=C[1], lw=1.8, zorder=5)
    ax.fill_between([Tst, Ts.max()], Yst * 1e3, Ys.max() * 1e3,
                    color=C[1], alpha=0.10, zorder=1)
    ax.annotate("the lens: hotter AND richer\n"
                r"$T>T^\ast$, $Y>Y^\ast$" "\n"
                "0.87% of volume\n37% of consumption",
                xy=(1480, 7.6), fontsize=8, color=C[1], ha="left", va="top")
    ax.annotate("hot and lean: burnt out, no reaction",
                xy=(2185, 0.12), fontsize=8, color=MUTED, ha="right", va="bottom")

    ax.set_xlabel("temperature [K]")
    ax.set_ylabel(r"$Y_{\rm H_2}$  [$\times 10^{-3}$]")
    ax.set_xlim(TU, 2210); ax.set_ylim(0, 8.0)
    fig.savefig(os.path.join(out, "fig_plane.pdf"))
    plt.close(fig)


# ---------------------------------------------------------------- figure 2
def fig_capture(src, out):
    c, _ = load(os.path.join(src, "fit_TY_L1_fine.dat"))
    T, Y, v = c["temp"], c["Y(H2)"], c["volume"]
    q = -c["I_R(H2)_mean"]
    C_ = v * q
    Vtot, Ctot = v.sum(), C_.sum()
    Ts, Ys, gi, gj, Vg, YL = centreline(T, Y, v)
    Cg = np.zeros_like(Vg); np.add.at(Cg, (gi, gj), C_)

    def quad(A):
        B = np.cumsum(A[::-1, ::-1], axis=0)
        return np.cumsum(B, axis=1)[::-1, ::-1]
    CL, VL = quad(Cg), quad(Vg)
    vq, cq = VL.ravel() / Vtot, CL.ravel() / Ctot

    budgets = np.logspace(np.log10(3e-4), np.log10(0.30), 60)

    def best_quad(b):
        m = vq <= b
        return cq[m].max() if m.any() else 0.0

    def curve(score):
        o = np.argsort(-score)
        return np.cumsum(v[o]) / Vtot, np.cumsum(C_[o]) / Ctot

    def at(vv, cc, b):
        m = vv <= b
        return cc[m][-1] if m.any() else 0.0

    series = [
        ("best possible region", [at(*curve(q), b) for b in budgets], C[0], "-"),
        ("the lens", [best_quad(b) for b in budgets], C[1], "-"),
        (r"$Y-Y_L(T)$", [at(*curve(Y - YL[gi]), b) for b in budgets], C[2], "-"),
        ("temperature alone", [at(*curve(T), b) for b in budgets], C[3], "-"),
    ]

    fig, ax = plt.subplots(figsize=(5.0, 3.3))
    style(ax)
    for name, ys, col, ls in series:
        ax.plot(budgets * 100, np.array(ys) * 100, color=col, lw=1.8, ls=ls,
                label=name, zorder=3)

    ys_t = series[-1][1]
    ax.annotate("temperature alone", xy=(12, ys_t[-12] * 100 + 4.5),
                color=C[3], fontsize=8, ha="center")
    ax.axvline(0.87, color=MUTED, lw=0.7, ls=":", zorder=2)
    ax.annotate("0.87%", xy=(0.9, 4), color=MUTED, fontsize=7.5, rotation=90)
    ax.set_xscale("log")
    ax.set_xlabel("volume budget [% of domain]")
    ax.set_ylabel("fuel consumption captured [%]")
    ax.set_xlim(3e-2, 30); ax.set_ylim(0, 104)
    ax.legend(loc="upper left", handlelength=1.6, borderpad=0.2,
              labelspacing=0.35)
    fig.savefig(os.path.join(out, "fig_capture.pdf"))
    plt.close(fig)


# ---------------------------------------------------------------- figure 3
def fig_lewis(src, out):
    c, _ = load(os.path.join(src, "lewis2d_L1.dat"))
    T, Y, v = c["temp"], c["Y(H2)"], c["volume"]
    gT, gY = c["gradMag_temp_mean"], c["gradMag_Y(H2)_mean"]
    xi, tau = YU - Y, T - TU
    with np.errstate(invalid="ignore", divide="ignore"):
        Le = (gY / gT) * tau / xi
    ok = (tau > 20) & (xi > 1e-5) & (gT > 0) & np.isfinite(Le) & (Le > 0) & (Le < 5)

    edges = np.arange(700, 1700, 40.0)
    mid, med, lo, hi = [], [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        s = ok & (T >= a) & (T < b)
        if s.sum() < 5:
            continue
        w, L = v[s], Le[s]
        o = np.argsort(L); cdf = np.cumsum(w[o]) / w.sum()
        mid.append(0.5 * (a + b))
        med.append(np.average(L, weights=w))
        lo.append(L[o][np.searchsorted(cdf, .25)])
        hi.append(L[o][np.searchsorted(cdf, .75)])

    fig, ax = plt.subplots(figsize=(5.0, 3.1))
    style(ax)
    ax.fill_between(mid, lo, hi, color=C[0], alpha=0.16, lw=0, zorder=2)
    ax.plot(mid, med, color=C[0], lw=1.9, zorder=4)
    ax.annotate("flux ratio, volume-weighted (IQR shaded)", xy=(905, 0.525),
                color=C[0], fontsize=8, ha="left")

    ax.axhline(0.438, color=C[1], lw=1.6, ls="--", zorder=3)
    ax.annotate("band curvature 0.44: biased high by 45%",
                xy=(818, 0.451), color=C[1], fontsize=8, ha="left")
    ax.axhline(0.30, color=MUTED, lw=1.0, ls=":", zorder=3)
    ax.annotate("literature, lean H$_2$–air", xy=(1420, 0.262),
                color=MUTED, fontsize=8, ha="left")
    ax.axvspan(700, 800, color="#e6e5e1", alpha=0.6, lw=0, zorder=1)
    ax.annotate("preheat\nLe = 0.313", xy=(706, 0.085), color=INK2, fontsize=8)

    ax.set_xlabel("temperature [K]")
    ax.set_ylabel("Lewis number")
    ax.set_xlim(700, 1680); ax.set_ylim(0.04, 0.58)
    fig.savefig(os.path.join(out, "fig_lewis.pdf"))
    plt.close(fig)


# ---------------------------------------------------------------- figure 4
def fig_curv(src, out):
    c, _ = load(os.path.join(src, "postulate_L1.dat"))
    T, K, v, Ym = c["temp"], c["meanCurv_temp"], c["volume"], c["Y(H2)_mean"]
    want = [898, 1088, 1246, 1404]
    Ts = np.unique(T)

    fig, ax = plt.subplots(figsize=(5.0, 3.3))
    style(ax)
    for i, tw in enumerate(want):
        Tv = Ts[np.argmin(np.abs(Ts - tw))]
        m = (T == Tv) & (v > 0)
        o = np.argsort(K[m])
        x, y = K[m][o] * 1e-4, Ym[m][o] * 1e3
        ax.plot(x, y, color=C[i], lw=1.8, zorder=3)
        ax.annotate(f"{Tv:.0f} K", xy=(x[0] + 0.3, y[0] + 0.18), color=C[i],
                    fontsize=8, ha="left", va="bottom")
    ax.axvline(0, color=MUTED, lw=0.7, zorder=2)
    ax.annotate("crest: convex toward the reactants", xy=(-0.5, 1.15),
                color=INK2, fontsize=8, ha="right")
    ax.annotate("cusp", xy=(0.5, 1.15), color=INK2, fontsize=8, ha="left")
    ax.set_xlabel(r"$\kappa = \nabla\!\cdot\mathbf{n}$  "
                  r"[$\times 10^{4}$ m$^{-1}$],  $\mathbf{n}$ toward products")
    ax.set_ylabel(r"$\langle Y_{\rm H_2}\rangle$ at fixed $T$  [$\times 10^{-3}$]")
    ax.set_xlim(-11, 11); ax.set_ylim(0.6, 9.9)
    fig.savefig(os.path.join(out, "fig_curvature.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    src = sys.argv[1] if len(sys.argv) > 1 else "."
    out = sys.argv[2] if len(sys.argv) > 2 else "."
    os.makedirs(out, exist_ok=True)
    for f in (fig_plane, fig_capture, fig_lewis, fig_curv):
        f(src, out)
        print("wrote", f.__name__)
