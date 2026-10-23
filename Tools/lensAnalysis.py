#!/usr/bin/env python3
"""
Post-processing for the thermodiffusive-lens analysis of lean premixed H2 flames.

Every number in Docs/ThermodiffusiveLens.tex comes from this file. It consumes
the ASCII tables written by amrJointFit (one row per occupied bin, carrying the
bin's volume, cell count, and the conditional mean and standard deviation of
each response) and nothing else, so a result can always be traced back to a
table on disk.

Columns are read from the header line, not by position:

    # <coord1> <coord2> ... volume cells <resp1>_mean <resp1>_sd ...

so the same command works whatever coords and resp were passed to amrJointFit.
Variable names may contain parentheses -- Y(H2) -- but not spaces.

Modes, each needing particular columns:

  band       coords (T, Y)                      the band centreline Y_L(T) and
                                                the local exponent of
                                                (Y_u - Y) ~ (T - T_u)^Le
  frontier   coords (T, Y), a response          the lens Pareto frontier, its
                                                efficiency against the optimal
                                                region, enrichment, and the
                                                retention of a fixed T threshold
  lewis      coords (T, Y), responses
             gradMag_<T> and gradMag_<Y>        the flux-ratio Lewis number,
                                                per bin and split by the sign of
                                                the excess-enthalpy scalar
  postulate  coords (T, curvature), response Y  enrichment at fixed temperature
                                                as a function of curvature
  r2         any coords, any responses          explained variance for every
                                                coordinate subset
  kinetics   coords (T, Y), species responses   the radical pool across the
                                                band: which side of chain
                                                branching each side of the
                                                band is on, and where it
                                                turns over

Usage:
    lensAnalysis.py band      tab.dat --Tu 693 --Yu 1.325792e-2
    lensAnalysis.py frontier  tab.dat --resp "I_R(H2)" --budgets 0.002,0.0087,0.05
    lensAnalysis.py lewis     tab.dat --Tu 693 --Yu 1.325792e-2
    lensAnalysis.py postulate tab.dat
    lensAnalysis.py r2        tab.dat
    lensAnalysis.py kinetics  tab.dat --Tu 693 --Yu 1.325792e-2 --crossover 1149

Sign conventions, which matter and are easy to get backwards:

  * meanCurv from amrGradStats is div(n) with n = grad(T)/|grad(T)|, pointing
    toward the PRODUCTS. The combustion literature and Src/curvature.cpp point
    the normal toward the reactants, so their curvature is the negative of
    this one. In the convention used here, kappa < 0 is a crest: convex toward
    the reactants.
  * I_R(H2) is negative where fuel is consumed. Consumption density is -I_R.
    Pass --negate for a response that should be flipped; it is applied
    automatically to any response whose name starts with I_R.
"""

import argparse
import sys

import numpy as np


# ----------------------------------------------------------------- table I/O

class Table:
    """An amrJointFit output table, addressed by column name."""

    def __init__(self, path):
        with open(path) as fh:
            header = fh.readline()
        if not header.startswith('#'):
            raise SystemExit(f"{path}: expected a '#' header line naming the columns")
        names = header[1:].split()
        try:
            iv = names.index('volume')
        except ValueError:
            raise SystemExit(f"{path}: header has no 'volume' column")
        if names[iv + 1] != 'cells':
            raise SystemExit(f"{path}: expected 'cells' to follow 'volume'")

        self.path = path
        self.coords = names[:iv]
        self.resp = [n[:-5] for n in names[iv + 2:] if n.endswith('_mean')]
        self.data = np.loadtxt(path, ndmin=2)
        if self.data.shape[1] != len(names):
            raise SystemExit(f"{path}: {self.data.shape[1]} columns but "
                             f"{len(names)} names in the header")
        self._col = {n: self.data[:, i] for i, n in enumerate(names)}
        self.volume = self._col['volume']
        self.cells = self._col['cells']

    def coord(self, name_or_index):
        if isinstance(name_or_index, int):
            return self._col[self.coords[name_or_index]]
        return self._col[name_or_index]

    def mean(self, resp):
        return self._col[resp + '_mean']

    def sd(self, resp):
        return self._col[resp + '_sd']

    def signed_mean(self, resp, negate=None):
        """Consumption-positive view of a response."""
        flip = resp.startswith('I_R') if negate is None else negate
        return -self.mean(resp) if flip else self.mean(resp)

    def require(self, names, mode):
        missing = [n for n in names if n not in self._col]
        if missing:
            raise SystemExit(f"{self.path}: mode '{mode}' needs column(s) "
                             f"{', '.join(missing)}; table has "
                             f"coords {self.coords} and resp {self.resp}")

    def describe(self):
        return (f"{self.path}: {len(self.data)} occupied bins, "
                f"coords {self.coords}, resp {self.resp}, "
                f"median {np.median(self.cells):.0f} cells/bin")


def grid(t, ix, iy, values):
    """Accumulate a per-bin quantity onto the (coord0, coord1) lattice."""
    xs = np.unique(t.coord(ix))
    ys = np.unique(t.coord(iy))
    gi = np.searchsorted(xs, t.coord(ix))
    gj = np.searchsorted(ys, t.coord(iy))
    out = np.zeros((len(xs), len(ys)))
    np.add.at(out, (gi, gj), values)
    return xs, ys, gi, gj, out


def guess_names(t):
    """Identify the temperature and fuel coordinates by name."""
    T = Y = None
    for c in t.coords:
        lc = c.lower()
        if T is None and ('temp' in lc or lc == 't'):
            T = c
        elif Y is None and lc.startswith('y('):
            Y = c
    return T, Y


# ----------------------------------------------------------------- band mode

def band_centreline(t, Tname, Yname):
    """Y_L(T): the volume-weighted conditional mean of Y at each T."""
    Ts, Ys, gi, gj, Vg = grid(t, Tname, Yname, t.volume)
    tot = Vg.sum(axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        YL = (Vg * Ys).sum(axis=1) / tot
    YL[tot <= 0] = np.nan
    return Ts, Ys, gi, gj, Vg, YL


def mode_band(t, args):
    Tname, Yname = args.T or guess_names(t)[0], args.Y or guess_names(t)[1]
    t.require([Tname, Yname], 'band')
    Ts, Ys, gi, gj, Vg, YL = band_centreline(t, Tname, Yname)

    Yu = args.Yu if args.Yu else np.nanmax(Ys)
    xi = Yu - YL
    tau = Ts - args.Tu
    ok = (tau > args.tau_min) & (xi > 1e-6) & np.isfinite(xi)
    if ok.sum() < 4:
        raise SystemExit("band: too few usable temperature bins; check --Tu")

    lx, ly = np.log(tau[ok]), np.log(xi[ok])
    slope = np.gradient(ly, lx)

    # Heat release, if present, marks where the no-reaction first integrals fail.
    hr = None
    for r in t.resp:
        if 'heatrelease' in r.lower():
            Hg = grid(t, Tname, Yname, t.volume * t.mean(r))[4]
            with np.errstate(invalid='ignore', divide='ignore'):
                hr = Hg.sum(axis=1) / Vg.sum(axis=1)
            hr = hr / np.nanmax(hr)
            break

    print(f"# band centreline and local exponent  (T_u = {args.Tu:g} K, "
          f"Y_u = {Yu:.6e})")
    print(f"# {'T':>8} {'T-Tu':>8} {'Y_L':>12} {'Yu-Y_L':>12} {'d lnXi/d lnTau':>15}"
          + (f" {'HR/HRmax':>9}" if hr is not None else ""))
    for k, i in enumerate(np.where(ok)[0]):
        line = (f"  {Ts[i]:>8.1f} {tau[i]:>8.1f} {YL[i]:>12.5e} "
                f"{xi[i]:>12.5e} {slope[k]:>15.4f}")
        if hr is not None:
            line += f" {hr[i]:>9.4f}"
        print(line)

    # The exponent is Le only where reaction is negligible AND we are below the
    # flame: heat release is near zero in the products too, so a bare threshold
    # on it selects both ends. Restrict to temperatures below the HR peak.
    if hr is not None:
        peak = np.nanargmax(hr)
        pre = ok & (np.arange(len(Ts)) < peak) & (hr < args.hr_max)
        label = f"HR < {args.hr_max:g} of peak and below the HR peak"
    else:
        pre = ok & (Ts < args.T_pre)
        label = f"T < {args.T_pre:g} K"
    if pre.sum() >= 3:
        sel = np.isin(np.where(ok)[0], np.where(pre)[0])
        fit = np.polyfit(lx[sel], ly[sel], 1)
        print(f"\n# preheat fit ({label}; {sel.sum()} bins, "
              f"T up to {Ts[pre].max():.0f} K)")
        print(f"#   exponent (band-curvature Lewis estimate) = {fit[0]:.4f}")
        print("#   NOTE this estimator is biased: Y_L(T) is an average over flame")
        print("#   elements of differing enrichment, so its slope is steeper than")
        print("#   the local flux ratio. Prefer the 'lewis' mode.")
    else:
        print("\n# too few preheat bins for a fit")


# ------------------------------------------------------------- frontier mode

def quadrant_cumulative(A, hot_high=True, rich_high=True):
    """Sum A over the quadrant {coord0 > t2, coord1 > t1} for every corner."""
    B = A[::-1, :] if hot_high else A
    B = np.cumsum(B, axis=0)
    B = B[::-1, :] if hot_high else B
    B = B[:, ::-1] if rich_high else B
    B = np.cumsum(B, axis=1)
    B = B[:, ::-1] if rich_high else B
    return B


def mode_frontier(t, args):
    Tname = args.T or guess_names(t)[0]
    Yname = args.Y or guess_names(t)[1]
    t.require([Tname, Yname], 'frontier')
    resp = args.resp or next((r for r in t.resp if r.startswith('I_R')), t.resp[0])
    t.require([resp + '_mean'], 'frontier')

    q = t.signed_mean(resp, args.negate)
    C = t.volume * q
    V, Ctot = t.volume.sum(), C.sum()
    Ts, Ys, gi, gj, Vg = grid(t, Tname, Yname, t.volume)
    Cg = grid(t, Tname, Yname, C)[4]

    # The lens is the HOT AND RICH quadrant. T and Y are anti-correlated in a
    # flame, so {T > t2, Y < t1} is the burnt corner and contains no reaction;
    # the lens is the off-band quadrant, where the gas is simultaneously hotter
    # and richer than the main T-Y band allows.
    CL = quadrant_cumulative(Cg)
    VL = quadrant_cumulative(Vg)
    vq, cq = VL.ravel() / V, CL.ravel() / Ctot

    # The optimal region of a given volume is a superlevel set of the response's
    # conditional mean. It is the upper bound on what ANY region shape achieves,
    # and its boundary is by definition a contour of consumption density.
    o = np.argsort(-q)
    Vo, Co = np.cumsum(t.volume[o]) / V, np.cumsum(C[o]) / Ctot

    # A one-parameter alternative: threshold the excess fuel above the band.
    YL = band_centreline(t, Tname, Yname)[5]
    dev = t.coord(Yname) - YL[gi]
    od = np.argsort(-dev)
    Vd, Cd = np.cumsum(t.volume[od]) / V, np.cumsum(C[od]) / Ctot

    def best_at(vv, cc, b):
        m = vv <= b
        return cc[m].max() if m.any() else 0.0

    def best_lens(b):
        m = vq <= b
        if not m.any():
            return 0.0, None
        j = int(np.argmax(np.where(m, cq, -np.inf)))
        i2, i1 = np.unravel_index(j, Cg.shape)
        return cq[j], (Ts[i2], Ys[i1], vq[j])

    print(f"# lens frontier for {resp} ({'consumption-positive' if q is not t.mean(resp) else 'as stored'})")
    print(f"# monotone under tightening t2: "
          f"{bool(np.all(np.diff(CL, axis=0) <= 1e-12 * abs(Ctot)))}")
    print(f"# monotone under tightening t1: "
          f"{bool(np.all(np.diff(CL, axis=1) <= 1e-12 * abs(Ctot)))}")
    print(f"# bins with negative density: {(q < 0).sum()} of {len(q)}")
    print(f"#\n# {'budget':>8} {'optimal':>9} {'lens':>9} {'Y-Y_L(T)':>10} "
          f"{'lens/opt':>9} {'t2':>9} {'t1':>12} {'enrich':>8}")
    for b in args.budgets:
        a = best_at(Vo, Co, b)
        c, th = best_lens(b)
        d = best_at(Vd, Cd, b)
        if th is None:
            continue
        print(f"  {100*b:>7.3f}% {100*a:>8.3f}% {100*c:>8.3f}% {100*d:>9.3f}% "
              f"{c/a if a > 0 else 0:>9.3f} {th[0]:>9.1f} {th[1]:>12.5e} "
              f"{c/th[2] if th[2] > 0 else 0:>7.1f}x")

    # How much of the two-parameter optimum survives freezing the T threshold.
    best = None
    for i2 in range(len(Ts)):
        if not (args.t2_lo <= Ts[i2] <= args.t2_hi):
            continue
        ratios = []
        for b in args.budgets:
            f = best_lens(b)[0]
            m = VL[i2, :] <= b * V
            ratios.append(CL[i2, :][m].max() / Ctot / f if (m.any() and f > 0) else 0.0)
        s = float(np.mean(ratios))
        if best is None or s > best[1]:
            best = (Ts[i2], s)
    if best:
        print(f"#\n# best single fixed t2 = {best[0]:.0f}, retaining "
              f"{best[1]:.3f} of the two-parameter optimum across these budgets")


# ---------------------------------------------------------------- lewis mode

def mode_lewis(t, args):
    Tname = args.T or guess_names(t)[0]
    Yname = args.Y or guess_names(t)[1]
    gT = args.gradT or f"gradMag_{Tname}"
    gY = args.gradY or f"gradMag_{Yname}"
    t.require([Tname, Yname, gT + '_mean', gY + '_mean'], 'lewis')

    T, Y = t.coord(Tname), t.coord(Yname)
    Yu = args.Yu if args.Yu else Y.max()
    xi, tau = Yu - Y, T - args.Tu

    # In the preheat zone the first integrals of the 1-D flame give
    #   rho D dY/dx = m (Y - Y_u)      lambda dT/dx = m c_p (T - T_u)
    # whose ratio is  |grad Y| / |grad T| = Le (Y_u - Y) / (T - T_u).
    # Using each bin's OWN Y -- not the T-conditional mean -- is what keeps this
    # free of the differential diffusion it is trying to measure.
    with np.errstate(invalid='ignore', divide='ignore'):
        Le = (t.mean(gY) / t.mean(gT)) * tau / xi
    ok = ((tau > args.tau_min) & (xi > 1e-5) & (t.mean(gT) > 0)
          & np.isfinite(Le) & (Le > 0) & (Le < args.Le_max))
    V = t.volume.sum()
    print(f"# flux-ratio Lewis number  (T_u = {args.Tu:g} K, Y_u = {Yu:.6e})")
    print(f"# usable bins {ok.sum()} of {len(t.data)}, "
          f"{100*t.volume[ok].sum()/V:.2f}% of the volume in range")
    print(f"#\n# {'T band':>15} {'vol':>9} {'Le vol-wt':>10} "
          f"{'p25':>8} {'p50':>8} {'p75':>8}")
    edges = args.bands
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = ok & (T >= lo) & (T < hi)
        if s.sum() < 3:
            continue
        w, L = t.volume[s], Le[s]
        o = np.argsort(L)
        cdf = np.cumsum(w[o]) / w.sum()
        p = [L[o][min(np.searchsorted(cdf, qq), len(L) - 1)] for qq in (.25, .5, .75)]
        print(f"  {lo:>6.0f}-{hi:<8.0f} {100*w.sum()/V:>8.4f}% "
              f"{np.average(L, weights=w):>10.4f} {p[0]:>8.4f} {p[1]:>8.4f} {p[2]:>8.4f}")

    # Apparent Le splits by the sign of the excess-enthalpy scalar. This is not
    # a transport coefficient changing -- the mixture is the same -- it is the
    # first integrals failing by exactly the differential-diffusion source, so
    # the split measures the imbalance with no reference to reaction rates.
    YL = band_centreline(t, Tname, Yname)[5]
    gi = np.searchsorted(np.unique(T), T)
    dev = Y - YL[gi]
    pre = ok & (T >= args.pre_lo) & (T < args.pre_hi)
    print(f"#\n# apparent Le by sign of S, over {args.pre_lo:g}-{args.pre_hi:g} K")
    for lab, s in ((" S > 0  (rich side, the lens)", pre & (dev > 0)),
                   (" S <= 0 (lean side, the cusps)", pre & (dev <= 0))):
        if s.sum() < 3:
            continue
        w = t.volume[s]
        print(f"# {lab}: vol {100*w.sum()/V:7.4f}%  "
              f"Le_apparent = {np.average(Le[s], weights=w):.4f}")


# ------------------------------------------------------------ postulate mode

def mode_postulate(t, args):
    Tname = args.T or guess_names(t)[0]
    Kname = args.K or next((c for c in t.coords if 'curv' in c.lower()), None)
    if Kname is None:
        raise SystemExit("postulate: no curvature coordinate found; pass --K")
    Yname = args.Y or next((r for r in t.resp if r.lower().startswith('y(')), None)
    if Yname is None:
        raise SystemExit("postulate: no Y(..) response found; pass --Y")
    t.require([Tname, Kname, Yname + '_mean'], 'postulate')

    T, K = t.coord(Tname), t.coord(Kname)
    Ym = t.mean(Yname)
    V = t.volume.sum()
    cons = next((r for r in t.resp if r.startswith('I_R')), None)

    # The prediction: at fixed T, the mean fuel mass fraction rises as the
    # curvature becomes more negative, because kappa < 0 here means convex
    # toward the reactants -- the crest, where a Le < 1 fuel out-diffuses the
    # heat that would otherwise balance it.
    print(f"# enrichment at fixed temperature versus curvature")
    print(f"# kappa = div(n), n = grad(T)/|grad(T)| toward products; "
          f"kappa < 0 is convex toward the reactants")
    bands = [(-np.inf, -3e4), (-3e4, -1e4), (-1e4, 1e4), (1e4, 3e4), (3e4, np.inf)]
    hdr = f"# {'T':>8} {'vol':>9}" + "".join(f"{f'{lo:.0e}':>11}" for lo, _ in bands)
    print(hdr + f" {'d<Y>/dk':>12}" + (f" {'cons ratio':>11}" if cons else ""))
    for Tv in np.unique(T):
        m = T == Tv
        if t.volume[m].sum() / V < args.vol_min or m.sum() < 5:
            continue
        line = f"  {Tv:>8.1f} {100*t.volume[m].sum()/V:>8.4f}%"
        for lo, hi in bands:
            s = m & (K >= lo) & (K < hi)
            if t.volume[s].sum() <= 0:
                line += f"{'-':>11}"
            else:
                line += f"{(t.volume[s]*Ym[s]).sum()/t.volume[s].sum():>11.4e}"
        slope = np.polyfit(K[m], Ym[m], 1, w=np.sqrt(t.volume[m]))[0]
        line += f" {slope:>12.4e}"
        if cons:
            a, b = m & (K < -1e4), m & (K > 1e4)
            if t.volume[a].sum() > 0 and t.volume[b].sum() > 0:
                qa = (t.volume[a] * -t.mean(cons)[a]).sum() / t.volume[a].sum()
                qb = (t.volume[b] * -t.mean(cons)[b]).sum() / t.volume[b].sum()
                line += f" {qa/qb if qb != 0 else np.nan:>11.2f}"
        print(line)

    # Partial slope of Y on kappa at fixed T, volume weighted.
    s = (T >= args.pre_lo) & (T < args.pre_hi) & (t.volume > 0)
    if s.sum() > 5:
        A = np.c_[np.ones(s.sum()), T[s], K[s]]
        w = np.sqrt(t.volume[s])
        coef = np.linalg.lstsq(A * w[:, None], Ym[s] * w, rcond=None)[0]
        print(f"#\n# over {args.pre_lo:g}-{args.pre_hi:g} K, "
              f"d<Y>/dkappa at fixed T = {coef[2]:.4e}")


# ----------------------------------------------------------------- r2 mode

def mode_r2(t, args):
    V = t.volume.sum()
    nc = len(t.coords)
    keys = []
    for c in t.coords:
        u = np.unique(t.coord(c))
        keys.append(np.searchsorted(u, t.coord(c)))
    keys = np.array(keys)

    for resp in (args.resp_list or t.resp):
        mm, ss = t.mean(resp), t.sd(resp)
        mu = (t.volume * mm).sum() / V
        tot = (t.volume * (ss**2 + (mm - mu)**2)).sum() / V
        print(f"\n# response {resp}   mean {mu:.6g}   sd {np.sqrt(tot):.6g}")
        print(f"# {'coords':>40} {'nOcc':>8} {'R^2':>10}")
        for mask in range(1, 1 << nc):
            use = [d for d in range(nc) if mask & (1 << d)]
            key = np.zeros(len(t.data), dtype=np.int64)
            for d in use:
                key = key * (keys[d].max() + 1) + keys[d]
            within = 0.0
            uniq, inv = np.unique(key, return_inverse=True)
            for b in range(len(uniq)):
                s = inv == b
                Vb = t.volume[s].sum()
                if Vb <= 0:
                    continue
                mb = (t.volume[s] * mm[s]).sum() / Vb
                within += (t.volume[s] * (ss[s]**2 + (mm[s] - mb)**2)).sum() / V
            name = ','.join(t.coords[d] for d in use)
            print(f"  {name:>40} {len(uniq):>8} {1.0 - within/tot:>10.6f}")


# ------------------------------------------------------------ kinetics mode

# The radical pool, split by which chain a species belongs to. H, O and OH
# carry chain branching; HO2 and H2O2 are what the competing termination
# channel makes. Overridable, because a mechanism may name them differently
# or carry only some of them.
BRANCH_DEFAULT = ['Y(H)', 'Y(O)', 'Y(OH)']
TERM_DEFAULT   = ['Y(HO2)', 'Y(H2O2)']


def pool_sum(t, names, sel):
    """Volume-weighted total of several species over a selection of bins."""
    w = t.volume[sel]
    if w.sum() <= 0:
        return np.nan
    return sum(np.average(t.mean(n)[sel], weights=w) for n in names)


def mode_kinetics(t, args):
    """Which chain the lens and the cusps are on, from the radical pool alone.

    In hydrogen the structure of a flame is set by a competition between
    chain branching, H + O2 -> OH + O, and termination, H + O2 + M -> HO2 + M.
    Which wins is a function of temperature and pressure, and the temperature
    where they balance -- crossover -- is near 1150 K for this mixture at
    8 atm, rising as water accumulates and raises the third-body efficiency.

    Rather than evaluate either rate, this reads the answer off the pool the
    two channels leave behind:

        Phi = (Y(HO2) + Y(H2O2)) / (Y(H) + Y(O) + Y(OH))

    large below crossover and small above it. Nothing here needs a rate
    constant, a mechanism file or a pressure, which is the point: the
    conclusion does not inherit the uncertainty of rate parameters it would
    otherwise be argued from.

    Two caveats, because on its own this is suggestive rather than decisive.
    Phi = 1 is a convention, not a rate equality -- HO2 is consumed as well
    as produced, so the absolute level of Phi carries mechanism detail. And a
    ratio of mass fractions is not a ratio of fluxes.

    What is subject to neither caveat is the comparison BETWEEN THE TWO SIDES
    of the T-Y band at the same temperature. Both sides are the same mixture
    at the same T and the same pressure, so a difference in Phi across the
    band is differential diffusion sorting the radical pool, measured with no
    reference to kinetics at all. That is the same logic as the apparent
    Lewis number split in the 'lewis' mode, and the two are read together.
    """
    Tname = args.T or guess_names(t)[0]
    Yname = args.Y or guess_names(t)[1]
    t.require([Tname, Yname], 'kinetics')

    have = set(t.resp)
    branch = [n for n in (args.branch or BRANCH_DEFAULT) if n in have]
    term   = [n for n in (args.term   or TERM_DEFAULT)   if n in have]
    if not branch or not term:
        raise SystemExit(
            "kinetics: need at least one branching and one termination "
            f"species among the responses.\n  looked for "
            f"{args.branch or BRANCH_DEFAULT} and {args.term or TERM_DEFAULT}\n"
            f"  table has {sorted(have)}\n"
            "  pass --branch/--term to name them differently")

    T = t.coord(Tname)
    Y = t.coord(Yname)
    V = t.volume.sum()

    # Side of the band, by the sign of the excess-enthalpy scalar: the same
    # construction the lewis mode splits on, so the two tables line up band
    # for band and can be read side by side.
    YL = band_centreline(t, Tname, Yname)[5]
    gi = np.searchsorted(np.unique(T), T)
    dev = Y - YL[gi]

    print("# radical pool across the T-Y band")
    print(f"#   branching  : {' + '.join(branch)}")
    print(f"#   termination: {' + '.join(term)}")
    print("#   Phi = termination / branching; Phi > 1 is the terminating side"
          " of crossover")
    if args.crossover:
        print(f"#   reference crossover, for comparison only: "
              f"{args.crossover:g} K")
    print(f"#\n# {'T band':>15} {'vol':>9} {'Phi all':>10} "
          f"{'Phi S>0':>10} {'Phi S<=0':>10} {'ratio':>8}")

    edges = args.bands
    prof = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        base = (T >= lo) & (T < hi)
        if t.volume[base].sum() <= 0:
            continue
        pa = pool_sum(t, term, base) / pool_sum(t, branch, base)
        rich = base & (dev > 0)
        lean = base & (dev <= 0)
        pr = pl = np.nan
        if t.volume[rich].sum() > 0:
            pr = pool_sum(t, term, rich) / pool_sum(t, branch, rich)
        if t.volume[lean].sum() > 0:
            pl = pool_sum(t, term, lean) / pool_sum(t, branch, lean)
        ratio = (pr/pl if (np.isfinite(pr) and np.isfinite(pl) and pl > 0)
                 else np.nan)
        print(f"  {lo:>6.0f}-{hi:<8.0f} {100*t.volume[base].sum()/V:>8.4f}% "
              f"{pa:>10.4f} {pr:>10.4f} {pl:>10.4f} {ratio:>8.3f}")
        prof.append((0.5*(lo+hi), pa, pr, pl))

    # Where Phi crosses one, by interpolation between the two bands that
    # straddle it. Phi falls over orders of magnitude, so the interpolation
    # is in log Phi; a linear one places the crossing systematically low.
    def crossing(idx, label):
        pts = [(p[0], p[idx]) for p in prof
               if np.isfinite(p[idx]) and p[idx] > 0]
        for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
            if (y0 - 1.0)*(y1 - 1.0) <= 0 and y0 != y1:
                f = (0.0 - np.log(y0))/(np.log(y1) - np.log(y0))
                print(f"#   {label:<28} Phi = 1 at {x0 + f*(x1 - x0):>7.1f} K")
                return
        print(f"#   {label:<28} no crossing within the bands given")

    print("#")
    crossing(1, "whole band")
    crossing(2, "rich side (the lens)")
    crossing(3, "lean side (the cusps)")

    # Net H production: an independent reading of the same competition that
    # does not pass through concentrations at all. Where branching wins the
    # radical pool is growing, so the conditional mean of I_R(H) is positive.
    hsrc = args.hrate if args.hrate else 'I_R(H)'
    if hsrc in have:
        print(f"#\n# net {hsrc} by temperature band, and where it turns")
        print(f"# {'T band':>15} {'all':>14} {'S>0':>14} {'S<=0':>14}")
        last = None
        for lo, hi in zip(edges[:-1], edges[1:]):
            base = (T >= lo) & (T < hi)
            if t.volume[base].sum() <= 0:
                continue
            vals = []
            for s in (base, base & (dev > 0), base & (dev <= 0)):
                vals.append(np.average(t.mean(hsrc)[s], weights=t.volume[s])
                            if t.volume[s].sum() > 0 else np.nan)
            print(f"  {lo:>6.0f}-{hi:<8.0f} {vals[0]:>14.5e} "
                  f"{vals[1]:>14.5e} {vals[2]:>14.5e}")
            if last is not None and np.isfinite(vals[0]) and last[1]*vals[0] < 0:
                print(f"#   net {hsrc} changes sign between {last[0]:.0f} and "
                      f"{0.5*(lo + hi):.0f} K")
            last = (0.5*(lo + hi), vals[0])
    else:
        print(f"#\n# no {hsrc} among the responses; carry it for the "
              "independent reading")


# --------------------------------------------------------------------- main

def main(argv=None):
    p = argparse.ArgumentParser(
        description="Thermodiffusive-lens analysis of amrJointFit tables.",
        epilog="See Docs/ThermodiffusiveLens.tex for what each mode measures.")
    p.add_argument('mode', choices=['band', 'frontier', 'lewis', 'postulate',
                                    'r2', 'kinetics'])
    p.add_argument('table')
    p.add_argument('--Tu', type=float, default=693.0,
                   help="unburnt temperature (K); the Lewis estimates are very "
                        "sensitive to this, so use the prescribed inflow value")
    p.add_argument('--Yu', type=float, default=None,
                   help="unburnt fuel mass fraction; default is the table maximum")
    p.add_argument('--T', default=None, help="name of the temperature coordinate")
    p.add_argument('--Y', default=None, help="name of the fuel coordinate")
    p.add_argument('--K', default=None, help="name of the curvature coordinate")
    p.add_argument('--gradT', default=None)
    p.add_argument('--gradY', default=None)
    p.add_argument('--resp', default=None, help="response for 'frontier'")
    p.add_argument('--resp-list', default=None, type=lambda s: s.split(','),
                   help="comma-separated responses for 'r2'")
    p.add_argument('--negate', action='store_true', default=None,
                   help="flip the sign of the response (automatic for I_R*)")
    p.add_argument('--budgets', type=lambda s: [float(x) for x in s.split(',')],
                   default=[0.002, 0.005, 0.0087, 0.02, 0.05],
                   help="volume budgets, as fractions")
    p.add_argument('--bands', type=lambda s: [float(x) for x in s.split(',')],
                   default=[713, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1750],
                   help="temperature band edges for 'lewis'")
    p.add_argument('--tau-min', type=float, default=20.0)
    p.add_argument('--hr-max', type=float, default=0.02)
    p.add_argument('--T-pre', type=float, default=900.0)
    p.add_argument('--Le-max', type=float, default=5.0)
    p.add_argument('--pre-lo', type=float, default=713.0)
    p.add_argument('--pre-hi', type=float, default=1000.0)
    p.add_argument('--t2-lo', type=float, default=1150.0)
    p.add_argument('--t2-hi', type=float, default=1750.0)
    p.add_argument('--vol-min', type=float, default=2e-5)
    p.add_argument('--branch', type=lambda s: s.split(','), default=None,
                   help="comma-separated branching-pool species for 'kinetics' "
                        f"(default {','.join(BRANCH_DEFAULT)})")
    p.add_argument('--term', type=lambda s: s.split(','), default=None,
                   help="comma-separated termination-pool species for 'kinetics' "
                        f"(default {','.join(TERM_DEFAULT)})")
    p.add_argument('--hrate', default=None,
                   help="net H production rate response for 'kinetics' "
                        "(default I_R(H))")
    p.add_argument('--crossover', type=float, default=None,
                   help="a reference crossover temperature to print alongside "
                        "the measured one; printed, never used in a calculation")
    args = p.parse_args(argv)

    t = Table(args.table)
    print(f"# {t.describe()}")
    {'band': mode_band, 'frontier': mode_frontier, 'lewis': mode_lewis,
     'postulate': mode_postulate, 'r2': mode_r2,
     'kinetics': mode_kinetics}[args.mode](t, args)
    return 0


if __name__ == '__main__':
    sys.exit(main())
