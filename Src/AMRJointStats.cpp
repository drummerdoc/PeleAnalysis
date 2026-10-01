#include "AMRJointStats.H"

#include <AMReX_ParallelDescriptor.H>
#include <algorithm>
#include <cmath>

using namespace amrex;

namespace amrcomplex {

AMRJointStats::AMRJointStats (const AMRComplex& cx, const Params& p)
    : m_cx(cx), m_p(p)
{
    m_nlev = cx.nLevels();
    const int nc = m_p.nCoord;
    const int nr = m_p.nResp;

    AMREX_ALWAYS_ASSERT(nc >= 1 && nr >= 0);
    AMREX_ALWAYS_ASSERT(nc <= 8);                 // mask is a bitfield
    AMREX_ALWAYS_ASSERT(cx.nComp() == nc + nr);
    AMREX_ALWAYS_ASSERT(int(m_p.clo.size())   == nc);
    AMREX_ALWAYS_ASSERT(int(m_p.chi.size())   == nc);
    AMREX_ALWAYS_ASSERT(int(m_p.cbins.size()) == nc);

    m_nbinTot = 1;
    for (int d = 0; d < nc; ++d) {
        AMREX_ALWAYS_ASSERT(m_p.cbins[d] >= 1);
        AMREX_ALWAYS_ASSERT(m_p.chi[d] > m_p.clo[d]);
        m_nbinTot *= long(m_p.cbins[d]);
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_nbinTot <= 20000000L,
        "AMRJointStats: too many bins; the table is dense, so the product of "
        "cbins has to stay within memory");

    m_cellVol.resize(m_nlev);
    for (int lev = 0; lev < m_nlev; ++lev) {
        const auto dx = cx.geom(lev).CellSize();
        Real v = 1.0;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { v *= dx[d]; }
        m_cellVol[lev] = v;
    }

    m_cnt.resize(m_nbinTot*long(m_nlev), 0L);
    m_s1.resize(m_nbinTot*long(std::max(nr,1)), Real(0.0));
    m_ss.resize(m_nbinTot*long(std::max(nr,1)), Real(0.0));

    // ---- pass one: counts, and the first moment of each response ----------
    //
    // Compensated per bin. A bin on a 300M-cell hierarchy can hold a million
    // cells, where a naive sum carries about 1e-10 of relative error -- enough
    // to show in the tenth digit of an R^2 that matters. Kahan takes that to a
    // couple of ulp whatever order the cells arrive in, which also means the
    // only remaining rank dependence is the cross-rank reduction itself.
    Vector<Real> comp(m_s1.size(), Real(0.0));
    Vector<Long> outCnt(m_nlev, 0L);
    cx.forEachNode([&] (const Node& n, const Real* v, Real cvol, bool /*shared*/)
    {
        const long b = binOf(v);
        if (b < 0) { ++outCnt[n.level]; return; }
        ++m_cnt[b*long(m_nlev) + long(n.level)];
        for (int r = 0; r < nr; ++r) {
            const long  i = b*long(nr)+long(r);
            const Real  y = cvol*v[nc+r] - comp[i];
            const Real  t = m_s1[i] + y;
            comp[i] = (t - m_s1[i]) - y;
            m_s1[i] = t;
        }
    });
    if (!m_cnt.empty()) {
        ParallelDescriptor::ReduceLongSum(m_cnt.dataPtr(), int(m_cnt.size()));
    }
    if (nr > 0) {
        ParallelDescriptor::ReduceRealSum(m_s1.dataPtr(), int(m_s1.size()));
    }
    ParallelDescriptor::ReduceLongSum(outCnt.dataPtr(), int(outCnt.size()));

    m_volIn = 0;
    for (long b = 0; b < m_nbinTot; ++b) { m_volIn += volumeOfBin(b); }
    m_volOut = 0;
    for (int lev = 0; lev < m_nlev; ++lev) {
        m_volOut += Real(outCnt[lev])*m_cellVol[lev];
    }

    // ---- pass two: squared deviation from THIS bin's mean -----------------
    //
    // Not S2 - S1*S1/V. That form subtracts two nearly equal numbers whenever
    // the response barely varies inside a bin, which is exactly the regime the
    // second-order convergence gate lives in.
    if (nr > 0) {
        Vector<Real> psi(m_nbinTot*long(nr), Real(0.0));
        for (long b = 0; b < m_nbinTot; ++b) {
            const Real vb = volumeOfBin(b);
            if (vb <= Real(0.0)) { continue; }
            for (int r = 0; r < nr; ++r) {
                psi[b*long(nr)+long(r)] = m_s1[b*long(nr)+long(r)]/vb;
            }
        }
        comp.assign(m_ss.size(), Real(0.0));
        cx.forEachNode([&] (const Node& /*n*/, const Real* v, Real cvol, bool)
        {
            const long b = binOf(v);
            if (b < 0) { return; }
            for (int r = 0; r < nr; ++r) {
                const long i = b*long(nr)+long(r);
                const Real d = v[nc+r] - psi[i];
                const Real y = cvol*d*d - comp[i];
                const Real t = m_ss[i] + y;
                comp[i] = (t - m_ss[i]) - y;
                m_ss[i] = t;
            }
        });
        ParallelDescriptor::ReduceRealSum(m_ss.dataPtr(), int(m_ss.size()));
    }

    if (m_p.verbose && ParallelDescriptor::IOProcessor()) {
        Print() << "AMRJointStats: " << m_nbinTot << " bins, volume in range "
                << m_volIn << ", dropped " << m_volOut << "\n";
    }
}

long AMRJointStats::binOf (const Real* v) const
{
    long idx = 0, stride = 1;
    for (int d = 0; d < m_p.nCoord; ++d) {
        const Real t = (v[d] - m_p.clo[d])/(m_p.chi[d] - m_p.clo[d]);
        if (t < Real(0.0) || t > Real(1.0)) { return -1; }
        int b = int(t*Real(m_p.cbins[d]));
        if (b >= m_p.cbins[d]) { b = m_p.cbins[d]-1; }   // closed upper edge
        if (b < 0)             { b = 0; }
        idx    += long(b)*stride;
        stride *= long(m_p.cbins[d]);
    }
    return idx;
}

void AMRJointStats::binIndex (long b, int* out) const
{
    for (int d = 0; d < m_p.nCoord; ++d) {
        out[d] = int(b % long(m_p.cbins[d]));
        b     /= long(m_p.cbins[d]);
    }
}

Real AMRJointStats::binCentre (long b, int d) const
{
    int idx[8];
    binIndex(b, idx);
    const Real w = (m_p.chi[d]-m_p.clo[d])/Real(m_p.cbins[d]);
    return m_p.clo[d] + (Real(idx[d])+Real(0.5))*w;
}

Real AMRJointStats::volumeOfBin (long b) const
{
    Real v = 0;
    for (int lev = 0; lev < m_nlev; ++lev) {
        v += Real(m_cnt[b*long(m_nlev)+long(lev)])*m_cellVol[lev];
    }
    return v;
}

Long AMRJointStats::cellsInBin (long b) const
{
    Long n = 0;
    for (int lev = 0; lev < m_nlev; ++lev) { n += m_cnt[b*long(m_nlev)+long(lev)]; }
    return n;
}

Real AMRJointStats::conditionalMean (long b, int r) const
{
    const Real vb = volumeOfBin(b);
    return (vb > Real(0.0)) ? m_s1[b*long(m_p.nResp)+long(r)]/vb : Real(0.0);
}

Real AMRJointStats::conditionalVar (long b, int r) const
{
    const Real vb = volumeOfBin(b);
    return (vb > Real(0.0)) ? m_ss[b*long(m_p.nResp)+long(r)]/vb : Real(0.0);
}

AMRJointStats::Fit AMRJointStats::fit (int r, unsigned mask, int coarsen) const
{
    AMREX_ALWAYS_ASSERT(r >= 0 && r < m_p.nResp);
    AMREX_ALWAYS_ASSERT(coarsen >= 1);

    const int nc = m_p.nCoord;
    const int nr = m_p.nResp;

    // Aggregated bin counts: a retained coordinate keeps ceil(nbins/coarsen)
    // groups, a dropped one collapses to a single group.
    int ab[8], cb[8];
    long nagg = 1;
    for (int d = 0; d < nc; ++d) {
        cb[d] = m_p.cbins[d];
        ab[d] = (mask & (1u<<unsigned(d)))
              ? (cb[d] + coarsen - 1)/coarsen
              : 1;
        nagg *= long(ab[d]);
    }

    Vector<Real> V(nagg, Real(0.0)), S(nagg, Real(0.0));

    int idx[8];
    auto aggOf = [&] (long b) -> long
    {
        binIndex(b, idx);
        long a = 0, stride = 1;
        for (int d = 0; d < nc; ++d) {
            const int g = (mask & (1u<<unsigned(d))) ? idx[d]/coarsen : 0;
            a      += long(g)*stride;
            stride *= long(ab[d]);
        }
        return a;
    };

    for (long b = 0; b < m_nbinTot; ++b) {
        const Real vb = volumeOfBin(b);
        if (vb <= Real(0.0)) { continue; }
        const long a = aggOf(b);
        V[a] += vb;
        S[a] += m_s1[b*long(nr)+long(r)];
    }

    Fit f;
    f.nAgg = nagg;
    for (long a = 0; a < nagg; ++a) {
        if (V[a] > Real(0.0)) { ++f.nOccupied; f.volume += V[a]; }
    }

    for (long b = 0; b < m_nbinTot; ++b) {
        const Real vb = volumeOfBin(b);
        if (vb <= Real(0.0)) { continue; }
        const long a    = aggOf(b);
        const Real psiF = m_s1[b*long(nr)+long(r)]/vb;
        const Real psiA = S[a]/V[a];
        const Real d    = psiF - psiA;
        f.ssWithin += m_ss[b*long(nr)+long(r)] + vb*d*d;
    }

    if (mask == 0u) {
        f.ssTotal = f.ssWithin;
        f.mean    = (f.volume > Real(0.0)) ? S[0]/V[0] : Real(0.0);
        f.r2      = 0;
    } else {
        const Fit tot = fit(r, 0u, 1);
        f.ssTotal = tot.ssWithin;
        f.mean    = tot.mean;
        f.r2      = (tot.ssWithin > Real(0.0))
                  ? Real(1.0) - f.ssWithin/tot.ssWithin
                  : Real(1.0);
    }
    return f;
}

} // namespace amrcomplex
