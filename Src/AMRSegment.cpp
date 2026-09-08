#include "AMRSegment.H"

#include "AMRParallel.H"

#include <AMReX_ParallelDescriptor.H>
#include <algorithm>
#include <map>
#include <set>
#include <vector>

using namespace amrex;

namespace amrcomplex {

long AMRSegment::findRoot (long g)
{
    long r = g;
    while (m_parent[r] != r) { r = m_parent[r]; }
    while (m_parent[g] != r) { const long nx = m_parent[g]; m_parent[g] = r; g = nx; }
    return r;
}

void AMRSegment::unite (long a, long b)
{
    long ra = findRoot(a), rb = findRoot(b);
    if (ra == rb) { return; }
    // Union by minimum gid, so the root of a component is the smallest gid it
    // contains -- a property of the component itself, not of the merge order.
    if (rb < ra) { std::swap(ra, rb); }
    m_parent[rb] = ra;
}

AMRSegment::AMRSegment (const AMRComplex& cx, int comp, Real threshold, bool above)
    : m_cx(cx), m_comp(comp), m_threshold(threshold), m_above(above)
{
    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < cx.nComp());
    localUnionFind();
    mergeAcrossRanks();
    tally();
}

void AMRSegment::localUnionFind ()
{
    const int nc = m_cx.nComp();
    const int ic = m_comp;
    const Real th = m_threshold;
    const bool ab = m_above;
    auto inSet = [=] (Real f) { return ab ? (f > th) : (f < th); };

    // Owned nodes: make a set for each, and record its cell volume. Doing this
    // before the element pass means an isolated node -- one no surviving
    // element touches -- is still a component with a volume, rather than
    // silently absent from the tally.
    std::set<long> owned;
    m_cx.forEachNode([&] (const Node& n, const Real* vals, Real, bool)
    {
        owned.insert(n.gid);
        if (inSet(vals[ic])) {
            m_parent[n.gid] = n.gid;
            m_lev[n.gid] = n.level;
        }
    });

    // Elements give the adjacency, including the periodic wrap and the
    // coarse-fine stitch.
    m_cx.forEachElement([&] (const Node (&cn)[Element::NNODE], const Real* vals)
    {
        int first = -1;
        for (int c = 0; c < Element::NNODE; ++c) {
            if (!inSet(vals[std::size_t(c)*std::size_t(nc)+std::size_t(ic)])) { continue; }
            const long g = cn[c].gid;
            if (m_parent.find(g) == m_parent.end()) { m_parent[g] = g; }  // a ghost
            if (first < 0) { first = c; } else { unite(cn[first].gid, g); }
        }
    });

    m_ownedSet.swap(owned);
}

void AMRSegment::mergeAcrossRanks ()
{
    // Every gid this rank has an opinion about, and the local root it assigns.
    // Two ranks emitting different roots for the same gid is exactly the
    // information needed to join their components; the gid itself is the
    // shared name, which is why identity had to be global from the start.
    Vector<Long> pairs;
    for (auto& kv : m_parent) {
        const long g = kv.first;
        pairs.push_back(Long(g));
        pairs.push_back(Long(findRoot(g)));
    }

    allGather(pairs);

    // Every rank replays the same union-find over the same gathered pairs, so
    // all ranks reach identical labels without a broadcast.
    for (std::size_t i = 0; i < pairs.size(); i += 2) {
        const long g = long(pairs[i]), r = long(pairs[i+1]);
        if (m_parent.find(g) == m_parent.end()) { m_parent[g] = g; }
        if (m_parent.find(r) == m_parent.end()) { m_parent[r] = r; }
        unite(g, r);
    }

    for (auto& kv : m_lev) { m_label[kv.first] = findRoot(kv.first); }
}

void AMRSegment::tally ()
{
    // Integer counts per (label, level), summed over OWNED nodes only so
    // nothing is double counted. Volume is formed from those counts at the
    // end, so it is an exact function of integers and identical at any rank
    // count.
    std::map<std::pair<long,int>,long> local;
    for (auto& kv : m_lev) { ++local[std::make_pair(m_label[kv.first], kv.second)]; }

    Vector<Long> labs, levs, cnts_;
    for (auto& kv : local) {
        labs.push_back(Long(kv.first.first));
        levs.push_back(Long(kv.first.second));
        cnts_.push_back(Long(kv.second));
    }

    allGather(labs);
    allGather(levs);
    allGather(cnts_);

    std::map<std::pair<long,int>,long> merged;
    for (std::size_t i = 0; i < labs.size(); ++i) {
        merged[std::make_pair(long(labs[i]), int(levs[i]))] += long(cnts_[i]);
    }

    Vector<Real> dV;
    dV.resize(std::size_t(m_cx.nLevels()), 0.0);
    for (int lev = 0; lev < m_cx.nLevels(); ++lev) {
        const auto dx = m_cx.geom(lev).CellSize();
        Real v = 1.0;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { v *= dx[d]; }
        dV[std::size_t(lev)] = v;
    }

    std::map<long,std::pair<Real,long>> byLabel;
    for (auto& kv : merged) {
        auto& e = byLabel[kv.first.first];
        e.first  += Real(kv.second) * dV[std::size_t(kv.first.second)];
        e.second += kv.second;
    }

    m_comps.clear();
    for (auto& kv : byLabel) {
        Component c;
        c.label  = kv.first;
        c.volume = kv.second.first;
        c.nNodes = kv.second.second;
        m_comps.push_back(c);
    }
}

Real AMRSegment::totalVolume () const noexcept
{
    Real v = 0;
    for (const auto& c : m_comps) { v += c.volume; }
    return v;
}

long AMRSegment::nNodesTotal () const noexcept
{
    long n = 0;
    for (const auto& c : m_comps) { n += c.nNodes; }
    return n;
}

} // namespace amrcomplex
