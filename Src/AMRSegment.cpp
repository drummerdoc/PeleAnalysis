#include "AMRSegment.H"

#include "AMRParallel.H"

#include <AMReX_ParallelDescriptor.H>
#include <algorithm>
#include <map>
#include <vector>

using namespace amrex;

namespace amrcomplex {

// The parent slot of a gid that is already in the union-find: the array for
// owned nodes, the map for everything else.
long& AMRSegment::parentRef (long g)
{
    const long i = m_cx.localIndex(g);
    if (i >= 0) { return m_parent[std::size_t(i)]; }
    return m_ghostParent.at(g);
}

long AMRSegment::findRoot (long g)
{
    long r = g;
    while (parentRef(r) != r) { r = parentRef(r); }
    while (parentRef(g) != r) { long& p = parentRef(g); const long nx = p; p = r; g = nx; }
    return r;
}

void AMRSegment::unite (long a, long b)
{
    long ra = findRoot(a), rb = findRoot(b);
    if (ra == rb) { return; }
    // Union by minimum gid, so the root of a component is the smallest gid it
    // contains -- a property of the component itself, not of the merge order.
    if (rb < ra) { std::swap(ra, rb); }
    parentRef(rb) = ra;
}

long AMRSegment::label (long gid) const noexcept
{
    const long i = m_cx.localIndex(gid);
    return (i >= 0) ? m_parent[std::size_t(i)] : -1;
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

    m_parent.assign(std::size_t(m_cx.nLocalCells()), -1L);

    // Owned nodes: make a set for each. Doing this before the element pass
    // means an isolated node -- one no surviving element touches -- is still
    // a component with a volume, rather than silently absent from the tally.
    m_cx.forEachNode([&] (const Node& n, const Real* vals, Real)
    {
        if (!inSet(vals[ic])) { return; }
        m_parent[std::size_t(m_cx.localIndex(n.gid))] = n.gid;
    });

    // Elements give the adjacency, including the periodic wrap and the
    // coarse-fine stitch. A corner owned here is already in the set, from the
    // node pass; any other corner is a ghost owned by another rank.
    m_cx.forEachElement([&] (const Node (&cn)[Element::NNODE], const Real* vals)
    {
        int first = -1;
        for (int c = 0; c < Element::NNODE; ++c) {
            if (!inSet(vals[std::size_t(c)*std::size_t(nc)+std::size_t(ic)])) { continue; }
            const long g = cn[c].gid;
            if (m_cx.localIndex(g) < 0) { m_ghostParent.emplace(g, g); }
            if (first < 0) { first = c; } else { unite(cn[first].gid, g); }
        }
    });
}

void AMRSegment::mergeAcrossRanks ()
{
    // Join the per-rank pieces into global components without replicating
    // anything node-sized.
    //
    // The first version all-gathered a (gid, root) pair for every node another
    // rank might see, and every rank replayed all of them. At the lowest
    // threshold that is about a quarter of all nodes, on every rank: on levels
    // 0-2 of the first real snapshot (261M nodes, 64 ranks) it did not fit in
    // 122 GB.
    //
    // Every cross-rank adjacency is an element with a corner another rank
    // owns -- a ghost here. So:
    //
    //   1. each ghost's (gid, local root) goes to the rank that OWNS the gid,
    //      and only there;
    //   2. the owner pairs that root with its own local root for the gid: one
    //      edge between two pieces. Many ghosts along the same stretch of
    //      boundary give the same edge, and duplicates are dropped;
    //   3. only those edges are gathered to every rank. Their number is set by
    //      how many component pieces meet at rank boundaries, not by nodes.
    //
    // Every rank then unites the same edge set by minimum gid, so labels are
    // still the minimum gid in each component and identical on every rank.
    const int nproc = ParallelDescriptor::NProcs();

    Vector<Vector<Long>> byRank(nproc);
    for (const auto& kv : m_ghostParent) {
        const long g = kv.first;
        auto& out = byRank[m_cx.gidOwner(g)];
        out.push_back(Long(g));
        out.push_back(Long(findRoot(g)));
    }
    Vector<Long> recv = allToAll(byRank);
    Vector<Vector<Long>>().swap(byRank);

    std::vector<std::pair<long,long>> edges;
    for (std::size_t i = 0; i < recv.size(); i += 2) {
        const long a = findRoot(long(recv[i]));   // owned here, so in the array
        const long b = long(recv[i+1]);
        if (a != b) { edges.emplace_back(std::min(a,b), std::max(a,b)); }
    }
    Vector<Long>().swap(recv);
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

    Vector<Long> flat;
    flat.reserve(2*edges.size());
    for (const auto& e : edges) { flat.push_back(Long(e.first)); flat.push_back(Long(e.second)); }
    std::vector<std::pair<long,long>>().swap(edges);

    allGather(flat);
    m_nExchanged = long(flat.size() / 2);

    // Union-find over piece roots alone: a map, but only as large as the edge
    // set.
    std::unordered_map<long,long> rp;
    auto gfind = [&] (long g) {
        long r = g;
        for (auto it = rp.find(r); it != rp.end() && it->second != r; it = rp.find(r)) { r = it->second; }
        while (g != r) { long& p = rp[g]; const long nx = p; p = r; g = nx; }
        return r;
    };
    for (std::size_t i = 0; i < flat.size(); i += 2) {
        const long a = long(flat[i]), b = long(flat[i+1]);
        rp.emplace(a, a);
        rp.emplace(b, b);
        long ra = gfind(a), rb = gfind(b);
        if (ra == rb) { continue; }
        if (rb < ra) { std::swap(ra, rb); }
        rp[rb] = ra;
    }

    // Point every owned node straight at its component's root: its label.
    // Two passes, and the order matters. Writing a global label into the
    // array while other nodes' paths still run through it sends those paths
    // to a gid this rank may hold in neither store. So first collapse every
    // path to its LOCAL root, which is always present, and only then map
    // local roots to global labels.
    for (std::size_t i = 0; i < m_parent.size(); ++i) {
        if (m_parent[i] >= 0) { m_parent[i] = findRoot(m_cx.localGid(long(i))); }
    }
    for (std::size_t i = 0; i < m_parent.size(); ++i) {
        if (m_parent[i] >= 0) { m_parent[i] = gfind(m_parent[i]); }
    }
    std::unordered_map<long,long>().swap(m_ghostParent);
}

void AMRSegment::tally ()
{
    // Integer counts per (label, level), summed over OWNED nodes only so
    // nothing is double counted. Volume is formed from those counts at the
    // end, so it is an exact function of integers and identical at any rank
    // count.
    std::map<std::pair<long,int>,long> local;
    for (std::size_t i = 0; i < m_parent.size(); ++i) {
        if (m_parent[i] < 0) { continue; }
        const int lev = m_cx.gidLevel(m_cx.localGid(long(i)));
        ++local[std::make_pair(m_parent[i], lev)];
    }

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
