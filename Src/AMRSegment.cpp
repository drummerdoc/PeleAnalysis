#include "AMRSegment.H"

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
    m_cx.forEachNode([&] (const Node& n, const Real* vals, Real cvol, bool)
    {
        owned.insert(n.gid);
        if (inSet(vals[ic])) {
            m_parent[n.gid] = n.gid;
            m_vol[n.gid] = cvol;
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

#ifdef AMREX_USE_MPI
    if (ParallelDescriptor::NProcs() > 1) {
        const int nproc  = ParallelDescriptor::NProcs();
        int cnt = int(pairs.size());
        Vector<int> cnts(nproc,0), offs(nproc,0);
        MPI_Allgather(&cnt,1,MPI_INT,cnts.dataPtr(),1,MPI_INT,
                      ParallelDescriptor::Communicator());
        for (int i=1;i<nproc;++i) { offs[i]=offs[i-1]+cnts[i-1]; }
        Vector<Long> all(std::size_t(offs[nproc-1]+cnts[nproc-1]));
        MPI_Allgatherv(pairs.dataPtr(),cnt,MPI_LONG_LONG,all.dataPtr(),
                       cnts.dataPtr(),offs.dataPtr(),MPI_LONG_LONG,
                       ParallelDescriptor::Communicator());
        pairs.swap(all);
    }
#endif

    // Every rank replays the same union-find over the same gathered pairs, so
    // all ranks reach identical labels without a broadcast.
    for (std::size_t i = 0; i < pairs.size(); i += 2) {
        const long g = long(pairs[i]), r = long(pairs[i+1]);
        if (m_parent.find(g) == m_parent.end()) { m_parent[g] = g; }
        if (m_parent.find(r) == m_parent.end()) { m_parent[r] = r; }
        unite(g, r);
    }

    for (auto& kv : m_vol) { m_label[kv.first] = findRoot(kv.first); }
}

void AMRSegment::tally ()
{
    // Per-label volume and node count, summed over OWNED nodes only so
    // nothing is double counted, then merged across ranks.
    std::map<long,std::pair<Real,long>> local;
    for (auto& kv : m_vol) {
        const long lab = m_label[kv.first];
        auto& e = local[lab];
        e.first  += kv.second;
        e.second += 1;
    }

    Vector<Long>  labs;
    Vector<Real>  vols;
    Vector<Long>  cnts_;
    for (auto& kv : local) {
        labs.push_back(Long(kv.first));
        vols.push_back(kv.second.first);
        cnts_.push_back(Long(kv.second.second));
    }

#ifdef AMREX_USE_MPI
    if (ParallelDescriptor::NProcs() > 1) {
        const int nproc = ParallelDescriptor::NProcs();
        int n = int(labs.size());
        Vector<int> cs(nproc,0), os(nproc,0);
        MPI_Allgather(&n,1,MPI_INT,cs.dataPtr(),1,MPI_INT,
                      ParallelDescriptor::Communicator());
        for (int i=1;i<nproc;++i) { os[i]=os[i-1]+cs[i-1]; }
        const int tot = os[nproc-1]+cs[nproc-1];
        Vector<Long> aL(std::size_t(tot)), aC(std::size_t(tot));
        Vector<Real> aV(std::size_t(tot));
        MPI_Allgatherv(labs.dataPtr(),n,MPI_LONG_LONG,aL.dataPtr(),
                       cs.dataPtr(),os.dataPtr(),MPI_LONG_LONG,
                       ParallelDescriptor::Communicator());
        MPI_Allgatherv(cnts_.dataPtr(),n,MPI_LONG_LONG,aC.dataPtr(),
                       cs.dataPtr(),os.dataPtr(),MPI_LONG_LONG,
                       ParallelDescriptor::Communicator());
        MPI_Allgatherv(vols.dataPtr(),n,ParallelDescriptor::Mpi_typemap<Real>::type(),
                       aV.dataPtr(),cs.dataPtr(),os.dataPtr(),
                       ParallelDescriptor::Mpi_typemap<Real>::type(),
                       ParallelDescriptor::Communicator());
        labs.swap(aL); cnts_.swap(aC); vols.swap(aV);
    }
#endif

    std::map<long,std::pair<Real,long>> merged;
    for (std::size_t i = 0; i < labs.size(); ++i) {
        auto& e = merged[long(labs[i])];
        e.first  += vols[i];
        e.second += long(cnts_[i]);
    }

    m_comps.clear();
    for (auto& kv : merged) {
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
