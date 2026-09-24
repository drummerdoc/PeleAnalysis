#ifndef AMRSEGMENT_H_
#define AMRSEGMENT_H_

#include "AMRComplex.H"
#include <unordered_map>
#include <set>

namespace amrcomplex {

//
// Connected components of a threshold set on an AMRComplex: the regions a
// candidate flame-cell definition picks out.
//
// Two properties matter more than speed here.
//
// Periodicity is free. Components are found by union-find over the complex's
// own adjacency, and that adjacency already carries the periodic wrap and the
// coarse-fine stitching. No boundary plane enters, so a doubly or triply
// periodic domain needs nothing special and the arbitrary placement of the
// domain faces cannot affect the answer.
//
// Labels are rank-independent. A component is named by the MINIMUM gid it
// contains, which is a property of the component and not of who found it. The
// same run at any rank count produces identical labels, so segmentations are
// directly comparable -- across rank counts, across time snapshots, and
// across different candidate fields, which is what the comparison this tool
// exists for requires.
//
class AMRSegment
{
public:

    // comp indexes the complex's carried components; above selects
    // {f > threshold} rather than {f < threshold}.
    AMRSegment (const AMRComplex& cx, int comp, amrex::Real threshold,
                bool above = true);

    struct Component
    {
        long        label  = -1;   // minimum gid in the component
        amrex::Real volume = 0;    // exact: sum over levels of count * dV(level)
        long        nNodes = 0;
    };

    // Global, sorted by label. Populated on all ranks.
    const amrex::Vector<Component>& components () const noexcept { return m_comps; }

    long        nComponents () const noexcept { return long(m_comps.size()); }
    amrex::Real totalVolume () const noexcept;
    long        nNodesTotal () const noexcept;

    // Global number of (gid, root) pairs gathered in mergeAcrossRanks: the
    // replicated volume the cross-rank merge costs every rank.
    long        nExchanged () const noexcept { return m_nExchanged; }

    // Component label for each owned node, keyed by gid; -1 for nodes not in
    // the threshold set.
    const std::unordered_map<long,long>& nodeLabels () const noexcept { return m_label; }

private:

    void localUnionFind ();
    void mergeAcrossRanks ();
    void tally ();

    const AMRComplex& m_cx;
    int               m_comp;
    amrex::Real       m_threshold;
    bool              m_above;

    std::unordered_map<long,long> m_parent;   // union-find over gids
    std::unordered_map<long,long> m_label;    // owned gid -> component label
    // Owned gid -> the LEVEL of the cell it names, not that cell's volume.
    //
    // Cell volumes take only one value per level, so tallying integer counts
    // per level and multiplying at the end makes a component's volume an
    // exact function of integers. Accumulating reals instead would make the
    // total depend on the order the partial sums were formed, and so on the
    // rank count -- a scientific output that moves in its last bits when you
    // change the decomposition is a bad output, and here it is avoidable.
    std::unordered_map<long,int> m_lev;
    std::set<long>                       m_ownedSet;

    // Nodes another rank might also see: ghosts here, and owned nodes on a
    // box boundary or beside a COVERED cell. Only these are exchanged.
    std::set<long>                       m_shared;

    amrex::Vector<Component> m_comps;
    long                     m_nExchanged = 0;

    long findRoot (long g);
    void unite (long a, long b);
};

} // namespace amrcomplex

#endif
