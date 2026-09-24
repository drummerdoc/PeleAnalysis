#ifndef AMRSEGMENT_H_
#define AMRSEGMENT_H_

#include "AMRComplex.H"
#include <unordered_map>
#include <vector>

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

    // Global number of piece-to-piece edges gathered in mergeAcrossRanks: the
    // one thing the cross-rank merge replicates on every rank.
    long        nExchanged () const noexcept { return m_nExchanged; }

    // Component label of a node this rank owns; -1 if it is not in the
    // threshold set or another rank owns it.
    long label (long gid) const noexcept;

private:

    void localUnionFind ();
    void mergeAcrossRanks ();
    void tally ();

    const AMRComplex& m_cx;
    int               m_comp;
    amrex::Real       m_threshold;
    bool              m_above;

    // Union-find over gids, in two stores.
    //
    // Nodes this rank owns live in a plain array over the complex's dense
    // local numbering: 8 bytes per local cell, -1 for "not in the set". This
    // replaced hash maps and sets keyed by gid that cost well over 100 bytes
    // per in-set node -- tens of GB at production size. After the merge every
    // owned entry points straight at its root, which is the node's label.
    //
    // Nodes another rank owns -- ghosts, and whatever the merge brings in --
    // live in a hash map. That set is surface-sized, not volume-sized.
    //
    // A node's level, needed for the exact volume tally, is not stored at
    // all: the gid determines it.
    std::vector<long>             m_parent;
    std::unordered_map<long,long> m_ghostParent;

    amrex::Vector<Component> m_comps;
    long                     m_nExchanged = 0;

    long& parentRef (long g);
    long findRoot (long g);
    void unite (long a, long b);
};

} // namespace amrcomplex

#endif
