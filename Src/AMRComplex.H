#ifndef AMRCOMPLEX_H_
#define AMRCOMPLEX_H_

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_PlotFileUtil.H>

#include <array>
#include <string>
#include <algorithm>

namespace amrcomplex {

//
// A vertex of the complex.
//
// Identity is (level, iv) -- topological, never positional. Two elements share
// a vertex because they name the same cell, not because two coordinates
// compared close. That is what makes the complex watertight across coarse-fine
// interfaces, and it is what lets a periodic FillBoundary perform the periodic
// identification for free: FillBoundary copies the Node *value*, and cannot
// shift an IntVect stored inside it, so the ghost beyond the high face ends up
// literally holding the low-face node.
//
// Deliberately a plain aggregate with no user-provided constructors and no
// default member initialisers, so it is trivially default-constructible,
// trivially copyable and trivially destructible -- the requirements for living
// inside a BaseFab and riding FillBoundary/ParallelCopy.
//
struct Node
{
    enum Type : int { INIT = 0, COVERED = 1, VALID = 2 };

    int            level;
    int            type;
    long           gid;    // global index, assigned at fill time
    amrex::IntVect iv;

    // Present only to satisfy template instantiation. FabArray::ParallelCopy
    // compiles its ADD branch whatever operation is actually requested, and
    // that branch needs operator+= on the value type. Only COPY is ever used
    // here, so this is never invoked; adding it does not affect triviality of
    // copy or destruction, which is what BaseFab requires.
    Node& operator+= (const Node&) noexcept { return *this; }
};

static_assert(std::is_trivially_copyable<Node>::value,
              "Node must be trivially copyable to live in a BaseFab");
static_assert(std::is_trivially_destructible<Node>::value,
              "Node must be trivially destructible to live in a BaseFab");

AMREX_FORCE_INLINE
Node makeNode (int lev, const amrex::IntVect& iv, int type) noexcept
{
    Node n;
    n.level = lev;
    n.type  = type;
    n.gid   = -1;
    n.iv    = iv;
    return n;
}

// Strict weak ordering. IntVect::operator< is lexicographic, so this is a
// genuine total order -- unlike an epsilon-based comparison on coordinates,
// which is not transitive and makes std::set behaviour insertion-order
// dependent.
struct NodeLT
{
    bool operator() (const Node& a, const Node& b) const noexcept
    {
        if (a.level != b.level) { return a.level < b.level; }
        return a.iv < b.iv;
    }
};

AMREX_FORCE_INLINE
bool sameNode (const Node& a, const Node& b) noexcept
{
    return a.level == b.level && a.iv == b.iv;
}

//
// One cell of the complex: the dual cell whose corners are the centres of
// 2^DIM neighbouring cells. Corner order follows VTK_QUAD / VTK_HEXAHEDRON:
//
//   n0 = iv                n4 = iv           + e2
//   n1 = iv + e0           n5 = iv + e0      + e2
//   n2 = iv + e0 + e1      n6 = iv + e0 + e1 + e2
//   n3 = iv      + e1      n7 = iv      + e1 + e2
//
// Elements spanning a coarse-fine interface carry repeated gids, because a
// block of fine ghost positions all resolve to the one coarse node. Such an
// element is a degenerate hexahedron -- a pyramid, prism or tet written in
// brick form -- and that degeneracy is the stitching. It needs no case table
// and imposes no restriction on refinement ratio.
//
struct Element
{
    static constexpr int NNODE = (AMREX_SPACEDIM == 2 ? 4 : 8);
    long n[NNODE];
};

// Corner offsets of the dual cell, VTK_QUAD / VTK_HEXAHEDRON order.
AMREX_FORCE_INLINE
amrex::IntVect cornerOffset (int c) noexcept
{
#if (AMREX_SPACEDIM == 2)
    const int ox[4] = {0,1,1,0};
    const int oy[4] = {0,0,1,1};
    return amrex::IntVect(ox[c], oy[c]);
#else
    const int ox[8] = {0,1,1,0,0,1,1,0};
    const int oy[8] = {0,0,1,1,0,0,1,1};
    const int oz[8] = {0,0,0,0,1,1,1,1};
    return amrex::IntVect(ox[c], oy[c], oz[c]);
#endif
}

//
// A conforming, optionally periodic cell complex over an AMReX plotfile
// hierarchy, distributed in the AMReX decomposition throughout. Nothing is
// ever assembled on a single rank.
//
class AMRComplex
{
public:

    struct Params
    {
        std::string                        infile;
        amrex::Vector<std::string>         varnames;           // empty => none
        int                                finestLevel = -1;   // -1 => all
        amrex::Array<int,AMREX_SPACEDIM>   is_per      = {{AMREX_D_DECL(0,0,0)}};
        int                                verbose     = 0;
    };

    // Report from verify(). This is the gate described in section 9.1(b) of
    // Docs/reeb-graph-assessment.md: nothing downstream is trustworthy until
    // the volume closes.
    struct VerifyReport
    {
        amrex::Real domainVolume  = 0.0;
        amrex::Real elementVolume = 0.0;
        amrex::Real relError      = 0.0;

        long nNodes    = 0;   // VALID nodes, global
        long nElements = 0;   // global

        // shape[k] = number of elements with exactly k distinct corner nodes.
        // 8 hex, 6 prism, 5 pyramid, 4 tet in 3D; 4 quad, 3 triangle in 2D.
        // Anything else is a surprise worth understanding, and this histogram
        // is what enumerates the cases the tet subdivision must handle.
        std::array<long, Element::NNODE + 1> shape{};

        // Face-incidence census. The complex tiling the domain by volume does
        // NOT imply it is face-conforming: volume is additive over any
        // decomposition, so a mesh with T-junctions at coarse-fine interfaces
        // closes exactly while still failing to share faces. Marching tets
        // needs face-conformity, so it has to be measured separately.
        // NOTE these four are RANK-LOCAL. Each rank sees only its own
        // elements, so a face shared across a rank boundary is counted once
        // on each side and looks like two boundary faces. They are exact in
        // serial and useful for debugging; they are not an MPI invariant and
        // must not be compared across rank counts.
        long facesInterior = 0;   // incidence 2 -- correct
        long facesBoundary = 0;   // incidence 1 -- domain wall, expected
        long facesBad      = 0;   // incidence >2 -- non-manifold
        long facesDegenerate = 0; // collapsed to a point or an edge

        // The rank-independent face invariant: XOR of a hash of every face
        // instance over the whole complex. A face shared by two elements
        // contributes its hash twice and cancels, wherever those elements
        // live, so on a closed complex -- every direction periodic -- this is
        // exactly zero. With walls it is the XOR of the boundary faces alone:
        // nonzero, but still the same value at any rank count.
        unsigned long long faceParity = 0;

        // VALID nodes still carrying the fill sentinel: the value plumbing
        // has a hole the identity plumbing does not. Must be zero.
        long nodesUnfilled = 0;

        // Order-independent checksums over the node set and the element set.
        //
        // Matching element *counts* across rank counts is weak evidence: two
        // different meshes can have the same count. These are sums of per-item
        // hashes, so they are independent of the order items are visited in
        // and of which rank holds what, but sensitive to any difference in
        // membership, connectivity, or -- for nodes -- the exact bits of the
        // propagated field values. They are the actual MPI invariant.
        long nodeChecksum = 0;
        long eltChecksum  = 0;

        bool passed (amrex::Real tol = 1.0e-12) const noexcept
            { return std::abs(relError) < tol; }
        bool faceConforming () const noexcept { return facesBad == 0; }
        bool faceClosed () const noexcept { return faceParity == 0; }
        bool fieldsFilled () const noexcept { return nodesUnfilled == 0; }
    };

    explicit AMRComplex (const Params& params);

    // Build stages, in the order they must run. The constructor calls all
    // three; they are public so each can be exercised on its own in tests.
    void classifyNodes (amrex::PlotFileData& pf);  // VALID/COVERED, periodic id, coarse stamping, field fill
    void numberNodes ();     // BoxArray-index numbering; rank-count independent
    void buildElements ();   // base cell restricted to the valid box

    VerifyReport verify () const;
    void reportVerify (const VerifyReport& r) const;

    int nLevels () const noexcept { return m_nlev; }
    int nComp () const noexcept { return m_ncomp; }

    // Field values carried on the nodes, with one ghost, propagated by
    // exactly the same route as the node identity -- see stampCoarseParents.
    const amrex::MultiFab& field (int lev) const noexcept { return m_field[lev]; }
    long nNodesGlobal () const noexcept { return m_nNodesGlobal; }
    const amrex::Vector<Element>& elements () const noexcept { return m_elements; }
    const amrex::Geometry& geom (int lev) const noexcept { return m_geom[lev]; }

    // Centre of the cell naming this node, at that node's own level.
    amrex::RealVect nodePosition (const Node& n) const noexcept;

    // The canonical element traversal. Every consumer goes through this, so
    // that ownership, corner ordering and the VALID test are stated once. Two
    // copies of this loop is two chances to disagree about which elements
    // exist -- and the geometry would still verify while the surface came out
    // wrong.
    //
    // func(const Node (&corners)[NNODE], const Real* values)
    //   values is nComp() per corner, corner-major: values[c*nComp()+n].
    template <class F>
    void forEachElement (F&& func) const
    {
        amrex::Vector<amrex::Real> vals(std::size_t(Element::NNODE)*std::size_t(std::max(m_ncomp,1)));
        for (int lev = 0; lev < m_nlev; ++lev) {
            for (amrex::MFIter mfi(m_nodes[lev]); mfi.isValid(); ++mfi) {
                auto const& a = m_nodes[lev].const_array(mfi);
                amrex::Array4<const amrex::Real> fv;
                if (m_ncomp > 0) { fv = m_field[lev].const_array(mfi); }
                const amrex::BoxList owned = ownedBaseRegions(lev, mfi.index());
                for (const amrex::Box& bbx : owned) {
                    AMREX_LOOP_3D(bbx, i, j, k,
                    {
                        const amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                        Node cn[Element::NNODE];
                        bool ok = true;
                        for (int c = 0; c < Element::NNODE; ++c) {
                            const amrex::IntVect jv = iv + cornerOffset(c);
                            const Node& n = a(AMREX_D_DECL(jv[0],jv[1],jv[2]));
                            if (n.type != Node::VALID) { ok = false; break; }
                            cn[c] = n;
                            for (int m = 0; m < m_ncomp; ++m) {
                                vals[std::size_t(c)*std::size_t(m_ncomp)+std::size_t(m)]
                                    = fv(AMREX_D_DECL(jv[0],jv[1],jv[2]),m);
                            }
                        }
                        if (ok) { func(cn, vals.dataPtr()); }
                    });
                }
            }
        }
    }

private:

    void stampCoarseParents (int lev);

    // The base cells this box owns, as regions. See AMRComplex.cpp; the
    // answer is not simply "this box's valid region", because the element
    // bridging a coarse-fine interface has its base in a ghost cell stamped
    // from the coarser level.
    amrex::BoxList ownedBaseRegions (int lev, int boxIdx) const;

    // Volume the complex should occupy: the full domain in periodic
    // directions, and the domain less one level-0 cell in non-periodic ones,
    // where cell-centred data leaves the complex inset by dx/2 at each wall.
    amrex::Real expectedVolume () const;

    Params                         m_params;
    int                            m_nlev = 0;
    long                           m_nNodesGlobal = 0;

    amrex::Vector<amrex::BoxArray>            m_grids;
    amrex::Vector<amrex::DistributionMapping> m_dmap;
    amrex::Vector<amrex::Geometry>            m_geom;
    amrex::Vector<amrex::IntVect>             m_refRatio;   // size nlev-1

    amrex::Vector<amrex::FabArray<amrex::BaseFab<Node>>> m_nodes;
    amrex::Vector<amrex::MultiFab>                       m_field;
    int                                                  m_ncomp = 0;

    // levelOffset[lev] + boxOffset[lev][i] + local index = gid
    amrex::Vector<long>                       m_levelOffset;
    amrex::Vector<amrex::Vector<long>>        m_boxOffset;

    amrex::Vector<Element>                    m_elements;   // this rank's share
};

} // namespace amrcomplex

#endif
