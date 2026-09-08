#ifndef AMRCOMPLEX_H_
#define AMRCOMPLEX_H_

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_PlotFileUtil.H>

#include <array>
#include <string>

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

        bool passed (amrex::Real tol = 1.0e-12) const noexcept
            { return std::abs(relError) < tol; }
    };

    explicit AMRComplex (const Params& params);

    // Build stages, in the order they must run. The constructor calls all
    // three; they are public so each can be exercised on its own in tests.
    void classifyNodes ();   // VALID/COVERED, periodic identification, coarse stamping
    void numberNodes ();     // BoxArray-index numbering; rank-count independent
    void buildElements ();   // base cell restricted to the valid box

    VerifyReport verify () const;
    void reportVerify (const VerifyReport& r) const;

    int nLevels () const noexcept { return m_nlev; }
    long nNodesGlobal () const noexcept { return m_nNodesGlobal; }
    const amrex::Vector<Element>& elements () const noexcept { return m_elements; }
    const amrex::Geometry& geom (int lev) const noexcept { return m_geom[lev]; }

    // Centre of the cell naming this node, at that node's own level.
    amrex::RealVect nodePosition (const Node& n) const noexcept;

private:

    void stampCoarseParents (int lev);

    // Which box owns the element whose base cell is iv. See AMRComplex.cpp;
    // the answer is not simply "the box whose valid region contains iv",
    // because the element that bridges a coarse-fine interface has its base
    // in a ghost cell stamped from the coarser level.
    bool ownsBaseCell (int lev, const amrex::IntVect& iv, int boxIdx) const;

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

    // levelOffset[lev] + boxOffset[lev][i] + local index = gid
    amrex::Vector<long>                       m_levelOffset;
    amrex::Vector<amrex::Vector<long>>        m_boxOffset;

    amrex::Vector<Element>                    m_elements;   // this rank's share
};

} // namespace amrcomplex

#endif
