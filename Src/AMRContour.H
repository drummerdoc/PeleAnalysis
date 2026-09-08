#ifndef AMRCONTOUR_H_
#define AMRCONTOUR_H_

#include "AMRComplex.H"

#include <map>
#include <utility>

namespace amrcomplex {

//
// An isosurface extracted from an AMRComplex.
//
// Two choices carry the whole design.
//
// 1. The complex is subdivided into tetrahedra and contoured with marching
//    TETRAHEDRA, not marching cubes. Marching tets has three distinct cases
//    and no ambiguous faces, so the result is topologically consistent by
//    construction -- there is no asymptotic-decider question to get wrong,
//    and none of the ~400 lines of case tables in isosurface.cpp.
//
// 2. A surface vertex is identified by the EDGE it lies on, and that edge is
//    named by the unordered pair of node gids at its ends. Two elements
//    sharing an edge therefore find the same vertex by identity, at any
//    level and across the periodic wrap -- no coordinate welding, no epsilon,
//    and none of the crack and mis-placement failures that follow from it.
//
struct Triangle { long v[3]; };

class AMRContour
{
public:

    // minGidDiagonal=false is a deliberately WRONG subdivision, kept as a
    // negative control: it splits every face on its local 0-2 diagonal
    // regardless of gids, so two cells sharing a face disagree whenever they
    // list its corners differently. If the closed-surface check cannot tell
    // that apart from the correct rule, the check is worthless. Never use it
    // for anything but testing.
    AMRContour (const AMRComplex& cx, int comp, amrex::Real isoVal,
                bool minGidDiagonal = true);

    struct Report
    {
        amrex::Real area = 0;

        // RANK-LOCAL: a surface vertex on a node-edge shared between two
        // ranks is created independently on each, so this over-counts under
        // MPI. Exact at one rank. Triangles are not affected -- each belongs
        // to exactly one element, which belongs to exactly one rank.
        long nVertices  = 0;

        long nTriangles = 0;
        long nEdges     = 0;
        long eulerChar  = 0;     // V - E + F

        // Rank-local, like the face census in AMRComplex: each rank sees only
        // its own triangles, so an edge shared across a rank boundary looks
        // like two boundary edges. Exact at one rank.
        long edgesInterior = 0;  // incidence 2 -- closed surface
        long edgesBoundary = 0;  // incidence 1
        long edgesBad      = 0;  // incidence >2 -- non-manifold

        // Rank-independent: XOR of a hash of every triangle-edge instance.
        // An edge shared by two triangles cancels wherever they live, so a
        // closed surface gives exactly zero at any rank count.
        //
        // This is the check that catches an inconsistent face diagonal in the
        // tet subdivision. The volume gate in AMRComplex cannot: tets of a
        // hexahedron sum to its volume whichever diagonals they use.
        unsigned long long edgeParity = 0;

        // Modular sum over triangles of a hash of their three vertex KEYS.
        // Additive and order-independent, and rank-independent because each
        // triangle is produced exactly once and the keys are global.
        long triChecksum = 0;

        bool closed () const noexcept { return edgeParity == 0; }
        bool manifold () const noexcept { return edgesBad == 0; }
    };

    Report verify () const;
    void reportVerify (const Report& r) const;

    const amrex::Vector<amrex::RealVect>& vertices () const noexcept { return m_verts; }
    const amrex::Vector<Triangle>&        triangles () const noexcept { return m_tris; }

private:

    // Split a (possibly degenerate) dual cell into tetrahedra by coning from
    // its minimum-gid corner. Returns the number of non-degenerate tets, with
    // corner indices 0..NNODE-1 written into tet[][4].
    //
    // Why this is conforming. Every face of the cell ends up split by the
    // diagonal through its own minimum-gid corner: for the three faces
    // containing the cone apex that is automatic, since a corner minimal over
    // the cell is minimal over any face containing it; for the three opposite
    // faces it is imposed explicitly. Two cells sharing a face see the same
    // four gids and so choose the same diagonal -- across levels too, because
    // gid is global, not per-level. Nothing depends on index parity, which is
    // what would otherwise fail exactly where two levels meet and the two
    // sides index differently.
    static int tetrahedralize (const Node (&cn)[Element::NNODE],
                               int tet[6][4], bool minGidDiagonal);

    long edgeVertex (long ga, long gb,
                     const amrex::RealVect& pa, const amrex::RealVect& pb,
                     amrex::Real fa, amrex::Real fb,
                     const amrex::RealVect& origin);

    const AMRComplex& m_cx;
    int               m_comp;
    amrex::Real       m_isoVal;
    bool              m_minGidDiag;

    std::map<std::pair<long,long>, long> m_vertCache;  // edge (gid,gid) -> vertex
    amrex::Vector<amrex::RealVect>       m_verts;      // canonical, wrapped

    // The GLOBAL name of each surface vertex: the node-edge it lies on. The
    // index into m_verts is rank-local and means nothing to another rank, so
    // every checksum and every incidence test must be phrased in these, not
    // in indices. Getting that wrong is what made the first MPI run report an
    // open surface on a closed one.
    amrex::Vector<std::pair<long,long>>  m_vertKey;
    amrex::Vector<Triangle>              m_tris;
    amrex::Real                          m_area = 0;   // summed in element frames
};

} // namespace amrcomplex

#endif
