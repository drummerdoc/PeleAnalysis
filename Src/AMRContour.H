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
// One element of the extracted surface: a triangle in 3D, a segment in 2D.
// Its boundary is correspondingly its 3 edges or its 2 endpoints, and the
// closure test is the same statement in both -- every boundary item shared by
// exactly two surface elements.
struct SurfElement
{
    static constexpr int NVERT = AMREX_SPACEDIM;
    long v[NVERT];
};

class AMRContour
{
public:

    // minGidDiagonal=false is a deliberately WRONG subdivision, kept as a
    // negative control: it splits every face on its local 0-2 diagonal
    // regardless of gids, so two cells sharing a face disagree whenever they
    // list its corners differently. If the closed-surface check cannot tell
    // that apart from the correct rule, the check is worthless. Never use it
    // for anything but testing.
    // checkDeterminism recomputes every cached vertex and compares bitwise.
    // It must be set at construction because the constructor does all the
    // work. See the note on m_checkDeterminism.
    // orientElements=false leaves the winding as the marching case happens
    // to produce it. That is the state of the code before orientation was
    // imposed, kept as a negative control: if the consistency check cannot
    // tell it apart from the oriented result, the check is worthless.
    AMRContour (const AMRComplex& cx, int comp, amrex::Real isoVal,
                bool minGidDiagonal = true, bool checkDeterminism = false,
                bool orientElements = true);

    struct Report
    {
        amrex::Real area = 0;

        // RANK-LOCAL: a surface vertex on a node-edge shared between two
        // ranks is created independently on each, so this over-counts under
        // MPI. Exact at one rank. Triangles are not affected -- each belongs
        // to exactly one element, which belongs to exactly one rank.
        long nVertices  = 0;

        long nElements  = 0;     // triangles in 3D, segments in 2D
        long nEdges     = 0;     // boundary items: edges in 3D, vertices in 2D
        long eulerChar  = 0;     // V-E+F in 3D (2 for a sphere); V-E in 2D (0)

        // Rank-local, like the face census in AMRComplex: each rank sees only
        // its own triangles, so an edge shared across a rank boundary looks
        // like two boundary edges. Exact at one rank.
        long edgesInterior = 0;  // incidence 2 -- closed surface/curve
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

        // Orientation.
        //
        // A closed surface is consistently oriented exactly when every
        // interior edge is traversed in OPPOSITE directions by the two
        // elements sharing it. That is combinatorial and exact -- no normals,
        // no tolerance, no reference direction needed after the fact.
        //
        // orientBad is the rank-local count of boundary items that fail it.
        // orientParity is the rank-independent form: each traversal
        // contributes +hash or -hash by direction, so a correctly oriented
        // pair cancels wherever the two elements live, and the global sum is
        // exactly zero.
        long orientBad = 0;
        long orientParity = 0;

        // Enclosed volume by the divergence theorem, with the orientation as
        // built (normals toward increasing f). Meaningful only for a surface
        // that does not cross a periodic face; diagnostic, not an invariant.
        amrex::Real signedVolume = 0;

        bool oriented () const noexcept { return orientParity == 0; }

        // Modular sum over triangles of a hash of their three vertex KEYS.
        // Additive and order-independent, and rank-independent because each
        // triangle is produced exactly once and the keys are global.
        long triChecksum = 0;

        bool closed () const noexcept { return edgeParity == 0; }
        bool manifold () const noexcept { return edgesBad == 0; }
    };

    Report verify () const;
    void reportVerify (const Report& r) const;

    // Write the surface in MEF, the format the *MEF toolchain and
    // surfMEFtoDAT read: a title line, a variable-name line, "nElts
    // nodesPerElt", a node blob written as an FArrayBox in NODE-major order,
    // and 1-based connectivity as raw ints.
    //
    // Under MPI this gathers to the IO processor, as isosurface.cpp does --
    // the same scaling wall, and the reason a .pvtu path exists on the
    // roadmap. It differs in one respect that matters: duplicate vertices are
    // merged by their node-edge KEY, not by comparing coordinates, so the
    // merge is exact and the global node ordering is a sort on those keys and
    // therefore identical at any rank count.
    void writeMEF (const std::string& outfile,
                   const std::string& label,
                   const amrex::Vector<std::string>& varnames) const;

    const amrex::Vector<amrex::RealVect>& vertices () const noexcept { return m_verts; }

    // Interpolated field values carried on the surface, nComp per vertex,
    // vertex-major. Every component is interpolated at the SAME parameter as
    // the position -- the one determined by the isosurface component -- which
    // is the convention isosurface.cpp uses.
    const amrex::Vector<amrex::Real>& vertexData () const noexcept { return m_vertData; }
    const amrex::Vector<SurfElement>&    surfElements () const noexcept { return m_tris; }

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
    // Split a (possibly degenerate) dual cell into simplices by coning from
    // its minimum-gid corner: tetrahedra in 3D, triangles in 2D. Returns the
    // number of non-degenerate simplices.
    static int simplexify (const Node (&cn)[Element::NNODE],
                           int simp[6][AMREX_SPACEDIM+1], bool minGidDiagonal);

    // Everything a surface vertex is depends on the EDGE alone: the two end
    // nodes, their values, and their carried data. Nothing about the element
    // that happened to reach the edge first may enter, or two elements -- and
    // so two ranks -- will produce copies differing in the last bits.
    //
    // da/db are the full component vectors at the two ends; every one is
    // interpolated at the parameter the isosurface component determines.
    long edgeVertex (const Node& na, const Node& nb,
                     amrex::Real fa, amrex::Real fb,
                     const amrex::Real* da, const amrex::Real* db);

    // Recompute a cached vertex and compare bitwise. A serial test for a
    // parallel defect: within one rank the cache hides a frame-dependent
    // computation entirely, and it only surfaces as two ranks disagreeing.
    bool m_checkDeterminism = false;
    long m_determinismMismatches = 0;

public:
    long determinismMismatches () const noexcept { return m_determinismMismatches; }
private:

    const AMRComplex& m_cx;
    int               m_comp;
    amrex::Real       m_isoVal;
    bool              m_minGidDiag;
    bool              m_orient = true;

    std::map<std::pair<long,long>, long> m_vertCache;  // edge (gid,gid) -> vertex
    amrex::Vector<amrex::RealVect>       m_verts;      // canonical, wrapped

    // The GLOBAL name of each surface vertex: the node-edge it lies on. The
    // index into m_verts is rank-local and means nothing to another rank, so
    // every checksum and every incidence test must be phrased in these, not
    // in indices. Getting that wrong is what made the first MPI run report an
    // open surface on a closed one.
    amrex::Vector<std::pair<long,long>>  m_vertKey;
    amrex::Vector<amrex::Real>           m_vertData;   // nComp per vertex
    amrex::Vector<SurfElement>           m_tris;
    amrex::Real                          m_area = 0;   // summed in element frames
};

} // namespace amrcomplex

#endif
