#include "AMRContour.H"

#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <map>
#include <array>

using namespace amrex;

namespace amrcomplex {

namespace {

#if (AMREX_SPACEDIM == 3)
// Faces in VTK_HEXAHEDRON order, each listed cyclically.
const int HFACE[6][4] = {{0,1,2,3},{4,5,6,7},{0,1,5,4},
                         {1,2,6,5},{2,3,7,6},{3,0,4,7}};
#endif

AMREX_FORCE_INLINE
void mixIn (unsigned long long& h, unsigned long long v) noexcept
{
    h ^= v;
    h *= 1099511628211ULL;
    h ^= (h >> 29);
}

constexpr long CHECKSUM_MOD = 2147483647L;   // 2^31 - 1

AMREX_FORCE_INLINE
long foldChecksum (unsigned long long h) noexcept
{
    return long(h % (unsigned long long)(CHECKSUM_MOD));
}

// Hash a triangle edge by the GLOBAL names of its two endpoints -- each a
// node-edge (gid,gid) -- rather than by local vertex indices, which differ
// from rank to rank and so make a shared edge look like two.
AMREX_FORCE_INLINE
unsigned long long edgeHash (std::pair<long,long> ka,
                             std::pair<long,long> kb) noexcept
{
    if (kb < ka) { std::swap(ka,kb); }
    unsigned long long h = 1469598103934665603ULL;
    mixIn(h, (unsigned long long)(ka.first));
    mixIn(h, (unsigned long long)(ka.second));
    mixIn(h, (unsigned long long)(kb.first));
    mixIn(h, (unsigned long long)(kb.second));
    return h;
}

} // anonymous namespace


int AMRContour::simplexify (const Node (&cn)[Element::NNODE],
                            int simp[6][AMREX_SPACEDIM+1], bool minGidDiagonal)
{
    // Cone from the minimum-gid corner. In 3D that means the three faces not
    // containing the apex, each split on its own min-gid diagonal; in 2D the
    // two edges not containing it. Either way every face of the cell ends up
    // split by the diagonal through its minimum-gid corner, so two cells
    // sharing a face split it identically -- across levels too, since gid is
    // global and nothing here depends on index parity.
    int imin = 0;
    for (int c = 1; c < Element::NNODE; ++c) {
        if (cn[c].gid < cn[imin].gid) { imin = c; }
    }

    int n = 0;

    auto emit = [&] (const int* v) -> void
    {
        // A stitch cell has repeated corners, so some simplices collapse.
        // Dropping them here is what turns the degenerate cell into the
        // pyramid, prism or triangle it actually is.
        for (int a = 0; a < AMREX_SPACEDIM+1; ++a) {
            for (int b = a+1; b < AMREX_SPACEDIM+1; ++b) {
                if (cn[v[a]].gid == cn[v[b]].gid) { return; }
            }
        }
        for (int a = 0; a < AMREX_SPACEDIM+1; ++a) { simp[n][a] = v[a]; }
        ++n;
    };

#if (AMREX_SPACEDIM == 2)
    // The four edges of the quad, cyclically.
    const int QEDGE[4][2] = {{0,1},{1,2},{2,3},{3,0}};
    for (int e = 0; e < 4; ++e) {
        if (QEDGE[e][0] == imin || QEDGE[e][1] == imin) { continue; }
        const int v[3] = {imin, QEDGE[e][0], QEDGE[e][1]};
        emit(v);
    }
#else
    for (int f = 0; f < 6; ++f) {

        bool hasApex = false;
        for (int c = 0; c < 4; ++c) { if (HFACE[f][c] == imin) { hasApex = true; } }
        if (hasApex) { continue; }   // coning to it would give zero volume

        int m = 0;
        if (minGidDiagonal) {
            for (int c = 1; c < 4; ++c) {
                if (cn[HFACE[f][c]].gid < cn[HFACE[f][m]].gid) { m = c; }
            }
        }   // else: always local 0-2, the negative control

        int t0[3], t1[3];
        if (m == 0 || m == 2) {             // diagonal 0-2
            t0[0]=HFACE[f][0]; t0[1]=HFACE[f][1]; t0[2]=HFACE[f][2];
            t1[0]=HFACE[f][0]; t1[1]=HFACE[f][2]; t1[2]=HFACE[f][3];
        } else {                            // diagonal 1-3
            t0[0]=HFACE[f][1]; t0[1]=HFACE[f][2]; t0[2]=HFACE[f][3];
            t1[0]=HFACE[f][1]; t1[1]=HFACE[f][3]; t1[2]=HFACE[f][0];
        }
        const int c0[4] = {imin,t0[0],t0[1],t0[2]};
        const int c1[4] = {imin,t1[0],t1[1],t1[2]};
        emit(c0);
        emit(c1);
    }
#endif
#if (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(minGidDiagonal);
#endif
    return n;
}

long AMRContour::edgeVertex (const Node& na, const Node& nb,
                             Real fa, Real fb,
                             const Real* da, const Real* db)
{
    // Canonicalise on the edge itself, so the answer cannot depend on which
    // element -- or which rank -- got here first.
    const bool aFirst = (na.gid < nb.gid);
    const Node& n0 = aFirst ? na : nb;
    const Node& n1 = aFirst ? nb : na;
    const Real  v0 = aFirst ? fa : fb;
    const Real  v1 = aFirst ? fb : fa;
    const Real* c0 = aFirst ? da : db;
    const Real* c1 = aFirst ? db : da;

    const std::pair<long,long> key(n0.gid, n1.gid);

    const auto& g   = m_cx.geom(0);
    const auto  plo = g.ProbLo();
    const auto  phi = g.ProbHi();
    const int   nc  = m_cx.nComp();

    const Real den = v1 - v0;
    const Real t   = (den == Real(0.0)) ? Real(0.5) : (m_isoVal - v0)/den;

    // Absolute positions of the two end nodes, the far one brought to its
    // minimum image of the near one. This is a function of the edge only --
    // the previous version added the enclosing element's corner-0 offset,
    // which differs between the elements sharing the edge and so produced
    // copies differing by an ulp.
    const RealVect q0 = m_cx.nodePosition(n0);
    RealVect       q1 = m_cx.nodePosition(n1);
    RealVect pw;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const Real L = phi[d] - plo[d];
        if (g.isPeriodic(d)) {
            if (q1[d] - q0[d] >  Real(0.5)*L) { q1[d] -= L; }
            if (q1[d] - q0[d] < -Real(0.5)*L) { q1[d] += L; }
        }
        Real x = q0[d] + t*(q1[d] - q0[d]);
        if (g.isPeriodic(d)) {
            while (x <  plo[d]) { x += L; }
            while (x >= phi[d]) { x -= L; }
        }
        pw[d] = x;
    }

    auto it = m_vertCache.find(key);
    if (it != m_vertCache.end()) {
        if (m_checkDeterminism) {
            const long j = it->second;
            bool same = true;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                if (m_verts[j][d] != pw[d]) { same = false; }
            }
            for (int n = 0; n < nc && same; ++n) {
                const Real want = c0[n] + t*(c1[n] - c0[n]);
                if (m_vertData[std::size_t(j)*std::size_t(nc)+std::size_t(n)] != want) {
                    same = false;
                }
            }
            if (!same) { ++m_determinismMismatches; }
        }
        return it->second;
    }

    for (int n = 0; n < nc; ++n) { m_vertData.push_back(c0[n] + t*(c1[n] - c0[n])); }

    const long idx = long(m_verts.size());
    m_verts.push_back(pw);
    m_vertKey.push_back(key);
    m_vertCache[key] = idx;
    return idx;
}

AMRContour::AMRContour (const AMRComplex& cx, int comp, Real isoVal,
                        bool minGidDiagonal, bool checkDeterminism,
                        bool orientElements)
    : m_cx(cx), m_comp(comp), m_isoVal(isoVal), m_minGidDiag(minGidDiagonal),
      m_orient(orientElements), m_checkDeterminism(checkDeterminism)
{
    const int nc = cx.nComp();
    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < nc);

    const auto& g0  = cx.geom(0);
    const auto  plo = g0.ProbLo();
    const auto  phi = g0.ProbHi();
    RealVect L;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { L[d] = phi[d] - plo[d]; }

    cx.forEachElement([&] (const Node (&cn)[Element::NNODE], const Real* vals)
    {
        // Element-local frame: corner 0 at the origin, every other corner
        // reached by its minimum image. A wrap element is unbroken here, which
        // is why lengths and areas taken in this frame are correct while the
        // same computation on raw coordinates is the classic periodic defect.
        const RealVect o = m_cx.nodePosition(cn[0]);
        RealVect q[Element::NNODE];
        Real     fv[Element::NNODE];
        for (int c = 0; c < Element::NNODE; ++c) {
            const RealVect p = m_cx.nodePosition(cn[c]);
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                Real sd = p[d] - o[d];
                if (g0.isPeriodic(d)) {
                    if (sd >  Real(0.5)*L[d]) { sd -= L[d]; }
                    if (sd < -Real(0.5)*L[d]) { sd += L[d]; }
                }
                q[c][d] = sd;
            }
            fv[c] = vals[std::size_t(c)*std::size_t(nc)+std::size_t(comp)];
        }

        int simp[6][AMREX_SPACEDIM+1];
        const int ns = simplexify(cn, simp, m_minGidDiag);

        for (int t = 0; t < ns; ++t) {

            const int NS = AMREX_SPACEDIM+1;
            int above[NS], below[NS];
            int na = 0, nb = 0;
            for (int a = 0; a < NS; ++a) {
                const int c = simp[t][a];
                if (fv[c] > m_isoVal) { above[na++] = c; } else { below[nb++] = c; }
            }
            if (na == 0 || nb == 0) { continue; }

            long v[4];
            int nv = 0;

#if (AMREX_SPACEDIM == 2)
            // A triangle cut by a level set always yields exactly two
            // crossings: the apex is whichever side has one vertex.
            {
                const int apex  = (na == 1) ? above[0] : below[0];
                const int* rest = (na == 1) ? below : above;
                for (int k = 0; k < 2; ++k) {
                    v[nv++] = edgeVertex(cn[apex], cn[rest[k]],
                                         fv[apex], fv[rest[k]],
                                         &vals[std::size_t(apex)*std::size_t(nc)],
                                         &vals[std::size_t(rest[k])*std::size_t(nc)]);
                }
            }
#else
            if (na == 1 || nb == 1) {
                const int apex  = (na == 1) ? above[0] : below[0];
                const int* rest = (na == 1) ? below : above;
                for (int k = 0; k < 3; ++k) {
                    v[nv++] = edgeVertex(cn[apex], cn[rest[k]],
                                         fv[apex], fv[rest[k]],
                                         &vals[std::size_t(apex)*std::size_t(nc)],
                                         &vals[std::size_t(rest[k])*std::size_t(nc)]);
                }
            } else {
                // Two-two split: the four crossings form a quad, listed
                // cyclically so the two triangles do not overlap.
                const int a0 = above[0], a1 = above[1];
                const int b0 = below[0], b1 = below[1];
                const int ea[4] = {a0,a0,a1,a1};
                const int eb[4] = {b0,b1,b1,b0};
                for (int k = 0; k < 4; ++k) {
                    v[nv++] = edgeVertex(cn[ea[k]], cn[eb[k]],
                                         fv[ea[k]], fv[eb[k]],
                                         &vals[std::size_t(ea[k])*std::size_t(nc)],
                                         &vals[std::size_t(eb[k])*std::size_t(nc)]);
                }
            }
#endif

#if (AMREX_SPACEDIM == 2)
            const int NSUB = 1;
            const int SUB[1][2] = {{0,1}};
#else
            const int NSUB = (nv == 3) ? 1 : 2;
            const int SUB[2][3] = {{0,1,2},{0,2,3}};
#endif
            for (int k = 0; k < NSUB; ++k) {
                SurfElement se;
                bool degen = false;
                for (int a = 0; a < SurfElement::NVERT; ++a) {
                    se.v[a] = v[SUB[k][a]];
                    for (int b = 0; b < a; ++b) {
                        if (se.v[b] == se.v[a]) { degen = true; }
                    }
                }
                if (degen) { continue; }

                // Measure in the element frame, from the same interpolated
                // positions, before any wrapping: area in 3D, length in 2D.
                RealVect w[SurfElement::NVERT];
                for (int a = 0; a < SurfElement::NVERT; ++a) {
                    const RealVect& pw = m_verts[se.v[a]];
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                        Real sd = pw[d] - o[d];
                        if (g0.isPeriodic(d)) {
                            if (sd >  Real(0.5)*L[d]) { sd -= L[d]; }
                            if (sd < -Real(0.5)*L[d]) { sd += L[d]; }
                        }
                        w[a][d] = sd;
                    }
                }
                // Orient the element so its normal points toward INCREASING
                // f -- the high side of the contour, the same side for every
                // element on the surface.
                //
                // Nothing upstream fixes this: the vertex order comes out of
                // the marching case, which follows the tet's corner order,
                // which follows the cone construction. That is arbitrary
                // relative to the geometry, so without this step neighbouring
                // elements wind opposite ways at random. Normals, curvature,
                // rendering and any enclosed-volume integral all depend on
                // it, and isosurface.cpp's CheckSurfaceNormal is an empty
                // stub, so this has never actually been guaranteed.
                //
                // The reference direction is from the below-value corners to
                // the above-value corners of this simplex -- exact, local, and
                // independent of the tet's handedness.
                RealVect ref;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    Real sa = 0, sb = 0;
                    for (int a = 0; a < na; ++a) { sa += q[above[a]][d]; }
                    for (int b = 0; b < nb; ++b) { sb += q[below[b]][d]; }
                    ref[d] = sa/Real(na) - sb/Real(nb);
                }
#if (AMREX_SPACEDIM == 2)
                {
                    // In 2D the element is a segment; its normal is the
                    // left-hand normal of v0->v1. Reverse the segment when
                    // that points down-gradient.
                    const Real tx = w[1][0]-w[0][0], ty = w[1][1]-w[0][1];
                    if (m_orient && (ty*ref[0] - tx*ref[1]) < Real(0.0)) {
                        std::swap(se.v[0], se.v[1]);
                        std::swap(w[0], w[1]);
                    }
                }
#else
                {
                    const Real a1x = w[1][0]-w[0][0], a1y = w[1][1]-w[0][1], a1z = w[1][2]-w[0][2];
                    const Real a2x = w[2][0]-w[0][0], a2y = w[2][1]-w[0][1], a2z = w[2][2]-w[0][2];
                    const Real nx = a1y*a2z - a1z*a2y;
                    const Real ny = a1z*a2x - a1x*a2z;
                    const Real nz = a1x*a2y - a1y*a2x;
                    if (m_orient && (nx*ref[0] + ny*ref[1] + nz*ref[2]) < Real(0.0)) {
                        std::swap(se.v[1], se.v[2]);
                        std::swap(w[1], w[2]);
                    }
                }
#endif
                m_tris.push_back(se);

#if (AMREX_SPACEDIM == 2)
                const Real dx0 = w[1][0]-w[0][0];
                const Real dy0 = w[1][1]-w[0][1];
                m_area += std::sqrt(dx0*dx0 + dy0*dy0);
#else
                const Real e1x = w[1][0]-w[0][0], e1y = w[1][1]-w[0][1], e1z = w[1][2]-w[0][2];
                const Real e2x = w[2][0]-w[0][0], e2y = w[2][1]-w[0][1], e2z = w[2][2]-w[0][2];
                const Real cx0 = e1y*e2z - e1z*e2y;
                const Real cy0 = e1z*e2x - e1x*e2z;
                const Real cz0 = e1x*e2y - e1y*e2x;
                m_area += Real(0.5)*std::sqrt(cx0*cx0 + cy0*cy0 + cz0*cz0);
#endif
            }
        }
    });
}

AMRContour::Report AMRContour::verify () const
{
    Report r;

    Real area = m_area;
    Long nv = Long(m_verts.size());
    Long nt = Long(m_tris.size());

    std::map<std::pair<long,long>,int> edgeCount;   // local indices; serial only
    unsigned long long par = 0ULL;
    Long triSum = 0;

    // The boundary of a surface element: its 3 edges in 3D, its 2 endpoints
    // in 2D. Closure is the same statement either way -- every boundary item
    // shared by exactly two surface elements.
    const int NV = SurfElement::NVERT;
    std::map<std::pair<long,long>,int> orientNet;   // local indices; serial only
    Long orientSum = 0;
    Real svol = 0;

    for (const auto& t : m_tris) {
#if (AMREX_SPACEDIM == 3)
        {   // Divergence theorem on the surface as oriented.
            const RealVect& p0 = m_verts[t.v[0]];
            const RealVect& p1 = m_verts[t.v[1]];
            const RealVect& p2 = m_verts[t.v[2]];
            svol += (p0[0]*(p1[1]*p2[2]-p1[2]*p2[1])
                   - p0[1]*(p1[0]*p2[2]-p1[2]*p2[0])
                   + p0[2]*(p1[0]*p2[1]-p1[1]*p2[0]))/Real(6.0);
        }
#endif
        // Directed boundary items: in 3D the three directed edges; in 2D the
        // tail and the head of the segment.
        for (int k = 0; k < NV; ++k) {
#if (AMREX_SPACEDIM == 2)
            const std::pair<long,long> ka = m_vertKey[t.v[k]];
            const std::pair<long,long> kb = ka;
            const int sgn = (k == 0) ? -1 : +1;          // tail, then head
            const long ia = t.v[k], ib = t.v[k];
#else
            const long ia = t.v[k], ib = t.v[(k+1)%NV];
            const std::pair<long,long> ka = m_vertKey[ia];
            const std::pair<long,long> kb = m_vertKey[ib];
            const int sgn = (ka < kb) ? +1 : -1;
#endif
            long lo = ia, hi = ib;
            if (lo > hi) { std::swap(lo,hi); }
            orientNet[std::make_pair(lo,hi)] += sgn;

            std::pair<long,long> k0 = ka, k1 = kb;
            if (k1 < k0) { std::swap(k0,k1); }
            const long hh = foldChecksum(edgeHash(k0,k1));
            orientSum = ((orientSum + sgn*hh) % CHECKSUM_MOD + CHECKSUM_MOD) % CHECKSUM_MOD;
        }

        for (int k = 0; k < NV; ++k) {
#if (AMREX_SPACEDIM == 2)
            long a = t.v[k], b = t.v[k];
#else
            long a = t.v[k], b = t.v[(k+1)%NV];
#endif
            if (a > b) { std::swap(a,b); }
            ++edgeCount[std::make_pair(a,b)];
            // Parity over the global names, so a boundary item shared between
            // two ranks cancels exactly as one shared within a rank does.
#if (AMREX_SPACEDIM == 2)
            par ^= edgeHash(m_vertKey[t.v[k]], m_vertKey[t.v[k]]);
#else
            par ^= edgeHash(m_vertKey[t.v[k]], m_vertKey[t.v[(k+1)%NV]]);
#endif
        }

        std::pair<long,long> kk[NV];
        for (int k = 0; k < NV; ++k) { kk[k] = m_vertKey[t.v[k]]; }
        std::sort(kk, kk+NV);
        unsigned long long th = 1469598103934665603ULL;
        for (int k = 0; k < NV; ++k) {
            mixIn(th, (unsigned long long)(kk[k].first));
            mixIn(th, (unsigned long long)(kk[k].second));
        }
        triSum = (triSum + foldChecksum(th)) % CHECKSUM_MOD;
    }

    for (const auto& on : orientNet) { if (on.second != 0) { ++r.orientBad; } }

    Long ne = Long(edgeCount.size());
    for (const auto& ec : edgeCount) {
        if      (ec.second == 1) { ++r.edgesBoundary; }
        else if (ec.second == 2) { ++r.edgesInterior; }
        else                     { ++r.edgesBad; }
    }

    ParallelDescriptor::ReduceRealSum(area);
    ParallelDescriptor::ReduceLongSum(nv);
    ParallelDescriptor::ReduceLongSum(nt);
    ParallelDescriptor::ReduceLongSum(ne);
    ParallelDescriptor::ReduceLongSum(triSum);
    ParallelDescriptor::ReduceLongSum(orientSum);
    ParallelDescriptor::ReduceRealSum(svol);
#ifdef AMREX_USE_MPI
    if (ParallelDescriptor::NProcs() > 1) {
        unsigned long long tmp = par;
        MPI_Allreduce(&tmp, &par, 1, MPI_UNSIGNED_LONG_LONG, MPI_BXOR,
                      ParallelDescriptor::Communicator());
    }
#endif

    r.area       = area;
    r.nVertices  = nv;
    r.nElements  = nt;
    r.nEdges     = ne;
#if (AMREX_SPACEDIM == 2)
    r.eulerChar  = nv - ne;              // 0 for a closed curve
#else
    r.eulerChar  = nv - ne + nt;         // 2 for a sphere
#endif
    r.edgeParity  = par;
    r.triChecksum  = triSum % CHECKSUM_MOD;
    r.orientParity = ((orientSum % CHECKSUM_MOD) + CHECKSUM_MOD) % CHECKSUM_MOD;
    r.signedVolume = svol;
    return r;
}

void AMRContour::reportVerify (const Report& r) const
{
    Print() << "\n=== AMRContour verify ===\n";
    Print() << "  isoValue       : " << m_isoVal << "\n";
    Print() << "~ vertices       : " << r.nVertices
            << "  (rank-local: shared vertices counted once per rank)\n";
    Print() << "  surf elements  : " << r.nElements << "\n";
    Print() << "~ measure        : " << r.area
            << ((AMREX_SPACEDIM == 2) ? "  (arc length)\n" : "  (area)\n");
    Print() << "  checksum tris  : " << r.triChecksum << "\n";
    Print() << "  checksum edges : " << r.edgeParity << "\n";
    Print() << "  orientation    : " << (r.oriented() ? "consistent" : "INCONSISTENT")
            << "  (normals toward increasing f; every interior edge traversed"
            << " both ways)\n";
    if (ParallelDescriptor::NProcs() == 1) {
        Print() << "~ misoriented boundary items : " << r.orientBad << "\n";
    }
#if (AMREX_SPACEDIM == 3)
    Print() << "~ enclosed volume: " << r.signedVolume
            << "  (divergence theorem; only meaningful if the surface does not"
            << " cross a periodic face)\n";
#endif
    Print() << "  closed         : " << (r.closed() ? "yes" : "NO")
            << "  (every boundary item shared by exactly two elements)\n";
    if (ParallelDescriptor::NProcs() == 1) {
        // Marked "~" like every other non-invariant: these need the
        // deduplicated vertex and boundary-item counts, which are only exact
        // at one rank, so they are PRESENT at one rank and ABSENT above it.
        // A harness comparing the block for equality must drop them, and a
        // line that appears conditionally is exactly as bad as one whose
        // value drifts.
        Print() << "~ boundary items : " << r.nEdges << "\n";
        Print() << "~ Euler          : " << r.eulerChar << "\n";
        Print() << "~ edge census (rank-local, exact only at one rank):\n";
        Print() << "~   interior (2)  : " << r.edgesInterior << "\n";
        Print() << "~   boundary (1)  : " << r.edgesBoundary << "\n";
        Print() << "~   non-manifold  : " << r.edgesBad << "\n";
    }
    Print() << std::endl;
}

//
// MEF output.
//
// Format, matching what surfMEFtoDAT and the *MEF tools read:
//   line   title
//   line   "X Y [Z] var1 var2 ..."
//   line   "<nElts> <nodesPerElt>"
//   binary FArrayBox over Box(0..nPts-1,0,0) with SPACEDIM+nComp components,
//          laid out NODE-major -- the fab is a container for the blob, not an
//          AMReX-ordered array
//   binary nElts*nodesPerElt ints, 1-BASED
//
void AMRContour::writeMEF (const std::string& outfile,
                           const std::string& label,
                           const Vector<std::string>& varnames) const
{
    const int nc  = m_cx.nComp();
    const int nnd = AMREX_SPACEDIM + nc;          // per-node reals
    const int NV  = SurfElement::NVERT;

    // Pack this rank's share: each vertex as (key, coords, values), each
    // element as the KEYS of its vertices. Keys, never indices -- an index is
    // meaningless on another rank, which is the mistake this file has already
    // made once.
    Vector<Long> vkey; vkey.reserve(std::size_t(2)*m_verts.size());
    Vector<Real> vdat; vdat.reserve(std::size_t(nnd)*m_verts.size());
    for (std::size_t i = 0; i < m_verts.size(); ++i) {
        vkey.push_back(Long(m_vertKey[i].first));
        vkey.push_back(Long(m_vertKey[i].second));
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { vdat.push_back(m_verts[i][d]); }
        for (int n = 0; n < nc; ++n) {
            vdat.push_back(m_vertData[std::size_t(i)*std::size_t(nc)+std::size_t(n)]);
        }
    }
    Vector<Long> ekey; ekey.reserve(std::size_t(2)*NV*m_tris.size());
    for (const auto& t : m_tris) {
        for (int k = 0; k < NV; ++k) {
            ekey.push_back(Long(m_vertKey[t.v[k]].first));
            ekey.push_back(Long(m_vertKey[t.v[k]].second));
        }
    }

#ifdef AMREX_USE_MPI
    if (ParallelDescriptor::NProcs() > 1) {
        auto gatherL = [] (Vector<Long>& v) {
            const int nproc = ParallelDescriptor::NProcs();
            const int ioproc = ParallelDescriptor::IOProcessorNumber();
            int cnt = int(v.size());
            Vector<int> cnts(nproc,0), offs(nproc,0);
            MPI_Gather(&cnt,1,MPI_INT,cnts.dataPtr(),1,MPI_INT,ioproc,
                       ParallelDescriptor::Communicator());
            Vector<Long> recv;
            if (ParallelDescriptor::IOProcessor()) {
                for (int i=1;i<nproc;++i) { offs[i]=offs[i-1]+cnts[i-1]; }
                recv.resize(std::size_t(offs[nproc-1]+cnts[nproc-1]));
            }
            MPI_Gatherv(v.dataPtr(),cnt,MPI_LONG_LONG,recv.dataPtr(),
                        cnts.dataPtr(),offs.dataPtr(),MPI_LONG_LONG,ioproc,
                        ParallelDescriptor::Communicator());
            v.swap(recv);
        };
        auto gatherR = [] (Vector<Real>& v) {
            const int nproc = ParallelDescriptor::NProcs();
            const int ioproc = ParallelDescriptor::IOProcessorNumber();
            int cnt = int(v.size());
            Vector<int> cnts(nproc,0), offs(nproc,0);
            MPI_Gather(&cnt,1,MPI_INT,cnts.dataPtr(),1,MPI_INT,ioproc,
                       ParallelDescriptor::Communicator());
            Vector<Real> recv;
            if (ParallelDescriptor::IOProcessor()) {
                for (int i=1;i<nproc;++i) { offs[i]=offs[i-1]+cnts[i-1]; }
                recv.resize(std::size_t(offs[nproc-1]+cnts[nproc-1]));
            }
            MPI_Gatherv(v.dataPtr(),cnt,ParallelDescriptor::Mpi_typemap<Real>::type(),
                        recv.dataPtr(),cnts.dataPtr(),offs.dataPtr(),
                        ParallelDescriptor::Mpi_typemap<Real>::type(),ioproc,
                        ParallelDescriptor::Communicator());
            v.swap(recv);
        };
        gatherL(vkey);
        gatherR(vdat);
        gatherL(ekey);
    }
#endif

    if (!ParallelDescriptor::IOProcessor()) { return; }

    // Merge duplicates by key and order the global node list by key. Both the
    // merge and the ordering are exact and rank-count independent: no
    // tolerance, and no dependence on who saw a vertex first.
    const std::size_t nrecv = vkey.size()/2;
    std::map<std::pair<Long,Long>, std::size_t> first;
    for (std::size_t i = 0; i < nrecv; ++i) {
        first.emplace(std::make_pair(vkey[2*i], vkey[2*i+1]), i);
    }
    std::map<std::pair<Long,Long>, long> renum;
    long idx = 0;
    for (const auto& kv : first) { renum[kv.first] = idx++; }
    const long nPts = idx;

    Vector<Real> nodes(std::size_t(nPts)*std::size_t(nnd));
    for (const auto& kv : first) {
        const long j = renum[kv.first];
        const std::size_t src = kv.second;
        for (int k = 0; k < nnd; ++k) {
            nodes[std::size_t(j)*std::size_t(nnd)+std::size_t(k)]
                = vdat[src*std::size_t(nnd)+std::size_t(k)];
        }
    }

    const long nElts = long(ekey.size()/(2*NV));
    // Note: elts(std::size_t(nElts)) here would be a function declaration.
    Vector<std::array<int,SurfElement::NVERT>> elts;
    elts.resize(std::size_t(nElts));
    for (long e = 0; e < nElts; ++e) {
        for (int k = 0; k < NV; ++k) {
            const std::size_t o = std::size_t(e)*std::size_t(2*NV)+std::size_t(2*k);
            const auto it = renum.find(std::make_pair(ekey[o], ekey[o+1]));
            AMREX_ALWAYS_ASSERT(it != renum.end());
            elts[std::size_t(e)][k] = int(it->second);
        }
    }

    // Canonicalise the ELEMENT order too. Nodes were already sorted by key,
    // but the element list arrives in per-rank traversal order and was
    // written in gather order, so the connectivity block differed between
    // rank counts while the node block matched exactly.
    //
    // Sorting is on the index tuple AS WRITTEN, which leaves each element's
    // own vertex order -- and therefore its orientation -- untouched. Every
    // surface element is produced by exactly one cell, so the tuple is unique
    // and this is a total order.
    std::sort(elts.begin(), elts.end());

    Vector<int> conn(std::size_t(nElts)*std::size_t(NV));
    for (long e = 0; e < nElts; ++e) {
        for (int k = 0; k < NV; ++k) {
            conn[std::size_t(e)*std::size_t(NV)+std::size_t(k)]
                = elts[std::size_t(e)][k] + 1;      // 1-based
        }
    }

    std::string vars = "X Y";
#if (AMREX_SPACEDIM == 3)
    vars += " Z";
#endif
    for (const auto& v : varnames) { vars += " " + v; }

    // The fab is only a carrier for the blob; the reader indexes it
    // node-major, so it is filled node-major here too.
    FArrayBox fab(Box(IntVect::TheZeroVector(),
                      IntVect(AMREX_D_DECL(int(nPts)-1,0,0))), nnd);
    std::memcpy(fab.dataPtr(), nodes.dataPtr(), nodes.size()*sizeof(Real));

    std::ofstream ofs(outfile.c_str(),
                      std::ios::out|std::ios::trunc|std::ios::binary);
    ofs << label << std::endl;
    ofs << vars  << std::endl;
    ofs << nElts << " " << NV << std::endl;
    fab.writeOn(ofs);
    ofs.write((char*)conn.dataPtr(), sizeof(int)*conn.size());
    ofs.close();

    Print() << "  wrote " << outfile << ": " << nPts << " nodes, "
            << nElts << " elements, " << nnd << " components per node\n";
}

} // namespace amrcomplex
