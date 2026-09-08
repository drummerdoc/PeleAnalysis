#include "AMRContour.H"

#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <cmath>

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

long AMRContour::edgeVertex (long ga, long gb,
                             const RealVect& pa, const RealVect& pb,
                             Real fa, Real fb,
                             const RealVect& origin)
{
    std::pair<long,long> key = (ga < gb) ? std::make_pair(ga,gb)
                                         : std::make_pair(gb,ga);
    auto it = m_vertCache.find(key);
    if (it != m_vertCache.end()) { return it->second; }

    // Interpolate in the direction the key was canonicalised, so the position
    // does not depend on which element reached the edge first.
    const RealVect& p0 = (ga < gb) ? pa : pb;
    const RealVect& p1 = (ga < gb) ? pb : pa;
    const Real      v0 = (ga < gb) ? fa : fb;
    const Real      v1 = (ga < gb) ? fb : fa;

    const Real den = v1 - v0;
    const Real t   = (den == Real(0.0)) ? Real(0.5) : (m_isoVal - v0)/den;

    RealVect p;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { p[d] = p0[d] + t*(p1[d] - p0[d]); }

    // Stored positions are absolute and wrapped back into the domain, so a
    // vertex on a wrap edge has one canonical location. Areas are NOT computed
    // from these -- they are computed in the element's own frame, below, where
    // the geometry is unbroken.
    const auto& g   = m_cx.geom(0);
    const auto  plo = g.ProbLo();
    const auto  phi = g.ProbHi();
    RealVect pw;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        Real x = origin[d] + p[d];
        if (g.isPeriodic(d)) {
            const Real L = phi[d] - plo[d];
            while (x <  plo[d]) { x += L; }
            while (x >= phi[d]) { x -= L; }
        }
        pw[d] = x;
    }

    const long idx = long(m_verts.size());
    m_verts.push_back(pw);
    m_vertKey.push_back(key);
    m_vertCache[key] = idx;
    return idx;
}

AMRContour::AMRContour (const AMRComplex& cx, int comp, Real isoVal,
                        bool minGidDiagonal)
    : m_cx(cx), m_comp(comp), m_isoVal(isoVal), m_minGidDiag(minGidDiagonal)
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
                    v[nv++] = edgeVertex(cn[apex].gid, cn[rest[k]].gid,
                                         q[apex], q[rest[k]],
                                         fv[apex], fv[rest[k]], o);
                }
            }
#else
            if (na == 1 || nb == 1) {
                const int apex  = (na == 1) ? above[0] : below[0];
                const int* rest = (na == 1) ? below : above;
                for (int k = 0; k < 3; ++k) {
                    v[nv++] = edgeVertex(cn[apex].gid, cn[rest[k]].gid,
                                         q[apex], q[rest[k]],
                                         fv[apex], fv[rest[k]], o);
                }
            } else {
                // Two-two split: the four crossings form a quad, listed
                // cyclically so the two triangles do not overlap.
                const int a0 = above[0], a1 = above[1];
                const int b0 = below[0], b1 = below[1];
                const int ea[4] = {a0,a0,a1,a1};
                const int eb[4] = {b0,b1,b1,b0};
                for (int k = 0; k < 4; ++k) {
                    v[nv++] = edgeVertex(cn[ea[k]].gid, cn[eb[k]].gid,
                                         q[ea[k]], q[eb[k]],
                                         fv[ea[k]], fv[eb[k]], o);
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
                m_tris.push_back(se);

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
    for (const auto& t : m_tris) {
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
    r.triChecksum = triSum % CHECKSUM_MOD;
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

} // namespace amrcomplex
