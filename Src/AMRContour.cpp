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


int AMRContour::tetrahedralize (const Node (&cn)[Element::NNODE], int tet[6][4],
                                bool minGidDiagonal)
{
#if (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(cn, tet);
    return 0;
#else
    int imin = 0;
    for (int c = 1; c < 8; ++c) { if (cn[c].gid < cn[imin].gid) { imin = c; } }

    int n = 0;
    for (int f = 0; f < 6; ++f) {

        bool hasApex = false;
        for (int c = 0; c < 4; ++c) { if (HFACE[f][c] == imin) { hasApex = true; } }
        if (hasApex) { continue; }   // coning to it would give zero volume

        // Diagonal through the face's own minimum-gid corner.
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

        const int cand[2][4] = {{imin,t0[0],t0[1],t0[2]},
                                {imin,t1[0],t1[1],t1[2]}};
        for (int k = 0; k < 2; ++k) {
            // A stitch cell has repeated corners, so some of these tets
            // collapse. Dropping them here is what turns the degenerate
            // hexahedron into the pyramid or prism it actually is.
            bool degen = false;
            for (int a = 0; a < 4 && !degen; ++a) {
                for (int b = a+1; b < 4; ++b) {
                    if (cn[cand[k][a]].gid == cn[cand[k][b]].gid) { degen = true; break; }
                }
            }
            if (degen) { continue; }
            for (int a = 0; a < 4; ++a) { tet[n][a] = cand[k][a]; }
            ++n;
        }
    }
    return n;
#endif
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

    const auto& g0 = cx.geom(0);
    const auto  plo = g0.ProbLo();
    const auto  phi = g0.ProbHi();
    RealVect L;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { L[d] = phi[d] - plo[d]; }

    cx.forEachElement([&] (const Node (&cn)[Element::NNODE], const Real* vals)
    {
#if (AMREX_SPACEDIM == 3)
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
                Real s = p[d] - o[d];
                if (g0.isPeriodic(d)) {
                    if (s >  Real(0.5)*L[d]) { s -= L[d]; }
                    if (s < -Real(0.5)*L[d]) { s += L[d]; }
                }
                q[c][d] = s;
            }
            fv[c] = vals[std::size_t(c)*std::size_t(nc)+std::size_t(comp)];
        }

        int tet[6][4];
        const int nt = tetrahedralize(cn, tet, m_minGidDiag);

        for (int t = 0; t < nt; ++t) {
            int above[4], below[4];
            int na = 0, nb = 0;
            for (int a = 0; a < 4; ++a) {
                const int c = tet[t][a];
                if (fv[c] > m_isoVal) { above[na++] = c; } else { below[nb++] = c; }
            }
            if (na == 0 || nb == 0) { continue; }

            long v[4];
            int nv = 0;
            if (na == 1 || nb == 1) {
                const int apex  = (na == 1) ? above[0] : below[0];
                const int* rest = (na == 1) ? below : above;
                for (int k = 0; k < 3; ++k) {
                    v[nv++] = edgeVertex(cn[apex].gid, cn[rest[k]].gid,
                                         q[apex], q[rest[k]],
                                         fv[apex], fv[rest[k]], o);
                }
            } else {
                // Two-two split: the four crossings form a quad. The corner
                // order below is cyclic around it, so the two triangles do not
                // overlap.
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

            const int TRI[2][3] = {{0,1,2},{0,2,3}};
            const int ntri = (nv == 3) ? 1 : 2;
            for (int k = 0; k < ntri; ++k) {
                Triangle tri;
                for (int a = 0; a < 3; ++a) { tri.v[a] = v[TRI[k][a]]; }
                if (tri.v[0]==tri.v[1] || tri.v[1]==tri.v[2] || tri.v[0]==tri.v[2]) {
                    continue;
                }
                m_tris.push_back(tri);

                // Area in the element frame, from the same interpolated
                // positions, before any wrapping.
                RealVect w[3];
                for (int a = 0; a < 3; ++a) {
                    const RealVect& pw = m_verts[tri.v[a]];
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                        Real s = pw[d] - o[d];
                        if (g0.isPeriodic(d)) {
                            if (s >  Real(0.5)*L[d]) { s -= L[d]; }
                            if (s < -Real(0.5)*L[d]) { s += L[d]; }
                        }
                        w[a][d] = s;
                    }
                }
                const RealVect e1{AMREX_D_DECL(w[1][0]-w[0][0], w[1][1]-w[0][1], w[1][2]-w[0][2])};
                const RealVect e2{AMREX_D_DECL(w[2][0]-w[0][0], w[2][1]-w[0][1], w[2][2]-w[0][2])};
                const Real cx0 = e1[1]*e2[2] - e1[2]*e2[1];
                const Real cy0 = e1[2]*e2[0] - e1[0]*e2[2];
                const Real cz0 = e1[0]*e2[1] - e1[1]*e2[0];
                m_area += Real(0.5)*std::sqrt(cx0*cx0 + cy0*cy0 + cz0*cz0);
            }
        }
#else
        amrex::ignore_unused(cn, vals, L);
#endif
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

    for (const auto& t : m_tris) {
        for (int k = 0; k < 3; ++k) {
            long a = t.v[k], b = t.v[(k+1)%3];
            if (a > b) { std::swap(a,b); }
            ++edgeCount[std::make_pair(a,b)];
            // Parity over the global names, so an edge shared between two
            // ranks cancels exactly as one shared within a rank does.
            par ^= edgeHash(m_vertKey[t.v[k]], m_vertKey[t.v[(k+1)%3]]);
        }

        std::pair<long,long> kk[3] = {m_vertKey[t.v[0]],
                                      m_vertKey[t.v[1]],
                                      m_vertKey[t.v[2]]};
        std::sort(kk, kk+3);
        unsigned long long th = 1469598103934665603ULL;
        for (int k = 0; k < 3; ++k) {
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
    r.nTriangles = nt;
    r.nEdges     = ne;
    r.eulerChar  = nv - ne + nt;
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
    Print() << "  triangles      : " << r.nTriangles << "\n";
    Print() << "~ area           : " << r.area << "\n";
    Print() << "  checksum tris  : " << r.triChecksum << "\n";
    Print() << "  checksum edges : " << r.edgeParity << "\n";
    Print() << "  closed surface : " << (r.closed() ? "yes" : "NO")
            << "  (every edge shared by exactly two triangles)\n";
    if (ParallelDescriptor::NProcs() == 1) {
        Print() << "  edges          : " << r.nEdges << "\n";
        Print() << "  Euler V-E+F    : " << r.eulerChar << "\n";
        Print() << "~ edge census (rank-local, exact only at one rank):\n";
        Print() << "~   interior (2)  : " << r.edgesInterior << "\n";
        Print() << "~   boundary (1)  : " << r.edgesBoundary << "\n";
        Print() << "~   non-manifold  : " << r.edgesBad << "\n";
    }
    Print() << std::endl;
}

} // namespace amrcomplex
