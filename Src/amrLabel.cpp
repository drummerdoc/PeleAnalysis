//
// Segmentation written as a plotfile: one integer component id per cell, on
// exactly the input's grids, for selected thresholds.
//
// The complex's nodes ARE the input's cell centres, so a segmentation is
// cell-centred data on the input hierarchy with nothing to interpolate. Each
// rank fills its own boxes and the plotfile is written in parallel. Because
// the grids are identical, the output composes with the existing tools:
// combinePlts merges it with the source fields, conditionalMean conditions on
// it, and a viewer overlays it directly.
//
//   amrLabel3d.gnu.MPI.ex infile=plt isoCompName=HeatRelease is_per="1 1 0" \
//       thresholds="4.35e10 4.27e11 6.88e11" above="0 0 1" outfile=plt_seg
//
// Values written, per cut:
//    id >= 1  component id; 1 is the largest by volume, then descending, ties
//             broken by the minimum gid in the component. Deterministic and
//             independent of the rank count. Compact rather than the gid
//             itself: a gid reaches ~7e8 on a production hierarchy, beyond
//             the 2^24 a single-precision plotfile holds exactly.
//    0        not in the set. {f > v} for above=1, {f < v} for above=0, so a
//             cell exactly at v is in neither.
//    -1       a covered cell whose finer children disagree. A covered cell
//             whose children all agree carries their value; labels cannot be
//             averaged, so this replaces average-down.
//
// Beside the plotfile, components_<name>.txt per cut: id, the gid label,
// volume, cell count per level, and whether the component touches a
// non-periodic wall. Whether it wraps a periodic direction is not yet
// computed; that needs offsets carried through the union-find.
//
// Surfaces (surfaces=1): the isosurface at each threshold as VTK polydata,
// <surf_dir>/<name>.pvtp indexing one binary <name>/piece_NNNNN.vtp per rank,
// written in parallel with no gather. Every triangle carries the id of the
// component it bounds -- the same id as in the plotfile and the table -- and
// that component's volume. surf_skip (one value, or one per threshold) leaves
// out the N largest components, so surf_skip=1 on a flame sheet keeps only
// the satellites. The table gains each component's surface area. labels=0
// writes the surfaces and tables alone. Open the .pvtp in ParaView or VisIt.
//
//   amrLabel3d.gnu.MPI.ex infile=plt isoCompName=HeatRelease is_per="1 1 0" \
//       thresholds="8.7e9 8.7e10" labels=0 surfaces=1 surf_skip=1 \
//       outfile=plt_seg surf_dir=plt_surf
//
#include "AMRSegment.H"
#include "AMRContour.H"
#include "AMRParallel.H"

#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include <AMReX_Utility.H>

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <map>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <unordered_map>
#include <vector>

using namespace amrex;
using namespace amrcomplex;

namespace {

std::string fmtThreshold (Real v)
{
    std::ostringstream os;
    os << std::setprecision(4) << std::scientific << v;
    return os.str();
}

std::string fileSafe (std::string s)
{
    for (char& c : s) {
        if (!std::isalnum(static_cast<unsigned char>(c)) && c != '.' && c != '-' && c != '+') { c = '_'; }
    }
    return s;
}

Real elapsedMax (Real t0)
{
    Real t = ParallelDescriptor::second() - t0;
    ParallelDescriptor::ReduceRealMax(t);
    return t;
}

// One block of a VTK XML "appended raw" section: a UInt64 byte count, then
// the bytes.
template <class T>
void appendBlock (std::string& blob, const std::vector<T>& v)
{
    const std::uint64_t n = std::uint64_t(v.size() * sizeof(T));
    blob.append(reinterpret_cast<const char*>(&n), sizeof(n));
    if (n > 0) { blob.append(reinterpret_cast<const char*>(v.data()), std::size_t(n)); }
}

//
// The isosurface at iso, tagged by component, as one .vtp piece per rank and
// a .pvtp index. Triangles only ever exist on the rank that owns their
// element, so pieces never overlap and nothing is gathered.
//
// The tag is exact, not interpolated. The contour's rule is "above" means
// f > iso, the segmentation's too, so every surface vertex lies on an edge
// with exactly one endpoint in the set, and that endpoint's label is the
// component the triangle bounds. For a below cut (f < iso) an endpoint lying
// exactly at iso is in neither; such a triangle takes its tag from another
// vertex, and one with no tagged vertex at all is counted as untagged.
//
// Coordinates are unwrapped per triangle around its first vertex, so a
// triangle crossing a periodic seam is drawn whole just outside the domain
// rather than stretched across it.
//
void writeSurface (const AMRComplex& cx, const AMRSegment& seg,
                   const std::unordered_map<long,int>& idOf,
                   const std::vector<AMRSegment::Component>& byVol,
                   Real iso, int skip, int maxId, const std::string& dir, const std::string& name,
                   const Array<int,AMREX_SPACEDIM>& isper,
                   Vector<Real>& areaOf, long& nWritten, long& nUntagged)
{
    AMRContour ct(cx, 0, iso);
    const auto& verts = ct.vertices();
    const auto& keys  = ct.vertexKeys();
    const auto& tris  = ct.surfElements();
    constexpr int NV  = SurfElement::NVERT;

    const auto plo = cx.geom(0).ProbLo();
    const auto phi = cx.geom(0).ProbHi();
    RealVect L;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { L[d] = phi[d] - plo[d]; }

    std::vector<float>        pts;
    std::vector<std::int64_t> conn, offs;
    std::vector<std::int32_t> cid;
    std::vector<float>        cvol;
    std::map<std::array<long,4>, std::int64_t> vmap;   // (vertex, image shift) -> output index

    nWritten = 0;
    nUntagged = 0;
    for (const auto& t : tris) {
        long lab = -1;
        for (int v = 0; v < NV && lab < 0; ++v) {
            const auto& k = keys[std::size_t(t.v[v])];
            lab = seg.labelAny(k.first);
            if (lab < 0) { lab = seg.labelAny(k.second); }
        }
        if (lab < 0) { ++nUntagged; continue; }
        const int id = idOf.at(lab);

        // Unwrap around the first vertex.
        const RealVect& p0 = verts[std::size_t(t.v[0])];
        RealVect q[NV];
        std::array<long,4> key[NV];
        for (int v = 0; v < NV; ++v) {
            q[v] = verts[std::size_t(t.v[v])];
            key[v] = {long(t.v[v]), 0, 0, 0};
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                if (!isper[d]) { continue; }
                const long s = std::lround((q[v][d] - p0[d]) / L[d]);
                q[v][d] -= Real(s) * L[d];
                key[v][std::size_t(d+1)] = s;
            }
        }

#if (AMREX_SPACEDIM == 2)
        const Real a = (q[1] - q[0]).vectorLength();
#else
        const RealVect e1 = q[1] - q[0], e2 = q[2] - q[0];
        const Real a = Real(0.5) * e1.crossProduct(e2).vectorLength();
#endif
        areaOf[std::size_t(id)] += a;

        if (id <= skip || (maxId > 0 && id > maxId)) { continue; }
        for (int v = 0; v < NV; ++v) {
            auto it = vmap.find(key[v]);
            if (it == vmap.end()) {
                it = vmap.emplace(key[v], std::int64_t(pts.size() / 3)).first;
                pts.push_back(float(q[v][0]));
                pts.push_back(float(q[v][1]));
                pts.push_back(AMREX_SPACEDIM == 3 ? float(q[v][AMREX_SPACEDIM-1]) : 0.0f);
            }
            conn.push_back(it->second);
        }
        offs.push_back(std::int64_t(conn.size()));
        cid.push_back(std::int32_t(id));
        cvol.push_back(float(byVol[std::size_t(id-1)].volume));
        ++nWritten;
    }

    ParallelDescriptor::ReduceRealSum(areaOf.dataPtr(), areaOf.size());

    // One piece per rank that has triangles.
    const std::string pieceDir = dir + "/" + name;
    if (ParallelDescriptor::IOProcessor()) { UtilCreateCleanDirectory(pieceDir, false); }
    ParallelDescriptor::Barrier();

    const int me = ParallelDescriptor::MyProc();
    char pname[32];
    std::snprintf(pname, sizeof(pname), "piece_%05d.vtp", me);
    const char* cellTag = (AMREX_SPACEDIM == 3) ? "Polys" : "Lines";
    if (nWritten > 0) {
        std::string blob;
        std::ostringstream x;
        std::uint64_t off = 0;
        auto next = [&] (std::size_t bytes) { const std::uint64_t o = off; off += 8 + bytes; return o; };
        const std::uint64_t oPts  = next(pts.size()  * sizeof(float));
        const std::uint64_t oCid  = next(cid.size()  * sizeof(std::int32_t));
        const std::uint64_t oVol  = next(cvol.size() * sizeof(float));
        const std::uint64_t oConn = next(conn.size() * sizeof(std::int64_t));
        const std::uint64_t oOffs = next(offs.size() * sizeof(std::int64_t));
        appendBlock(blob, pts);
        appendBlock(blob, cid);
        appendBlock(blob, cvol);
        appendBlock(blob, conn);
        appendBlock(blob, offs);

        const long np = long(pts.size() / 3);
        x << "<?xml version=\"1.0\"?>\n"
          << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"
          << " <PolyData>\n"
          << "  <Piece NumberOfPoints=\"" << np << "\" NumberOfVerts=\"0\" NumberOfStrips=\"0\" "
          << (AMREX_SPACEDIM == 3 ? "NumberOfLines=\"0\" NumberOfPolys=\"" : "NumberOfPolys=\"0\" NumberOfLines=\"")
          << nWritten << "\">\n"
          << "   <Points><DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << oPts << "\"/></Points>\n"
          << "   <CellData Scalars=\"component\">\n"
          << "    <DataArray type=\"Int32\" Name=\"component\" format=\"appended\" offset=\"" << oCid << "\"/>\n"
          << "    <DataArray type=\"Float32\" Name=\"volume\" format=\"appended\" offset=\"" << oVol << "\"/>\n"
          << "   </CellData>\n"
          << "   <" << cellTag << ">\n"
          << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"appended\" offset=\"" << oConn << "\"/>\n"
          << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"appended\" offset=\"" << oOffs << "\"/>\n"
          << "   </" << cellTag << ">\n"
          << "  </Piece>\n"
          << " </PolyData>\n"
          << " <AppendedData encoding=\"raw\">\n_";
        std::ofstream o(pieceDir + "/" + pname, std::ios::binary);
        const std::string head = x.str();
        o.write(head.data(), std::streamsize(head.size()));
        o.write(blob.data(), std::streamsize(blob.size()));
        o << "\n </AppendedData>\n</VTKFile>\n";
    }

    Vector<Long> counts(1, Long(nWritten));
    allGather(counts);
    Long nu = nUntagged, nw = nWritten;
    ParallelDescriptor::ReduceLongSum(nu);
    ParallelDescriptor::ReduceLongSum(nw);
    nUntagged = long(nu);
    nWritten  = long(nw);

    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream o(dir + "/" + name + ".pvtp");
        o << "<?xml version=\"1.0\"?>\n"
          << "<VTKFile type=\"PPolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"
          << " <PPolyData GhostLevel=\"0\">\n"
          << "  <PPoints><PDataArray type=\"Float32\" NumberOfComponents=\"3\"/></PPoints>\n"
          << "  <PCellData Scalars=\"component\">\n"
          << "   <PDataArray type=\"Int32\" Name=\"component\"/>\n"
          << "   <PDataArray type=\"Float32\" Name=\"volume\"/>\n"
          << "  </PCellData>\n";
        for (int r = 0; r < int(counts.size()); ++r) {
            if (counts[r] == 0) { continue; }
            std::snprintf(pname, sizeof(pname), "piece_%05d.vtp", r);
            o << "  <Piece Source=\"" << name << "/" << pname << "\"/>\n";
        }
        o << " </PPolyData>\n</VTKFile>\n";
    }
}

} // namespace

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        ParmParse pp;

        AMRComplex::Params p;
        pp.get("infile", p.infile);
        pp.query("finestLevel", p.finestLevel);

        std::string var;
        pp.get("isoCompName", var);
        p.varnames = {var};

        Vector<int> isper(AMREX_SPACEDIM, 0);
        pp.queryarr("is_per", isper, 0, AMREX_SPACEDIM);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { p.is_per[d] = isper[d]; }

        const int ncut = pp.countval("thresholds");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ncut > 0, "give thresholds=\"v1 v2 ...\"");
        Vector<Real> thr;
        pp.getarr("thresholds", thr, 0, ncut);
        Vector<int> above(ncut, 1);
        if (pp.countval("above") > 0) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pp.countval("above") == ncut,
                                             "above needs one entry per threshold");
            pp.getarr("above", above, 0, ncut);
        }

        std::string outfile;
        pp.get("outfile", outfile);
        {
            // Never write over the input.
            auto strip = [] (std::string x) { while (x.size() > 1 && x.back() == '/') { x.pop_back(); } return x; };
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(strip(outfile) != strip(p.infile),
                                             "outfile must differ from infile");
        }

        // What to write. labels: the plotfile of ids. surfaces: the isosurface
        // at each threshold as a parallel VTK polydata set (.pvtp + one .vtp
        // per rank), each triangle tagged with the id of the component it
        // bounds. surf_skip=N leaves out the N largest components -- the
        // sheet, or the two bulks of a below cut -- so what remains is the
        // satellites.
        int doLabels = 1;   pp.query("labels", doLabels);
        int doSurf   = 0;   pp.query("surfaces", doSurf);
        // One value for all cuts, or one per threshold.
        Vector<int> surfSkip(ncut, 0);
        if (const int ns = pp.countval("surf_skip")) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ns == 1 || ns == ncut,
                                             "surf_skip takes one value or one per threshold");
            Vector<int> v;
            pp.getarr("surf_skip", v, 0, ns);
            for (int k = 0; k < ncut; ++k) { surfSkip[k] = v[ns == 1 ? 0 : k]; }
        }
        // surf_max=N writes only the N largest components (0 = all) -- with
        // finestLevel set low, a cheap context surface for the satellites.
        Vector<int> surfMax(ncut, 0);
        if (const int ns = pp.countval("surf_max")) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ns == 1 || ns == ncut,
                                             "surf_max takes one value or one per threshold");
            Vector<int> v;
            pp.getarr("surf_max", v, 0, ns);
            for (int k = 0; k < ncut; ++k) { surfMax[k] = v[ns == 1 ? 0 : k]; }
        }
        // Appended to every surface and table name from this run, so runs
        // sharing a surf_dir -- full resolution and context -- cannot collide.
        std::string surfSuffix;
        pp.query("surf_suffix", surfSuffix);
        std::string surfDir = outfile + "_surf";
        pp.query("surf_dir", surfDir);
        {
            // UtilCreateCleanDirectory renames whatever is already there, so
            // the surface directory must not be the input either.
            auto strip = [] (std::string x) { while (x.size() > 1 && x.back() == '/') { x.pop_back(); } return x; };
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(strip(surfDir) != strip(p.infile),
                                             "surf_dir must differ from infile");
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(doLabels || doSurf, "nothing to write");

        Real t0 = ParallelDescriptor::second();
        AMRComplex cx(p);
        const auto cr = cx.verify(false);
        if (!cr.passed()) { cx.reportVerify(cr); Abort("complex does not verify"); }
        Print() << "# complex: " << cx.nLevels() << " levels, " << cr.nElements
                << " elements, built and verified in " << elapsedMax(t0) << " s\n";

        const int nlev = cx.nLevels();
        Vector<MultiFab> out(nlev);
        if (doLabels) {
            for (int lev = 0; lev < nlev; ++lev) {
                out[lev].define(cx.boxArray(lev), cx.distributionMap(lev), ncut, 0);
                out[lev].setVal(0.0);
            }
        }
        // Create, never clean: runs are meant to share a surf_dir (a context
        // surface beside full-resolution satellites), and cleaning would move
        // the earlier run's output aside. Each surface's own piece directory
        // is cleaned when it is written.
        if (doSurf && ParallelDescriptor::IOProcessor()) {
            if (!UtilCreateDirectory(surfDir, 0755)) { CreateDirectoryFailed(surfDir); }
        }
        ParallelDescriptor::Barrier();

        Vector<std::string> names(ncut);
        std::vector<std::string> tables(ncut);

        for (int k = 0; k < ncut; ++k) {
            t0 = ParallelDescriptor::second();
            AMRSegment seg(cx, 0, thr[k], above[k] != 0, doSurf != 0);

            names[k] = "seg_" + var + "_" + (above[k] ? "above" : "below") + "_" + fmtThreshold(thr[k]);

            // Compact ids: 1 = largest by volume; ties broken by label.
            std::vector<AMRSegment::Component> byVol(seg.components().begin(), seg.components().end());
            std::sort(byVol.begin(), byVol.end(), [] (const auto& a, const auto& b) {
                return (a.volume != b.volume) ? (a.volume > b.volume) : (a.label < b.label);
            });
            std::unordered_map<long,int> idOf;
            idOf.reserve(byVol.size());
            for (std::size_t i = 0; i < byVol.size(); ++i) { idOf[byVol[i].label] = int(i) + 1; }

            const long ncomp = long(byVol.size());
            Vector<Long> perLev;               // cells per (id, level), global
            perLev.resize(std::size_t(ncomp*nlev), 0);
            Vector<Long> wallLabels;           // ids touching a non-periodic wall

            // Valid cells. Covered cells are marked and resolved below.
            constexpr Real COVERED_MARK = -2.0;
            for (int lev = 0; lev < nlev; ++lev) {
                const Box& dom = cx.geom(lev).Domain();
                for (MFIter mfi(cx.nodes(lev)); mfi.isValid(); ++mfi) {
                    const Box& vbx = mfi.validbox();
                    auto const& nd = cx.nodes(lev).const_array(mfi);
                    Array4<Real> o;
                    if (doLabels) { o = out[lev].array(mfi); }
                    AMREX_LOOP_3D(vbx, i, j, k2,
                    {
                        const Node& n = nd(i,j,k2);
                        if (n.type == Node::COVERED) {
                            if (doLabels) { o(i,j,k2,k) = COVERED_MARK; }
                            continue;
                        }
                        const long lab = seg.label(n.gid);
                        if (lab < 0) {
                            if (doLabels) { o(i,j,k2,k) = 0.0; }
                            continue;
                        }
                        const int id = idOf.at(lab);
                        if (doLabels) { o(i,j,k2,k) = Real(id); }
                        ++perLev[std::size_t((id-1)*nlev + lev)];
                        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                            if (p.is_per[d]) { continue; }
                            if (n.iv[d] == dom.smallEnd(d) || n.iv[d] == dom.bigEnd(d)) {
                                wallLabels.push_back(Long(id));
                                break;
                            }
                        }
                    });
                }
            }
            ParallelDescriptor::ReduceLongSum(perLev.dataPtr(), perLev.size());
            std::sort(wallLabels.begin(), wallLabels.end());
            wallLabels.erase(std::unique(wallLabels.begin(), wallLabels.end()), wallLabels.end());
            allGather(wallLabels);
            std::vector<char> atWall(std::size_t(ncomp) + 1, 0);
            for (Long id : wallLabels) { atWall[std::size_t(id)] = 1; }

            // Covered cells, finest first: the common value of the 2^DIM
            // children if they agree, else -1.
            for (int lev = nlev-2; lev >= 0 && doLabels; --lev) {
                const IntVect rr = cx.refRatio(lev);
                BoxArray cba = cx.boxArray(lev+1);
                cba.coarsen(rr);
                MultiFab mm(cba, cx.distributionMap(lev+1), 2, 0);
                for (MFIter mfi(mm); mfi.isValid(); ++mfi) {
                    const Box& cbx = mfi.validbox();
                    auto const& f = out[lev+1].const_array(mfi);
                    auto const& m = mm.array(mfi);
                    AMREX_LOOP_3D(cbx, ic, jc, kc,
                    {
                        Real lo =  1e300;
                        Real hi = -1e300;
                        const IntVect civ(AMREX_D_DECL(ic,jc,kc));
                        const Box kids = Box(civ, civ).refine(rr);
                        const auto kl = lbound(kids);
                        const auto kh = ubound(kids);
                        for (int k2 = kl.z; k2 <= kh.z; ++k2) {
                        for (int j  = kl.y; j  <= kh.y; ++j ) {
                        for (int i  = kl.x; i  <= kh.x; ++i ) {
                            lo = std::min(lo, f(i,j,k2,k));
                            hi = std::max(hi, f(i,j,k2,k));
                        }}}
                        m(ic,jc,kc,0) = lo;
                        m(ic,jc,kc,1) = hi;
                    });
                }
                MultiFab cm(cx.boxArray(lev), cx.distributionMap(lev), 2, 0);
                cm.setVal(0.0);
                cm.ParallelCopy(mm, 0, 0, 2);
                for (MFIter mfi(out[lev]); mfi.isValid(); ++mfi) {
                    const Box& vbx = mfi.validbox();
                    auto const& o = out[lev].array(mfi);
                    auto const& c = cm.const_array(mfi);
                    AMREX_LOOP_3D(vbx, i, j, k2,
                    {
                        if (o(i,j,k2,k) != COVERED_MARK) { continue; }
                        o(i,j,k2,k) = (c(i,j,k2,0) == c(i,j,k2,1)) ? c(i,j,k2,0) : Real(-1.0);
                    });
                }
            }

            // The isosurface, each triangle tagged with its component.
            Vector<Real> areaOf;
            areaOf.resize(std::size_t(ncomp) + 1, 0.0);
            if (doSurf) {
                const Real ts = ParallelDescriptor::second();
                long nWritten = 0, nUntagged = 0;
                writeSurface(cx, seg, idOf, byVol, thr[k], surfSkip[k], surfMax[k], surfDir,
                             names[k] + surfSuffix,
                             p.is_per, areaOf, nWritten, nUntagged);
                Print() << "#   surface: " << nWritten << " triangles written (components "
                        << surfSkip[k]+1 << ".." << (surfMax[k] > 0 ? std::min<long>(surfMax[k], ncomp) : ncomp)
                        << "), " << nUntagged
                        << " untagged, " << elapsedMax(ts) << " s\n";
            }

            // The per-component table.
            std::ostringstream tb;
            tb << "# " << names[k] << "\n";
            tb << "# source " << p.infile << ", " << var << (above[k] ? " > " : " < ")
               << std::setprecision(10) << thr[k] << ", is_per";
            for (int d = 0; d < AMREX_SPACEDIM; ++d) { tb << " " << p.is_per[d]; }
            tb << "\n# wraps a periodic direction: not computed\n";
            tb << "# id label volume" << (doSurf ? " area" : "") << " nCells";
            for (int lev = 0; lev < nlev; ++lev) { tb << " cellsL" << lev; }
            tb << " touchesWall\n" << std::setprecision(10);
            for (long i = 0; i < ncomp; ++i) {
                const auto& c = byVol[std::size_t(i)];
                tb << (i+1) << " " << c.label << " " << c.volume;
                if (doSurf) { tb << " " << areaOf[std::size_t(i+1)]; }
                tb << " " << c.nNodes;
                for (int lev = 0; lev < nlev; ++lev) { tb << " " << perLev[std::size_t(i*nlev + lev)]; }
                tb << " " << int(atWall[std::size_t(i+1)]) << "\n";
            }
            tables[k] = tb.str();

            long nWall = 0;
            for (long i = 1; i <= ncomp; ++i) { nWall += atWall[std::size_t(i)]; }
            Print() << "# " << names[k] << ": " << ncomp << " components, "
                    << nWall << " touching a wall, " << elapsedMax(t0) << " s\n";
        }

        // Tables go with the plotfile when there is one, else with the
        // surfaces.
        if (doLabels) {
            t0 = ParallelDescriptor::second();
            Vector<const MultiFab*> mfs(nlev);
            Vector<Geometry> geoms(nlev);
            Vector<int> steps(nlev, 0);
            Vector<IntVect> rr(std::max(nlev-1, 1), IntVect(1));
            for (int lev = 0; lev < nlev; ++lev) {
                mfs[lev]   = &out[lev];
                geoms[lev] = cx.geom(lev);
                if (lev < nlev-1) { rr[lev] = cx.refRatio(lev); }
            }
            WriteMultiLevelPlotfile(outfile, nlev, mfs, names, geoms, cx.time(), steps, rr);
            Print() << "# wrote " << outfile << " in " << elapsedMax(t0) << " s\n";
        }
        const std::string tableDir = doLabels ? outfile : surfDir;
        if (ParallelDescriptor::IOProcessor()) {
            for (int k = 0; k < ncut; ++k) {
                std::ofstream ofs(tableDir + "/components_" + fileSafe(names[k] + surfSuffix) + ".txt");
                ofs << tables[k];
            }
        }
    }
    amrex::Finalize();
    return 0;
}
