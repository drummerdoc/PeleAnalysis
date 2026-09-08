#include "AMRParallel.H"

#include <AMReX_ParallelDescriptor.H>

using namespace amrex;

namespace amrcomplex {

namespace {

// One implementation, two entry points per type. The plain C++ is here,
// unconditionally compiled; only the two MPI calls are guarded.
template <class T>
void gatherImpl (Vector<T>& v, bool toIOOnly)
{
    const int nproc = ParallelDescriptor::NProcs();
    if (nproc < 2) { return; }

    const int ioproc = ParallelDescriptor::IOProcessorNumber();
    const int cnt    = int(v.size());

    Vector<int> cnts, offs;
    cnts.resize(nproc, 0);
    offs.resize(nproc, 0);

    Vector<T> recv;

#ifdef AMREX_USE_MPI
    MPI_Datatype dt = ParallelDescriptor::Mpi_typemap<T>::type();

    if (toIOOnly) {
        MPI_Gather(const_cast<int*>(&cnt), 1, MPI_INT, cnts.dataPtr(), 1, MPI_INT,
                   ioproc, ParallelDescriptor::Communicator());
    } else {
        MPI_Allgather(const_cast<int*>(&cnt), 1, MPI_INT, cnts.dataPtr(), 1, MPI_INT,
                      ParallelDescriptor::Communicator());
    }

    const bool wantData = (!toIOOnly) || ParallelDescriptor::IOProcessor();
    if (wantData) {
        for (int i = 1; i < nproc; ++i) { offs[i] = offs[i-1] + cnts[i-1]; }
        recv.resize(std::size_t(offs[nproc-1] + cnts[nproc-1]));
    }

    if (toIOOnly) {
        MPI_Gatherv(v.dataPtr(), cnt, dt, recv.dataPtr(), cnts.dataPtr(),
                    offs.dataPtr(), dt, ioproc, ParallelDescriptor::Communicator());
    } else {
        MPI_Allgatherv(v.dataPtr(), cnt, dt, recv.dataPtr(), cnts.dataPtr(),
                       offs.dataPtr(), dt, ParallelDescriptor::Communicator());
    }

    v.swap(recv);
#else
    amrex::ignore_unused(ioproc, cnt, cnts, offs, recv, toIOOnly);
#endif
}

} // anonymous namespace

void allGather  (Vector<Long>& v) { gatherImpl(v, false); }
void allGather  (Vector<Real>& v) { gatherImpl(v, false); }
void gatherToIO (Vector<Long>& v) { gatherImpl(v, true ); }
void gatherToIO (Vector<Real>& v) { gatherImpl(v, true ); }

} // namespace amrcomplex
