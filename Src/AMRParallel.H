#ifndef AMRPARALLEL_H_
#define AMRPARALLEL_H_

#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

namespace amrcomplex {

//
// Collective helpers, in one place and with a specific structural purpose.
//
// Code inside "#ifdef AMREX_USE_MPI" is never compiled by a serial build, so
// a plain C++ mistake there survives every serial test and surfaces only when
// someone builds with MPI. That has already happened once here: a
// most-vexing-parse -- "Vector<Long> v(std::size_t(n))" declares a function --
// sat undetected in a gather because nothing serial ever parsed it.
//
// So the conditional region is made as small as it can be. Everything that is
// ordinary C++ -- the counts, the offsets, the resizes, the swap -- lives
// outside the guard and is compiled by every build. Only the MPI calls
// themselves are conditional, and they are written once rather than repeated
// at each call site.
//
// Rule that follows, and worth keeping: declare then resize. Never
// "Vector<T> v(std::size_t(n))".
//
void allGather (amrex::Vector<amrex::Long>& v);
void allGather (amrex::Vector<amrex::Real>& v);

// Personalised exchange: byRank[p] goes to rank p (byRank has NProcs
// entries). Returns everything received, concatenated in rank order.
amrex::Vector<amrex::Long> allToAll (const amrex::Vector<amrex::Vector<amrex::Long>>& byRank);

// Gathers to the IO processor; other ranks are left with an empty vector.
void gatherToIO (amrex::Vector<amrex::Long>& v);
void gatherToIO (amrex::Vector<amrex::Real>& v);

} // namespace amrcomplex

#endif
