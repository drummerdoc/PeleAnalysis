#ifndef AMRJOINTSTATS_H_
#define AMRJOINTSTATS_H_

#include "AMRComplex.H"

namespace amrcomplex {

//
// Volume-weighted joint statistics on an AMRComplex: how well does a chosen
// set of progress variables parametrise a response?
//
// The question this answers. Pick coordinates (T, Y_H2, ...) and a response q
// (heat release, a consumption rate). Among ALL functions of the coordinates,
// the one minimising mean-square error in q is the conditional mean
// psi = <q | coords>; there is no better map of that form, so no functional
// form has to be chosen or defended. The error is the conditional variance,
// and the figure of merit is the explained-variance fraction
//
//     R^2 = 1 - E[ var(q | coords) ] / var(q),
//
// the share of the variability in q that the coordinates account for. In
// statistics that is the correlation ratio; in variance-based sensitivity
// analysis the first-order Sobol index; in combustion it is the standard way
// to ask whether a set of variables adequately parametrises the thermochemical
// state.
//
// Why the measure is exact. The traversal is over VALID cells, which partition
// the domain across levels, each weighted by the volume of its own cell. So
// the volume in each bin is an integer count per level times that level's cell
// volume: bit-identical at any rank count, and independent of how the data is
// distributed. Nothing is interpolated and no quadrature enters.
//
// Why two passes. The one-pass identity var = S2/V - (S1/V)^2 cancels away
// about 2*log10(mean/sigma) digits, which for a smooth field on fine bins is
// most of them -- and it is precisely the smooth, fine-bin case that the
// convergence gate rests on. So pass one accumulates counts and S1, giving the
// per-bin mean, and pass two accumulates the squared deviation from THAT. The
// result is stable whatever the ratio of mean to standard deviation.
//
// Why one pass covers every subset and every bin width. Marginalising the
// table over an axis gives the lower-dimensional coordinate subset exactly;
// summing adjacent bins gives the coarser binning exactly. The within-bin sum
// of squares transforms correctly under both by the law of total variance,
//
//     SS_within(coarse) = sum_fine [ SS_within(fine) + V_fine (psi_fine - psi_coarse)^2 ],
//
// evaluated as written rather than as a difference of squares. So {T}, {Y} and
// {T,Y} at every bin width come out of a single traversal of the cells, and the
// bin-width dependence -- R^2 rises monotonically as bins shrink and reaches 1
// trivially when each bin holds one cell -- becomes a column of the table
// instead of a trap.
//
class AMRJointStats
{
public:

    struct Params
    {
        // The complex must carry exactly the coordinates followed by the
        // responses, in that order, as its field components.
        int                        nCoord = 0;
        int                        nResp  = 0;
        amrex::Vector<amrex::Real> clo, chi;   // per coordinate
        amrex::Vector<int>         cbins;      // per coordinate
        int                        verbose = 0;
    };

    AMRJointStats (const AMRComplex& cx, const Params& p);

    struct Fit
    {
        amrex::Real r2       = 0;
        amrex::Real ssWithin = 0;   // volume-weighted, within aggregated bins
        amrex::Real ssTotal  = 0;   // volume-weighted, about the global mean
        amrex::Real volume   = 0;
        amrex::Real mean     = 0;
        long        nOccupied = 0;  // aggregated bins actually holding volume
        long        nAgg      = 0;  // aggregated bins in total
    };

    // mask selects which coordinates are conditioned on, bit d for coordinate
    // d; coarsen divides the bin count of each retained coordinate. mask == 0
    // is the unconditional case and returns the total variance with r2 == 0.
    Fit fit (int resp, unsigned mask, int coarsen) const;

    amrex::Real volumeOfBin (long b) const;
    amrex::Long cellsInBin (long b) const;
    amrex::Real conditionalMean (long b, int resp) const;  // psi; 0 on an empty bin
    amrex::Real conditionalVar  (long b, int resp) const;

    // Bin centre of coordinate d for flattened bin b.
    amrex::Real binCentre (long b, int d) const;
    void        binIndex  (long b, int* out) const;
    long        binOf     (const amrex::Real* coordValues) const;  // -1 if outside

    long nBins () const noexcept { return m_nbinTot; }
    amrex::Real volumeTotal   () const noexcept { return m_volIn;  }
    amrex::Real volumeDropped () const noexcept { return m_volOut; }

private:

    const AMRComplex& m_cx;
    Params            m_p;
    int               m_nlev = 0;
    long              m_nbinTot = 0;

    amrex::Vector<amrex::Real> m_cellVol;      // per level
    amrex::Vector<amrex::Long> m_cnt;          // nbinTot * nlev
    amrex::Vector<amrex::Real> m_s1;           // nbinTot * nResp
    amrex::Vector<amrex::Real> m_ss;           // nbinTot * nResp
    amrex::Real                m_volIn  = 0;
    amrex::Real                m_volOut = 0;
};

} // namespace amrcomplex

#endif
