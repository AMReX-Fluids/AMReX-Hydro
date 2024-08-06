/**
 * \file hydro_bcs_K.H
 *
 * This header file contains the inlined __host__ __device__ functions required for
 * the boundary condition routines for Godunov and MOL. It also contains function declarations
 * for controlling host functions.
 *
 */

/** \addtogroup Utilities
 * @{
 */


#ifndef HYDRO_BCS_K_H
#define HYDRO_BCS_K_H

#include <AMReX_Gpu.H>
#include <AMReX_BCRec.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Array.H>
#include <iomanip>
#include <hydro_constants.H>

/**
 * \namespace HydroBC
 *
 */
namespace HydroBC{
//
// Choose between single BC per domain face or position dependent BC array
//
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::BCRec getBC (const int i, const int j, const int k, const int n,
                    const amrex::Box& m_domain, const amrex::BCRec* bcr,
                    amrex::Array4<int const> const& bca)
{
    if ( !bca ) { return bcr[n]; }

    int lo[AMREX_SPACEDIM];
    int hi[AMREX_SPACEDIM];
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
        int index[] = {i,j,k};

        for ( int dim = 0; dim < AMREX_SPACEDIM; dim++){
            if (index[dim] < m_domain.smallEnd(dim)) { index[dim] = m_domain.smallEnd(dim); }
            if (index[dim] > m_domain.bigEnd(dim))   { index[dim] = m_domain.bigEnd(dim); }
            if (dim == dir) { index[dim] = m_domain.smallEnd(dim)-1; }
        }
        lo[dir] = bca.contains(index[0], index[1], index[2]) ?
            bca(index[0], index[1], index[2], n) : 0;
// FIXME?? if we don't contain (i,j,k) then it doesn't matter what the bc is there, because we don't touch it

        index[dir] = m_domain.bigEnd(dir)+1;
        hi[dir] = bca.contains(index[0], index[1], index[2]) ?
            bca(index[0], index[1], index[2], n) : 0;
    }

    amrex::BCRec bc(lo, hi);
    return bc;
}

/**
 *
 * <A ID="SetXEdgeBCs"></A>
 *
 * \brief Boundary condition effects
 *
 *
 * For a detailed discussion see the AMReX-Hydro Guide
 * <A HREF="https://amrex-codes.github.io/amrex/hydro_html/EBMOL.html#sec:slopes">
 * section for EBMOL</A>.
 *
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SetXEdgeBCs( int i, int j, int k, int n,
                  const amrex::Array4<const amrex::Real> &s,
                  amrex::Real &lo, amrex::Real &hi,
                  amrex::Real velm, amrex::Real velp,
                  int bclo, int domlo,
                  int bchi, int domhi,
                  bool is_velocity )
{
    using namespace amrex;

    // This function is for setting BCs ON domain faces only.
    // NOTE: this is called in two ways:
    //  1) in the extrapolation of velocity components only --
    //     in this case velm and velp may be predicted
    //  2) for all components in the construction of edge states
    //     in this case velm and velp are identical and equal to
    //     the MAC velocity on the current edge
    AMREX_ASSERT( domlo <= i && i <= domhi+1 );

    if (i == domlo)
    {
        if ( (bclo == BCType::direction_dependent) && is_velocity &&
                (n == XVEL) && (s(domlo-1,j,k,n) >= 0.0) )
        {
            lo = s(domlo-1, j, k, n);
            hi = s(domlo-1, j, k, n);
        }
        else if ( (bclo == BCType::ext_dir) ||
                  (bclo == BCType::direction_dependent && velp >= 0.0) )
        {
            lo = s(domlo-1, j, k, n);
            // For turbulent inflow, there are times when the inflow face
            // may have a predicted outflowing velocity. Here, we preserve
            // the normal component of the Dirichlet BC, but allow the
            // tangential components to transport values from the interior.
            if ( n==XVEL && is_velocity ) hi=lo;
        }
        else if ( bclo == BCType::foextrap || bclo == BCType::hoextrap ||
                  bclo == BCType::reflect_even ||
                 (bclo == BCType::direction_dependent && velp < 0.0))
        {
            lo = hi;
        }
        else if (bclo == BCType::reflect_odd)
        {
            hi = Real(0.);
            lo = Real(0.);
        }
    }
    else if (i == domhi+1)
    {
        if ( (bchi == BCType::direction_dependent) && is_velocity &&
                (n == XVEL) && (s(domhi+1,j,k,n) <= 0.0) )
        {
            hi = s(domhi+1, j, k, n);
            lo = s(domhi+1, j, k, n);
        }
        else if ( (bchi == BCType::ext_dir) ||
             (bchi == BCType::direction_dependent && velm <= 0.0) )
        {
            hi = s(domhi+1, j, k, n);
            if ( n==XVEL && is_velocity ) lo = hi;
        }
        else if ( bchi == BCType::foextrap || bchi == BCType::hoextrap ||
                  bchi == BCType::reflect_even ||
                 (bchi == BCType::direction_dependent && velm > 0.0))
        {
            hi = lo;
        }
        else if (bchi == BCType::reflect_odd)
        {
            lo = Real(0.);
            hi = Real(0.);
        }
    }
    else
        return;
}

/*
 * \fn SetYEdgeBCs
 */
/**
 * \brief Boundary condition effects
 *
 * For a detailed discussion see the AMReX-Hydro Guide section EBMOL.
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SetYEdgeBCs (int i, int j, int k, int n,
                  const amrex::Array4<const amrex::Real> &s,
                  amrex::Real &lo, amrex::Real &hi,
                  amrex::Real velm, amrex::Real velp,
                  int bclo, int domlo,
                  int bchi, int domhi,
                  bool is_velocity )
{
    using namespace amrex;

    // This function is for setting BCs ON domain faces only.
    AMREX_ASSERT( domlo <= j && j <= domhi+1 );

    if (j == domlo)
    {
        if ( (bclo == BCType::direction_dependent) && is_velocity &&
                (n == YVEL) && (s(i,domlo-1,k,n) >= 0.0) )
        {
            lo = s(i, domlo-1, k, n);
            hi = s(i, domlo-1, k, n);
        }
        else if ( (bclo == BCType::ext_dir) ||
                  (bclo == BCType::direction_dependent && velp >= 0.0) )
        {
            lo = s(i, domlo-1, k, n);
            if ( n==YVEL && is_velocity ) hi=lo;
        }
        else if ( bclo == BCType::foextrap || bclo == BCType::hoextrap ||
                  bclo == BCType::reflect_even ||
                 (bclo == BCType::direction_dependent && velp < 0.0))
        {
            lo = hi;
        }
        else if (bclo == BCType::reflect_odd)
        {
            hi = Real(0.);
            lo = Real(0.);
        }
    }
    else  if (j == domhi+1)
    {
        if ( (bchi == BCType::direction_dependent) && is_velocity &&
                (n == YVEL) && (s(i,domhi+1,k,n) <= 0.0) )
        {
            hi = s(i, domhi+1, k, n);
            lo = s(i, domhi+1, k, n);
        }
        else if ( (bchi == BCType::ext_dir) ||
                  (bchi == BCType::direction_dependent && velm <= 0.0) )
        {
            hi = s(i, domhi+1, k, n);
            if ( n==YVEL && is_velocity ) lo=hi;
        }
        else if ( bchi == BCType::foextrap || bchi == BCType::hoextrap ||
                  bchi == BCType::reflect_even ||
                 (bchi == BCType::direction_dependent && velm > 0.0))
        {
            hi = lo;
        }
        else if (bchi == BCType::reflect_odd)
        {
            lo = Real(0.);
            hi = Real(0.);
        }
    }
    else
        return;
}


//                                                                             //
//   =============================== 3D ONLY ===============================   //
//                                                                             //
#if (AMREX_SPACEDIM==3)

/**
 *
 * \brief Boundary condition effects
 *
 * For a detailed discussion see the AMReX-Hydro Guide.
 * For 3D only.
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SetZEdgeBCs (int i, int j, int k, int n,
                  const amrex::Array4<const amrex::Real> &s,
                  amrex::Real &lo, amrex::Real &hi,
                  amrex::Real velm, amrex::Real velp,
                  int bclo, int domlo,
                  int bchi, int domhi,
                  bool is_velocity )
{
    using namespace amrex;

    // This function is for setting BCs ON domain faces only.
    AMREX_ASSERT( domlo <= k && k <= domhi+1 );

    if (k == domlo)
    {
        if ( (bclo == BCType::direction_dependent) && is_velocity &&
                (n == ZVEL) && (s(i,j,domlo-1,n) >= 0.0) )
        {
            lo = s(i, j, domlo-1, n);
            hi = s(i, j, domlo-1, n);
        }
        else if ( (bclo == BCType::ext_dir) ||
                  (bclo == BCType::direction_dependent && velp >= 0.0) )
        {
            lo = s(i, j, domlo-1, n);
            if ( n==ZVEL && is_velocity ) hi=lo;
        }
        else if ( bclo == BCType::foextrap || bclo == BCType::hoextrap ||
                  bclo == BCType::reflect_even ||
                 (bclo == BCType::direction_dependent && velp < 0.0))
        {
            lo = hi;
        }
        else if(bclo == BCType::reflect_odd)
        {
            hi = Real(0.);
            lo = Real(0.);
        }
    }
    else if (k == domhi+1)
    {
        if ( (bchi == BCType::direction_dependent) && is_velocity &&
                (n == ZVEL) && (s(i,j,domhi+1,n) <= 0.0) )
        {
            hi = s(i, j, domhi+1, n);
            lo = s(i, j, domhi+1, n);
        }
        else if ( (bchi == BCType::ext_dir) ||
                  (bchi == BCType::direction_dependent && velm <= 0.0) )
        {
            hi = s(i, j, domhi+1, n);
            if ( n==ZVEL && is_velocity ) lo=hi;
        }
        else if ( bchi == BCType::foextrap || bchi == BCType::hoextrap ||
                  bchi == BCType::reflect_even ||
                 (bchi == BCType::direction_dependent && velm > 0.0))
        {
            hi = lo;
        }
        else if(bchi == BCType::reflect_odd)
        {
            lo = Real(0.);
            hi = Real(0.);
        }
    }
    else
        return;
}
#endif
} // namespace
#endif
/** @}*/
