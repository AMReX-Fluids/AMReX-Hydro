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
                  int bclo, int domlo,
                  int bchi, int domhi,
                  bool is_velocity,
                  amrex::Array4<const int> const& bc_arr = {})
{

    using namespace amrex;

    // This function is for setting BCs ON domain faces only.
    AMREX_ASSERT( domlo <= i && i <= domhi+1 );

    if (i == domlo)
    {
        int bc = bc_arr ? bc_arr(domlo-1, j, k, n) : bclo;

        if (bc==BCType::ext_dir)
        {
            lo = s(domlo-1, j, k, n);
            // For turbulent inflow, there are times when the inflow face
            // may have a predicted outflowing velocity. Here, we preserve
            // the normal component of the Dirichlet BC, but allow the
            // tangential components to transport values from the interior.
            if ( n==XVEL && is_velocity ) hi=lo;
        }
        else if ( bc == BCType::foextrap || bc == BCType::hoextrap ||
                  bc == BCType::reflect_even )
        {
            lo = hi;
        }
        else if (bc == BCType::reflect_odd)
        {
            hi = 0.;
            lo = hi;
        }
    }
    else if (i == domhi+1)
    {
        int bc = bc_arr ? bc_arr(domhi+1, j, k, n) : bchi;

        if (bc==BCType::ext_dir)
        {
            hi = s(domhi+1, j, k, n);
            if ( n==XVEL && is_velocity ) lo = hi;
        }
        else if ( bc == BCType::foextrap || bc == BCType::hoextrap ||
                  bc == BCType::reflect_even )
        {
            hi = lo;
        }
        else if (bc == BCType::reflect_odd)
        {
            lo = 0.;
            hi = lo;
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
void SetYEdgeBCs ( int i, int j, int k, int n,
                   const amrex::Array4<const amrex::Real> &s,
                   amrex::Real &lo, amrex::Real &hi,
                   int bclo, int domlo,
                   int bchi, int domhi,
                   bool is_velocity,
                   amrex::Array4<const int> const& bc_arr = {})
{
    using namespace amrex;

    // This function is for setting BCs ON domain faces only.
    AMREX_ASSERT( domlo <= j && j <= domhi+1 );

    if (j == domlo)
    {
        int bc = bc_arr ? bc_arr(i, domlo-1, k, n) : bclo;

        if (bc==BCType::ext_dir)
        {
            lo = s(i, domlo-1, k, n);
            if ( n==YVEL && is_velocity ) hi=lo;
        }
        else if ( bc == BCType::foextrap || bc == BCType::hoextrap ||
                  bc == BCType::reflect_even )
        {
            lo = hi;
        }
        else if (bc == BCType::reflect_odd)
        {
            hi = 0.;
            lo = hi;
        }
    }
    else  if (j == domhi+1)
    {
        int bc = bc_arr ? bc_arr(i, domhi+1, k, n) : bchi;

        if (bc==BCType::ext_dir)
        {
            hi = s(i, domhi+1, k, n);
            if ( n==YVEL && is_velocity ) lo=hi;
        }
        else if ( bc == BCType::foextrap || bc == BCType::hoextrap ||
                  bc == BCType::reflect_even )
        {
            hi = lo;
        }
        else if(bc == BCType::reflect_odd)
        {
            lo = 0.;
            hi = lo;
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
void SetZEdgeBCs ( int i, int j, int k, int n,
                   const amrex::Array4<const amrex::Real> &s,
                   amrex::Real &lo, amrex::Real &hi,
                   int bclo, int domlo,
                   int bchi, int domhi,
                   bool is_velocity,
                   amrex::Array4<const int> const& bc_arr = {})
{
    using namespace amrex;

    // This function is for setting BCs ON domain faces only.
    AMREX_ASSERT( domlo <= k && k <= domhi+1 );

    if (k == domlo)
    {
        int bc = bc_arr ? bc_arr(i, j, domlo-1, n) : bclo;

        if (bc==BCType::ext_dir)
        {
            lo = s(i, j, domlo-1, n);
            if ( n==ZVEL && is_velocity ) hi=lo;
        }
        else if ( bc == BCType::foextrap || bc == BCType::hoextrap ||
                  bc == BCType::reflect_even )
        {
            lo = hi;
        }
        else if(bc == BCType::reflect_odd)
        {
            hi = 0.;
            lo = hi;
        }
    }
    else if (k == domhi+1)
    {
        int bc = bc_arr ? bc_arr(i, j, domhi+1, n) : bchi;

        if (bc==BCType::ext_dir)
        {
            hi = s(i,j,domhi+1, n);
            if ( n==ZVEL && is_velocity ) lo=hi;
        }
        else if ( bc == BCType::foextrap || bc == BCType::hoextrap ||
                  bc == BCType::reflect_even )
        {
            hi = lo;
        }
        else if(bc == BCType::reflect_odd)
        {
            lo = 0.;
            hi = lo;
        }
    }
    else
        return;
}
#endif
} // namespace
#endif
/** @}*/
