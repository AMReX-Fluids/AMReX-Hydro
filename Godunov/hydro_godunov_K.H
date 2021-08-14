/**
 * \file hydro_godunov_K.H
 *
 * \addtogroup Godunov
 *  @{
 */

#ifndef HYDRO_GODUNOV_K_H
#define HYDRO_GODUNOV_K_H

#include <AMReX_Gpu.H>
// #include <AMReX_FArrayBox.H>
#include <AMReX_BCRec.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Array.H>
#include <iomanip>
#include <hydro_constants.H>

/* This header file contains the inlined __host__ __device__ functions required for
   the scalar advection routines for 3D Godunov. It also contains function declarations
   for controlling host functions.  */

namespace GodunovTransBC {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SetTransTermXBCs (const int i, const int j, const int k, const int n,
                       const amrex::Array4<const amrex::Real> &s,
                       amrex::Real &lo,
                       amrex::Real &hi,
                       const int bclo, const int bchi,
                       const int domlo, const int domhi,
                       const bool is_velocity )
{
    using namespace amrex;


    // Low X
    if (i <= domlo)
    {
        if (bclo==BCType::ext_dir)
        {
            lo = s(domlo-1,j,k,n);
            if( n == XVEL && is_velocity ) hi=lo;
        }
        else if(bclo == BCType::foextrap || bclo == BCType::hoextrap || bclo == BCType::reflect_even)
        {
            lo = hi;
        }
        else if(bclo == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }
    // High X
    else if (i > domhi)
    {
        if (bchi==BCType::ext_dir)
        {
            hi = s(domhi+1,j,k,n) ;
            if( n ==XVEL && is_velocity ) lo=hi;
        }
        else if (bchi == BCType::foextrap || bchi == BCType::hoextrap || bchi == BCType::reflect_even)
        {
            hi = lo;
        }
        else if(bchi == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }
    else
        return;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SetTransTermYBCs ( const int i, const int j, const int k, const int n,
                        const amrex::Array4<const amrex::Real> &s,
                        amrex::Real &lo,
                        amrex::Real &hi,
                        const int bclo, const int bchi,
                        const int domlo, const int domhi,
                        const bool is_velocity )
{
    using namespace amrex;


    // Low Y
    if (j <= domlo)
    {
        if (bclo==BCType::ext_dir)
        {
            lo = s(i,domlo-1,k,n);
            if ( n == YVEL && is_velocity ) hi = lo;
        }
        else if (bclo == BCType::foextrap || bclo == BCType::hoextrap || bclo == BCType::reflect_even)
        {
            lo = hi;
        }
        else if(bclo == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }
    // High Y
    else if (j > domhi)
    {
        if (bchi==BCType::ext_dir)
        {
            hi = s(i,domhi+1,k,n);
            if( n == YVEL && is_velocity ) lo = hi ;
        }
        else if (bchi == BCType::foextrap || bchi == BCType::hoextrap || bchi == BCType::reflect_even)
        {
            hi = lo;
        }
        else if (bchi == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }
    else
        return;
}

#if (AMREX_SPACEDIM==3)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SetTransTermZBCs(const int i, const int j, const int k, const int n,
                      const amrex::Array4<const amrex::Real> &s,
                      amrex::Real &lo,
                      amrex::Real &hi,
                      const int bclo, const int bchi,
                      const int domlo, const int domhi,
                      const bool is_velocity)
{
    using namespace amrex;


    // Low Z
    if (k <= domlo)
    {
        if (bclo==BCType::ext_dir)
        {
            lo =s(i,j,domlo-1,n);
            if ( n == ZVEL && is_velocity ) hi = lo;
        }
        else if (bclo == BCType::foextrap || bclo == BCType::hoextrap || bclo == BCType::reflect_even)
        {
            lo = hi;
        }
        else if(bclo == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }
    // High Z
    else if (k > domhi)
    {
        if (bchi==BCType::ext_dir)
        {
            hi = s(i,j,domhi+1,n);
            if ( n == ZVEL && is_velocity ) lo = hi ;
        }
        else if (bchi == BCType::foextrap || bchi == BCType::hoextrap || bchi == BCType::reflect_even)
        {
            hi = lo;
        }
        else if (bchi == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }
    else
        return;
}
#endif
}
#endif
/** @} */
