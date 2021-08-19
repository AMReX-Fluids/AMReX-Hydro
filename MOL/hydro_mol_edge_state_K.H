/** \addtogroup MOL
 *  @{
 */

#ifndef HYDRO_MOL_EDGE_STATE_K_H_
#define HYDRO_MOL_EDGE_STATE_K_H_

#include <AMReX_Slopes_K.H>
#include <hydro_bcs_K.H>
#include <hydro_constants.H>

namespace MOL {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real hydro_mol_xedge_state_extdir ( int i, int j, int k, int n,
                                           amrex::Array4<amrex::Real const> const& q,
                                           amrex::Array4<amrex::Real const> const& umac,
                                           amrex::BCRec const* const d_bcrec,
                                           int domlo, int domhi, bool is_velocity) noexcept
{
    amrex::Real  qs;

    //slope order
    int order = 2;

    bool edlo = (d_bcrec[n].lo(0) == amrex::BCType::ext_dir);
    bool edhi = (d_bcrec[n].hi(0) == amrex::BCType::ext_dir);

    bool extdir_or_ho_lo = (d_bcrec[n].lo(0) == amrex::BCType::ext_dir) || d_bcrec[n].lo(0) == amrex::BCType::hoextrap;
    bool extdir_or_ho_hi = (d_bcrec[n].hi(0) == amrex::BCType::ext_dir) || d_bcrec[n].lo(0) == amrex::BCType::hoextrap;

    if (edlo && i <= domlo)
    {
        qs = q(domlo-1,j,k,n);
    }
    else if ( edhi && i >= domhi+1)
    {
        qs = q(domhi+1,j,k,n);
    }
    else
    {
        amrex::Real qpls = q(i  ,j,k,n) - 0.5 *
	  amrex_calc_xslope_extdir( i  , j, k, n, order, q, extdir_or_ho_lo, extdir_or_ho_hi, domlo, domhi );
        amrex::Real qmns = q(i-1,j,k,n) + 0.5 *
	  amrex_calc_xslope_extdir( i-1, j, k, n, order, q, extdir_or_ho_lo, extdir_or_ho_hi, domlo, domhi );

        HydroBC::SetXEdgeBCs(i,j,k,n,q,qmns,qpls,d_bcrec[n].lo(0),domlo,d_bcrec[n].hi(0),domhi,is_velocity);

        if ( (i==domlo) && (d_bcrec[n].lo(0) == amrex::BCType::foextrap || d_bcrec[n].lo(0) == amrex::BCType::hoextrap) )
        {
            if ( umac(i,j,k) >= 0. && n==XVEL && is_velocity )  qpls = amrex::min(qpls,0.);
            qmns = qpls;
        }
        if ( (i==domhi+1) && (d_bcrec[n].hi(0) == amrex::BCType::foextrap || d_bcrec[n].hi(0) == amrex::BCType::hoextrap) )
        {
            if ( umac(i,j,k) <= 0. && n==XVEL && is_velocity ) qmns = amrex::max(qmns,0.);
            qpls = qmns;
        }

        if ( umac(i,j,k) > small_vel)
        {
            qs = qmns;
        }
        else if ( umac(i,j,k) < -small_vel)
        {
            qs = qpls;
        }
        else
        {
            qs = 0.5*(qmns+qpls);
        }
    }

    return qs;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real hydro_mol_xedge_state ( int i, int j, int k, int n,
                                    amrex::Array4<amrex::Real const> const& q,
                                    amrex::Array4<amrex::Real const> const& umac,
                                    amrex::BCRec const* const d_bcrec,
                                    int domlo, int domhi, bool is_velocity) noexcept
{
    //slope order
    int order = 2;

    amrex::Real qs;
    amrex::Real qpls = q(i  ,j,k,n) - 0.5 * amrex_calc_xslope( i  , j, k, n, order, q );
    amrex::Real qmns = q(i-1,j,k,n) + 0.5 * amrex_calc_xslope( i-1, j, k, n, order, q );

    HydroBC::SetXEdgeBCs(i,j,k,n,q,qmns,qpls,d_bcrec[n].lo(0),domlo,d_bcrec[n].hi(0),domhi,is_velocity);

    if ( (i==domlo) && (d_bcrec[n].lo(0) == amrex::BCType::foextrap || d_bcrec[n].lo(0) == amrex::BCType::hoextrap) )
    {
        if ( umac(i,j,k) >= 0. && n==XVEL && is_velocity )  qpls = amrex::min(qpls,0.);
        qmns = qpls;
    }
    if ( (i==domhi+1) && (d_bcrec[n].hi(0) == amrex::BCType::foextrap || d_bcrec[n].hi(0) == amrex::BCType::hoextrap) )
    {
        if ( umac(i,j,k) <= 0. && n==XVEL && is_velocity ) qmns = amrex::max(qmns,0.);
         qpls = qmns;
    }

    if ( umac(i,j,k) > small_vel)
    {
        qs = qmns;
    }
    else if ( umac(i,j,k) < -small_vel)
    {
        qs = qpls;
    }
    else
    {
        qs = 0.5*(qmns+qpls);
    }

    return qs;
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real hydro_mol_yedge_state_extdir ( int i, int j, int k, int n,
                                           amrex::Array4<amrex::Real const> const& q,
                                           amrex::Array4<amrex::Real const> const& vmac,
                                           amrex::BCRec const* const d_bcrec,
                                           int domlo, int domhi, bool is_velocity) noexcept
{
    //slope order
    int order = 2;

    bool edlo = (d_bcrec[n].lo(1) == amrex::BCType::ext_dir);
    bool edhi = (d_bcrec[n].hi(1) == amrex::BCType::ext_dir);

    bool extdir_or_ho_lo = (d_bcrec[n].lo(1) == amrex::BCType::ext_dir) || d_bcrec[n].lo(1) == amrex::BCType::hoextrap;
    bool extdir_or_ho_hi = (d_bcrec[n].hi(1) == amrex::BCType::ext_dir) || d_bcrec[n].lo(1) == amrex::BCType::hoextrap;

    amrex::Real qs;

    if (edlo && j <= domlo)
    {
        qs = q(i,domlo-1,k,n);
    }
    else if ( edhi && j >= domhi+1)
    {
        qs = q(i,domhi+1,k,n);
    }
    else
    {
        amrex::Real qpls = q(i,j  ,k,n) - 0.5 *
	  amrex_calc_yslope_extdir( i, j  , k, n, order, q, extdir_or_ho_lo, extdir_or_ho_hi, domlo, domhi );
        amrex::Real qmns = q(i,j-1,k,n) + 0.5 *
	  amrex_calc_yslope_extdir( i, j-1, k, n, order, q, extdir_or_ho_lo, extdir_or_ho_hi, domlo, domhi );

        HydroBC::SetYEdgeBCs(i,j,k,n,q,qmns,qpls,d_bcrec[n].lo(1),domlo,d_bcrec[n].hi(1),domhi,is_velocity);

        if ( (j==domlo) && (d_bcrec[n].lo(1) == amrex::BCType::foextrap || d_bcrec[n].lo(1) == amrex::BCType::hoextrap) )
        {
            if ( vmac(i,j,k) >= 0. && n==YVEL && is_velocity )  qpls = amrex::min(qpls,0.);
            qmns = qpls;
        }
        if ( (j==domhi+1) && (d_bcrec[n].hi(1) == amrex::BCType::foextrap || d_bcrec[n].hi(1) == amrex::BCType::hoextrap) )
        {
            if ( vmac(i,j,k) <= 0. && n==YVEL && is_velocity ) qmns = amrex::max(qmns,0.);
             qpls = qmns;
        }

        if ( vmac(i,j,k) > small_vel)
        {
            qs = qmns;
        }
        else if ( vmac(i,j,k) < -small_vel)
        {
            qs = qpls;
        }
        else
        {
            qs = 0.5*(qmns+qpls);
        }
    }

    return qs;
}



AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real hydro_mol_yedge_state ( int i, int j, int k, int n,
                                    amrex::Array4<amrex::Real const> const& q,
                                    amrex::Array4<amrex::Real const> const& vmac,
                                    amrex::BCRec const* const d_bcrec,
                                    int domlo, int domhi, bool is_velocity) noexcept
{
    //slope order
    int order = 2;

    amrex::Real qs;
    amrex::Real qpls = q(i,j  ,k,n) - 0.5 * amrex_calc_yslope( i, j  , k, n, order, q );
    amrex::Real qmns = q(i,j-1,k,n) + 0.5 * amrex_calc_yslope( i, j-1, k, n, order, q );

    HydroBC::SetYEdgeBCs(i,j,k,n,q,qmns,qpls,d_bcrec[n].lo(1),domlo,d_bcrec[n].hi(1),domhi,is_velocity);

    if ( (j==domlo) && (d_bcrec[n].lo(1) == amrex::BCType::foextrap || d_bcrec[n].lo(1) == amrex::BCType::hoextrap) )
    {
        if ( vmac(i,j,k) >= 0. && n==YVEL && is_velocity )  qpls = amrex::min(qpls,0.);
        qmns = qpls;
    }
    if ( (j==domhi+1) && (d_bcrec[n].hi(1) == amrex::BCType::foextrap || d_bcrec[n].hi(1) == amrex::BCType::hoextrap) )
    {
        if ( vmac(i,j,k) <= 0. && n==YVEL && is_velocity ) qmns = amrex::max(qmns,0.);
         qpls = qmns;
    }

    if ( vmac(i,j,k) > small_vel)
    {
        qs = qmns;
    }
    else if ( vmac(i,j,k) < -small_vel)
    {
        qs = qpls;
    }
    else
    {
        qs = 0.5*(qmns+qpls);
    }

    return qs;

}

#if (AMREX_SPACEDIM==3)

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real hydro_mol_zedge_state_extdir ( int i, int j, int k, int n,
                                           amrex::Array4<amrex::Real const> const& q,
                                           amrex::Array4<amrex::Real const> const& wmac,
                                           amrex::BCRec const* const d_bcrec,
                                           int domlo, int domhi, bool is_velocity) noexcept
{
    //slope order
    int order = 2;

    amrex::Real qs;

    bool edlo = (d_bcrec[n].lo(2) == amrex::BCType::ext_dir);
    bool edhi = (d_bcrec[n].hi(2) == amrex::BCType::ext_dir);

    bool extdir_or_ho_lo = (d_bcrec[n].lo(2) == amrex::BCType::ext_dir) || d_bcrec[n].lo(2) == amrex::BCType::hoextrap;
    bool extdir_or_ho_hi = (d_bcrec[n].hi(2) == amrex::BCType::ext_dir) || d_bcrec[n].lo(2) == amrex::BCType::hoextrap;

    if (edlo && k <= domlo)
    {
        qs = q(i,j,domlo-1,n);
    }
    else if ( edhi && k >= domhi+1)
    {
        qs = q(i,j,domhi+1,n);
    }
    else
    {
        amrex::Real qpls = q(i,j,k  ,n) - 0.5 *
	  amrex_calc_zslope_extdir( i, j, k  , n, order, q, extdir_or_ho_lo, extdir_or_ho_hi, domlo, domhi );
        amrex::Real qmns = q(i,j,k-1,n) + 0.5 *
	  amrex_calc_zslope_extdir( i, j, k-1, n, order, q, extdir_or_ho_lo, extdir_or_ho_hi, domlo, domhi );

        HydroBC::SetZEdgeBCs(i,j,k,n,q,qmns,qpls,d_bcrec[n].lo(2),domlo,d_bcrec[n].hi(2),domhi,is_velocity);

        if ( (k==domlo) && (d_bcrec[n].lo(2) == amrex::BCType::foextrap || d_bcrec[n].lo(2) == amrex::BCType::hoextrap) )
        {
            if ( wmac(i,j,k) >= 0. && n==ZVEL && is_velocity )  qpls = amrex::min(qpls,0.);
            qmns = qpls;
        }
        if ( (k==domhi+1) && (d_bcrec[n].hi(2) == amrex::BCType::foextrap || d_bcrec[n].hi(2) == amrex::BCType::hoextrap) )
        {
            if ( wmac(i,j,k) <= 0. && n==ZVEL && is_velocity ) qmns = amrex::max(qmns,0.);
             qpls = qmns;
        }

        if ( wmac(i,j,k) > small_vel)
        {
            qs = qmns;
        }
        else if ( wmac(i,j,k) < -small_vel)
        {
            qs = qpls;
        }
        else
        {
            qs = 0.5*(qmns+qpls);
        }
    }

    return qs;
}



AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real hydro_mol_zedge_state ( int i, int j, int k, int n,
                                    amrex::Array4<amrex::Real const> const& q,
                                    amrex::Array4<amrex::Real const> const& wmac,
                                    amrex::BCRec const* const d_bcrec,
                                    int domlo, int domhi, bool is_velocity) noexcept
{
    //slope order
    int order = 2;

    amrex::Real qs;
    amrex::Real qpls = q(i,j,k  ,n) - 0.5 * amrex_calc_zslope( i, j, k  , n, order, q );
    amrex::Real qmns = q(i,j,k-1,n) + 0.5 * amrex_calc_zslope( i, j, k-1, n, order, q );

    HydroBC::SetZEdgeBCs(i,j,k,n,q,qmns,qpls,d_bcrec[n].lo(2),domlo,d_bcrec[n].hi(2),domhi,is_velocity);

    if ( (k==domlo) && (d_bcrec[n].lo(2) == amrex::BCType::foextrap || d_bcrec[n].lo(2) == amrex::BCType::hoextrap) )
    {
        if ( wmac(i,j,k) >= 0. && n==ZVEL && is_velocity )  qpls = amrex::min(qpls,0.);
        qmns = qpls;
    }
    if ( (k==domhi+1) && (d_bcrec[n].hi(2) == amrex::BCType::foextrap || d_bcrec[n].hi(2) == amrex::BCType::hoextrap) )
    {
        if ( wmac(i,j,k) <= 0. && n==ZVEL && is_velocity ) qmns = amrex::max(qmns,0.);
         qpls = qmns;
    }

    if ( wmac(i,j,k) > small_vel)
    {
        qs = qmns;
    }
    else if ( wmac(i,j,k) < -small_vel)
    {
        qs = qpls;
    }
    else
    {
        qs = 0.5*(qmns+qpls);
    }

    return qs;

}

#endif


}

#endif
/** @}*/
