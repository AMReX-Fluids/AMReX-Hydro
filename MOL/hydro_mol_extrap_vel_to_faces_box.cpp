/** \addtogroup MOL
 *  @{
 */

#include <hydro_constants.H>
#include <hydro_bcs_K.H>
#include <hydro_mol.H>
#include <AMReX_BCRec.H>
#include <AMReX_BC_TYPES.H>
#include <hydro_slopes_K.H>

using namespace amrex;

namespace {
    std::pair<bool,bool> has_extdir_or_ho (BCRec const* bcrec, int ncomp, int dir)
    {
        std::pair<bool,bool> r{false,false};
        for (int n = 0; n < ncomp; ++n) {
            r.first = r.first
                 || (bcrec[n].lo(dir) == BCType::ext_dir)
                 || (bcrec[n].lo(dir) == BCType::hoextrap);
            r.second = r.second
                 || (bcrec[n].hi(dir) == BCType::ext_dir)
                 || (bcrec[n].hi(dir) == BCType::hoextrap);
        }
        return r;
    }
}



void
MOL::ExtrapVelToFacesBox (  AMREX_D_DECL( Box const& ubx,
                                          Box const& vbx,
                                          Box const& wbx),
                            AMREX_D_DECL( Array4<Real> const& u,
                                          Array4<Real> const& v,
                                          Array4<Real> const& w),
                            Array4<Real const> const& vcc,
                            const Geometry&  geom,
                            Vector<BCRec> const& h_bcrec,
                            BCRec  const* d_bcrec)
{
    BL_PROFILE("MOL::ExtrapVelToFacesBox");

    int ncomp = AMREX_SPACEDIM; // This is only used because h_bcrec and d_bcrec hold the
                                // bc's for all three velocity components

    const Box& domain_box = geom.Domain();
    const int domain_ilo = domain_box.smallEnd(0);
    const int domain_ihi = domain_box.bigEnd(0);
    const int domain_jlo = domain_box.smallEnd(1);
    const int domain_jhi = domain_box.bigEnd(1);
#if (AMREX_SPACEDIM == 3)
    const int domain_klo = domain_box.smallEnd(2);
    const int domain_khi = domain_box.bigEnd(2);
#endif

    constexpr int average_not_upwind = 0;

    constexpr int order = 2;

    // At an ext_dir or hoextrap boundary,
    //    the boundary value is on the face, not cell center.
    auto extdir_lohi = has_extdir_or_ho(h_bcrec.data(), ncomp, static_cast<int>(Direction::x));
    bool has_extdir_or_ho_lo = extdir_lohi.first;
    bool has_extdir_or_ho_hi = extdir_lohi.second;

    if ((has_extdir_or_ho_lo && domain_ilo >= ubx.smallEnd(0)-1) ||
        (has_extdir_or_ho_hi && domain_ihi <= ubx.bigEnd(0)))
    {
        amrex::ParallelFor(ubx, [vcc,domain_ilo,domain_ihi,u,d_bcrec]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            bool extdir_or_ho_ilo = (d_bcrec[0].lo(0) == BCType::ext_dir) ||
                                    (d_bcrec[0].lo(0) == BCType::hoextrap);
            bool extdir_or_ho_ihi = (d_bcrec[0].hi(0) == BCType::ext_dir) ||
                                    (d_bcrec[0].hi(0) == BCType::hoextrap);

            const Real vcc_pls = vcc(i,j,k,0);
            const Real vcc_mns = vcc(i-1,j,k,0);

            constexpr int     n = 0;

            Real upls = vcc_pls - 0.5 * amrex_calc_xslope_extdir(
                 i  ,j,k,0,order,vcc,extdir_or_ho_ilo, extdir_or_ho_ihi, domain_ilo, domain_ihi);

            Real umns = vcc_mns + 0.5 * amrex_calc_xslope_extdir(
                 i-1,j,k,0,order,vcc,extdir_or_ho_ilo, extdir_or_ho_ihi, domain_ilo, domain_ihi);

            HydroBC::SetXEdgeBCs(i, j, k, n, vcc, umns, upls, d_bcrec[0].lo(0), domain_ilo, d_bcrec[0].hi(0), domain_ihi, true);

            if ( (i==domain_ilo) && (d_bcrec[0].lo(0) == BCType::foextrap || d_bcrec[0].lo(0) == BCType::hoextrap) )
            {
                upls = amrex::min(upls,0.0_rt);
                umns = upls;
            }
            if ( (i==domain_ihi+1) && (d_bcrec[0].hi(0) == BCType::foextrap || d_bcrec[0].hi(0) == BCType::hoextrap) )
            {
                 umns = amrex::max(umns,0.0_rt);
                 upls = umns;
            }

            Real u_val(0);
            if (average_not_upwind) {
                u_val = 0.5 * (upls + umns);
            } else if (umns >= 0.0 || upls <= 0.0) {

                Real avg = 0.5 * (upls + umns);

                if (avg >= small_vel) {
                    u_val = umns;
                }
                else if (avg <= -small_vel){
                    u_val = upls;
                }
            }

            if (i == domain_ilo && (d_bcrec[0].lo(0) == BCType::ext_dir)) {
                u_val = vcc_mns;
            } else if (i == domain_ihi+1 && (d_bcrec[0].hi(0) == BCType::ext_dir)) {
                u_val = vcc_pls;
            }

            u(i,j,k) = u_val;
        });
    }
    else
    {
        amrex::ParallelFor(ubx, [vcc,domain_ilo,domain_ihi,u,d_bcrec]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            constexpr int     n = 0;

            Real upls = vcc(i  ,j,k,0) - 0.5 * amrex_calc_xslope(i  ,j,k,0,order,vcc);
            Real umns = vcc(i-1,j,k,0) + 0.5 * amrex_calc_xslope(i-1,j,k,0,order,vcc);

            HydroBC::SetXEdgeBCs(i, j, k, n, vcc, umns, upls, d_bcrec[0].lo(0), domain_ilo, d_bcrec[0].hi(0), domain_ihi, true);

            if ( (i==domain_ilo) && (d_bcrec[0].lo(0) == BCType::foextrap || d_bcrec[0].lo(0) == BCType::hoextrap) )
            {
                upls = amrex::min(upls,0.0_rt);
                umns = upls;
            }
            if ( (i==domain_ihi+1) && (d_bcrec[0].hi(0) == BCType::foextrap || d_bcrec[0].hi(0) == BCType::hoextrap) )
            {
                 umns = amrex::max(umns,0.0_rt);
                 upls = umns;
            }

            Real u_val(0);

            if (average_not_upwind) {
                u_val = 0.5 * (upls + umns);
            } else if (umns >= 0.0 || upls <= 0.0) {

                Real avg = 0.5 * (upls + umns);

                if (avg >= small_vel) {
                    u_val = umns;
                }
                else if (avg <= -small_vel){
                    u_val = upls;
                }
            }

            u(i,j,k) = u_val;
        });
    }

    // At an ext_dir or hoextrap boundary,
    //    the boundary value is on the face, not cell center.
    extdir_lohi = has_extdir_or_ho(h_bcrec.data(), ncomp, static_cast<int>(Direction::y));
    has_extdir_or_ho_lo = extdir_lohi.first;
    has_extdir_or_ho_hi = extdir_lohi.second;

    if ((has_extdir_or_ho_lo && domain_jlo >= vbx.smallEnd(1)-1) ||
        (has_extdir_or_ho_hi && domain_jhi <= vbx.bigEnd(1)))
    {
        amrex::ParallelFor(vbx, [vcc,domain_jlo,domain_jhi,v,d_bcrec]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            bool extdir_or_ho_jlo = (d_bcrec[1].lo(1) == BCType::ext_dir) ||
                                    (d_bcrec[1].lo(1) == BCType::hoextrap);
            bool extdir_or_ho_jhi = (d_bcrec[1].hi(1) == BCType::ext_dir) ||
                                    (d_bcrec[1].hi(1) == BCType::hoextrap);

            const Real vcc_pls = vcc(i,j,k,1);
            const Real vcc_mns = vcc(i,j-1,k,1);

            constexpr int     n = 1;

            Real vpls = vcc_pls - 0.5 * amrex_calc_yslope_extdir(
                 i,j,k,1,order,vcc,extdir_or_ho_jlo,extdir_or_ho_jhi,domain_jlo,domain_jhi);
            Real vmns = vcc_mns + 0.5 * amrex_calc_yslope_extdir(
                 i,j-1,k,1,order,vcc,extdir_or_ho_jlo,extdir_or_ho_jhi,domain_jlo,domain_jhi);

            HydroBC::SetYEdgeBCs(i, j, k, n, vcc, vmns, vpls, d_bcrec[1].lo(1), domain_jlo, d_bcrec[1].hi(1), domain_jhi, true);

            if ( (j==domain_jlo) && (d_bcrec[1].lo(1) == BCType::foextrap || d_bcrec[1].lo(1) == BCType::hoextrap) )
            {
                vpls = amrex::min(vpls,0.0_rt);
                vmns = vpls;
            }
            if ( (j==domain_jhi+1) && (d_bcrec[1].hi(1) == BCType::foextrap || d_bcrec[1].hi(1) == BCType::hoextrap) )
            {
                 vmns = amrex::max(vmns,0.0_rt);
                 vpls = vmns;
            }

            Real v_val(0);

            if (average_not_upwind) {
                v_val = 0.5 * (vpls + vmns);
            } else if (vmns >= 0.0 || vpls <= 0.0) {
                Real avg = 0.5 * (vpls + vmns);

                if (avg >= small_vel) {
                    v_val = vmns;
                }
                else if (avg <= -small_vel){
                    v_val = vpls;
                }
            }

            if (j == domain_jlo && (d_bcrec[1].lo(1) == BCType::ext_dir)) {
                v_val = vcc_mns;
            } else if (j == domain_jhi+1 && (d_bcrec[1].hi(1) == BCType::ext_dir)) {
                v_val = vcc_pls;
            }

            v(i,j,k) = v_val;
        });
    }
    else
    {
        amrex::ParallelFor(vbx, [vcc,domain_jlo,domain_jhi,v,d_bcrec]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            constexpr int     n = 1;

            Real vpls = vcc(i,j  ,k,1) - 0.5 * amrex_calc_yslope(i,j  ,k,1,order,vcc);
            Real vmns = vcc(i,j-1,k,1) + 0.5 * amrex_calc_yslope(i,j-1,k,1,order,vcc);

            HydroBC::SetYEdgeBCs(i, j, k, n, vcc, vmns, vpls, d_bcrec[1].lo(1), domain_jlo, d_bcrec[1].hi(1), domain_jhi, true);

            if ( (j==domain_jlo) && (d_bcrec[1].lo(1) == BCType::foextrap || d_bcrec[1].lo(1) == BCType::hoextrap) )
            {
                vpls = amrex::min(vpls,0.0_rt);
                vmns = vpls;
            }
            if ( (j==domain_jhi+1) && (d_bcrec[1].hi(1) == BCType::foextrap || d_bcrec[1].hi(1) == BCType::hoextrap) )
            {
                 vmns = amrex::max(vmns,0.0_rt);
                 vpls = vmns;
            }

            Real v_val(0);

            if (average_not_upwind) {
                v_val = 0.5 * (vpls + vmns);
            } else if (vmns >= 0.0 || vpls <= 0.0) {
                Real avg = 0.5 * (vpls + vmns);

                if (avg >= small_vel) {
                    v_val = vmns;
                }
                else if (avg <= -small_vel) {
                    v_val = vpls;
                }
            }

            v(i,j,k) = v_val;
        });
    }

#if (AMREX_SPACEDIM == 3)
    // At an ext_dir or hoextrap boundary,
    //    the boundary value is on the face, not cell center.
    extdir_lohi = has_extdir_or_ho(h_bcrec.data(), ncomp, static_cast<int>(Direction::z));
    has_extdir_or_ho_lo = extdir_lohi.first;
    has_extdir_or_ho_hi = extdir_lohi.second;

    if ((has_extdir_or_ho_lo && domain_klo >= wbx.smallEnd(2)-1) ||
        (has_extdir_or_ho_hi && domain_khi <= wbx.bigEnd(2)))
    {
        amrex::ParallelFor(wbx, [vcc,domain_klo,domain_khi,w,d_bcrec]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            bool extdir_or_ho_klo = (d_bcrec[2].lo(2) == BCType::ext_dir) ||
                                    (d_bcrec[2].lo(2) == BCType::hoextrap);
            bool extdir_or_ho_khi = (d_bcrec[2].hi(2) == BCType::ext_dir) ||
                                    (d_bcrec[2].hi(2) == BCType::hoextrap);

            const Real vcc_pls = vcc(i,j,k,2);
            const Real vcc_mns = vcc(i,j,k-1,2);

            constexpr int     n = 2;

            Real wpls = vcc_pls - 0.5 * amrex_calc_zslope_extdir(
                 i,j,k  ,2,order,vcc,extdir_or_ho_klo,extdir_or_ho_khi,domain_klo,domain_khi);
            Real wmns = vcc_mns + 0.5 * amrex_calc_zslope_extdir(
                 i,j,k-1,2,order,vcc,extdir_or_ho_klo,extdir_or_ho_khi,domain_klo,domain_khi);

            HydroBC::SetZEdgeBCs(i, j, k, n, vcc, wmns, wpls, d_bcrec[2].lo(2), domain_klo, d_bcrec[2].hi(2), domain_khi, true);

            if ( (k==domain_klo) && (d_bcrec[2].lo(2) == BCType::foextrap || d_bcrec[2].lo(2) == BCType::hoextrap) )
            {
                wpls = amrex::min(wpls,0.0_rt);
                wmns = wpls;
            }
            if ( (k==domain_khi+1) && (d_bcrec[2].hi(2) == BCType::foextrap || d_bcrec[2].hi(2) == BCType::hoextrap) )
            {
                 wmns = amrex::max(wmns,0.0_rt);
                 wpls = wmns;
            }

            Real w_val(0);

            if (average_not_upwind) {
                w_val = 0.5 * (wpls + wmns);
            } else if (wmns >= 0.0 || wpls <= 0.0) {
                Real avg = 0.5 * (wpls + wmns);

                if (avg >= small_vel) {
                    w_val = wmns;
                }
                else if (avg <= -small_vel) {
                    w_val = wpls;
                }
            }

            if (k == domain_klo && (d_bcrec[2].lo(2) == BCType::ext_dir)) {
                w_val = vcc_mns;
            } else if (k == domain_khi+1 && (d_bcrec[2].hi(2) == BCType::ext_dir)) {
                w_val = vcc_pls;
            }

            w(i,j,k) = w_val;
        });
    }
    else
    {
        amrex::ParallelFor(wbx, [vcc,domain_klo,domain_khi,w,d_bcrec]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            constexpr int     n = 2;

            Real wpls = vcc(i,j,k  ,2) - 0.5 * amrex_calc_zslope(i,j,k  ,2,order,vcc);
            Real wmns = vcc(i,j,k-1,2) + 0.5 * amrex_calc_zslope(i,j,k-1,2,order,vcc);

            HydroBC::SetZEdgeBCs(i, j, k, n, vcc, wmns, wpls, d_bcrec[2].lo(2), domain_klo, d_bcrec[2].hi(2), domain_khi, true);

            if ( (k==domain_klo) && (d_bcrec[2].lo(2) == BCType::foextrap || d_bcrec[2].lo(2) == BCType::hoextrap) )
            {
                wpls = amrex::min(wpls,0.0_rt);
                wmns = wpls;
            }
            if ( (k==domain_khi+1) && (d_bcrec[2].hi(2) == BCType::foextrap || d_bcrec[2].hi(2) == BCType::hoextrap) )
            {
                 wmns = amrex::max(wmns,0.0_rt);
                 wpls = wmns;
            }

            Real w_val(0);

            if (average_not_upwind) {
                w_val = 0.5 * (wpls + wmns);
            } else if (wmns >= 0.0 || wpls <= 0.0) {
                Real avg = 0.5 * (wpls + wmns);

                if (avg >= small_vel) {
                    w_val = wmns;
                }
                else if (avg <= -small_vel) {
                    w_val = wpls;
                }
            }

            w(i,j,k) = w_val;
        });
    }
#endif
}
/** @}*/
