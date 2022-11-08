/**
 * \file hydro_ebmol_extrap_vel_to_faces_box.cpp
 * \addtogroup EBMOL
 * @{
 *
 */

#include <hydro_constants.H>
#include <hydro_bcs_K.H>
#include <hydro_ebmol.H>
#if (AMREX_SPACEDIM == 2)
#include <hydro_eb_slopes_2D_K.H>
#elif (AMREX_SPACEDIM == 3)
#include <hydro_eb_slopes_3D_K.H>
#endif

using namespace amrex;

namespace {
    std::pair<bool,bool> has_extdir_or_ho (BCRec const* bcrec, int n, int dir)
    {
        std::pair<bool,bool> r{false,false};
        r.first = r.first
             || (bcrec[n].lo(dir) == BCType::ext_dir)
             || (bcrec[n].lo(dir) == BCType::hoextrap);
        r.second = r.second
             || (bcrec[n].hi(dir) == BCType::ext_dir)
             || (bcrec[n].hi(dir) == BCType::hoextrap);
        return r;
    }
}

void
EBMOL::ExtrapVelToFacesBox ( AMREX_D_DECL( Box const& ubx,
                                           Box const& vbx,
                                           Box const& wbx ),
                             AMREX_D_DECL( Array4<Real> const& u,
                                           Array4<Real> const& v,
                                           Array4<Real> const& w ),
                             Array4<Real const> const& vcc,
                             Array4<EBCellFlag const> const& flag,
                             AMREX_D_DECL( Array4<Real const> const& fcx,
                                           Array4<Real const> const& fcy,
                                           Array4<Real const> const& fcz ),
                             Array4<Real const> const& ccc,
                             Array4<Real const> const& vfrac,
                             const Geometry&  geom,
                             Vector<BCRec> const& h_bcrec,
                             const BCRec* d_bcrec )
{

    const Box& domain_box = geom.Domain();
    AMREX_D_TERM(
        const int domain_ilo = domain_box.smallEnd(0);
        const int domain_ihi = domain_box.bigEnd(0);,
        const int domain_jlo = domain_box.smallEnd(1);
        const int domain_jhi = domain_box.bigEnd(1);,
        const int domain_klo = domain_box.smallEnd(2);
        const int domain_khi = domain_box.bigEnd(2););

    int order = 2;

    // ****************************************************************************
    // Decide whether the stencil at each cell might need to see values that
    //     live on face centroids rather than cell centroids, i.e.
    //     are at a domain boundary with ext_dir or hoextrap boundary conditions
    // ****************************************************************************


    int n_for_xbc = 0;

    auto extdir_lohi_x_for_u = has_extdir_or_ho(h_bcrec.data(), n_for_xbc, static_cast<int>(Direction::x));
    bool has_extdir_or_ho_lo_x_for_u = extdir_lohi_x_for_u.first;
    bool has_extdir_or_ho_hi_x_for_u = extdir_lohi_x_for_u.second;

    auto extdir_lohi_y_for_u = has_extdir_or_ho(h_bcrec.data(), n_for_xbc, static_cast<int>(Direction::y));
    bool has_extdir_or_ho_lo_y_for_u = extdir_lohi_y_for_u.first;
    bool has_extdir_or_ho_hi_y_for_u = extdir_lohi_y_for_u.second;

#if (AMREX_SPACEDIM == 3)
    auto extdir_lohi_z_for_u = has_extdir_or_ho(h_bcrec.data(), n_for_xbc, static_cast<int>(Direction::z));
    bool has_extdir_or_ho_lo_z_for_u = extdir_lohi_z_for_u.first;
    bool has_extdir_or_ho_hi_z_for_u = extdir_lohi_z_for_u.second;
#endif



    // ****************************************************************************
    // Predict to x-faces
    // ****************************************************************************
    if ((has_extdir_or_ho_lo_x_for_u && domain_ilo >= ubx.smallEnd(0)-1) ||
        (has_extdir_or_ho_hi_x_for_u && domain_ihi <= ubx.bigEnd(0)    ) ||
        (has_extdir_or_ho_lo_y_for_u && domain_jlo >= vbx.smallEnd(1)-1) ||
        (has_extdir_or_ho_hi_y_for_u && domain_jhi <= vbx.bigEnd(1)    )
#if (AMREX_SPACEDIM == 2)
        )
#elif (AMREX_SPACEDIM == 3)
        ||
        (has_extdir_or_ho_lo_z_for_u && domain_klo >= wbx.smallEnd(2)-1) ||
        (has_extdir_or_ho_hi_z_for_u && domain_khi <= wbx.bigEnd(2)    ) )
#endif
    {
        amrex::ParallelFor(Box(ubx),
        [u,vcc,flag,ccc,vfrac,d_bcrec,
         AMREX_D_DECL(fcx,fcy,fcz),
         AMREX_D_DECL(domain_ilo,domain_jlo,domain_klo),
         AMREX_D_DECL(domain_ihi,domain_jhi,domain_khi),
         order]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real u_val(0);

            AMREX_D_TERM(bool extdir_or_ho_ilo = (d_bcrec[0].lo(0) == BCType::ext_dir) ||
                                                 (d_bcrec[0].lo(0) == BCType::hoextrap);,
                         bool extdir_or_ho_jlo = (d_bcrec[0].lo(1) == BCType::ext_dir) ||
                                                 (d_bcrec[0].lo(1) == BCType::hoextrap);,
                         bool extdir_or_ho_klo = (d_bcrec[0].lo(2) == BCType::ext_dir) ||
                                                 (d_bcrec[0].lo(2) == BCType::hoextrap););

            AMREX_D_TERM(bool extdir_or_ho_ihi = (d_bcrec[0].hi(0) == BCType::ext_dir) ||
                                                 (d_bcrec[0].hi(0) == BCType::hoextrap);,
                         bool extdir_or_ho_jhi = (d_bcrec[0].hi(1) == BCType::ext_dir) ||
                                                 (d_bcrec[0].hi(1) == BCType::hoextrap);,
                         bool extdir_or_ho_khi = (d_bcrec[0].hi(2) == BCType::ext_dir) ||
                                                 (d_bcrec[0].hi(2) == BCType::hoextrap););

            if (flag(i,j,k).isConnected(-1,0,0))
            {
               Real yf = fcx(i,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
               Real zf = fcx(i,j,k,1);
#endif
               AMREX_D_TERM(Real delta_x = Real(0.5) + ccc(i,j,k,0);,
                            Real delta_y = yf  - ccc(i,j,k,1);,
                            Real delta_z = zf  - ccc(i,j,k,2););

               Real vcc_mns = vcc(i-1,j,k,0);
               Real vcc_pls = vcc(i,j,k,0);

               Real cc_umax = amrex::max(vcc_pls, vcc_mns);
               Real cc_umin = amrex::min(vcc_pls, vcc_mns);

               // Compute slopes of component "0" of vcc
               const auto& slopes_eb_hi = amrex_lim_slopes_extdir_eb(i,j,k,0,vcc,ccc,vfrac,
                                          AMREX_D_DECL(fcx,fcy,fcz), flag,
                                          AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                          AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                          AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                          AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                          order);

#if (AMREX_SPACEDIM == 3)
               Real upls = vcc_pls - delta_x * slopes_eb_hi[0]
                                   + delta_y * slopes_eb_hi[1]
                                   + delta_z * slopes_eb_hi[2];
#else
               Real upls = vcc_pls - delta_x * slopes_eb_hi[0]
                                   + delta_y * slopes_eb_hi[1];
#endif
               upls = amrex::max(amrex::min(upls, cc_umax), cc_umin);

               AMREX_D_TERM(delta_x = Real(0.5) - ccc(i-1,j,k,0);,
                            delta_y = yf  - ccc(i-1,j,k,1);,
                            delta_z = zf  - ccc(i-1,j,k,2););

               // Compute slopes of component "0" of vcc
               const auto& slopes_eb_lo = amrex_lim_slopes_extdir_eb(i-1,j,k,0,vcc,ccc,vfrac,
                                          AMREX_D_DECL(fcx,fcy,fcz), flag,
                                          AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                          AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                          AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                          AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                          order);

#if (AMREX_SPACEDIM == 3)
               Real umns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1]
                                   + delta_z * slopes_eb_lo[2];
#else
               Real umns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1];
#endif
               umns = amrex::max(amrex::min(umns, cc_umax), cc_umin);

               HydroBC::SetXEdgeBCs(i, j, k, 0, vcc, umns, upls, d_bcrec[0].lo(0), domain_ilo, d_bcrec[0].hi(0), domain_ihi, true);

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

               if ( umns >= 0.0 || upls <= 0.0 ) {
                  Real avg = Real(0.5) * ( upls + umns );

                  if (avg >= small_vel) {
                    u_val = umns;
                  }
                  else if (avg <= -small_vel) {
                    u_val = upls;
                  }
               }

               if (i == domain_ilo && (d_bcrec[0].lo(0) == BCType::ext_dir)) {
                   u_val = vcc_mns;
               } else if (i == domain_ihi+1 && (d_bcrec[0].hi(0) == BCType::ext_dir)) {
                   u_val = vcc_pls;
               }
            }

            u(i,j,k) = u_val;
        });
    }
    else
    {
        amrex::ParallelFor(Box(ubx),
        [u,vcc,flag,AMREX_D_DECL(fcx,fcy,fcz),ccc,vfrac,order,d_bcrec,domain_ilo,domain_ihi]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real u_val(0);

            if (flag(i,j,k).isConnected(-1,0,0))
            {
               Real yf = fcx(i,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
               Real zf = fcx(i,j,k,1);
#endif

               AMREX_D_TERM(Real delta_x = Real(0.5) + ccc(i,j,k,0);,
                            Real delta_y = yf  - ccc(i,j,k,1);,
                            Real delta_z = zf  - ccc(i,j,k,2););

               const Real vcc_mns = vcc(i-1,j,k,0);
               const Real vcc_pls = vcc(i,j,k,0);

               Real cc_umax = amrex::max(vcc_pls, vcc_mns);
               Real cc_umin = amrex::min(vcc_pls, vcc_mns);

               // Compute slopes of component "0" of vcc
               const auto slopes_eb_hi = amrex_lim_slopes_eb(i,j,k,0,vcc,ccc,vfrac,AMREX_D_DECL(fcx,fcy,fcz),flag,order);

#if (AMREX_SPACEDIM == 3)
               Real upls = vcc_pls - delta_x * slopes_eb_hi[0]
                                   + delta_y * slopes_eb_hi[1]
                                   + delta_z * slopes_eb_hi[2];
#else
               Real upls = vcc_pls - delta_x * slopes_eb_hi[0]
                                   + delta_y * slopes_eb_hi[1];
#endif
               upls = amrex::max(amrex::min(upls, cc_umax), cc_umin);

               AMREX_D_TERM(delta_x = Real(0.5) - ccc(i-1,j,k,0);,
                            delta_y = yf  - ccc(i-1,j,k,1);,
                            delta_z = zf  - ccc(i-1,j,k,2););

               // Compute slopes of component "0" of vcc
               const auto& slopes_eb_lo = amrex_lim_slopes_eb(i-1,j,k,0,vcc,ccc,vfrac,AMREX_D_DECL(fcx,fcy,fcz),flag,order);

#if (AMREX_SPACEDIM == 3)
               Real umns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1]
                                   + delta_z * slopes_eb_lo[2];
#else
               Real umns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1];
#endif
               umns = amrex::max(amrex::min(umns, cc_umax), cc_umin);

               HydroBC::SetXEdgeBCs(i, j, k, 0, vcc, umns, upls, d_bcrec[0].lo(0), domain_ilo, d_bcrec[0].hi(0), domain_ihi, true);

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

               if ( umns >= 0.0 || upls <= 0.0 ) {
                  Real avg = Real(0.5) * ( upls + umns );

                  if (avg >= small_vel) {
                    u_val = umns;
                  }
                  else if (avg <= -small_vel) {
                    u_val = upls;
                  }
               }
            }

            u(i,j,k) = u_val;
        });
    }

    // ****************************************************************************
    // Predict to y-faces
    // ****************************************************************************

    int n_for_ybc = 1;

    auto extdir_lohi_x_for_v = has_extdir_or_ho(h_bcrec.data(), n_for_ybc, static_cast<int>(Direction::x));
    bool has_extdir_or_ho_lo_x_for_v = extdir_lohi_x_for_v.first;
    bool has_extdir_or_ho_hi_x_for_v = extdir_lohi_x_for_v.second;

    auto extdir_lohi_y_for_v = has_extdir_or_ho(h_bcrec.data(), n_for_ybc, static_cast<int>(Direction::y));
    bool has_extdir_or_ho_lo_y_for_v = extdir_lohi_y_for_v.first;
    bool has_extdir_or_ho_hi_y_for_v = extdir_lohi_y_for_v.second;

#if (AMREX_SPACEDIM == 3)
    auto extdir_lohi_z_for_v = has_extdir_or_ho(h_bcrec.data(), n_for_ybc, static_cast<int>(Direction::z));
    bool has_extdir_or_ho_lo_z_for_v = extdir_lohi_z_for_v.first;
    bool has_extdir_or_ho_hi_z_for_v = extdir_lohi_z_for_v.second;
#endif

    if ((has_extdir_or_ho_lo_x_for_v && domain_ilo >= ubx.smallEnd(0)-1) ||
        (has_extdir_or_ho_hi_x_for_v && domain_ihi <= ubx.bigEnd(0)    ) ||
        (has_extdir_or_ho_lo_y_for_v && domain_jlo >= vbx.smallEnd(1)-1) ||
        (has_extdir_or_ho_hi_y_for_v && domain_jhi <= vbx.bigEnd(1)    )
#if (AMREX_SPACEDIM == 2)
        )
#elif (AMREX_SPACEDIM == 3)
        ||
        (has_extdir_or_ho_lo_z_for_v && domain_klo >= wbx.smallEnd(2)-1) ||
        (has_extdir_or_ho_hi_z_for_v && domain_khi <= wbx.bigEnd(2)    ) )
#endif
    {
        amrex::ParallelFor(Box(vbx),
        [v,vcc,flag,ccc,vfrac,d_bcrec,
         AMREX_D_DECL(fcx,fcy,fcz),
         AMREX_D_DECL(domain_ilo,domain_jlo,domain_klo),
         AMREX_D_DECL(domain_ihi,domain_jhi,domain_khi),
         order]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real v_val(0);

            AMREX_D_TERM(bool extdir_or_ho_ilo = (d_bcrec[1].lo(0) == BCType::ext_dir) ||
                                                 (d_bcrec[1].lo(0) == BCType::hoextrap);,
                         bool extdir_or_ho_jlo = (d_bcrec[1].lo(1) == BCType::ext_dir) ||
                                                 (d_bcrec[1].lo(1) == BCType::hoextrap);,
                         bool extdir_or_ho_klo = (d_bcrec[1].lo(2) == BCType::ext_dir) ||
                                                 (d_bcrec[1].lo(2) == BCType::hoextrap););

            AMREX_D_TERM(bool extdir_or_ho_ihi = (d_bcrec[1].hi(0) == BCType::ext_dir) ||
                                                 (d_bcrec[1].hi(0) == BCType::hoextrap);,
                         bool extdir_or_ho_jhi = (d_bcrec[1].hi(1) == BCType::ext_dir) ||
                                                 (d_bcrec[1].hi(1) == BCType::hoextrap);,
                         bool extdir_or_ho_khi = (d_bcrec[1].hi(2) == BCType::ext_dir) ||
                                                 (d_bcrec[1].hi(2) == BCType::hoextrap););

            if (flag(i,j,k).isConnected(0,-1,0))
            {
               Real xf = fcy(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
               Real zf = fcy(i,j,k,1);
#endif

               AMREX_D_TERM(Real delta_x = xf  - ccc(i,j,k,0);,
                            Real delta_y = Real(0.5) + ccc(i,j,k,1);,
                            Real delta_z = zf  - ccc(i,j,k,2););

               const Real vcc_mns = vcc(i,j-1,k,1);
               const Real vcc_pls = vcc(i,j,k,1);

               Real cc_vmax = amrex::max(vcc_pls, vcc_mns);
               Real cc_vmin = amrex::min(vcc_pls, vcc_mns);

               // Compute slopes of component "1" of vcc
               const auto& slopes_eb_hi = amrex_lim_slopes_extdir_eb(i,j,k,1,vcc,ccc,vfrac,
                                          AMREX_D_DECL(fcx,fcy,fcz), flag,
                                          AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                          AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                          AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                          AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                          order);

#if (AMREX_SPACEDIM == 3)
               Real vpls = vcc_pls + delta_x * slopes_eb_hi[0]
                                   - delta_y * slopes_eb_hi[1]
                                   + delta_z * slopes_eb_hi[2];
#else
               Real vpls = vcc_pls + delta_x * slopes_eb_hi[0]
                                   - delta_y * slopes_eb_hi[1];
#endif

               vpls = amrex::max(amrex::min(vpls, cc_vmax), cc_vmin);

               AMREX_D_TERM(delta_x = xf  - ccc(i,j-1,k,0);,
                            delta_y = Real(0.5) - ccc(i,j-1,k,1);,
                            delta_z = zf  - ccc(i,j-1,k,2););

               // Compute slopes of component "1" of vcc
               const auto& slopes_eb_lo = amrex_lim_slopes_extdir_eb(i,j-1,k,1,vcc,ccc,vfrac,
                                          AMREX_D_DECL(fcx,fcy,fcz), flag,
                                          AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                          AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                          AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                          AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                          order);

#if (AMREX_SPACEDIM == 3)
               Real vmns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1]
                                   + delta_z * slopes_eb_lo[2];
#else
               Real vmns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1];
#endif

               vmns = amrex::max(amrex::min(vmns, cc_vmax), cc_vmin);

               HydroBC::SetYEdgeBCs(i, j, k, 1, vcc, vmns, vpls, d_bcrec[1].lo(1), domain_jlo, d_bcrec[1].hi(1), domain_jhi, true);

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

               if ( vmns >= 0.0 || vpls <= 0.0_rt ) {
                  Real avg = Real(0.5) * ( vpls + vmns );

                  if (avg >= small_vel) {
                    v_val = vmns;
                  }
                  else if (avg <= -small_vel) {
                    v_val = vpls;
                  }
               }

               if (j == domain_jlo && (d_bcrec[1].lo(1) == BCType::ext_dir)) {
                   v_val = vcc_mns;
               }
               else if (j == domain_jhi+1 && (d_bcrec[1].hi(1) == BCType::ext_dir)) {
                   v_val = vcc_pls;
               }
            }

            v(i,j,k) = v_val;
        });
    }
    else
    {
        amrex::ParallelFor(Box(vbx),
        [v,vcc,flag,AMREX_D_DECL(fcx,fcy,fcz),ccc,vfrac,order,d_bcrec,domain_jlo,domain_jhi]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real v_val(0);

            if (flag(i,j,k).isConnected(0,-1,0))
            {
               Real xf = fcy(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
               Real zf = fcy(i,j,k,1);
#endif

               AMREX_D_TERM(Real delta_x = xf  - ccc(i,j,k,0);,
                            Real delta_y = Real(0.5) + ccc(i,j,k,1);,
                            Real delta_z = zf  - ccc(i,j,k,2););

               const Real vcc_mns = vcc(i,j-1,k,1);
               const Real vcc_pls = vcc(i,j,k,1);

               Real cc_vmax = amrex::max(vcc_pls, vcc_mns);
               Real cc_vmin = amrex::min(vcc_pls, vcc_mns);

               // Compute slopes of component "1" of vcc
               const auto slopes_eb_hi = amrex_lim_slopes_eb(i,j,k,1,vcc,ccc,vfrac,AMREX_D_DECL(fcx,fcy,fcz),flag,order);

#if (AMREX_SPACEDIM == 3)
               Real vpls = vcc_pls + delta_x * slopes_eb_hi[0]
                                   - delta_y * slopes_eb_hi[1]
                                   + delta_z * slopes_eb_hi[2];
#else
               Real vpls = vcc_pls + delta_x * slopes_eb_hi[0]
                                   - delta_y * slopes_eb_hi[1];
#endif

               vpls = amrex::max(amrex::min(vpls, cc_vmax), cc_vmin);

               AMREX_D_TERM(delta_x = xf  - ccc(i,j-1,k,0);,
                            delta_y = Real(0.5) - ccc(i,j-1,k,1);,
                            delta_z = zf  - ccc(i,j-1,k,2););

               // Compute slopes of component "1" of vcc
               const auto& slopes_eb_lo = amrex_lim_slopes_eb(i,j-1,k,1,vcc,ccc,vfrac,AMREX_D_DECL(fcx,fcy,fcz),flag,order);

#if (AMREX_SPACEDIM == 3)
               Real vmns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1]
                                   + delta_z * slopes_eb_lo[2];
#else
               Real vmns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1];
#endif

               vmns = amrex::max(amrex::min(vmns, cc_vmax), cc_vmin);

               HydroBC::SetYEdgeBCs(i, j, k, 1, vcc, vmns, vpls, d_bcrec[1].lo(1), domain_jlo, d_bcrec[1].hi(1), domain_jhi, true);

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

               if ( vmns >= 0.0 || vpls <= 0.0_rt ) {
                  Real avg = Real(0.5) * ( vpls + vmns );

                  if (avg >= small_vel) {
                    v_val = vmns;
                  }
                  else if (avg <= -small_vel) {
                    v_val = vpls;
                  }
               }
            }

            v(i,j,k) = v_val;
        });
    }

#if (AMREX_SPACEDIM == 3)
    // ****************************************************************************
    // Predict to z-faces
    // ****************************************************************************

    int n_for_zbc = 2;

    auto extdir_lohi_x_for_w = has_extdir_or_ho(h_bcrec.data(), n_for_zbc, static_cast<int>(Direction::x));
    bool has_extdir_or_ho_lo_x_for_w = extdir_lohi_x_for_w.first;
    bool has_extdir_or_ho_hi_x_for_w = extdir_lohi_x_for_w.second;

    auto extdir_lohi_y_for_w = has_extdir_or_ho(h_bcrec.data(), n_for_zbc, static_cast<int>(Direction::y));
    bool has_extdir_or_ho_lo_y_for_w = extdir_lohi_y_for_w.first;
    bool has_extdir_or_ho_hi_y_for_w = extdir_lohi_y_for_w.second;

#if (AMREX_SPACEDIM == 3)
    auto extdir_lohi_z_for_w = has_extdir_or_ho(h_bcrec.data(), n_for_zbc, static_cast<int>(Direction::z));
    bool has_extdir_or_ho_lo_z_for_w = extdir_lohi_z_for_w.first;
    bool has_extdir_or_ho_hi_z_for_w = extdir_lohi_z_for_w.second;
#endif

    if ((has_extdir_or_ho_lo_x_for_w && domain_ilo >= ubx.smallEnd(0)-1) ||
        (has_extdir_or_ho_hi_x_for_w && domain_ihi <= ubx.bigEnd(0)    ) ||
        (has_extdir_or_ho_lo_y_for_w && domain_jlo >= vbx.smallEnd(1)-1) ||
        (has_extdir_or_ho_hi_y_for_w && domain_jhi <= vbx.bigEnd(1)    )
#if (AMREX_SPACEDIM == 2)
        )
#elif (AMREX_SPACEDIM == 3)
        ||
        (has_extdir_or_ho_lo_z_for_w && domain_klo >= wbx.smallEnd(2)-1) ||
        (has_extdir_or_ho_hi_z_for_w && domain_khi <= wbx.bigEnd(2)    ) )
#endif
    {
        amrex::ParallelFor(Box(wbx),
        [w,vcc,flag,ccc,vfrac,d_bcrec,
         AMREX_D_DECL(fcx,fcy,fcz),
         domain_ilo,domain_ihi,domain_jlo,domain_jhi,domain_klo,domain_khi,
         order]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real w_val(0);

            bool extdir_or_ho_ilo = (d_bcrec[2].lo(0) == BCType::ext_dir) ||
                                    (d_bcrec[2].lo(0) == BCType::hoextrap);
            bool extdir_or_ho_ihi = (d_bcrec[2].hi(0) == BCType::ext_dir) ||
                                    (d_bcrec[2].hi(0) == BCType::hoextrap);

            bool extdir_or_ho_jlo = (d_bcrec[2].lo(1) == BCType::ext_dir) ||
                                    (d_bcrec[2].lo(1) == BCType::hoextrap);
            bool extdir_or_ho_jhi = (d_bcrec[2].hi(1) == BCType::ext_dir) ||
                                    (d_bcrec[2].hi(1) == BCType::hoextrap);

            bool extdir_or_ho_klo = (d_bcrec[2].lo(2) == BCType::ext_dir) ||
                                    (d_bcrec[2].lo(2) == BCType::hoextrap);
            bool extdir_or_ho_khi = (d_bcrec[2].hi(2) == BCType::ext_dir) ||
                                    (d_bcrec[2].hi(2) == BCType::hoextrap);

            if (flag(i,j,k).isConnected(0,0,-1))
            {
               Real xf = fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
               Real yf = fcz(i,j,k,1);

               Real delta_x = xf  - ccc(i,j,k,0);
               Real delta_y = yf  - ccc(i,j,k,1);
               Real delta_z = Real(0.5) + ccc(i,j,k,2);

               const Real vcc_mns = vcc(i,j,k-1,2);
               const Real vcc_pls = vcc(i,j,k,2);

               Real cc_wmax = amrex::max(vcc_pls, vcc_mns);
               Real cc_wmin = amrex::min(vcc_pls, vcc_mns);

               // Compute slopes of component "2" of vcc
               const auto& slopes_eb_hi = amrex_lim_slopes_extdir_eb(i,j,k,2,vcc,ccc,vfrac,
                                          AMREX_D_DECL(fcx,fcy,fcz), flag,
                                          AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                          AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                          AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                          AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                          order);

               Real wpls = vcc_pls + delta_x * slopes_eb_hi[0]
                                   + delta_y * slopes_eb_hi[1]
                                   - delta_z * slopes_eb_hi[2];

               wpls = amrex::max(amrex::min(wpls, cc_wmax), cc_wmin);

               delta_x = xf  - ccc(i,j,k-1,0);
               delta_y = yf  - ccc(i,j,k-1,1);
               delta_z = Real(0.5) - ccc(i,j,k-1,2);

               // Compute slopes of component "2" of vcc
               const auto& slopes_eb_lo = amrex_lim_slopes_extdir_eb(i,j,k-1,2,vcc,ccc,vfrac,
                                          AMREX_D_DECL(fcx,fcy,fcz), flag,
                                          AMREX_D_DECL(extdir_or_ho_ilo, extdir_or_ho_jlo, extdir_or_ho_klo),
                                          AMREX_D_DECL(extdir_or_ho_ihi, extdir_or_ho_jhi, extdir_or_ho_khi),
                                          AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                          AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                          order);

               Real wmns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1]
                                   + delta_z * slopes_eb_lo[2];

               wmns = amrex::max(amrex::min(wmns, cc_wmax), cc_wmin);

               HydroBC::SetZEdgeBCs(i, j, k, 2, vcc, wmns, wpls, d_bcrec[2].lo(2), domain_klo, d_bcrec[2].hi(2), domain_khi, true);

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

               if ( wmns >= 0.0 || wpls <= 0.0_rt ) {
                  Real avg = Real(0.5) * ( wpls + wmns );

                  if (avg >= small_vel) {
                    w_val = wmns;
                  }
                  else if (avg <= -small_vel) {
                    w_val = wpls;
                  }
               }

                if (k == domain_klo && (d_bcrec[2].lo(2) == BCType::ext_dir)) {
                    w_val = vcc_mns;
                }
                else if (k == domain_khi+1 && (d_bcrec[2].hi(2) == BCType::ext_dir)) {
                    w_val = vcc_pls;
                }
            }

            w(i,j,k) = w_val;
        });
    }
    else
    {
        amrex::ParallelFor(Box(wbx),
        [w,vcc,flag,AMREX_D_DECL(fcx,fcy,fcz),ccc,vfrac,order,d_bcrec,domain_klo,domain_khi]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real w_val(0);

            if (flag(i,j,k).isConnected(0,0,-1))
            {
               Real xf = fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
               Real yf = fcz(i,j,k,1);

               Real delta_x = xf  - ccc(i,j,k,0);
               Real delta_y = yf  - ccc(i,j,k,1);
               Real delta_z = Real(0.5) + ccc(i,j,k,2);

               const Real vcc_mns = vcc(i,j,k-1,2);
               const Real vcc_pls = vcc(i,j,k,2);

               Real cc_wmax = amrex::max(vcc_pls, vcc_mns);
               Real cc_wmin = amrex::min(vcc_pls, vcc_mns);

               // Compute slopes of component "2" of vcc
               const auto slopes_eb_hi = amrex_lim_slopes_eb(i,j,k,2,vcc,ccc,vfrac,AMREX_D_DECL(fcx,fcy,fcz),flag,order);

               Real wpls = vcc_pls + delta_x * slopes_eb_hi[0]
                                   + delta_y * slopes_eb_hi[1]
                                   - delta_z * slopes_eb_hi[2];

               wpls = amrex::max(amrex::min(wpls, cc_wmax), cc_wmin);

               delta_x = xf  - ccc(i,j,k-1,0);
               delta_y = yf  - ccc(i,j,k-1,1);
               delta_z = Real(0.5) - ccc(i,j,k-1,2);

               // Compute slopes of component "2" of vcc
               const auto& slopes_eb_lo = amrex_lim_slopes_eb(i,j,k-1,2,vcc,ccc,vfrac,AMREX_D_DECL(fcx,fcy,fcz),flag,order);

               Real wmns = vcc_mns + delta_x * slopes_eb_lo[0]
                                   + delta_y * slopes_eb_lo[1]
                                   + delta_z * slopes_eb_lo[2];

               wmns = amrex::max(amrex::min(wmns, cc_wmax), cc_wmin);

               HydroBC::SetZEdgeBCs(i, j, k, 2, vcc, wmns, wpls, d_bcrec[2].lo(2), domain_klo, d_bcrec[2].hi(2), domain_khi, true);

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

               if ( wmns >= 0.0 || wpls <= 0.0_rt ) {
                  Real avg = Real(0.5) * ( wpls + wmns );

                  if (avg >= small_vel) {
                    w_val = wmns;
                  }
                  else if (avg <= -small_vel) {
                    w_val = wpls;
                  }
               }
            }

            w(i,j,k) = w_val;
        });
    }
#endif
}
/** @} */
