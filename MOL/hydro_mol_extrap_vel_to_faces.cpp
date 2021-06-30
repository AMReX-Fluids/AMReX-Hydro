#include <hydro_mol.H>

using namespace amrex;

/**
 * \fn void MOL::ExtrapVelToFaces ( const MultiFab&  a_vel,
 *                                   AMREX_D_DECL( MultiFab& a_umac,
 *                                                 MultiFab& a_vmac,
 *                                                 MultiFab& a_wmac ),
 *                                   const Geometry&  a_geom,
 *           			const Vector<BCRec>& h_bcrec,
 *                                   BCRec  const* d_bcrec)
 * \brief Compute upwinded FC velocities by extrapolating CC values in SPACE ONLY
 *
 * \namespace MOL
 *
 * \param a_vel   The cell-centered velocity field
 * \param a_umac  x-component of velocity on centroid of x-face
 * \param a_vmac  y-component of velocity on centroid of y-face
 * \param a_wmac  z-component of velocity on centroid of z-face
 * \param a_geom  Geometry object at this level 
 * \param h_bcrec Host version of BCRec
 * \param d_bcrec Device version of BCRec
 *
 * Compute upwinded FC velocities by extrapolating CC values in SPACE ONLY
 * This is NOT a Godunov type extrapolation: there is NO dependence on time!
 * The resulting FC velocities are computed at the CENTROID of the face.
 *
 */
void 
MOL::ExtrapVelToFaces ( const MultiFab&  a_vel,
                        AMREX_D_DECL( MultiFab& a_umac,
                                      MultiFab& a_vmac,
                                      MultiFab& a_wmac ),
                        const Geometry&  a_geom,
			const Vector<BCRec>& h_bcrec,
                        BCRec  const* d_bcrec)
{
    BL_PROFILE("MOL::ExtrapVelToFaces");

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(a_vel, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            AMREX_D_TERM( Box const& ubx = mfi.nodaltilebox(0);,
                          Box const& vbx = mfi.nodaltilebox(1);,
                          Box const& wbx = mfi.nodaltilebox(2););

            AMREX_D_TERM( Array4<Real> const& u = a_umac.array(mfi);,
                          Array4<Real> const& v = a_vmac.array(mfi);,
                          Array4<Real> const& w = a_wmac.array(mfi););

            Array4<Real const> const& vcc = a_vel.const_array(mfi);
            ExtrapVelToFacesBox( AMREX_D_DECL(ubx,vbx,wbx),
                                 AMREX_D_DECL(u,v,w),
                                 vcc,a_geom,h_bcrec, d_bcrec);
        }

    }

}
