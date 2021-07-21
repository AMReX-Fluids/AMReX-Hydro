/** \addtogroup MOL
 *  @{
 */

#include <hydro_mol.H>

using namespace amrex;


/*
 * Doxygen's documentation system wants the description where the component is 
 * first declared. I take this to mean that the workaround `fix`, was to move
 * the docs for this function to the file hydro_mol.H  
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
/** @}*/
