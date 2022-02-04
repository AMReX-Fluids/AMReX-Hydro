/**
 * \file hydro_bds_edge_state_2D.cpp
 *
 * \addtogroup BDS
 *  @{
 */

#include <hydro_bds.H>

using namespace amrex;

constexpr amrex::Real eps = 1.0e-8;

/**
 * Uses the Bell-Dawson-Shubin (BDS) algorithm, a higher order Godunov
 * method for scalar conservation laws in two dimensions, to compute
 * edge states.
 *
 * \param [in]     s_mf       MultiFab of state
 * \param [in]     state_comp The component of the state MultiFab.
 * \param [in]     geom       Box geometry.
 * \param [in,out] xedge      MultiFab containing x-edges.
 * \param [in,out] yedge      MultiFab containing y-edges.
 * \param [in]     edge_comp  The component of the edge MultiFabs.
 * \param [in]     umac       Face velocities.
 * \param [in]     vmac       Face velocities.
 * \param [in]     fq         Multifab for forces.
 * \param [in]     fq_comp    Component for Multifab for forces.
 * \param [in]     iconserv   Indicates conservative dimensions.
 * \param [in]     dt         Time step.
 *
 */

void
BDS::ComputeEdgeState ( const MultiFab& s_mf,
                               const int state_comp,
                               const Geometry& geom,
                               MultiFab& xedge,
                               MultiFab& yedge,
                               const int edge_comp,
                               MultiFab const& umac,
                               MultiFab const& vmac,
                               MultiFab const& fq,
                               const int fq_comp,
                               const int is_conservative,
                               const Real dt)
{
    BoxArray ba = s_mf.boxArray();
    DistributionMapping dmap = s_mf.DistributionMap();

    MultiFab slope_mf(ba,dmap,3,1);

    BDS::ComputeSlopes(s_mf,geom,slope_mf,state_comp);

    BDS::ComputeConc(s_mf, state_comp,
                         geom,
                         xedge, yedge, edge_comp,
                         slope_mf,
                         umac,
                         vmac,
                         fq, fq_comp,
                         is_conservative,
                         dt);
}

/**
 * Compute bilinear slopes for BDS algorithm.
 *
 * \param [in]  s_mf MultiFab of state.
 * \param [in]  geom Box geometry.
 * \param [out] slope_mf MuliFab to store slope information.
 * \param [in]  comp The component of the MultiFab.
 *
 */

void
BDS::ComputeSlopes (MultiFab const& s_mf,
                        const Geometry& geom,
                        MultiFab& slope_mf,
                        const int state_comp)
{
    constexpr bool limit_slopes = true;


    BoxArray ba = s_mf.boxArray();
    DistributionMapping dmap = s_mf.DistributionMap();
    GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    // local variables
    MultiFab sint_mf(convert(ba,IntVect(AMREX_D_DECL(1,1,1))), dmap, 1, 1);

    Real hx = dx[0];
    Real hy = dx[1];

    for ( MFIter mfi(sint_mf); mfi.isValid(); ++mfi){

        const Box& bx = mfi.growntilebox(1);
        Array4<const Real> const& s    = s_mf.array(mfi, state_comp);
        Array4<      Real> const& sint = sint_mf.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

            // bicubic interpolation to corner points
            // (i,j,k) refers to lower corner of cell
            // Added k index -- placeholder for 2d
            sint(i,j,k) = (s(i-2,j-2,k) + s(i-2,j+1,k) + s(i+1,j-2,k) + s(i+1,j+1,k)
                    - 7.0*(s(i-2,j-1,k) + s(i-2,j  ,k) + s(i-1,j-2,k) + s(i  ,j-2,k) +
                           s(i-1,j+1,k) + s(i  ,j+1,k) + s(i+1,j-1,k) + s(i+1,j  ,k))
                   + 49.0*(s(i-1,j-1,k) + s(i  ,j-1,k) + s(i-1,j  ,k) + s(i  ,j  ,k)) ) / 144.0;
        });

    }

    for ( MFIter mfi(s_mf); mfi.isValid(); ++mfi){

        const Box& bx = mfi.growntilebox(1);

        Array4<const Real> const& s     = s_mf.array(mfi, state_comp);
        Array4<const Real> const& sint  = sint_mf.array(mfi);
        Array4<      Real> const& slope = slope_mf.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

            // compute initial estimates of slopes from unlimited corner points

            // local variables
            Real sumloc, redfac, redmax, div, kdp, sumdif, sgndif;

            Array1D<Real, 1, 4> diff;
            Array1D<Real, 1, 4> smin;
            Array1D<Real, 1, 4> smax;
            Array1D<Real, 1, 4> sc;

            // sx
            slope(i,j,k,0) = 0.5*(sint(i+1,j+1,k) + sint(i+1,j,k) - sint(i,j+1,k) - sint(i,j,k)) / hx;
            // sy
            slope(i,j,k,1) = 0.5*(sint(i+1,j+1,k) - sint(i+1,j,k) + sint(i,j+1,k) - sint(i,j,k)) / hy;
            // sxy
            slope(i,j,k,2) =     (sint(i+1,j+1,k) - sint(i+1,j,k) - sint(i,j+1,k) + sint(i,j,k)) / (hx*hy);

            if (limit_slopes) {

                // ++ / sint(i+1,j+1)
                sc(4) = s(i,j,k) + 0.5*(hx*slope(i,j,k,0) + hy*slope(i,j,k,1)) + 0.25*hx*hy*slope(i,j,k,2);

                // +- / sint(i+1,j  )
                sc(3) = s(i,j,k) + 0.5*(hx*slope(i,j,k,0) - hy*slope(i,j,k,1)) - 0.25*hx*hy*slope(i,j,k,2);

                // -+ / sint(i  ,j+1)
                sc(2) = s(i,j,k) - 0.5*(hx*slope(i,j,k,0) - hy*slope(i,j,k,1)) - 0.25*hx*hy*slope(i,j,k,2);

                // -- / sint(i  ,j  )
                sc(1) = s(i,j,k) - 0.5*(hx*slope(i,j,k,0) + hy*slope(i,j,k,1)) + 0.25*hx*hy*slope(i,j,k,2);

                // enforce max/min bounds
                smin(4) = min(s(i,j,k), s(i+1,j,k), s(i,j+1,k), s(i+1,j+1,k));
                smax(4) = max(s(i,j,k), s(i+1,j,k), s(i,j+1,k), s(i+1,j+1,k));

                smin(3) = min(s(i,j,k), s(i+1,j,k), s(i,j-1,k), s(i+1,j-1,k));
                smax(3) = max(s(i,j,k), s(i+1,j,k), s(i,j-1,k), s(i+1,j-1,k));

                smin(2) = min(s(i,j,k), s(i-1,j,k), s(i,j+1,k), s(i-1,j+1,k));
                smax(2) = max(s(i,j,k), s(i-1,j,k), s(i,j+1,k), s(i-1,j+1,k));

                smin(1) = min(s(i,j,k), s(i-1,j,k), s(i,j-1,k), s(i-1,j-1,k));
                smax(1) = max(s(i,j,k), s(i-1,j,k), s(i,j-1,k), s(i-1,j-1,k));

                for(int mm=1; mm<=4; ++mm){
                   sc(mm) = max(min(sc(mm), smax(mm)), smin(mm));
                }

                // iterative loop
                for(int ll=1; ll<=3; ++ll){
                   sumloc = 0.25*(sc(4) + sc(3) + sc(2) + sc(1));
                   sumdif = (sumloc - s(i,j,k))*4.0;
                   sgndif = std::copysign(1.0,sumdif);

                   for(int mm=1; mm<=4; ++mm){
                      diff(mm) = (sc(mm) - s(i,j,k))*sgndif;
                   }

                   kdp = 0;

                   for(int mm=1; mm<=4; ++mm){
                      if (diff(mm) > eps) {
                         kdp = kdp+1;
                      }
                   }

                   for(int mm=1; mm<=4; ++mm){
                      if (kdp<1) {
                         div = 1.0;
                      } else {
                         div = kdp;
                      }

                      if (diff(mm)>eps) {
                         redfac = sumdif*sgndif/div;
                         kdp = kdp-1;
                      } else {
                         redfac = 0.0;
                      }

                      if (sgndif > 0.0) {
                         redmax = sc(mm) - smin(mm);
                      } else {
                         redmax = smax(mm) - sc(mm);
                      }

                      redfac = min(redfac,redmax);
                      sumdif = sumdif - redfac*sgndif;
                      sc(mm) = sc(mm) - redfac*sgndif;
                   }
                }

                // final slopes
                // sx
                slope(i,j,k,0) = 0.5*( sc(4) + sc(3) -sc(1) - sc(2) )/hx;
                // sy
                slope(i,j,k,1) = 0.5*( sc(4) + sc(2) -sc(1) - sc(3) )/hy;
                // sxy
                slope(i,j,k,2) =     ( sc(1) + sc(4) -sc(2) - sc(3) )/(hx*hy);

            }
        });
    }
}

/**
 * Returns updated edge value.
 *
 * \param [in] s
 * \param [in] slope
 * \param [in] del
 *
 *
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real eval (const Real s,
           Array1D<Real,1,3>& slope,
           Array1D<Real,1,2>& del)
{
    Real val = s + del(1)*slope(1) + del(2)*slope(2) + del(1)*del(2)*slope(3);

    return val;
}

/**
 * Compute Conc for BDS algorithm.
 *
 * \param [in] s_mf MultiFab of state.
 * \param [in] state_comp Component of the MultiFab of state.
 * \param [in] geom Box geometry.
 * \param [in,out] xedge MuliFab containing x-edges.
 * \param [in,out] yedge MuliFab containing y-edges.
 * \param [in] edge_comp The component of the edge MultiFab.
 * \param [in] slope_mf MuliFab containing slope information.
 * \param [in] umac MuliFab for u-face velocity.
 * \param [in] vmac MuliFab for v-face velocity.
 * \param [in] fq Multifab for forces.
 * \param [in] fq_comp Component for Multifab for forces.
 * \param [in] dt Time step.
 *
 *
 */

void
BDS::ComputeConc (const MultiFab& s_mf,
                      const int state_comp,
                      const Geometry& geom,
                      MultiFab& xedge,
                      MultiFab& yedge,
                      const int edge_comp,
                      const MultiFab& slope_mf,
                      MultiFab const& umac,
                      MultiFab const& vmac,
                      MultiFab const& fq,
                      const int fq_comp,
                      const int is_conservative,
                      const Real dt)
{

    BoxArray ba = s_mf.boxArray();
    DistributionMapping dmap = s_mf.DistributionMap();
    GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    // local variables
    int Nghost = 1;

    MultiFab   ux_mf(ba, dmap, 1, Nghost);
    MultiFab   vy_mf(ba, dmap, 1, Nghost);
    MultiFab divu_mf(ba, dmap, 1, Nghost);
    
    Real hx = dx[0];
    Real hy = dx[1];

    Real dt2 = dt/2.0;
    Real dt3 = dt/3.0;

    constexpr Real half = 0.5;

    // compute cell-centered ux, vy, and divu
    for ( MFIter mfi(ux_mf); mfi.isValid(); ++mfi){

        const Box& bx = mfi.growntilebox(1);

        Array4<const Real> const& uadv  = umac.array(mfi);
        Array4<const Real> const& vadv  = vmac.array(mfi);
        Array4<      Real> const& ux    = ux_mf.array(mfi);
        Array4<      Real> const& vy    = vy_mf.array(mfi);
        Array4<      Real> const& divu  = divu_mf.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

             ux(i,j,k) = (uadv(i+1,j,k) - uadv(i,j,k)) / hx;
             vy(i,j,k) = (vadv(i,j+1,k) - vadv(i,j,k)) / hy;
             divu(i,j,k) = ux(i,j,k) + vy(i,j,k);

       });
    }

    // compute sedgex on x-faces
    for ( MFIter mfi(umac); mfi.isValid(); ++mfi){

        const Box& bx = mfi.tilebox();

        Array4<const Real> const& s      = s_mf.array(mfi, state_comp);
        Array4<const Real> const& slope  = slope_mf.array(mfi);
        Array4<const Real> const& uadv   = umac.array(mfi);
        Array4<const Real> const& vadv   = vmac.array(mfi);
        Array4<const Real> const& force  = fq.array(mfi,fq_comp);
        Array4<const Real> const& divu  = divu_mf.array(mfi);

        Array4<      Real> const& ux     = ux_mf.array(mfi);
        Array4<      Real> const& vy     = vy_mf.array(mfi);
        Array4<      Real> const& sedgex = xedge.array(mfi,edge_comp);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

            //local variables
            Array1D<Real, 1, 2> del;
            Array1D<Real, 1, 2> p1;
            Array1D<Real, 1, 2> p2;
            Array1D<Real, 1, 2> p3;

            Array1D<Real, 1, 3> slope_tmp;

            int ioff, joff;
            Real isign,jsign;
            Real val1,val2,val3;
            Real u, v;
            Real gamma;

            ///////////////////////////////////////////////
            // compute sedgex without transverse corrections
            ///////////////////////////////////////////////

            if (uadv(i,j,k) > 0) {
                isign = 1.;
                ioff = -1;
            } else {
                isign = -1.;
                ioff = 0;
            }

            for(int n=1; n<=3; ++n){
                slope_tmp(n) = slope(i+ioff,j,k,n-1);
            }
            
            // centroid of rectangular volume
            del(1) = isign*0.5*hx - 0.5*uadv(i,j,k)*dt;
            del(2) = 0.;
            sedgex(i,j,k) = eval(s(i+ioff,j,k),slope_tmp,del);

            // source term
            if (is_conservative) {
                sedgex(i,j,k) = sedgex(i,j,k)*(1. - dt2*ux(i+ioff,j,k)) + dt2*force(i+ioff,j,k);
            } else {
                sedgex(i,j,k) = sedgex(i,j,k)*(1. + dt2*vy(i+ioff,j,k)) + dt2*force(i+ioff,j,k);
            }

            ///////////////////////////////////////////////
            // compute \Gamma^{y+}
            ///////////////////////////////////////////////

            if (vadv(i+ioff,j+1,k) > 0.) {
                jsign = 1.;
                joff = 0;
            } else {
                jsign = -1.;
                joff = 1;
            }

            u = 0.;
            if (uadv(i,j,k)*uadv(i,j+joff,k) > 0.) {
                u = uadv(i,j+joff,k);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = jsign*0.5*hy;

            p2(1) = isign*0.5*hx - uadv(i,j,k)*dt;
            p2(2) = jsign*0.5*hy;

            p3(1) = isign*0.5*hx - u*dt;
            p3(2) = jsign*0.5*hy - vadv(i+ioff,j+1,k)*dt;

            for(int n=1; n<=3; ++n){
                slope_tmp(n) = slope(i+ioff,j+joff,k,n-1);
            }

            for (int ll=1; ll<=2; ++ll) {
                del(ll) = (p2(ll)+p3(ll))/2.;
            }
            val1 = eval(s(i+ioff,j+joff,k),slope_tmp,del);

            for (int ll=1; ll<=2; ++ll) {
                del(ll) = (p1(ll)+p3(ll))/2.;
            }
            val2 = eval(s(i+ioff,j+joff,k),slope_tmp,del);

            for (int ll=1; ll<=2; ++ll) {
                del(ll) = (p1(ll)+p2(ll))/2.;
            }
            val3 = eval(s(i+ioff,j+joff,k),slope_tmp,del);

            // average these centroid values to get the average value
            gamma = (val1+val2+val3)/3.;

            // source term
            if (is_conservative) {
                gamma = gamma*(1. - dt3*divu(i+ioff,j+joff,k));
            }

            ///////////////////////////////////////////////
            // correct sedgex with \Gamma^{y+}
            ///////////////////////////////////////////////

            gamma = gamma * vadv(i+ioff,j+1,k);
            sedgex(i,j,k) = sedgex(i,j,k) - dt*gamma/(2.*hy);
          
            ///////////////////////////////////////////////
            // compute \Gamma^{y-}
            ///////////////////////////////////////////////
          
            if (vadv(i+ioff,j,k) > 0.) {
                jsign = 1.;
                joff = -1;
            } else {
                jsign = -1.;
                joff = 0;
                    }

            u = 0.;
            if (uadv(i,j,k)*uadv(i,j+joff,k) > 0.) {
                u = uadv(i,j+joff,k);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = jsign*0.5*hy;

            p2(1) = isign*0.5*hx - uadv(i,j,k)*dt;
            p2(2) = jsign*0.5*hy;

            p3(1) = isign*0.5*hx - u*dt;
            p3(2) = jsign*0.5*hy - vadv(i+ioff,j,k)*dt;

            for(int n=1; n<=3; ++n){
                slope_tmp(n) = slope(i+ioff,j+joff,k,n-1);
            }
             
            for (int ll=1; ll<=2; ++ll) {
                del(ll) = (p2(ll)+p3(ll))/2.;
            }
            val1 = eval(s(i+ioff,j+joff,k),slope_tmp,del);

            for (int ll=1; ll<=2; ++ll) {
                del(ll) = (p1(ll)+p3(ll))/2.;
            }
            val2 = eval(s(i+ioff,j+joff,k),slope_tmp,del);

            for (int ll=1; ll<=2; ++ll) {
                del(ll) = (p1(ll)+p2(ll))/2.;
            }
            val3 = eval(s(i+ioff,j+joff,k),slope_tmp,del);

            // average these centroid values to get the average value
            gamma = (val1+val2+val3)/3.;

            // source term
            if (is_conservative) {
                gamma = gamma*(1. - dt3*divu(i+ioff,j+joff,k));
            }

            ///////////////////////////////////////////////
            // correct sedgex with \Gamma^{y-}
            ///////////////////////////////////////////////

            gamma = gamma * vadv(i+ioff,j,k);
            sedgex(i,j,k) = sedgex(i,j,k) + dt*gamma/(2.*hy);
        });
    }

    // compute sedgey on y-faces
    for ( MFIter mfi(vmac); mfi.isValid(); ++mfi){

        const Box& bx = mfi.tilebox();

        Array4<const Real> const& s      = s_mf.array(mfi, state_comp);
        Array4<const Real> const& slope  = slope_mf.array(mfi);
        Array4<const Real> const& uadv   = umac.array(mfi);
        Array4<const Real> const& vadv   = vmac.array(mfi);
        Array4<const Real> const& force  = fq.array(mfi,fq_comp);
        Array4<const Real> const& divu  = divu_mf.array(mfi);

        Array4<      Real> const& ux     = ux_mf.array(mfi);
        Array4<      Real> const& vy     = vy_mf.array(mfi);
        Array4<      Real> const& sedgey = yedge.array(mfi,edge_comp);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

            //local variables
            Array1D<Real, 1, 2> del;
            Array1D<Real, 1, 2> p1;
            Array1D<Real, 1, 2> p2;
            Array1D<Real, 1, 2> p3;

            Array1D<Real, 1, 3> slope_tmp;

            int ioff, joff;
            Real isign,jsign;
            Real val1,val2,val3;
            Real u, v;
            Real gamma;

            ///////////////////////////////////////////////
            // compute sedgey without transverse corrections
            ///////////////////////////////////////////////

            // centroid of rectangular volume
            if (vadv(i,j,k) > 0.) {
                jsign = 1.;
                joff = -1;
            } else {
                jsign = -1.;
                joff = 0;
            }

            for(int n=1; n<=3; ++n){
                slope_tmp(n) = slope(i,j+joff,k,n-1);
            }

            del(1) = 0.;
            del(2) = jsign*0.5*hy - 0.5*vadv(i,j,k)*dt;
            sedgey(i,j,k) = eval(s(i,j+joff,k),slope_tmp,del);

            // source term
            if (is_conservative) {
                sedgey(i,j,k) = sedgey(i,j,k)*(1. - dt2*vy(i,j+joff,k)) + dt2*force(i,j+joff,k);
            } else {
                sedgey(i,j,k) = sedgey(i,j,k)*(1. + dt2*ux(i,j+joff,k)) + dt2*force(i,j+joff,k);
            }

            ///////////////////////////////////////////////
            // compute \Gamma^{x+}
            ///////////////////////////////////////////////

            if (uadv(i+1,j+joff,k) > 0.) {
                isign = 1.;
                ioff = 0;
            } else {
                isign = -1.;
                ioff = 1;
            }

            v = 0.;
            if (vadv(i,j,k)*vadv(i+ioff,j,k) > 0.) {
                v = vadv(i+ioff,j,k);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = jsign*0.5*hy;

            p2(1) = isign*0.5*hx;
            p2(2) = jsign*0.5*hy - vadv(i,j,k)*dt;

            p3(1) = isign*0.5*hx - uadv(i+1,j+joff,k)*dt;
            p3(2) = jsign*0.5*hy - v*dt;

            for(int n=1; n<=3; ++n){
                slope_tmp(n) = slope(i+ioff,j+joff,k,n-1);
            }
             
            for (int ll=1; ll<=2; ++ll) {
                del(ll) = (p2(ll)+p3(ll))/2.;
            }
            val1 = eval(s(i+ioff,j+joff,k),slope_tmp,del);

            for (int ll=1; ll<=2; ++ll) {
                del(ll) = (p1(ll)+p3(ll))/2.;
            }
            val2 = eval(s(i+ioff,j+joff,k),slope_tmp,del);

            for (int ll=1; ll<=2; ++ll) {
                del(ll) = (p1(ll)+p2(ll))/2.;
            }
            val3 = eval(s(i+ioff,j+joff,k),slope_tmp,del);

            // average these centroid values to get the average value
            gamma = (val1+val2+val3)/3.;

            // source term
            if (is_conservative) {
                gamma = gamma*(1. - dt3*divu(i+ioff,j+joff,k));
            }

            ///////////////////////////////////////////////
            // correct sedgey with \Gamma^{x+}
            ///////////////////////////////////////////////
             
            gamma = gamma * uadv(i+1,j+joff,k);
            sedgey(i,j,k) = sedgey(i,j,k) - dt*gamma/(2.*hx);

            ///////////////////////////////////////////////
            // compute \Gamma^{x-}
            ///////////////////////////////////////////////
          
            if (uadv(i,j+joff,k) > 0.) {
                isign = 1.;
                ioff = -1;
            } else {
                isign = -1.;
                ioff = 0;
            }

            v = 0.;
            if (vadv(i,j,k)*vadv(i+ioff,j,k) > 0.) {
                v = vadv(i+ioff,j,k);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = jsign*0.5*hy;

            p2(1) = isign*0.5*hx;
            p2(2) = jsign*0.5*hy - vadv(i,j,k)*dt;

            p3(1) = isign*0.5*hx - uadv(i,j+joff,k)*dt;
            p3(2) = jsign*0.5*hy - v*dt;

            for(int n=1; n<=3; ++n){
                slope_tmp(n) = slope(i+ioff,j+joff,k,n-1);
            }

            for (int ll=1; ll<=2; ++ll) {
                del(ll) = (p2(ll)+p3(ll))/2.;
            }
            val1 = eval(s(i+ioff,j+joff,k),slope_tmp,del);

            for (int ll=1; ll<=2; ++ll) {
                del(ll) = (p1(ll)+p3(ll))/2.;
            }
            val2 = eval(s(i+ioff,j+joff,k),slope_tmp,del);

            for (int ll=1; ll<=2; ++ll) {
                del(ll) = (p1(ll)+p2(ll))/2.;
            }
            val3 = eval(s(i+ioff,j+joff,k),slope_tmp,del);

            // average these centroid values to get the average value
            gamma = (val1+val2+val3)/3.;

            // source term
            if (is_conservative) {
                gamma = gamma*(1. - dt3*divu(i+ioff,j+joff,k));
            }

            ///////////////////////////////////////////////
            // correct sedgey with \Gamma^{x-}
            ///////////////////////////////////////////////

            gamma = gamma * uadv(i,j+joff,k);
            sedgey(i,j,k) = sedgey(i,j,k) + dt*gamma/(2.*hx);
        });
    }
}

/** @} */
