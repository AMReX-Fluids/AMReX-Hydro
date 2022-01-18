/**
 * \file hydro_bds_edge_state_2D.cpp
 *
 * \addtogroup Bds
 *  @{
 */

#include <hydro_godunov_plm.H>
#include <hydro_godunov_ppm.H>
#include <hydro_godunov.H>
#include <hydro_godunov_K.H>
#include <hydro_bcs_K.H>

/*
 * Check Understanding:
 * s_mf is the multifab state being worked on.
 * edge states are returned and computed from this.
 * macs are inputs.
 *
 *
 */



using namespace amrex;




/**
 * Uses the BDS algorithm to compute edge states in 2D.
 *
 * \param s_mf [in] MultiFab of state
 * \param geom [in] Box geometry.
 * \param edges [out] Array of MultiFabs containing one MultiFab for each, x-edge and y-edge.
 * \param macs [out] Face velocities.
 * \param dt [in] Time step.
 * \param comp [in] The component of the MultiFab.
 *
 */

void
Bds::ComputeEdgeState ( const MultiFab& s_mf,  //input multifab s
                        const Geometry& geom,
                        std::array<MultiFab, AMREX_SPACEDIM>& edges, //MultFab& sn_mf,
                        std::array<MultiFab, AMREX_SPACEDIM>& macs, //umac, vmac, wmac
                        Real dt,
                        int comp,
                        int is_conserv)

{
    BoxArray ba = s_mf.boxArray();
    DistributionMapping dmap = s_mf.DistributionMap();

    MultiFab slope_mf(ba,dmap,3,1);

    Bds::ComputeSlopes(s_mf,geom,slope_mf,comp);

    Bds::ComputeConc(s_mf, geom, edges, slope_mf, macs, dt, comp, is_convserv);


}

/**
 * Compute bilinear slopes for BDS algorithm.
 *
 * \param s_mf [in] MultiFab of state.
 * \param geom [in] Box geometry.
 * \param slope_mf [out] MuliFab containing slope information.
 * \param comp [in] The component of the MultiFab.
 *
 * No changes from bds.cpp
 */


void
Bds::ComputeSlopes( MultiFab const& s_mf,
                    const Geometry& geom,
                    MultiFab& slope_mf,
                    int comp)
{
    BoxArray ba = s_mf.boxArray();
    DistributionMapping dmap = s_mf.DistributionMap();
    GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    // local variables
    MultiFab sint_mf(convert(ba,IntVect(AMREX_D_DECL(1,1,1))), dmap, 1, 1);

    Real hx = dx[0];
    Real hy = dx[1];

    for ( MFIter mfi(sint_mf); mfi.isValid(); ++mfi){

        const Box& bx = mfi.growntilebox(1);
        Array4<const Real> const& s    = s_mf.array(mfi, comp);
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

        Array4<const Real> const& s     = s_mf.array(mfi, comp);
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
 * Compute Concs??? for BDS algorithm.
 *
 * \param s_mf [in] MultiFab of state.
 * \param geom [in] Box geometry.
 * \param edges [out] Array of MuliFabs containing edge states, i.e. x-edge, y-edge.
 * \param slope_mf [in] MuliFab containing slope information.
 * \param macs [in] Array of MuliFabs containing edge states, i.e. x-edge, y-edge.
 * \param dt [in] Time step.
 * \param comp [in] The component of the MultiFab.
 *
 *
 */

void
Bds::ComputeConc (const MultiFab& s_mf,
                  const Geometry& geom,
                  std::array<MultiFab, AMREX_SPACEDIM>& edges,
                  const MultiFab& slope_mf,
                  const std::array<MultiFab, AMREX_SPACEDIM>& macs,
                  const Real dt,
                  int comp,
                  int is_conserv)
{

    BoxArray ba = s_mf.boxArray();
    DistributionMapping dmap = s_mf.DistributionMap();
    GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    // local variables
    int Nghost = 0;

    // These have been replaces by `edges`
    //MultiFab siphj_mf(convert(ba,IntVect::TheDimensionVector(0)), dmap, 1, Nghost);
    //MultiFab sijph_mf(convert(ba,IntVect::TheDimensionVector(1)), dmap, 1, Nghost);

    Real hx = dx[0];
    Real hy = dx[1];

    // calculate Gamma plus for flux F
    for ( MFIter mfi(macs[0]); mfi.isValid(); ++mfi){

        const Box& bx = mfi.tilebox();

        Array4<const Real> const& s      = s_mf.array(mfi, comp);
        Array4<const Real> const& slope  = slope_mf.array(mfi);
        Array4<const Real> const& uadv  = macs[0].array(mfi);
        Array4<const Real> const& vadv  = macs[1].array(mfi);

        //local variables
        //Array4<      Real> const& siphj = siphj_mf.array(mfi);
        Array4<      Real> const& siphj = edges[0].array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

            i--; //adjust indices

            //local variables

            Real hxs,hys;
            Real gamp,gamm;
            Real vtrans,stem,vaddif,vdif;
            Real u1,u2;
            Real vv;

            int iup,jup;
            Real isign, jsign;
            Real divu;


            if (uadv(i+1,j,k) > 0.0) {
               iup   = i;
               isign = 1.0;
            } else {
               iup   = i+1;
               isign = -1.0;
            }

            vtrans = vadv(iup,j+1,k);
            u1 = uadv(i+1,j,k);
            if (vtrans > 0.0) {
               jup   = j;
               jsign = 1.0;
               u2 = uadv(i+1,j,k);
            } else {
               jup   = j+1;
               jsign = -1.0;
               u2 = 0.0;
               if (uadv(i+1,j,k)*uadv(i+1,j+1,k) > 0.0) {
                  u2 = uadv(i+1,j+1,k);
               }
            }

            vv = vadv(iup,j+1,k);

            hxs = hx*isign;
            hys = hy*jsign;

            gamp = s(iup,jup,k)+
                 (hxs*.5 - (u1+u2)*dt/3.0)*slope(iup,jup,k,0) +
                 (hys*.5 -    vv*dt/3.0)  *slope(iup,jup,k,1) +
                 (3.*hxs*hys-2.*(u1+u2)*dt*hys-2.*vv*hxs*dt+
                 vv*(2.*u2+u1)*dt*dt)     *slope(iup,jup,k,2)/12.0;

            // end of calculation of Gamma plus for flux F
            // ****************************************

            // *****************************************
            // calculate Gamma minus for flux F

            if (uadv(i+1,j,k) > 0.0) {
               iup   = i;
               isign = 1.0;
            } else {
               iup   = i+1;
               isign = -1.0;
            }

            vtrans = vadv(iup,j,k);
            u1 = uadv(i+1,j,k);
            if (vtrans > 0.0) {
               jup   = j-1;
               jsign = 1.0;
               u2 = 0.0;
               if (uadv(i+1,j,k)*uadv(i+1,j-1,k) > 0.0) {
                  u2 = uadv(i+1,j-1,k);
               }
            } else {
               jup   = j;
               jsign = -1.0;
               u2 = uadv(i+1,j,k);
            }

            vv = vadv(iup,j,k);

            hxs = hx*isign;
            hys = hy*jsign;

            gamm = s(iup,jup,k)+
                 (hxs*0.5 - (u1+u2)*dt/3.0)*slope(iup,jup,k,0) +
                 (hys*0.5 -      vv*dt/3.0)*slope(iup,jup,k,1) +
                 (3.0*hxs*hys-2.0*(u1+u2)*dt*hys-2.0*vv*hxs*dt +
                      vv*(2.0*u2+u1)*dt*dt)*slope(iup,jup,k,2)/12.0;

            // end of calculation of Gamma minus for flux F
            // ****************************************

            // *********************************
            // calculate siphj

            if (uadv(i+1,j,k) > 0.0) {
               iup   = i;
               isign = 1.0;
            } else {
               iup   = i+1;
               isign = -1.0;
            }

            vdif = 0.5*dt*(vadv(iup,j+1,k)*gamp -
                 vadv(iup,j,k)*gamm ) / hy;
            stem = s(iup,j,k) + (isign*hx - uadv(i+1,j,k)*dt)*0.5*slope(iup,j,k,0);
            vaddif = stem*0.5*dt*(
                    uadv(iup+1,j,k)-uadv(iup,j,k))/hx;
            divu = (uadv(iup+1,j,k)-uadv(iup,j,k))/hx +
                   (vadv(iup,j+1,k)-vadv(iup,j,k))/hy;
            siphj(i+1,j,k) = stem - vdif - vaddif + 0.5*dt*stem*divu;


        });
    } // end of calculation of siphj}


    // calculate Gamma plus for flux G
    for ( MFIter mfi(macs[1]); mfi.isValid(); ++mfi){

        const Box& bx = mfi.tilebox();

        Array4<const Real> const& s      = s_mf.array(mfi, comp);
        Array4<const Real> const& slope  = slope_mf.array(mfi);
        Array4<const Real> const& uadv  = macs[0].array(mfi);
        Array4<const Real> const& vadv  = macs[1].array(mfi);

        //local variables
        Array4<      Real> const& sijph = edges[1].array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

            j--; //adjust indices

            //local variables
            Real hxs,hys;
            Real gamp,gamm;
            Real vtrans,stem,vaddif,vdif;
            Real v1,v2;
            Real uu;

            int iup,jup;
            Real isign, jsign;
            Real divu;

            if (vadv(i,j+1,k) > 0.0) {
               jup   = j;
               jsign = 1.0;
            } else {
               jup   = j+1;
               jsign = -1.0;
            }


            vtrans = uadv(i+1,jup,k);
            v1 = vadv(i,j+1,k);
            if (vtrans > 0.0) {
               iup   = i;
               isign = 1.0;
               v2 = vadv(i,j+1,k);
            } else {
               iup   = i+1;
               isign = -1.0;
               v2 = 0.0;
               if (vadv(i,j+1,k)*vadv(i+1,j+1,k) > 0.0) {
                  v2 = vadv(i+1,j+1,k);
               }
            }

            uu = uadv(i+1,jup,k);

            hxs = hx*isign;
            hys = hy*jsign;

            gamp = s(iup,jup,k)+
                 (hys*0.5 - (v1+v2)*dt/3.0)*slope(iup,jup,k,1) +
                 (hxs*0.5 - uu*dt/3.0)     *slope(iup,jup,k,0) +
                 (3.0*hxs*hys-2.0*(v1+v2)*dt*hxs-2.*uu*hys*dt+
                 (2.0*v2+v1)*uu*dt*dt)     *slope(iup,jup,k,2)/12.0;

            // end of calculation of Gamma plus for flux G
            // ****************************************

            // *****************************************
            // calculate Gamma minus for flux G

            if (vadv(i,j+1,k) > 0.0) {
               jup   = j;
               jsign = 1.0;
            } else {
               jup   = j+1;
               jsign = -1.0;
            }

            vtrans = uadv(i,jup,k);
            v1 = vadv(i,j+1,k);
            if (vtrans > 0.0) {
               iup   = i-1;
               isign = 1.0;
               v2 = 0.0;
               if (vadv(i,j+1,k)*vadv(i-1,j+1,k) > 0) {
                  v2 = vadv(i-1,j+1,k);
               }
            } else {
               iup   = i;
               isign = -1.0;
               v2 = vadv(i,j+1,k);
            }

            uu = uadv(i,jup,k);

            hxs = hx*isign;
            hys = hy*jsign;

            gamm = s(iup,jup,k) +
                 (hys*.5 - (v1+v2)*dt/3.)*slope(iup,jup,k,1) +
                 (hxs*.5 - uu*dt/3.)     *slope(iup,jup,k,0) +
                 (3.*hxs*hys-2.*(v1+v2)*dt*hxs-2.*uu*hys*dt+
                 (2.*v2+v1)*uu*dt*dt)    *slope(iup,jup,k,2)/12.0;

            // end of calculation of Gamma minus for flux G
            // ****************************************

            // *********************************
            // calculate sijph

            if (vadv(i,j+1,k) > 0) {
               jup   = j;
               jsign = 1.0;
            } else {
               jup   = j+1;
               jsign = -1.0;
            }

            vdif = 0.5*dt*
                 (uadv(i+1,jup,k)*gamp-uadv(i,jup,k)*gamm)/hx;
            stem = s(i,jup,k) + (jsign*hy - vadv(i,j+1,k)*dt)*0.5*slope(i,jup,k,1);
            vaddif = stem*0.5*dt*(vadv(i,jup+1,k) - vadv(i,jup,k))/hy;
            divu =  (uadv(i+1,jup,k)-uadv(i,jup,k))/hx +
                 (vadv(i,jup+1,k)-vadv(i,jup,k))/hy;
            sijph(i,j+1,k) = stem - vdif - vaddif + 0.5*dt*stem*divu;

            // end of calculation of sijph
            // *************************************

        });
    }
}

/** @} */

