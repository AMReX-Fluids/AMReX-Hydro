/**
 * \file hydro_bds_edge_state_3D.cpp
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
Godunov::ComputeEdgeStateBDS ( const MultiFab& s_mf,  //input multifab s
                               const Geometry& geom,
                               std::array<MultiFab, AMREX_SPACEDIM>& edges, //MultFab& sn_mf,
                               std::array<MultiFab, AMREX_SPACEDIM>& macs, //umac, vmac, wmac
                               Real dt,
                               int comp)

{
    BoxArray ba = s_mf.boxArray();
    DistributionMapping dmap = s_mf.DistributionMap();

    MultiFab slope_mf(ba,dmap,7,1);

    Godunov::ComputeSlopes(s_mf,geom,slope_mf,comp);

    Godunov::ComputeConc(s_mf, geom, edges, slope_mf, macs, dt, comp, is_convserv);


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
Godunov::ComputeSlopes( MultiFab const& s_mf,
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
    Real hz = dx[2];

    Real c1 = (343.0/1728.0);
    Real c2 = (49.0 /1728.0);
    Real c3 = (7.0  /1728.0);
    Real c4 = (1.0  /1728.0);

    for ( MFIter mfi(sint_mf); mfi.isValid(); ++mfi){

        const Box& bx = mfi.growntilebox(1);
        Array4<const Real> const& s    = s_mf.array(mfi, comp);
        Array4<      Real> const& sint = sint_mf.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

            // tricubic interpolation to corner points
            // (i,j,k) refers to lower corner of cell
            sint(i,j,k) = c1*( s(i  ,j  ,k  ) + s(i-1,j  ,k  ) + s(i  ,j-1,k  )
                              +s(i  ,j  ,k-1) + s(i-1,j-1,k  ) + s(i-1,j  ,k-1)
                              +s(i  ,j-1,k-1) + s(i-1,j-1,k-1) )
                         -c2*( s(i-1,j  ,k+1) + s(i  ,j  ,k+1) + s(i-1,j-1,k+1)
                              +s(i  ,j-1,k+1) + s(i-1,j+1,k  ) + s(i  ,j+1,k  )
                              +s(i-2,j  ,k  ) + s(i+1,j  ,k  ) + s(i-2,j-1,k  )
                              +s(i+1,j-1,k  ) + s(i-1,j-2,k  ) + s(i  ,j-2,k  )
                              +s(i-1,j+1,k-1) + s(i  ,j+1,k-1) + s(i-2,j  ,k-1)
                              +s(i+1,j  ,k-1) + s(i-2,j-1,k-1) + s(i+1,j-1,k-1)
                              +s(i-1,j-2,k-1) + s(i  ,j-2,k-1) + s(i-1,j  ,k-2)
                              +s(i  ,j  ,k-2) + s(i-1,j-1,k-2) + s(i  ,j-1,k-2) )
                         +c3*( s(i-1,j+1,k+1) + s(i  ,j+1,k+1) + s(i-2,j  ,k+1)
                              +s(i+1,j  ,k+1) + s(i-2,j-1,k+1) + s(i+1,j-1,k+1)
                              +s(i-1,j-2,k+1) + s(i  ,j-2,k+1) + s(i-2,j+1,k  )
                              +s(i+1,j+1,k  ) + s(i-2,j-2,k  ) + s(i+1,j-2,k  )
                              +s(i-2,j+1,k-1) + s(i+1,j+1,k-1) + s(i-2,j-2,k-1)
                              +s(i+1,j-2,k-1) + s(i-1,j+1,k-2) + s(i  ,j+1,k-2)
                              +s(i-2,j  ,k-2) + s(i+1,j  ,k-2) + s(i-2,j-1,k-2)
                              +s(i+1,j-1,k-2) + s(i-1,j-2,k-2) + s(i  ,j-2,k-2) )
                         -c4*( s(i-2,j+1,k+1) + s(i+1,j+1,k+1) + s(i-2,j-2,k+1)
                              +s(i+1,j-2,k+1) + s(i-2,j+1,k-2) + s(i+1,j+1,k-2)
                              +s(i-2,j-2,k-2) + s(i+1,j-2,k-2) );
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

            // Variables local to this loop
            Array1D<Real, 1, 8> diff;
            Array1D<Real, 1, 8> smin;
            Array1D<Real, 1, 8> smax;
            Array1D<Real, 1, 8> sc;

             // compute initial estimates of slopes from unlimited corner points
             // sx
             slope(i,j,k,0) = 0.25*(( sint(i+1,j  ,k  ) + sint(i+1,j+1,k  )
                                     +sint(i+1,j  ,k+1) + sint(i+1,j+1,k+1) )
                                   -( sint(i  ,j  ,k  ) + sint(i  ,j+1,k  )
                                     +sint(i  ,j  ,k+1) + sint(i  ,j+1,k+1) )) / hx;
             // sy
             slope(i,j,k,1) = 0.25*(( sint(i  ,j+1,k  ) + sint(i+1,j+1,k  )
                                     +sint(i  ,j+1,k+1) + sint(i+1,j+1,k+1) )
                                   -( sint(i  ,j  ,k  ) + sint(i+1,j  ,k  )
                                     +sint(i  ,j  ,k+1) + sint(i+1,j  ,k+1) )) / hy;

             // sz
             slope(i,j,k,2) = 0.25*(( sint(i  ,j  ,k+1) + sint(i+1,j  ,k+1)
                                     +sint(i  ,j+1,k+1) + sint(i+1,j+1,k+1) )
                                   -( sint(i  ,j  ,k  ) + sint(i+1,j  ,k  )
                                     +sint(i  ,j+1,k  ) + sint(i+1,j+1,k  ) )) / hz;

             // sxy
             slope(i,j,k,3) = 0.5*( ( sint(i  ,j  ,k  ) + sint(i  ,j  ,k+1)
                                     +sint(i+1,j+1,k  ) + sint(i+1,j+1,k+1) )
                                   -( sint(i+1,j  ,k  ) + sint(i+1,j  ,k+1)
                                     +sint(i  ,j+1,k  ) + sint(i  ,j+1,k+1) )) / (hx*hy);

             // sxz
             slope(i,j,k,4) = 0.5*( ( sint(i  ,j  ,k  ) + sint(i  ,j+1,k  )
                                     +sint(i+1,j  ,k+1) + sint(i+1,j+1,k+1) )
                                   -( sint(i+1,j  ,k  ) + sint(i+1,j+1,k  )
                                     +sint(i  ,j  ,k+1) + sint(i  ,j+1,k+1) )) / (hx*hz);

             // syz
             slope(i,j,k,5) = 0.5*( ( sint(i  ,j  ,k  ) + sint(i+1,j  ,k  )
                                     +sint(i  ,j+1,k+1) + sint(i+1,j+1,k+1) )
                                   -( sint(i  ,j  ,k+1) + sint(i+1,j  ,k+1)
                                     +sint(i  ,j+1,k  ) + sint(i+1,j+1,k  ) )) / (hy*hz);

             // sxyz
             slope(i,j,k,6) =       (-sint(i  ,j  ,k  ) + sint(i+1,j  ,k  ) + sint(i  ,j+1,k  )
                                     +sint(i  ,j  ,k+1) - sint(i+1,j+1,k  ) - sint(i+1,j  ,k+1)
                                     -sint(i  ,j+1,k+1) + sint(i+1,j+1,k+1) ) / (hx*hy*hz);

             if (limit_slopes) {

                 // +++ / sint(i+1,j+1,k+1)
                 sc(8) = s(i,j,k)
                      +0.5  *(     hx*slope(i,j,k,0)+   hy*slope(i,j,k,1)+   hz*slope(i,j,k,2))
                      +0.25 *(  hx*hy*slope(i,j,k,3)+hx*hz*slope(i,j,k,4)+hy*hz*slope(i,j,k,5))
                      +0.125*hx*hy*hz*slope(i,j,k,6);

                 // ++- / sint(i+1,j+1,k  )
                 sc(7) = s(i,j,k)
                      +0.5  *(     hx*slope(i,j,k,0)+   hy*slope(i,j,k,1)-   hz*slope(i,j,k,2))
                      +0.25 *(  hx*hy*slope(i,j,k,3)-hx*hz*slope(i,j,k,4)-hy*hz*slope(i,j,k,5))
                      -0.125*hx*hy*hz*slope(i,j,k,6);

                 // +-+ / sint(i+1,j  ,k+1)
                 sc(6) = s(i,j,k)
                      +0.5  *(     hx*slope(i,j,k,0)-   hy*slope(i,j,k,1)+   hz*slope(i,j,k,2))
                      +0.25 *( -hx*hy*slope(i,j,k,3)+hx*hz*slope(i,j,k,4)-hy*hz*slope(i,j,k,5))
                      -0.125*hx*hy*hz*slope(i,j,k,6);

                 // +-- / sint(i+1,j  ,k  )
                 sc(5) = s(i,j,k)
                      +0.5  *(     hx*slope(i,j,k,0)-   hy*slope(i,j,k,1)-   hz*slope(i,j,k,2))
                      +0.25 *( -hx*hy*slope(i,j,k,3)-hx*hz*slope(i,j,k,4)+hy*hz*slope(i,j,k,5))
                      +0.125*hx*hy*hz*slope(i,j,k,6);

                 // -++ / sint(i  ,j+1,k+1)
                 sc(4) = s(i,j,k)
                      +0.5  *(    -hx*slope(i,j,k,0)+   hy*slope(i,j,k,1)+   hz*slope(i,j,k,2))
                      +0.25 *( -hx*hy*slope(i,j,k,3)-hx*hz*slope(i,j,k,4)+hy*hz*slope(i,j,k,5))
                      -0.125*hx*hy*hz*slope(i,j,k,6);

                 // -+- / sint(i  ,j+1,k  )
                 sc(3) = s(i,j,k)
                      +0.5  *(    -hx*slope(i,j,k,0)+   hy*slope(i,j,k,1)-   hz*slope(i,j,k,2))
                      +0.25 *( -hx*hy*slope(i,j,k,3)+hx*hz*slope(i,j,k,4)-hy*hz*slope(i,j,k,5))
                      +0.125*hx*hy*hz*slope(i,j,k,6);

                 // --+ / sint(i  ,j  ,k+1)
                 sc(2) = s(i,j,k)
                      +0.5  *(    -hx*slope(i,j,k,0)-   hy*slope(i,j,k,1)+   hz*slope(i,j,k,2))
                      +0.25 *(  hx*hy*slope(i,j,k,3)-hx*hz*slope(i,j,k,4)-hy*hz*slope(i,j,k,5))
                      +0.125*hx*hy*hz*slope(i,j,k,6);

                 // ---/ sint(i  ,j  ,k  )
                 sc(1) = s(i,j,k)
                      +0.5  *(    -hx*slope(i,j,k,0)-   hy*slope(i,j,k,1)-   hz*slope(i,j,k,2))
                      +0.25 *(  hx*hy*slope(i,j,k,3)+hx*hz*slope(i,j,k,4)+hy*hz*slope(i,j,k,5))
                      -0.125*hx*hy*hz*slope(i,j,k,6);

                 // enforce max/min bounds
                 smin(8) = min(s(i  ,j  ,k  ),s(i+1,j  ,k  ),s(i  ,j+1,k  ),s(i  ,j  ,k+1),
                               s(i+1,j+1,k  ),s(i+1,j  ,k+1),s(i  ,j+1,k+1),s(i+1,j+1,k+1));
                 smax(8) = max(s(i  ,j  ,k  ),s(i+1,j  ,k  ),s(i  ,j+1,k  ),s(i  ,j  ,k+1),
                               s(i+1,j+1,k  ),s(i+1,j  ,k+1),s(i  ,j+1,k+1),s(i+1,j+1,k+1));

                 smin(7) = min(s(i  ,j  ,k-1),s(i+1,j  ,k-1),s(i  ,j+1,k-1),s(i  ,j  ,k  ),
                               s(i+1,j+1,k-1),s(i+1,j  ,k  ),s(i  ,j+1,k  ),s(i+1,j+1,k  ));
                 smax(7) = max(s(i  ,j  ,k-1),s(i+1,j  ,k-1),s(i  ,j+1,k-1),s(i  ,j  ,k  ),
                               s(i+1,j+1,k-1),s(i+1,j  ,k  ),s(i  ,j+1,k  ),s(i+1,j+1,k  ));

                 smin(6) = min(s(i  ,j-1,k  ),s(i+1,j-1,k  ),s(i  ,j  ,k  ),s(i  ,j-1,k+1),
                               s(i+1,j  ,k  ),s(i+1,j-1,k+1),s(i  ,j  ,k+1),s(i+1,j  ,k+1));
                 smax(6) = max(s(i  ,j-1,k  ),s(i+1,j-1,k  ),s(i  ,j  ,k  ),s(i  ,j-1,k+1),
                               s(i+1,j  ,k  ),s(i+1,j-1,k+1),s(i  ,j  ,k+1),s(i+1,j  ,k+1));

                 smin(5) = min(s(i  ,j-1,k-1),s(i+1,j-1,k-1),s(i  ,j  ,k-1),s(i  ,j-1,k  ),
                               s(i+1,j  ,k-1),s(i+1,j-1,k  ),s(i  ,j  ,k  ),s(i+1,j  ,k  ));
                 smax(5) = max(s(i  ,j-1,k-1),s(i+1,j-1,k-1),s(i  ,j  ,k-1),s(i  ,j-1,k  ),
                               s(i+1,j  ,k-1),s(i+1,j-1,k  ),s(i  ,j  ,k  ),s(i+1,j  ,k  ));

                 smin(4) = min(s(i-1,j  ,k  ),s(i  ,j  ,k  ),s(i-1,j+1,k  ),s(i-1,j  ,k+1),
                               s(i  ,j+1,k  ),s(i  ,j  ,k+1),s(i-1,j+1,k+1),s(i  ,j+1,k+1));
                 smax(4) = max(s(i-1,j  ,k  ),s(i  ,j  ,k  ),s(i-1,j+1,k  ),s(i-1,j  ,k+1),
                               s(i  ,j+1,k  ),s(i  ,j  ,k+1),s(i-1,j+1,k+1),s(i  ,j+1,k+1));

                 smin(3) = min(s(i-1,j  ,k-1),s(i  ,j  ,k-1),s(i-1,j+1,k-1),s(i-1,j  ,k  ),
                               s(i  ,j+1,k-1),s(i  ,j  ,k  ),s(i-1,j+1,k  ),s(i  ,j+1,k  ));
                 smax(3) = max(s(i-1,j  ,k-1),s(i  ,j  ,k-1),s(i-1,j+1,k-1),s(i-1,j  ,k  ),
                               s(i  ,j+1,k-1),s(i  ,j  ,k  ),s(i-1,j+1,k  ),s(i  ,j+1,k  ));

                 smin(2) = min(s(i-1,j-1,k  ),s(i  ,j-1,k  ),s(i-1,j  ,k  ),s(i-1,j-1,k+1),
                               s(i  ,j  ,k  ),s(i  ,j-1,k+1),s(i-1,j  ,k+1),s(i  ,j  ,k+1));
                 smax(2) = max(s(i-1,j-1,k  ),s(i  ,j-1,k  ),s(i-1,j  ,k  ),s(i-1,j-1,k+1),
                               s(i  ,j  ,k  ),s(i  ,j-1,k+1),s(i-1,j  ,k+1),s(i  ,j  ,k+1));

                 smin(1) = min(s(i-1,j-1,k-1),s(i  ,j-1,k-1),s(i-1,j  ,k-1),s(i-1,j-1,k  ),
                               s(i  ,j  ,k-1),s(i  ,j-1,k  ),s(i-1,j  ,k  ),s(i  ,j  ,k  ));
                 smax(1) = max(s(i-1,j-1,k-1),s(i  ,j-1,k-1),s(i-1,j  ,k-1),s(i-1,j-1,k  ),
                               s(i  ,j  ,k-1),s(i  ,j-1,k  ),s(i-1,j  ,k  ),s(i  ,j  ,k  ));

                 for(int mm=1; mm<=8; ++mm){
                    sc(mm) = max(min(sc(mm), smax(mm)), smin(mm));
                 }

                 // iterative loop
                 for(int ll = 1; ll<=6; ++ll){
                    sumloc = 0.125*(sc(1)+sc(2)+sc(3)+sc(4)+sc(5)+sc(6)+sc(7)+sc(8));
                    sumdif = (sumloc - s(i,j,k))*8.0;
                    sgndif = std::copysign(1.0,sumdif);

                    for(int mm=1; mm<=8; ++mm){
                       diff(mm) = (sc(mm) - s(i,j,k))*sgndif;
                    }

                    kdp = 0;

                    for(int mm=1; mm<=8; ++mm){
                       if (diff(mm) > eps) {
                          kdp = kdp+1;
                       }
                    }

                    for(int mm=1; mm<=8; ++mm){
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
                 slope(i,j,k,0) = 0.25*( ( sc(5) + sc(7)
                                          +sc(6) + sc(8))
                                        -( sc(1) + sc(3)
                                          +sc(2) + sc(4)) ) / hx;

                 // sy
                 slope(i,j,k,1) = 0.25*( ( sc(3) + sc(7)
                                          +sc(4) + sc(8))
                                        -( sc(1) + sc(5)
                                          +sc(2) + sc(6)) ) / hy;

                 // sz
                 slope(i,j,k,2) = 0.25*( ( sc(2) + sc(6)
                                          +sc(4) + sc(8))
                                        -( sc(1) + sc(5)
                                          +sc(3) + sc(7)) ) / hz;

                 // sxy
                 slope(i,j,k,3) = 0.5*( ( sc(1) + sc(2)
                                         +sc(7) + sc(8))
                                       -( sc(5) + sc(6)
                                         +sc(3) + sc(4)) ) / (hx*hy);

                 // sxz
                 slope(i,j,k,4) = 0.5*( ( sc(1) + sc(3)
                                         +sc(6) + sc(8))
                                       -( sc(5) + sc(7)
                                         +sc(2) + sc(4)) ) / (hx*hz);

                 // syz
                 slope(i,j,k,5) = 0.5*( ( sc(1) + sc(5)
                                         +sc(4) + sc(8))
                                       -( sc(2) + sc(6)
                                         +sc(3) + sc(7)) ) / (hy*hz);

                 // sxyz
                 slope(i,j,k,6) = (-sc(1) + sc(5) + sc(3)
                                   +sc(2) - sc(7) - sc(6)
                                   -sc(4) + sc(8) ) / (hx*hy*hz);

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
Godunov::ComputeConc (const MultiFab& s_mf,
                  const Geometry& geom,
                  std::array<MultiFab, AMREX_SPACEDIM>& edges,
                  const MultiFab& slope_mf,
                  const std::array<MultiFab, AMREX_SPACEDIM>& macs,
                  const Real dt,
                  int comp)
{

    BoxArray ba = s_mf.boxArray();
    DistributionMapping dmap = s_mf.DistributionMap();
    GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    // local variables
    int Nghost = 1;

    // Replaced by edges
    //MultiFab sedgex_mf(ba, dmap, 1, Nghost);
    //MultiFab sedgey_mf(ba, dmap, 1, Nghost);
    //MultiFab sedgez_mf(ba, dmap, 1, Nghost);

    MultiFab ux_mf(ba, dmap, 1, Nghost);
    MultiFab vy_mf(ba, dmap, 1, Nghost);
    MultiFab wz_mf(ba, dmap, 1, Nghost);

    Real hx = dx[0];
    Real hy = dx[1];
    Real hz = dx[2];

    Real dt2 = dt/2.0;
    Real dt3 = dt/3.0;
    Real dt4 = dt/4.0;

    constexpr Real half = 0.5;
    constexpr Real sixth = 1.0/6.0;

    // compute cell-centered ux, vy, and wz
    for ( MFIter mfi(ux_mf); mfi.isValid(); ++mfi){

        const Box& bx = mfi.growntilebox(1);

        Array4<const Real> const& uadv  = umac_mf[0].array(mfi);
        Array4<const Real> const& vadv  = umac_mf[1].array(mfi);
        Array4<const Real> const& wadv  = umac_mf[2].array(mfi);
        Array4<      Real> const& ux    = ux_mf.array(mfi);
        Array4<      Real> const& vy    = vy_mf.array(mfi);
        Array4<      Real> const& wz    = wz_mf.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){


             ux(i,j,k) = (uadv(i+1,j,k) - uadv(i,j,k)) / hx;
             vy(i,j,k) = (vadv(i,j+1,k) - vadv(i,j,k)) / hy;
             wz(i,j,k) = (wadv(i,j,k+1) - wadv(i,j,k)) / hz;

       });
    }

    // compute sedgex on x-faces
    for ( MFIter mfi(umac_mf[0]); mfi.isValid(); ++mfi){

        const Box& bx = mfi.tilebox();

        Array4<const Real> const& s      = s_mf.array(mfi, comp);
        Array4<const Real> const& slope  = slope_mf.array(mfi);
        Array4<const Real> const& uadv  = umac_mf[0].array(mfi);
        Array4<const Real> const& vadv  = umac_mf[1].array(mfi);
        Array4<const Real> const& wadv  = umac_mf[2].array(mfi);

        //local variables
        Array4<      Real> const& ux    = ux_mf.array(mfi);
        Array4<      Real> const& vy    = vy_mf.array(mfi);
        Array4<      Real> const& wz    = wz_mf.array(mfi);
        Array4<      Real> const& sedgex = edges[0].array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

            //local variables
            Array1D<Real, 1, 3> del;
            Array1D<Real, 1, 3> p1;
            Array1D<Real, 1, 3> p2;
            Array1D<Real, 1, 3> p3;
            Array1D<Real, 1, 3> p4;

            Array1D<Real, 1, 7> slope_tmp;

            int ioff, joff, koff;

            Real isign,jsign,ksign;
            Real val1,val2,val3,val4,val5;
            Real u;
            Real uu,vv,ww;
            Real gamma,gamma2;

            ////////////////////////////////////////////////
            // compute sedgex without transverse corrections
            ////////////////////////////////////////////////

            if (uadv(i,j,k) > 0.0) {
               isign = 1.0;
               ioff = -1;
            } else {
               isign = -1.0;
               ioff = 0;
            }

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i+ioff,j,k,n-1);
            }


            // centroid of rectangular volume
            del(1) = isign*0.5*hx - 0.5*uadv(i,j,k)*dt;
            del(2) = 0.0;
            del(3) = 0.0;
            sedgex(i,j,k) = eval(s(i+ioff,j,k),slope_tmp,del);

            // source term
            sedgex(i,j,k) = sedgex(i,j,k) - dt2*sedgex(i,j,k)*ux(i+ioff,j,k);

            ////////////////////////////////////////////////
            // compute \Gamma^{y+} without corner corrections
            ////////////////////////////////////////////////

            if (vadv(i+ioff,j+1,k) > 0.0) {
               jsign = 1.0;
               joff = 0;
            } else {
               jsign = -1.0;
               joff = 1;
            }


            u = 0.0;
            if (uadv(i,j,k)*uadv(i,j+joff,k) > 0.0) {
               u = uadv(i,j+joff,k);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = jsign*0.5*hy;
            p1(3) = 0.0;

            p2(1) = isign*0.5*hx - uadv(i,j,k)*dt;
            p2(2) = jsign*0.5*hy;
            p2(3) = 0.0;

            p3(1) = isign*0.5*hx - u*dt;
            p3(2) = jsign*0.5*hy - vadv(i+ioff,j+1,k)*dt;
            p3(3) = 0.0;

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i+ioff,j+joff,k,n-1);
            }

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p2(ll)+p3(ll))/2.0;
            }
            val1 = eval(s(i+ioff,j+joff,k),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p3(ll))/2.0;
            }
            val2 = eval(s(i+ioff,j+joff,k),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p2(ll))/2.0;
            }
            val3 = eval(s(i+ioff,j+joff,k),slope_tmp,del);

            // average these centroid values to get the average value
            gamma = (val1+val2+val3)/3.0;

            // source term
            gamma = gamma - dt3*(gamma*ux(i+ioff,j+joff,k) + gamma*vy(i+ioff,j+joff,k));

            ////////////////////////////////////////////////
            // correct \Gamma^{y+} with \Gamma^{y+,z+}
            ////////////////////////////////////////////////

            if (wadv(i+ioff,j+joff,k+1) > 0.0) {
               ksign = 1.0;
               koff = 0;
            } else {
               ksign = -1.0;
               koff = 1;
            }

            uu = 0.0;
            if (uadv(i,j,k)*uadv(i,j+joff,k+koff) > 0.0) {
               uu = uadv(i,j+joff,k+koff);
            }

            vv = 0.0;
            if (vadv(i+ioff,j+1,k)*vadv(i+ioff,j+1,k+koff) > 0.0) {
               vv = vadv(i+ioff,j+1,k+koff);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = jsign*0.5*hy;
            p1(3) = ksign*0.5*hz;

            p2(1) = isign*0.5*hx - uadv(i,j,k)*dt;
            p2(2) = jsign*0.5*hy;
            p2(3) = ksign*0.5*hz;

            p3(1) = isign*0.5*hx - uadv(i,j,k)*dt;
            p3(2) = jsign*0.5*hy - vadv(i+ioff,j+1,k)*dt;
            p3(3) = ksign*0.5*hz;

            p4(1) = isign*0.5*hx - uu*dt;
            p4(2) = jsign*0.5*hy - vv*dt;
            p4(3) = ksign*0.5*hz - wadv(i+ioff,j+joff,k+1)*dt;

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i+ioff,j+joff,k+koff,n-1);
            }

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
            }
            val1 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
            }
            val2 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
            }
            val3 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
            }
            val4 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
            }
            val5 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

            // source term
            gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff)
                                     +gamma2*vy(i+ioff,j+joff,k+koff)
                                     +gamma2*wz(i+ioff,j+joff,k+koff));

            gamma2 = gamma2 * wadv(i+ioff,j+joff,k+1);

            gamma = gamma - dt*gamma2/(3.0*hz);

            ////////////////////////////////////////////////
            // correct \Gamma^{y+} with \Gamma^{y+,z-}
            ////////////////////////////////////////////////

            if (wadv(i+ioff,j+joff,k) > 0.0) {
               ksign = 1.0;
               koff = -1;
            } else {
               ksign = -1.0;
               koff = 0;
            }

            uu = 0.0;
            if (uadv(i,j,k)*uadv(i,j+joff,k+koff) > 0.0) {
               uu = uadv(i,j+joff,k+koff);
            }

            vv = 0.0;
            if (vadv(i+ioff,j+1,k)*vadv(i+ioff,j+1,k+koff) > 0.0) {
               vv = vadv(i+ioff,j+1,k+koff);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = jsign*0.5*hy;
            p1(3) = ksign*0.5*hz;

            p2(1) = isign*0.5*hx - uadv(i,j,k)*dt;
            p2(2) = jsign*0.5*hy;
            p2(3) = ksign*0.5*hz;

            p3(1) = isign*0.5*hx - uadv(i,j,k)*dt;
            p3(2) = jsign*0.5*hy - vadv(i+ioff,j+1,k)*dt;
            p3(3) = ksign*0.5*hz;

            p4(1) = isign*0.5*hx - uu*dt;
            p4(2) = jsign*0.5*hy - vv*dt;
            p4(3) = ksign*0.5*hz - wadv(i+ioff,j+joff,k)*dt;

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i+ioff,j+joff,k+koff,n-1);
            }

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
            }
            val1 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
            }
            val2 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
            }
            val3 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
            }
            val4 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
            }
            val5 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

            // source term
            gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff)
                                     +gamma2*vy(i+ioff,j+joff,k+koff)
                                     +gamma2*wz(i+ioff,j+joff,k+koff));

            gamma2 = gamma2 * wadv(i+ioff,j+joff,k);

            gamma = gamma + dt*gamma2/(3.0*hz);

            ////////////////////////////////////////////////
            // correct sedgex with \Gamma^{y+}
            ////////////////////////////////////////////////

            gamma = gamma * vadv(i+ioff,j+1,k);
            sedgex(i,j,k) = sedgex(i,j,k) - dt*gamma/(2.0*hy);

            ////////////////////////////////////////////////
            // compute \Gamma^{y-} without corner corrections
            ////////////////////////////////////////////////

            if (vadv(i+ioff,j,k) > 0.0) {
               jsign = 1.0;
               joff = -1;
            } else {
               jsign = -1.0;
               joff = 0;
            }

            u = 0.0;
            if (uadv(i,j,k)*uadv(i,j+joff,k) > 0.0) {
               u = uadv(i,j+joff,k);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = jsign*0.5*hy;
            p1(3) = 0.0;

            p2(1) = isign*0.5*hx - uadv(i,j,k)*dt;
            p2(2) = jsign*0.5*hy;
            p2(3) = 0.0;

            p3(1) = isign*0.5*hx - u*dt;
            p3(2) = jsign*0.5*hy - vadv(i+ioff,j,k)*dt;
            p3(3) = 0.0;

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i+ioff,j+joff,k,n-1);
            }

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p2(ll)+p3(ll))/2.0;
            }
            val1 = eval(s(i+ioff,j+joff,k),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p3(ll))/2.0;
            }
            val2 = eval(s(i+ioff,j+joff,k),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p2(ll))/2.0;
            }
            val3 = eval(s(i+ioff,j+joff,k),slope_tmp,del);

            // average these centroid values to get the average value
            gamma = (val1+val2+val3)/3.0;

            // source term
            gamma = gamma - dt3*(gamma*ux(i+ioff,j+joff,k) + gamma*vy(i+ioff,j+joff,k));

            ////////////////////////////////////////////////
            // correct \Gamma^{y-} with \Gamma^{y-,z+}
            ////////////////////////////////////////////////

            if (wadv(i+ioff,j+joff,k+1) > 0.0) {
               ksign = 1.0;
               koff = 0;
            } else {
               ksign = -1.0;
               koff = 1;
            }

            uu = 0.0;
            if (uadv(i,j,k)*uadv(i,j+joff,k+koff) > 0.0) {
               uu = uadv(i,j+joff,k+koff);
            }

            vv = 0.0;
            if (vadv(i+ioff,j,k)*vadv(i+ioff,j,k+koff) > 0.0) {
               vv = vadv(i+ioff,j,k+koff);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = jsign*0.5*hy;
            p1(3) = ksign*0.5*hz;

            p2(1) = isign*0.5*hx - uadv(i,j,k)*dt;
            p2(2) = jsign*0.5*hy;
            p2(3) = ksign*0.5*hz;

            p3(1) = isign*0.5*hx - uadv(i,j,k)*dt;
            p3(2) = jsign*0.5*hy - vadv(i+ioff,j,k)*dt;
            p3(3) = ksign*0.5*hz;

            p4(1) = isign*0.5*hx - uu*dt;
            p4(2) = jsign*0.5*hy - vv*dt;
            p4(3) = ksign*0.5*hz - wadv(i+ioff,j+joff,k+1)*dt;

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i+ioff,j+joff,k+koff,n-1);
            }

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
            }
            val1 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
            }
            val2 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
            }
            val3 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
            }
            val4 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
            }
            val5 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

            // source term
            gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff)
                                     +gamma2*vy(i+ioff,j+joff,k+koff)
                                     +gamma2*wz(i+ioff,j+joff,k+koff));

            gamma2 = gamma2 * wadv(i+ioff,j+joff,k+1);

            gamma = gamma - dt*gamma2/(3.0*hz);

            ////////////////////////////////////////////////
            // correct \Gamma^{y-} with \Gamma^{y-,z-}
            ////////////////////////////////////////////////

            if (wadv(i+ioff,j+joff,k) > 0.0) {
               ksign = 1.0;
               koff = -1;
            } else {
               ksign = -1.0;
               koff = 0;
            }

            uu = 0.0;
            if (uadv(i,j,k)*uadv(i,j+joff,k+koff) > 0.0) {
               uu = uadv(i,j+joff,k+koff);
            }

            vv = 0.0;
            if (vadv(i+ioff,j,k)*vadv(i+ioff,j,k+koff) > 0.0) {
               vv = vadv(i+ioff,j,k+koff);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = jsign*0.5*hy;
            p1(3) = ksign*0.5*hz;

            p2(1) = isign*0.5*hx - uadv(i,j,k)*dt;
            p2(2) = jsign*0.5*hy;
            p2(3) = ksign*0.5*hz;

            p3(1) = isign*0.5*hx - uadv(i,j,k)*dt;
            p3(2) = jsign*0.5*hy - vadv(i+ioff,j,k)*dt;
            p3(3) = ksign*0.5*hz;

            p4(1) = isign*0.5*hx - uu*dt;
            p4(2) = jsign*0.5*hy - vv*dt;
            p4(3) = ksign*0.5*hz - wadv(i+ioff,j+joff,k)*dt;

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i+ioff,j+joff,k+koff,n-1);
            }

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
            }
            val1 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
            }
            val2 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
            }
            val3 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
            }
            val4 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
            }
            val5 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

            // source term
            gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff)
                                     +gamma2*vy(i+ioff,j+joff,k+koff)
                                     +gamma2*wz(i+ioff,j+joff,k+koff));

            gamma2 = gamma2 * wadv(i+ioff,j+joff,k);

            gamma = gamma + dt*gamma2/(3.0*hz);

            ////////////////////////////////////////////////
            // correct sedgex with \Gamma^{y-}
            ////////////////////////////////////////////////

            gamma = gamma * vadv(i+ioff,j,k);
            sedgex(i,j,k) = sedgex(i,j,k) + dt*gamma/(2.0*hy);

            ////////////////////////////////////////////////
            // compute \Gamma^{z+} without corner corrections
            ////////////////////////////////////////////////

            if (wadv(i+ioff,j,k+1) > 0.0) {
               ksign = 1.0;
               koff = 0;
            } else {
               ksign = -1.0;
               koff = 1;
            }

            u = 0.0;
            if (uadv(i,j,k)*uadv(i,j,k+koff) > 0.0) {
               u = uadv(i,j,k+koff);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = 0.0;
            p1(3) = ksign*0.5*hz;

            p2(1) = isign*0.5*hx - uadv(i,j,k)*dt;
            p2(2) = 0.0;
            p2(3) = ksign*0.5*hz;

            p3(1) = isign*0.5*hx - u*dt;
            p3(2) = 0.0;
            p3(3) = ksign*0.5*hz - wadv(i+ioff,j,k+1)*dt;

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i+ioff,j,k+koff,n-1);
            }

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p2(ll)+p3(ll))/2.0;
            }
            val1 = eval(s(i+ioff,j,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p3(ll))/2.0;
            }
            val2 = eval(s(i+ioff,j,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p2(ll))/2.0;
            }
            val3 = eval(s(i+ioff,j,k+koff),slope_tmp,del);

            // average these centroid values to get the average value
            gamma = (val1+val2+val3)/3.0;

            // source term
            gamma = gamma - dt3*(gamma*ux(i+ioff,j,k+koff) + gamma*wz(i+ioff,j,k+koff));

            ////////////////////////////////////////////////
            // correct \Gamma^{z+} with \Gamma^{z+,y+}
            ////////////////////////////////////////////////

            if (vadv(i+ioff,j+1,k+koff) > 0.0) {
               jsign = 1.0;
               joff = 0;
            } else {
               jsign = -1.0;
               joff = 1;
            }

            uu = 0.0;
            if (uadv(i,j,k)*uadv(i,j+joff,k+koff) > 0.0) {
               uu = uadv(i,j+joff,k+koff);
            }

            ww = 0.0;
            if (wadv(i+ioff,j,k+1)*wadv(i+ioff,j+joff,k+1) > 0.0) {
               ww = wadv(i+ioff,j+joff,k+1);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = jsign*0.5*hy;
            p1(3) = ksign*0.5*hz;

            p2(1) = isign*0.5*hx - uadv(i,j,k)*dt;
            p2(2) = jsign*0.5*hy;
            p2(3) = ksign*0.5*hz;

            p3(1) = isign*0.5*hx - uadv(i,j,k)*dt;
            p3(2) = jsign*0.5*hy;
            p3(3) = ksign*0.5*hz - wadv(i+ioff,j,k+1)*dt;

            p4(1) = isign*0.5*hx - uu*dt;
            p4(2) = jsign*0.5*hy - vadv(i+ioff,j+1,k+koff)*dt;
            p4(3) = ksign*0.5*hz - ww*dt;

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i+ioff,j+joff,k+koff,n-1);
            }

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
            }
            val1 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
            }
            val2 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
            }
            val3 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
            }
            val4 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
            }
            val5 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

            // source term
            gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff)
                                     +gamma2*vy(i+ioff,j+joff,k+koff)
                                     +gamma2*wz(i+ioff,j+joff,k+koff));

            gamma2 = gamma2 * vadv(i+ioff,j+1,k+koff);

            gamma = gamma - dt*gamma2/(3.0*hy);

            ////////////////////////////////////////////////
            // correct \Gamma^{z+} with \Gamma^{z+,y-}
            ////////////////////////////////////////////////

            if (vadv(i+ioff,j,k+koff) > 0.0) {
               jsign = 1.0;
               joff = -1;
            } else {
               jsign = -1.0;
               joff = 0;
            }

            uu = 0.0;
            if (uadv(i,j,k)*uadv(i,j+joff,k+koff) > 0.0) {
               uu = uadv(i,j+joff,k+koff);
            }

            ww = 0.0;
            if (wadv(i+ioff,j,k+1)*wadv(i+ioff,j+joff,k+1) > 0.0) {
               ww = wadv(i+ioff,j+joff,k+1);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = jsign*0.5*hy;
            p1(3) = ksign*0.5*hz;

            p2(1) = isign*0.5*hx - uadv(i,j,k)*dt;
            p2(2) = jsign*0.5*hy;
            p2(3) = ksign*0.5*hz;

            p3(1) = isign*0.5*hx - uadv(i,j,k)*dt;
            p3(2) = jsign*0.5*hy;
            p3(3) = ksign*0.5*hz - wadv(i+ioff,j,k+1)*dt;

            p4(1) = isign*0.5*hx - uu*dt;
            p4(2) = jsign*0.5*hy - vadv(i+ioff,j,k+koff)*dt;
            p4(3) = ksign*0.5*hz - ww*dt;

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i+ioff,j+joff,k+koff,n-1);
            }

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
            }
            val1 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
            }
            val2 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
            }
            val3 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
            }
            val4 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
            }
            val5 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

            // source term
            gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff)
                                     +gamma2*vy(i+ioff,j+joff,k+koff)
                                     +gamma2*wz(i+ioff,j+joff,k+koff));

            gamma2 = gamma2 * vadv(i+ioff,j,k+koff);

            gamma = gamma + dt*gamma2/(3.0*hy);

            ////////////////////////////////////////////////
            // correct sedgex with \Gamma^{z+}
            ////////////////////////////////////////////////

            gamma = gamma * wadv(i+ioff,j,k+1);
            sedgex(i,j,k) = sedgex(i,j,k) - dt*gamma/(2.0*hz);

            ////////////////////////////////////////////////
            // compute \Gamma^{z-} without corner corrections
            ////////////////////////////////////////////////

            if (wadv(i+ioff,j,k) > 0.0) {
               ksign = 1.0;
               koff = -1;
            } else {
               ksign = -1.0;
               koff = 0;
            }

            u = 0.0;
            if (uadv(i,j,k)*uadv(i,j,k+koff) > 0.0) {
               u = uadv(i,j,k+koff);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = 0.0;
            p1(3) = ksign*0.5*hz;

            p2(1) = isign*0.5*hx - uadv(i,j,k)*dt;
            p2(2) = 0.0;
            p2(3) = ksign*0.5*hz;

            p3(1) = isign*0.5*hx - u*dt;
            p3(2) = 0.0;
            p3(3) = ksign*0.5*hz - wadv(i+ioff,j,k)*dt;

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i+ioff,j,k+koff,n-1);
            }

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p2(ll)+p3(ll))/2.0;
            }
            val1 = eval(s(i+ioff,j,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p3(ll))/2.0;
            }
            val2 = eval(s(i+ioff,j,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p2(ll))/2.0;
            }
            val3 = eval(s(i+ioff,j,k+koff),slope_tmp,del);

            // average these centroid values to get the average value
            gamma = (val1+val2+val3)/3.0;

            // source term
            gamma = gamma - dt3*(gamma*ux(i+ioff,j,k+koff) + gamma*wz(i+ioff,j,k+koff));

            ////////////////////////////////////////////////
            // correct \Gamma^{z-} with \Gamma^{z-,y+}
            ////////////////////////////////////////////////

            if (vadv(i+ioff,j+1,k+koff) > 0.0) {
               jsign = 1.0;
               joff = 0;
            } else {
               jsign = -1.0;
               joff = 1;
            }

            uu = 0.0;
            if (uadv(i,j,k)*uadv(i,j+joff,k+koff) > 0.0) {
               uu = uadv(i,j+joff,k+koff);
            }

            ww = 0.0;
            if (wadv(i+ioff,j,k)*wadv(i+ioff,j+joff,k) > 0.0) {
               ww = wadv(i+ioff,j+joff,k);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = jsign*0.5*hy;
            p1(3) = ksign*0.5*hz;

            p2(1) = isign*0.5*hx - uadv(i,j,k)*dt;
            p2(2) = jsign*0.5*hy;
            p2(3) = ksign*0.5*hz;

            p3(1) = isign*0.5*hx - uadv(i,j,k)*dt;
            p3(2) = jsign*0.5*hy;
            p3(3) = ksign*0.5*hz - wadv(i+ioff,j,k)*dt;

            p4(1) = isign*0.5*hx - uu*dt;
            p4(2) = jsign*0.5*hy - vadv(i+ioff,j+1,k+koff)*dt;
            p4(3) = ksign*0.5*hz - ww*dt;

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i+ioff,j+joff,k+koff,n-1);
            }

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
            }
            val1 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
            }
            val2 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
            }
            val3 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
            }
            val4 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
            }
            val5 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

            // source term
            gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff)
                                     +gamma2*vy(i+ioff,j+joff,k+koff)
                                     +gamma2*wz(i+ioff,j+joff,k+koff));

            gamma2 = gamma2 * vadv(i+ioff,j+1,k+koff);

            gamma = gamma - dt*gamma2/(3.0*hy);

            ////////////////////////////////////////////////
            // correct \Gamma^{z-} with \Gamma^{z-,y-}
            ////////////////////////////////////////////////

            if (vadv(i+ioff,j,k+koff) > 0.0) {
               jsign = 1.0;
               joff = -1;
            } else {
               jsign = -1.0;
               joff = 0;
            }

            uu = 0.0;
            if (uadv(i,j,k)*uadv(i,j+joff,k+koff) > 0.0) {
               uu = uadv(i,j+joff,k+koff);
            }

            ww = 0.0;
            if (wadv(i+ioff,j,k)*wadv(i+ioff,j+joff,k) > 0.0) {
               ww = wadv(i+ioff,j+joff,k);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = jsign*0.5*hy;
            p1(3) = ksign*0.5*hz;

            p2(1) = isign*0.5*hx - uadv(i,j,k)*dt;
            p2(2) = jsign*0.5*hy;
            p2(3) = ksign*0.5*hz;

            p3(1) = isign*0.5*hx - uadv(i,j,k)*dt;
            p3(2) = jsign*0.5*hy;
            p3(3) = ksign*0.5*hz - wadv(i+ioff,j,k)*dt;

            p4(1) = isign*0.5*hx - uu*dt;
            p4(2) = jsign*0.5*hy - vadv(i+ioff,j,k+koff)*dt;
            p4(3) = ksign*0.5*hz - ww*dt;

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i+ioff,j+joff,k+koff,n-1);
            }

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
            }
            val1 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
            }
            val2 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
            }
            val3 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
            }
            val4 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
            }
            val5 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

            // source term
            gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff)
                                     +gamma2*vy(i+ioff,j+joff,k+koff)
                                     +gamma2*wz(i+ioff,j+joff,k+koff));

            gamma2 = gamma2 * vadv(i+ioff,j,k+koff);

            gamma = gamma + dt*gamma2/(3.0*hy);

            ////////////////////////////////////////////////
            // correct sedgex with \Gamma^{z-}
            ////////////////////////////////////////////////

            gamma = gamma * wadv(i+ioff,j,k);
            sedgex(i,j,k) = sedgex(i,j,k) + dt*gamma/(2.0*hz);

        });
    }



    // compute sedgey on y-faces
    for ( MFIter mfi(umac_mf[1]); mfi.isValid(); ++mfi){

        const Box& bx = mfi.tilebox();

        Array4<const Real> const& s      = s_mf.array(mfi, comp);
        Array4<const Real> const& slope  = slope_mf.array(mfi);
        Array4<const Real> const& uadv  = umac_mf[0].array(mfi);
        Array4<const Real> const& vadv  = umac_mf[1].array(mfi);
        Array4<const Real> const& wadv  = umac_mf[2].array(mfi);

        //local variables
        Array4<      Real> const& ux    = ux_mf.array(mfi);
        Array4<      Real> const& vy    = vy_mf.array(mfi);
        Array4<      Real> const& wz    = wz_mf.array(mfi);
        Array4<      Real> const& sedgey = edges[1].array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

            //local variables
            Array1D<Real, 1, 3> del;
            Array1D<Real, 1, 3> p1;
            Array1D<Real, 1, 3> p2;
            Array1D<Real, 1, 3> p3;
            Array1D<Real, 1, 3> p4;

            Array1D<Real, 1, 7> slope_tmp;

            int ioff, joff, koff;

            Real isign,jsign,ksign;
            Real val1,val2,val3,val4,val5;
            Real v;
            Real uu,vv,ww;
            Real gamma,gamma2;

            ////////////////////////////////////////////////
            // compute sedgey without transverse corrections
            ////////////////////////////////////////////////

            // centroid of rectangular volume
            if (vadv(i,j,k) > 0.0) {
               jsign = 1.0;
               joff = -1;
            } else {
               jsign = -1.0;
               joff = 0;
            }

            del(1) = 0.0;
            del(2) = jsign*0.5*hy - 0.5*vadv(i,j,k)*dt;
            del(3) = 0.0;

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i,j+joff,k,n-1);
            }

            sedgey(i,j,k) = eval(s(i,j+joff,k),slope_tmp,del);

            // source term
            sedgey(i,j,k) = sedgey(i,j,k) - dt2*sedgey(i,j,k)*vy(i,j+joff,k);

            ////////////////////////////////////////////////
            // compute \Gamma^{x+} without corner corrections
            ////////////////////////////////////////////////

            if (uadv(i+1,j+joff,k) > 0.0) {
               isign = 1.0;
               ioff = 0;
            } else {
               isign = -1.0;
               ioff = 1;
            }

            v = 0.0;
            if (vadv(i,j,k)*vadv(i+ioff,j,k) > 0.0) {
               v = vadv(i+ioff,j,k);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = jsign*0.5*hy;
            p1(3) = 0.0;

            p2(1) = isign*0.5*hx;
            p2(2) = jsign*0.5*hy - vadv(i,j,k)*dt;
            p2(3) = 0.0;

            p3(1) = isign*0.5*hx - uadv(i+1,j+joff,k)*dt;
            p3(2) = jsign*0.5*hy - v*dt;
            p3(3) = 0.0;

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i+ioff,j+joff,k,n-1);
            }

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p2(ll)+p3(ll))/2.0;
            }
            val1 = eval(s(i+ioff,j+joff,k),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p3(ll))/2.0;
            }
            val2 = eval(s(i+ioff,j+joff,k),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p2(ll))/2.0;
            }
            val3 = eval(s(i+ioff,j+joff,k),slope_tmp,del);

            // average these centroid values to get the average value
            gamma = (val1+val2+val3)/3.0;

            // source term
            gamma = gamma - dt3*(gamma*vy(i+ioff,j+joff,k) + gamma*ux(i+ioff,j+joff,k));

            ////////////////////////////////////////////////
            // correct \Gamma^{x+} with \Gamma^{x+,z+}
            ////////////////////////////////////////////////

            if (wadv(i+ioff,j+joff,k+1) > 0.0) {
               ksign = 1.0;
               koff = 0;
            } else {
               ksign = -1.0;
               koff = 1;
            }

            vv = 0.0;
            if (vadv(i,j,k)*vadv(i+ioff,j,k+koff) > 0.0) {
               vv = vadv(i+ioff,j,k+koff);
            }

            uu = 0.0;
            if (uadv(i+1,j+joff,k)*uadv(i+1,j+joff,k+koff) > 0.0) {
               uu = uadv(i+1,j+joff,k+koff);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = jsign*0.5*hy;
            p1(3) = ksign*0.5*hz;

            p2(1) = isign*0.5*hx;
            p2(2) = jsign*0.5*hy - vadv(i,j,k)*dt;
            p2(3) = ksign*0.5*hz;

            p3(1) = isign*0.5*hx - uadv(i+1,j+joff,k)*dt;
            p3(2) = jsign*0.5*hy - vadv(i,j,k)*dt;
            p3(3) = ksign*0.5*hz;

            p4(1) = isign*0.5*hx - uu*dt;
            p4(2) = jsign*0.5*hy - vv*dt;
            p4(3) = ksign*0.5*hz - wadv(i+ioff,j+joff,k+1)*dt;

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i+ioff,j+joff,k+koff,n-1);
            }

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
            }
            val1 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
            }
            val2 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
            }
            val3 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
            }
            val4 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
            }
            val5 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

            // source term
            gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff)
                                     +gamma2*vy(i+ioff,j+joff,k+koff)
                                     +gamma2*wz(i+ioff,j+joff,k+koff));

            gamma2 = gamma2 * wadv(i+ioff,j+joff,k+1);

            gamma = gamma - dt*gamma2/(3.0*hz);

            ////////////////////////////////////////////////
            // correct \Gamma^{x+} with \Gamma^{x+,z-}
            ////////////////////////////////////////////////

            if (wadv(i+ioff,j+joff,k) > 0.0) {
               ksign = 1.0;
               koff = -1;
            } else {
               ksign = -1.0;
               koff = 0;
            }

            vv = 0.0;
            if (vadv(i,j,k)*vadv(i+ioff,j,k+koff) > 0.0) {
               vv = vadv(i+ioff,j,k+koff);
            }

            uu = 0.0;
            if (uadv(i+1,j+joff,k)*uadv(i+1,j+joff,k+koff) > 0.0) {
               uu = uadv(i+1,j+joff,k+koff);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = jsign*0.5*hy;
            p1(3) = ksign*0.5*hz;

            p2(1) = isign*0.5*hx;
            p2(2) = jsign*0.5*hy - vadv(i,j,k)*dt;
            p2(3) = ksign*0.5*hz;

            p3(1) = isign*0.5*hx - uadv(i+1,j+joff,k)*dt;
            p3(2) = jsign*0.5*hy - vadv(i,j,k)*dt;
            p3(3) = ksign*0.5*hz;

            p4(1) = isign*0.5*hx - uu*dt;
            p4(2) = jsign*0.5*hy - vv*dt;
            p4(3) = ksign*0.5*hz - wadv(i+ioff,j+joff,k)*dt;

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i+ioff,j+joff,k+koff,n-1);
            }

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
            }
            val1 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
            }
            val2 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
            }
            val3 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
            }
            val4 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
            }
            val5 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

            // source term
            gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff)
                                     +gamma2*vy(i+ioff,j+joff,k+koff)
                                     +gamma2*wz(i+ioff,j+joff,k+koff));

            gamma2 = gamma2 * wadv(i+ioff,j+joff,k);

            gamma = gamma + dt*gamma2/(3.0*hz);

            ////////////////////////////////////////////////
            // correct sedgey with \Gamma^{x+}
            ////////////////////////////////////////////////

            gamma = gamma * uadv(i+1,j+joff,k);
            sedgey(i,j,k) = sedgey(i,j,k) - dt*gamma/(2.0*hx);

            ////////////////////////////////////////////////
            // compute \Gamma^{x-} without corner corrections
            ////////////////////////////////////////////////

            if (uadv(i,j+joff,k) > 0.0) {
               isign = 1.0;
               ioff = -1;
            } else {
               isign = -1.0;
               ioff = 0;
            }

            v = 0.0;
            if (vadv(i,j,k)*vadv(i+ioff,j,k) > 0.0) {
               v = vadv(i+ioff,j,k);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = jsign*0.5*hy;
            p1(3) = 0.0;

            p2(1) = isign*0.5*hx;
            p2(2) = jsign*0.5*hy - vadv(i,j,k)*dt;
            p2(3) = 0.0;

            p3(1) = isign*0.5*hx - uadv(i,j+joff,k)*dt;
            p3(2) = jsign*0.5*hy - v*dt;
            p3(3) = 0.0;

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i+ioff,j+joff,k,n-1);
            }

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p2(ll)+p3(ll))/2.0;
            }
            val1 = eval(s(i+ioff,j+joff,k),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p3(ll))/2.0;
            }
            val2 = eval(s(i+ioff,j+joff,k),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p2(ll))/2.0;
            }
            val3 = eval(s(i+ioff,j+joff,k),slope_tmp,del);

            // average these centroid values to get the average value
            gamma = (val1+val2+val3)/3.0;

            // source term
            gamma = gamma - dt3*(gamma*vy(i+ioff,j+joff,k) + gamma*ux(i+ioff,j+joff,k));

            ////////////////////////////////////////////////
            // correct \Gamma^{x-} with \Gamma^{x-,z+}
            ////////////////////////////////////////////////

            if (wadv(i+ioff,j+joff,k+1) > 0.0) {
               ksign = 1.0;
               koff = 0;
            } else {
               ksign = -1.0;
               koff = 1;
            }

            vv = 0.0;
            if (vadv(i,j,k)*vadv(i+ioff,j,k+koff) > 0.0) {
               vv = vadv(i+ioff,j,k+koff);
            }

            uu = 0.0;
            if (uadv(i,j+joff,k)*uadv(i,j+joff,k+koff) > 0.0) {
               uu = uadv(i,j+joff,k+koff);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = jsign*0.5*hy;
            p1(3) = ksign*0.5*hz;

            p2(1) = isign*0.5*hx;
            p2(2) = jsign*0.5*hy - vadv(i,j,k)*dt;
            p2(3) = ksign*0.5*hz;

            p3(1) = isign*0.5*hx - uadv(i,j+joff,k)*dt;
            p3(2) = jsign*0.5*hy - vadv(i,j,k)*dt;
            p3(3) = ksign*0.5*hz;

            p4(1) = isign*0.5*hx - uu*dt;
            p4(2) = jsign*0.5*hy - vv*dt;
            p4(3) = ksign*0.5*hz - wadv(i+ioff,j+joff,k+1)*dt;

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i+ioff,j+joff,k+koff,n-1);
            }

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
            }
            val1 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
            }
            val2 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
            }
            val3 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
            }
            val4 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
            }
            val5 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

            // source term
            gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff)
                                     +gamma2*vy(i+ioff,j+joff,k+koff)
                                     +gamma2*wz(i+ioff,j+joff,k+koff));

            gamma2 = gamma2 * wadv(i+ioff,j+joff,k+1);

            gamma = gamma - dt*gamma2/(3.0*hz);

            ////////////////////////////////////////////////
            // correct \Gamma^{x-} with \Gamma^{x-,z-}
            ////////////////////////////////////////////////

            if (wadv(i+ioff,j+joff,k) > 0.0) {
               ksign = 1.0;
               koff = -1;
            } else {
               ksign = -1.0;
               koff = 0;
            }

            vv = 0.0;
            if (vadv(i,j,k)*vadv(i+ioff,j,k+koff) > 0.0) {
               vv = vadv(i+ioff,j,k+koff);
            }

            uu = 0.0;
            if (uadv(i,j+joff,k)*uadv(i,j+joff,k+koff) > 0.0) {
               uu = uadv(i,j+joff,k+koff);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = jsign*0.5*hy;
            p1(3) = ksign*0.5*hz;

            p2(1) = isign*0.5*hx;
            p2(2) = jsign*0.5*hy - vadv(i,j,k)*dt;
            p2(3) = ksign*0.5*hz;

            p3(1) = isign*0.5*hx - uadv(i,j+joff,k)*dt;
            p3(2) = jsign*0.5*hy - vadv(i,j,k)*dt;
            p3(3) = ksign*0.5*hz;

            p4(1) = isign*0.5*hx - uu*dt;
            p4(2) = jsign*0.5*hy - vv*dt;
            p4(3) = ksign*0.5*hz - wadv(i+ioff,j+joff,k)*dt;

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i+ioff,j+joff,k+koff,n-1);
            }

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
            }
            val1 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
            }
            val2 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
            }
            val3 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
            }
            val4 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
            }
            val5 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

            // source term
            gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff)
                                     +gamma2*vy(i+ioff,j+joff,k+koff)
                                     +gamma2*wz(i+ioff,j+joff,k+koff));

            gamma2 = gamma2 * wadv(i+ioff,j+joff,k);

            gamma = gamma + dt*gamma2/(3.0*hz);

            ////////////////////////////////////////////////
            // correct sedgey with \Gamma^{x-}
            ////////////////////////////////////////////////

            gamma = gamma * uadv(i,j+joff,k);
            sedgey(i,j,k) = sedgey(i,j,k) + dt*gamma/(2.0*hx);

            ////////////////////////////////////////////////
            // compute \Gamma^{z+} without corner corrections
            ////////////////////////////////////////////////

            if (wadv(i,j+joff,k+1) > 0.0) {
               ksign = 1.0;
               koff = 0;
            } else {
               ksign = -1.0;
               koff = 1;
            }

            v = 0.0;
            if (vadv(i,j,k)*vadv(i,j,k+koff) > 0.0) {
               v = vadv(i,j,k+koff);
            }

            p1(1) = 0.0;
            p1(2) = jsign*0.5*hy;
            p1(3) = ksign*0.5*hz;

            p2(1) = 0.0;
            p2(2) = jsign*0.5*hy - vadv(i,j,k)*dt;
            p2(3) = ksign*0.5*hz;

            p3(1) = 0.0;
            p3(2) = jsign*0.5*hy - v*dt;
            p3(3) = ksign*0.5*hz - wadv(i,j+joff,k+1)*dt;

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i,j+joff,k+koff,n-1);
            }

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p2(ll)+p3(ll))/2.0;
            }
            val1 = eval(s(i,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p3(ll))/2.0;
            }
            val2 = eval(s(i,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p2(ll))/2.0;
            }
            val3 = eval(s(i,j+joff,k+koff),slope_tmp,del);

            // average these centroid values to get the average value
            gamma = (val1+val2+val3)/3.0;

            // source term
            gamma = gamma - dt3*(gamma*vy(i,j+joff,k+koff) + gamma*wz(i,j+joff,k+koff));

            ////////////////////////////////////////////////
            // correct \Gamma^{z+} with \Gamma^{z+,x+}
            ////////////////////////////////////////////////

            if (uadv(i+1,j+joff,k+koff) > 0.0) {
               isign = 1.0;
               ioff = 0;
            } else {
               isign = -1.0;
               ioff = 1;
            }

            vv = 0.0;
            if (vadv(i,j,k)*vadv(i+ioff,j,k+koff) > 0.0) {
               vv = vadv(i+ioff,j,k+koff);
            }

            ww = 0.0;
            if (wadv(i,j+joff,k+1)*wadv(i+ioff,j+joff,k+1) > 0.0) {
               ww = wadv(i+ioff,j+joff,k+1);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = jsign*0.5*hy;
            p1(3) = ksign*0.5*hz;

            p2(1) = isign*0.5*hx;
            p2(2) = jsign*0.5*hy - vadv(i,j,k)*dt;
            p2(3) = ksign*0.5*hz;

            p3(1) = isign*0.5*hx;
            p3(2) = jsign*0.5*hy - vadv(i,j,k)*dt;
            p3(3) = ksign*0.5*hz - wadv(i,j+joff,k+1)*dt;

            p4(1) = isign*0.5*hx - uadv(i+1,j+joff,k+koff)*dt;
            p4(2) = jsign*0.5*hy - vv*dt;
            p4(3) = ksign*0.5*hz - ww*dt;

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i+ioff,j+joff,k+koff,n-1);
            }

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
            }
            val1 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
            }
            val2 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
            }
            val3 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
            }
            val4 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
            }
            val5 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

            // source term
            gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff)
                                     +gamma2*vy(i+ioff,j+joff,k+koff)
                                     +gamma2*wz(i+ioff,j+joff,k+koff));

            gamma2 = gamma2 * uadv(i+1,j+joff,k+koff);

            gamma = gamma - dt*gamma2/(3.0*hx);

            ////////////////////////////////////////////////
            // correct \Gamma^{z+} with \Gamma^{z+,x-}
            ////////////////////////////////////////////////

            if (uadv(i,j+joff,k+koff) > 0.0) {
               isign = 1.0;
               ioff = -1;
            } else {
               isign = -1.0;
               ioff = 0;
            }

            vv = 0.0;
            if (vadv(i,j,k)*vadv(i+ioff,j,k+koff) > 0.0) {
               vv = vadv(i+ioff,j,k+koff);
            }

            ww = 0.0;
            if (wadv(i,j+joff,k+1)*wadv(i+ioff,j+joff,k+1) > 0.0) {
               ww = wadv(i+ioff,j+joff,k+1);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = jsign*0.5*hy;
            p1(3) = ksign*0.5*hz;

            p2(1) = isign*0.5*hx;
            p2(2) = jsign*0.5*hy - vadv(i,j,k)*dt;
            p2(3) = ksign*0.5*hz;

            p3(1) = isign*0.5*hx;
            p3(2) = jsign*0.5*hy - vadv(i,j,k)*dt;
            p3(3) = ksign*0.5*hz - wadv(i,j+joff,k+1)*dt;

            p4(1) = isign*0.5*hx - uadv(i,j+joff,k+koff)*dt;
            p4(2) = jsign*0.5*hy - vv*dt;
            p4(3) = ksign*0.5*hz - ww*dt;

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i+ioff,j+joff,k+koff,n-1);
            }

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
            }
            val1 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
            }
            val2 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
            }
            val3 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
            }
            val4 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
            }
            val5 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

            // source term
            gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff)
                                     +gamma2*vy(i+ioff,j+joff,k+koff)
                                     +gamma2*wz(i+ioff,j+joff,k+koff));

            gamma2 = gamma2 * uadv(i,j+joff,k+koff);

            gamma = gamma + dt*gamma2/(3.0*hx);

            ////////////////////////////////////////////////
            // correct sedgey with \Gamma^{z+}
            ////////////////////////////////////////////////

            gamma = gamma * wadv(i,j+joff,k+1);
            sedgey(i,j,k) = sedgey(i,j,k) - dt*gamma/(2.0*hz);

            ////////////////////////////////////////////////
            // compute \Gamma^{z-} without corner corrections
            ////////////////////////////////////////////////

            if (wadv(i,j+joff,k) > 0.0) {
               ksign = 1.0;
               koff = -1;
            } else {
               ksign = -1.0;
               koff = 0;
            }

            v = 0.0;
            if (vadv(i,j,k)*vadv(i,j,k+koff) > 0.0) {
               v = vadv(i,j,k+koff);
            }

            p1(1) = 0.0;
            p1(2) = jsign*0.5*hy;
            p1(3) = ksign*0.5*hz;

            p2(1) = 0.0;
            p2(2) = jsign*0.5*hy - vadv(i,j,k)*dt;
            p2(3) = ksign*0.5*hz;

            p3(1) = 0.0;
            p3(2) = jsign*0.5*hy - v*dt;
            p3(3) = ksign*0.5*hz - wadv(i,j+joff,k)*dt;

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i,j+joff,k+koff,n-1);
            }

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p2(ll)+p3(ll))/2.0;
            }
            val1 = eval(s(i,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p3(ll))/2.0;
            }
            val2 = eval(s(i,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p2(ll))/2.0;
            }
            val3 = eval(s(i,j+joff,k+koff),slope_tmp,del);

            // average these centroid values to get the average value
            gamma = (val1+val2+val3)/3.0;

            // source term
            gamma = gamma - dt3*(gamma*vy(i,j+joff,k+koff) + gamma*wz(i,j+joff,k+koff));

            ////////////////////////////////////////////////
            // correct \Gamma^{z-} with \Gamma^{z-,x+}
            ////////////////////////////////////////////////

            if (uadv(i+1,j+joff,k+koff) > 0.0) {
               isign = 1.0;
               ioff = 0;
            } else {
               isign = -1.0;
               ioff = 1;
            }

            vv = 0.0;
            if (vadv(i,j,k)*vadv(i+ioff,j,k+koff) > 0.0) {
               vv = vadv(i+ioff,j,k+koff);
            }

            ww = 0.0;
            if (wadv(i,j+joff,k)*wadv(i+ioff,j+joff,k) > 0.0) {
               ww = wadv(i+ioff,j+joff,k);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = jsign*0.5*hy;
            p1(3) = ksign*0.5*hz;

            p2(1) = isign*0.5*hx;
            p2(2) = jsign*0.5*hy - vadv(i,j,k)*dt;
            p2(3) = ksign*0.5*hz;

            p3(1) = isign*0.5*hx;
            p3(2) = jsign*0.5*hy - vadv(i,j,k)*dt;
            p3(3) = ksign*0.5*hz - wadv(i,j+joff,k)*dt;

            p4(1) = isign*0.5*hx - uadv(i+1,j+joff,k+koff)*dt;
            p4(2) = jsign*0.5*hy - vv*dt;
            p4(3) = ksign*0.5*hz - ww*dt;

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i+ioff,j+joff,k+koff,n-1);
            }

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
            }
            val1 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
            }
            val2 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
            }
            val3 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
            }
            val4 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
            }
            val5 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

            // source term
            gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff)
                                     +gamma2*vy(i+ioff,j+joff,k+koff)
                                     +gamma2*wz(i+ioff,j+joff,k+koff));

            gamma2 = gamma2 * uadv(i+1,j+joff,k+koff);

            gamma = gamma - dt*gamma2/(3.0*hx);

            ////////////////////////////////////////////////
            // correct \Gamma^{z-} with \Gamma^{z-,x-}
            ////////////////////////////////////////////////

            if (uadv(i,j+joff,k+koff) > 0.0) {
               isign = 1.0;
               ioff = -1;
            } else {
               isign = -1.0;
               ioff = 0;
            }

            vv = 0.0;
            if (vadv(i,j,k)*vadv(i+ioff,j,k+koff) > 0.0) {
               vv = vadv(i+ioff,j,k+koff);
            }

            ww = 0.0;
            if (wadv(i,j+joff,k)*wadv(i+ioff,j+joff,k) > 0.0) {
               ww = wadv(i+ioff,j+joff,k);
            }

            p1(1) = isign*0.5*hx;
            p1(2) = jsign*0.5*hy;
            p1(3) = ksign*0.5*hz;

            p2(1) = isign*0.5*hx;
            p2(2) = jsign*0.5*hy - vadv(i,j,k)*dt;
            p2(3) = ksign*0.5*hz;

            p3(1) = isign*0.5*hx;
            p3(2) = jsign*0.5*hy - vadv(i,j,k)*dt;
            p3(3) = ksign*0.5*hz - wadv(i,j+joff,k)*dt;

            p4(1) = isign*0.5*hx - uadv(i,j+joff,k+koff)*dt;
            p4(2) = jsign*0.5*hy - vv*dt;
            p4(3) = ksign*0.5*hz - ww*dt;

            for(int n=1; n<=7; ++n){
                slope_tmp(n) = slope(i+ioff,j+joff,k+koff,n-1);
            }

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
            }
            val1 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
            }
            val2 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
            }
            val3 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
            }
            val4 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            for(int ll=1; ll<=3; ++ll ){
               del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
            }
            val5 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

            gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

            // source term
            gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff)
                                     +gamma2*vy(i+ioff,j+joff,k+koff)
                                     +gamma2*wz(i+ioff,j+joff,k+koff));

            gamma2 = gamma2 * uadv(i,j+joff,k+koff);

            gamma = gamma + dt*gamma2/(3.0*hx);

            ////////////////////////////////////////////////
            // correct sedgey with \Gamma^{z-}
            ////////////////////////////////////////////////

            gamma = gamma * wadv(i,j+joff,k);
            sedgey(i,j,k) = sedgey(i,j,k) + dt*gamma/(2.0*hz);

        });
    }

    // compute sedgez on z-faces
    for ( MFIter mfi(umac_mf[2]); mfi.isValid(); ++mfi){

        const Box& bx = mfi.tilebox();

        Array4<const Real> const& s      = s_mf.array(mfi, comp);
        Array4<const Real> const& slope  = slope_mf.array(mfi);
        Array4<const Real> const& uadv  = umac_mf[0].array(mfi);
        Array4<const Real> const& vadv  = umac_mf[1].array(mfi);
        Array4<const Real> const& wadv  = umac_mf[2].array(mfi);

        //local variables
        Array4<      Real> const& ux    = ux_mf.array(mfi);
        Array4<      Real> const& vy    = vy_mf.array(mfi);
        Array4<      Real> const& wz    = wz_mf.array(mfi);
        Array4<      Real> const& sedgez = edges[2].array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

            //local variables
            Array1D<Real, 1, 3> del;
            Array1D<Real, 1, 3> p1;
            Array1D<Real, 1, 3> p2;
            Array1D<Real, 1, 3> p3;
            Array1D<Real, 1, 3> p4;

            Array1D<Real, 1, 7> slope_tmp;

            int ioff, joff, koff;

            Real isign,jsign,ksign;
            Real val1,val2,val3,val4,val5;
            Real w;
            Real uu,vv,ww;
            Real gamma,gamma2;

             ////////////////////////////////////////////////
             // compute sedgez without transverse corrections
             ////////////////////////////////////////////////

             // centroid of rectangular volume
             if (wadv(i,j,k) > 0.0) {
                ksign = 1.0;
                koff = -1;
             } else {
                ksign = -1.0;
                koff = 0;
             }

             for(int n=1; n<=7; ++n){
                 slope_tmp(n) = slope(i,j,k+koff,n-1);
             }

             del(1) = 0.0;
             del(2) = 0.0;
             del(3) = ksign*0.5*hz - 0.5*wadv(i,j,k)*dt;
             sedgez(i,j,k) = eval(s(i,j,k+koff),slope_tmp,del);

             // source term
             sedgez(i,j,k) = sedgez(i,j,k) - dt2*sedgez(i,j,k)*wz(i,j,k+koff);

             ////////////////////////////////////////////////
             // compute \Gamma^{x+} without corner corrections
             ////////////////////////////////////////////////

             if (uadv(i+1,j,k+koff) > 0.0) {
                isign = 1.0;
                ioff = 0;
             } else {
                isign = -1.0;
                ioff = 1;
             }

             w = 0.0;
             if (wadv(i,j,k)*wadv(i+ioff,j,k) > 0.0) {
                w = wadv(i+ioff,j,k);
             }

             p1(1) = isign*0.5*hx;
             p1(2) = 0.0;
             p1(3) = ksign*0.5*hz;

             p2(1) = isign*0.5*hx;
             p2(2) = 0.0;
             p2(3) = ksign*0.5*hz - wadv(i,j,k)*dt;

             p3(1) = isign*0.5*hx - uadv(i+1,j,k+koff)*dt;
             p3(2) = 0.0;
             p3(3) = ksign*0.5*hz - w*dt;

             for(int n=1; n<=7; ++n){
                 slope_tmp(n) = slope(i+ioff,j,k+koff,n-1);
             }

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = (p2(ll)+p3(ll))/2.0;
             }
             val1 = eval(s(i+ioff,j,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = (p1(ll)+p3(ll))/2.0;
             }
             val2 = eval(s(i+ioff,j,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = (p1(ll)+p2(ll))/2.0;
             }
             val3 = eval(s(i+ioff,j,k+koff),slope_tmp,del);

             // average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.0;

             // source term
             gamma = gamma - dt3*(gamma*wz(i+ioff,j,k+koff) + gamma*ux(i+ioff,j,k+koff));

             ////////////////////////////////////////////////
             // correct \Gamma^{x+} with \Gamma^{x+,y+}
             ////////////////////////////////////////////////

             if (vadv(i+ioff,j+1,k+koff) > 0.0) {
                jsign = 1.0;
                joff = 0;
             } else {
                jsign = -1.0;
                joff = 1;
             }

             ww = 0.0;
             if (wadv(i,j,k)*wadv(i+ioff,j+joff,k) > 0.0) {
                ww = wadv(i+ioff,j+joff,k);
             }

             uu = 0.0;
             if (uadv(i+1,j,k+koff)*uadv(i+1,j+joff,k+koff) > 0.0) {
                uu = uadv(i+1,j+joff,k+koff);
             }

             p1(1) = isign*0.5*hx;
             p1(2) = jsign*0.5*hy;
             p1(3) = ksign*0.5*hz;

             p2(1) = isign*0.5*hx;
             p2(2) = jsign*0.5*hy;
             p2(3) = ksign*0.5*hz - wadv(i,j,k)*dt;

             p3(1) = isign*0.5*hx - uadv(i+1,j+joff,k)*dt;
             p3(2) = jsign*0.5*hy;
             p3(3) = ksign*0.5*hz - wadv(i,j,k)*dt;

             p4(1) = isign*0.5*hx - uu*dt;
             p4(2) = jsign*0.5*hy - vadv(i+ioff,j+1,k+koff)*dt;
             p4(3) = ksign*0.5*hz - ww*dt;

             for(int n=1; n<=7; ++n){
                 slope_tmp(n) = slope(i+ioff,j+joff,k+koff,n-1);
             }

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
             }
             val1 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
             }
             val2 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
             }
             val3 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
             }
             val4 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
             }
             val5 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

             // source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff)
                                      +gamma2*vy(i+ioff,j+joff,k+koff)
                                      +gamma2*wz(i+ioff,j+joff,k+koff));

             gamma2 = gamma2 * vadv(i+ioff,j+1,k+koff);

             gamma = gamma - dt*gamma2/(3.0*hy);

             ////////////////////////////////////////////////
             // correct \Gamma^{x+} with \Gamma^{x+,y-}
             ////////////////////////////////////////////////

             if (vadv(i+ioff,j,k+koff) > 0.0) {
                jsign = 1.0;
                joff = -1;
             } else {
                jsign = -1.0;
                joff = 0;
             }

             ww = 0.0;
             if (wadv(i,j,k)*wadv(i+ioff,j+joff,k) > 0.0) {
                ww = wadv(i+ioff,j+joff,k);
             }

             uu = 0.0;
             if (uadv(i+1,j,k+koff)*uadv(i+1,j+joff,k+koff) > 0) {
                uu = uadv(i+1,j+joff,k+koff);
             }

             p1(1) = isign*0.5*hx;
             p1(2) = jsign*0.5*hy;
             p1(3) = ksign*0.5*hz;

             p2(1) = isign*0.5*hx;
             p2(2) = jsign*0.5*hy;
             p2(3) = ksign*0.5*hz - wadv(i,j,k)*dt;

             p3(1) = isign*0.5*hx - uadv(i+1,j+joff,k)*dt;
             p3(2) = jsign*0.5*hy;
             p3(3) = ksign*0.5*hz - wadv(i,j,k)*dt;

             p4(1) = isign*0.5*hx - uu*dt;
             p4(2) = jsign*0.5*hy - vadv(i+ioff,j,k+koff)*dt;
             p4(3) = ksign*0.5*hz - ww*dt;

             for(int n=1; n<=7; ++n){
                 slope_tmp(n) = slope(i+ioff,j+joff,k+koff,n-1);
             }

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
             }
             val1 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
             }
             val2 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
             }
             val3 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
             }
             val4 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
             }
             val5 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

             // source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff)
                                      +gamma2*vy(i+ioff,j+joff,k+koff)
                                      +gamma2*wz(i+ioff,j+joff,k+koff));

             gamma2 = gamma2 * vadv(i+ioff,j,k+koff);

             gamma = gamma + dt*gamma2/(3.0*hy);

             ////////////////////////////////////////////////
             // correct sedgez with \Gamma^{x+}
             ////////////////////////////////////////////////

             gamma = gamma * uadv(i+1,j,k+koff);
             sedgez(i,j,k) = sedgez(i,j,k) - dt*gamma/(2.0*hx);

             ////////////////////////////////////////////////
             // compute \Gamma^{x-} without corner corrections
             ////////////////////////////////////////////////

             if (uadv(i,j,k+koff) > 0.0) {
                isign = 1.0;
                ioff = -1;
             } else {
                isign = -1.0;
                ioff = 0;
             }

             w = 0.0;
             if (wadv(i,j,k)*wadv(i+ioff,j,k) > 0.0) {
                w = wadv(i+ioff,j,k);
             }

             p1(1) = isign*0.5*hx;
             p1(2) = 0.0;
             p1(3) = ksign*0.5*hz;

             p2(1) = isign*0.5*hx;
             p2(2) = 0.0;
             p2(3) = ksign*0.5*hz - wadv(i,j,k)*dt;

             p3(1) = isign*0.5*hx - uadv(i,j,k+koff)*dt;
             p3(2) = 0.0;
             p3(3) = ksign*0.5*hz - w*dt;

             for(int n=1; n<=7; ++n){
                 slope_tmp(n) = slope(i+ioff,j,k+koff,n-1);
             }

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = (p2(ll)+p3(ll))/2.0;
             }
             val1 = eval(s(i+ioff,j,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = (p1(ll)+p3(ll))/2.0;
             }
             val2 = eval(s(i+ioff,j,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = (p1(ll)+p2(ll))/2.0;
             }
             val3 = eval(s(i+ioff,j,k+koff),slope_tmp,del);

             // average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.0;

             // source term
             gamma = gamma - dt3*(gamma*wz(i+ioff,j,k+koff) + gamma*ux(i+ioff,j,k+koff));

             ////////////////////////////////////////////////
             // correct \Gamma^{x-} with \Gamma^{x-,y+}
             ////////////////////////////////////////////////

             if (vadv(i+ioff,j+1,k+koff) > 0.0) {
                jsign = 1.0;
                joff = 0;
             } else {
                jsign = -1.0;
                joff = 1;
             }

             ww = 0.0;
             if (wadv(i,j,k)*wadv(i+ioff,j+joff,k) > 0.0) {
                ww = wadv(i+ioff,j+joff,k);
             }

             uu = 0.0;
             if (uadv(i,j,k+koff)*uadv(i,j+joff,k+koff) > 0.0) {
                uu = uadv(i,j+joff,k+koff);
             }

             p1(1) = isign*0.5*hx;
             p1(2) = jsign*0.5*hy;
             p1(3) = ksign*0.5*hz;

             p2(1) = isign*0.5*hx;
             p2(2) = jsign*0.5*hy;
             p2(3) = ksign*0.5*hz - wadv(i,j,k)*dt;

             p3(1) = isign*0.5*hx - uadv(i,j+joff,k)*dt;
             p3(2) = jsign*0.5*hy;
             p3(3) = ksign*0.5*hz - wadv(i,j,k)*dt;

             p4(1) = isign*0.5*hx - uu*dt;
             p4(2) = jsign*0.5*hy - vadv(i+ioff,j+1,k+koff)*dt;
             p4(3) = ksign*0.5*hz - ww*dt;

             for(int n=1; n<=7; ++n){
                 slope_tmp(n) = slope(i+ioff,j+joff,k+koff,n-1);
             }

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
             }
             val1 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
             }
             val2 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
             }
             val3 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
             }
             val4 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
             }
             val5 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

             // source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff)
                                      +gamma2*vy(i+ioff,j+joff,k+koff)
                                      +gamma2*wz(i+ioff,j+joff,k+koff));

             gamma2 = gamma2 * vadv(i+ioff,j+1,k+koff);

             gamma = gamma - dt*gamma2/(3.0*hy);

             ////////////////////////////////////////////////
             // correct \Gamma^{x-} with \Gamma^{x-,y-}
             ////////////////////////////////////////////////

             if (vadv(i+ioff,j,k+koff) > 0.0) {
                jsign = 1.0;
                joff = -1;
             } else {
                jsign = -1.0;
                joff = 0;
             }

             ww = 0.0;
             if (wadv(i,j,k)*wadv(i+ioff,j+joff,k) > 0.0) {
                ww = wadv(i+ioff,j+joff,k);
             }

             uu = 0.0;
             if (uadv(i,j,k+koff)*uadv(i,j+joff,k+koff) > 0.0) {
                uu = uadv(i,j+joff,k+koff);
             }

             p1(1) = isign*0.5*hx;
             p1(2) = jsign*0.5*hy;
             p1(3) = ksign*0.5*hz;

             p2(1) = isign*0.5*hx;
             p2(2) = jsign*0.5*hy;
             p2(3) = ksign*0.5*hz - wadv(i,j,k)*dt;

             p3(1) = isign*0.5*hx - uadv(i,j+joff,k)*dt;
             p3(2) = jsign*0.5*hy;
             p3(3) = ksign*0.5*hz - wadv(i,j,k)*dt;

             p4(1) = isign*0.5*hx - uu*dt;
             p4(2) = jsign*0.5*hy - vadv(i+ioff,j,k+koff)*dt;
             p4(3) = ksign*0.5*hz - ww*dt;

             for(int n=1; n<=7; ++n){
                 slope_tmp(n) = slope(i+ioff,j+joff,k+koff,n-1);
             }

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
             }
             val1 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
             }
             val2 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
             }
             val3 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
             }
             val4 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
             }
             val5 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

             // source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff)
                                      +gamma2*vy(i+ioff,j+joff,k+koff)
                                      +gamma2*wz(i+ioff,j+joff,k+koff));

             gamma2 = gamma2 * vadv(i+ioff,j,k+koff);

             gamma = gamma + dt*gamma2/(3.0*hy);

             ////////////////////////////////////////////////
             // correct sedgez with \Gamma^{x-}
             ////////////////////////////////////////////////

             gamma = gamma * uadv(i,j,k+koff);
             sedgez(i,j,k) = sedgez(i,j,k) + dt*gamma/(2.0*hx);

             ////////////////////////////////////////////////
             // compute \Gamma^{y+} without corner corrections
             ////////////////////////////////////////////////

             if (vadv(i,j+1,k+koff) > 0.0) {
                jsign = 1.0;
                joff = 0;
             } else {
                jsign = -1.0;
                joff = 1;
             }

             w = 0.0;
             if (wadv(i,j,k)*wadv(i,j+joff,k) > 0.0) {
                w = wadv(i,j+joff,k);
             }

             p1(1) = 0.0;
             p1(2) = jsign*0.5*hy;
             p1(3) = ksign*0.5*hz;

             p2(1) = 0.0;
             p2(2) = jsign*0.5*hy;
             p2(3) = ksign*0.5*hz - wadv(i,j,k)*dt;

             p3(1) = 0.0;
             p3(2) = jsign*0.5*hy - vadv(i,j+1,k+koff)*dt;
             p3(3) = ksign*0.5*hz - w*dt;

             for(int n=1; n<=7; ++n){
                 slope_tmp(n) = slope(i,j+joff,k+koff,n-1);
             }

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = (p2(ll)+p3(ll))/2.0;
             }
             val1 = eval(s(i,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = (p1(ll)+p3(ll))/2.0;
             }
             val2 = eval(s(i,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = (p1(ll)+p2(ll))/2.0;
             }
             val3 = eval(s(i,j+joff,k+koff),slope_tmp,del);

             // average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.0;

             // source term
             gamma = gamma - dt3*(gamma*wz(i,j+joff,k+koff) + gamma*vy(i,j+joff,k+koff));

             ////////////////////////////////////////////////
             // correct \Gamma^{y+} with \Gamma^{y+,x+}
             ////////////////////////////////////////////////

             if (uadv(i+1,j+joff,k+koff) > 0.0) {
                isign = 1.0;
                ioff = 0;
             } else {
                isign = -1.0;
                ioff = 1;
             }

             ww = 0.0;
             if (wadv(i,j,k)*wadv(i+ioff,j+joff,k) > 0.0) {
                ww = wadv(i+ioff,j+joff,k);
             }

             vv = 0.0;
             if (vadv(i,j+1,k+koff)*vadv(i+ioff,j+1,k+koff) > 0.0) {
                vv = vadv(i+ioff,j+1,k+koff);
             }

             p1(1) = isign*0.5*hx;
             p1(2) = jsign*0.5*hy;
             p1(3) = ksign*0.5*hz;

             p2(1) = isign*0.5*hx;
             p2(2) = jsign*0.5*hy;
             p2(3) = ksign*0.5*hz - wadv(i,j,k)*dt;

             p3(1) = isign*0.5*hx;
             p3(2) = jsign*0.5*hy - vadv(i+ioff,j+1,k)*dt;
             p3(3) = ksign*0.5*hz - wadv(i,j,k)*dt;

             p4(1) = isign*0.5*hx - uadv(i+1,j+joff,k+koff)*dt;
             p4(2) = jsign*0.5*hy - vv*dt;
             p4(3) = ksign*0.5*hz - ww*dt;

             for(int n=1; n<=7; ++n){
                 slope_tmp(n) = slope(i+ioff,j+joff,k+koff,n-1);
             }

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
             }
             val1 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
             }
             val2 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
             }
             val3 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
             }
             val4 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
             }
             val5 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

             // source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff)
                                      +gamma2*vy(i+ioff,j+joff,k+koff)
                                      +gamma2*wz(i+ioff,j+joff,k+koff));

             gamma2 = gamma2 * uadv(i+1,j+joff,k+koff);

             gamma = gamma - dt*gamma2/(3.0*hx);

             ////////////////////////////////////////////////
             // correct \Gamma^{y+} with \Gamma^{y+,x-}
             ////////////////////////////////////////////////

             if (uadv(i,j+joff,k+koff) > 0.0) {
                isign = 1.0;
                ioff = -1;
             } else {
                isign = -1.0;
                ioff = 0;
             }

             ww = 0.0;
             if (wadv(i,j,k)*wadv(i+ioff,j+joff,k) > 0.0) {
                ww = wadv(i+ioff,j+joff,k);
             }

             vv = 0.0;
             if (vadv(i,j+1,k+koff)*vadv(i+ioff,j+1,k+koff) > 0.0) {
                vv = vadv(i+ioff,j+1,k+koff);
             }

             p1(1) = isign*0.5*hx;
             p1(2) = jsign*0.5*hy;
             p1(3) = ksign*0.5*hz;

             p2(1) = isign*0.5*hx;
             p2(2) = jsign*0.5*hy;
             p2(3) = ksign*0.5*hz - wadv(i,j,k)*dt;

             p3(1) = isign*0.5*hx;
             p3(2) = jsign*0.5*hy - vadv(i+ioff,j+1,k)*dt;
             p3(3) = ksign*0.5*hz - wadv(i,j,k)*dt;

             p4(1) = isign*0.5*hx - uadv(i,j+joff,k+koff)*dt;
             p4(2) = jsign*0.5*hy - vv*dt;
             p4(3) = ksign*0.5*hz - ww*dt;

             for(int n=1; n<=7; ++n){
                 slope_tmp(n) = slope(i+ioff,j+joff,k+koff,n-1);
             }

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
             }
             val1 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
             }
             val2 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
             }
             val3 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
             }
             val4 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
             }
             val5 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

             // source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff)
                                      +gamma2*vy(i+ioff,j+joff,k+koff)
                                      +gamma2*wz(i+ioff,j+joff,k+koff));

             gamma2 = gamma2 * uadv(i,j+joff,k+koff);

             gamma = gamma + dt*gamma2/(3.0*hx);

             ////////////////////////////////////////////////
             // correct sedgez with \Gamma^{y+}
             ////////////////////////////////////////////////

             gamma = gamma * vadv(i,j+1,k+koff);
             sedgez(i,j,k) = sedgez(i,j,k) - dt*gamma/(2.0*hy);

             ////////////////////////////////////////////////
             // compute \Gamma^{y-} without corner corrections
             ////////////////////////////////////////////////

             if (vadv(i,j,k+koff) > 0.0) {
                jsign = 1.0;
                joff = -1;
             } else {
                jsign = -1.0;
                joff = 0;
             }

             w = 0.0;
             if (wadv(i,j,k)*wadv(i,j+joff,k) > 0.0) {
                w = wadv(i,j+joff,k);
             }

             p1(1) = 0.0;
             p1(2) = jsign*0.5*hy;
             p1(3) = ksign*0.5*hz;

             p2(1) = 0.0;
             p2(2) = jsign*0.5*hy;
             p2(3) = ksign*0.5*hz - wadv(i,j,k)*dt;

             p3(1) = 0.0;
             p3(2) = jsign*0.5*hy - vadv(i,j,k+koff)*dt;
             p3(3) = ksign*0.5*hz - w*dt;

             for(int n=1; n<=7; ++n){
                 slope_tmp(n) = slope(i,j+joff,k+koff,n-1);
             }

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = (p2(ll)+p3(ll))/2.0;
             }
             val1 = eval(s(i,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = (p1(ll)+p3(ll))/2.0;
             }
             val2 = eval(s(i,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = (p1(ll)+p2(ll))/2.0;
             }
             val3 = eval(s(i,j+joff,k+koff),slope_tmp,del);

             // average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.0;

             // source term
             gamma = gamma - dt3*(gamma*wz(i,j+joff,k+koff) + gamma*vy(i,j+joff,k+koff));

             ////////////////////////////////////////////////
             // correct \Gamma^{y-} with \Gamma^{y-,x+};
             ////////////////////////////////////////////////

             if (uadv(i+1,j+joff,k+koff) > 0.0) {
                isign = 1.0;
                ioff = 0;
             } else {
                isign = -1.0;
                ioff = 1;
             }

             ww = 0.0;
             if (wadv(i,j,k)*wadv(i+ioff,j+joff,k) > 0.0) {
                ww = wadv(i+ioff,j+joff,k);
             }

             vv = 0.0;
             if (vadv(i,j,k+koff)*vadv(i+ioff,j,k+koff) > 0.0) {
                vv = vadv(i+ioff,j,k+koff);
             }

             p1(1) = isign*0.5*hx;
             p1(2) = jsign*0.5*hy;
             p1(3) = ksign*0.5*hz;

             p2(1) = isign*0.5*hx;
             p2(2) = jsign*0.5*hy;
             p2(3) = ksign*0.5*hz - wadv(i,j,k)*dt;

             p3(1) = isign*0.5*hx;
             p3(2) = jsign*0.5*hy - vadv(i+ioff,j,k)*dt;
             p3(3) = ksign*0.5*hz - wadv(i,j,k)*dt;

             p4(1) = isign*0.5*hx - uadv(i+1,j+joff,k+koff)*dt;
             p4(2) = jsign*0.5*hy - vv*dt;
             p4(3) = ksign*0.5*hz - ww*dt;

             for(int n=1; n<=7; ++n){
                 slope_tmp(n) = slope(i+ioff,j+joff,k+koff,n-1);
             }

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
             }
             val1 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
             }
             val2 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
             }
             val3 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
             }
             val4 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
             }
             val5 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

             // source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff)
                                      +gamma2*vy(i+ioff,j+joff,k+koff)
                                      +gamma2*wz(i+ioff,j+joff,k+koff));

             gamma2 = gamma2 * uadv(i+1,j+joff,k+koff);

             gamma = gamma - dt*gamma2/(3.0*hx);

             ////////////////////////////////////////////////
             // correct \Gamma^{y-} with \Gamma^{y-,x-}
             ////////////////////////////////////////////////

             if (uadv(i,j+joff,k+koff) > 0.0) {
                isign = 1.0;
                ioff = -1;
             } else {
                isign = -1.0;
                ioff = 0;
             }

             ww = 0.0;
             if (wadv(i,j,k)*wadv(i+ioff,j+joff,k) > 0.0) {
                ww = wadv(i+ioff,j+joff,k);
             }

             vv = 0.0;
             if (vadv(i,j,k+koff)*vadv(i+ioff,j,k+koff) > 0.0) {
                vv = vadv(i+ioff,j,k+koff);
             }

             p1(1) = isign*0.5*hx;
             p1(2) = jsign*0.5*hy;
             p1(3) = ksign*0.5*hz;

             p2(1) = isign*0.5*hx;
             p2(2) = jsign*0.5*hy;
             p2(3) = ksign*0.5*hz - wadv(i,j,k)*dt;

             p3(1) = isign*0.5*hx;
             p3(2) = jsign*0.5*hy - vadv(i+ioff,j,k)*dt;
             p3(3) = ksign*0.5*hz - wadv(i,j,k)*dt;

             p4(1) = isign*0.5*hx - uadv(i,j+joff,k+koff)*dt;
             p4(2) = jsign*0.5*hy - vv*dt;
             p4(3) = ksign*0.5*hz - ww*dt;

             for(int n=1; n<=7; ++n){
                 slope_tmp(n) = slope(i+ioff,j+joff,k+koff,n-1);
             }

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
             }
             val1 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
             }
             val2 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
             }
             val3 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
             }
             val4 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             for(int ll=1; ll<=3; ++ll ){
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
             }
             val5 = eval(s(i+ioff,j+joff,k+koff),slope_tmp,del);

             gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

             // source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff)
                                      +gamma2*vy(i+ioff,j+joff,k+koff)
                                      +gamma2*wz(i+ioff,j+joff,k+koff));

             gamma2 = gamma2 * uadv(i,j+joff,k+koff);

             gamma = gamma + dt*gamma2/(3.0*hx);

             ////////////////////////////////////////////////
             // correct sedgez with \Gamma^{y-}
             ////////////////////////////////////////////////

             gamma = gamma * vadv(i,j,k+koff);
             sedgez(i,j,k) = sedgez(i,j,k) + dt*gamma/(2.0*hy);

          });
    }
}

/** @} */

