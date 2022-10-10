/**
 * \file hydro_bds_edge_state_2D.cpp
 *
 * \addtogroup BDS
 *  @{
 */

#include <hydro_bds.H>
#include <hydro_constants.H>

using namespace amrex;

constexpr amrex::Real eps = 1.0e-8;

/**
 * Uses the Bell-Dawson-Shubin (BDS) algorithm, a higher order Godunov
 * method for scalar conservation laws in two dimensions, to compute
 * edge states.
 *
 * \param [in]     bx          Current grid patch
 * \param [in]     ncomp       Number of components to work on
 * \param [in]     q           Array4 of state, starting at component of interest
 * \param [in,out] xedge       Array4 containing x-edges, starting at component of interest
 * \param [in,out] yedge       Array4 containing y-edges, starting at component of interest
 * \param [in]     umac        x-Face velocities.
 * \param [in]     vmac        y-Face velocities.
 * \param [in]     fq          Array4 for forces, starting at component of interest
 * \param [in]     geom        Level geometry.
 * \param [in]     l_dt        Time step.
 * \param [in]     iconserv    Indicates conservative dimensions.
 * \param [in]     is_velocity Indicates a component is velocity so boundary conditions can
 *                             be properly addressed. The header hydro_constants.H
 *                             defines the component positon by [XY]VEL macro.
 *
 */

void
BDS::ComputeEdgeState ( Box const& bx, int ncomp,
                        Array4<Real const> const& q,
                        Array4<Real      > const& xedge,
                        Array4<Real      > const& yedge,
                        Array4<Real const> const& umac,
                        Array4<Real const> const& vmac,
                        Array4<Real const> const& divu,
                        Array4<Real const> const& fq,
                        Geometry geom,
                        Real l_dt,
                        BCRec const* pbc, int const* iconserv,
                        const bool is_velocity)
{
    // For now, loop on components here
    for( int icomp = 0; icomp < ncomp; ++icomp)
    {
        // Temporary slope container per components
        Box const& bxg1 = amrex::grow(bx,1);
        FArrayBox slopefab(bxg1,3,The_Async_Arena());

        BDS::ComputeSlopes(bx, geom, icomp,
                           q, slopefab.array(),
                           pbc);

        BDS::ComputeConc(bx, geom, icomp,
                         q, xedge, yedge, slopefab.array(),
                         umac, vmac, divu, fq,
                         iconserv,
                         l_dt, pbc, is_velocity);
    }
}

/**
 * Compute bilinear slopes for BDS algorithm.
 *
 * \param [in]  bx      Current grid patch
 * \param [in]  geom    Level geometry.
 * \param [in]  icomp   Component of the state Array4.
 * \param [in]  s       Array4<const> of state vector.
 * \param [out] slopes  Array4 to store slope information.
 *
 */

void
BDS::ComputeSlopes ( Box const& bx,
                     const Geometry& geom,
                     int icomp,
                     Array4<Real const> const& s,
                     Array4<Real      > const& slopes,
                     BCRec const* pbc)
{
    constexpr bool limit_slopes = true;

    // Define container for the nodal interpolated state
    Box const& ngbx = amrex::grow(amrex::convert(bx,IntVect(AMREX_D_DECL(1,1,1))),1);
    FArrayBox tmpnodefab(ngbx,1,The_Async_Arena());
    auto const& sint = tmpnodefab.array();

    Box const& gbx = amrex::grow(bx,1);
    GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    Real hx = dx[0];
    Real hy = dx[1];

    Box const& domain = geom.Domain();
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    auto bc = pbc[icomp];

    // Abort for cell-centered BC types
    if ( bc.lo(0) == BCType::reflect_even || bc.lo(0) == BCType::reflect_odd || bc.lo(0) == BCType::hoextrapcc ||
         bc.hi(0) == BCType::reflect_even || bc.hi(0) == BCType::reflect_odd || bc.hi(0) == BCType::hoextrapcc ||
         bc.lo(1) == BCType::reflect_even || bc.lo(1) == BCType::reflect_odd || bc.lo(1) == BCType::hoextrapcc ||
         bc.hi(1) == BCType::reflect_even || bc.hi(1) == BCType::reflect_odd || bc.hi(1) == BCType::hoextrapcc )
        amrex::Abort("BDS::Slopes: Unsupported BC type. Supported types are int_dir, ext_dir, foextrap, and hoextrap");

    bool lo_x_physbc = (bc.lo(0) == BCType::foextrap || bc.lo(0) == BCType::hoextrap || bc.lo(0) == BCType::ext_dir) ? true : false;
    bool hi_x_physbc = (bc.hi(0) == BCType::foextrap || bc.hi(0) == BCType::hoextrap || bc.hi(0) == BCType::ext_dir) ? true : false;
    bool lo_y_physbc = (bc.lo(1) == BCType::foextrap || bc.lo(1) == BCType::hoextrap || bc.lo(1) == BCType::ext_dir) ? true : false;
    bool hi_y_physbc = (bc.hi(1) == BCType::foextrap || bc.hi(1) == BCType::hoextrap || bc.hi(1) == BCType::ext_dir) ? true : false;

    // bicubic interpolation to corner points
    // (i,j,k) refers to lower corner of cell
    // Added k index -- placeholder for 2d
    ParallelFor(ngbx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

        // set node values equal to the average of the ghost cell values since they store the physical condition on the boundary
        if ( i<=dlo.x && lo_x_physbc ) {
            sint(i,j,k) = 0.5*(s(dlo.x-1,j,k,icomp) + s(dlo.x-1,j-1,k,icomp));
            return;
        }
        if ( i>=dhi.x+1 && hi_x_physbc ) {
            sint(i,j,k) = 0.5*(s(dhi.x+1,j,k,icomp) + s(dhi.x+1,j-1,k,icomp));
            return;
        }
        if ( j<=dlo.y && lo_y_physbc ) {
            sint(i,j,k) = 0.5*(s(i,dlo.y-1,k,icomp) + s(i-1,dlo.y-1,k,icomp));
            return;
        }
        if ( j>=dhi.y+1 && hi_y_physbc ) {
            sint(i,j,k) = 0.5*(s(i,dhi.y+1,k,icomp) + s(i-1,dhi.y+1,k,icomp));
            return;
        }

        // one cell inward from any physical boundary, revert to 4-point average
        if ( (i==dlo.x+1 && lo_x_physbc) ||
             (i==dhi.x   && hi_x_physbc) ||
             (j==dlo.y+1 && lo_y_physbc) ||
             (j==dhi.y   && hi_y_physbc) ) {

            sint(i,j,k) = 0.25* (s(i,j,k,icomp) + s(i-1,j,k,icomp) + s(i,j-1,k,icomp) + s(i-1,j-1,k,icomp));
            return;
        }

        sint(i,j,k) = (s(i-2,j-2,k,icomp) + s(i-2,j+1,k,icomp) + s(i+1,j-2,k,icomp) + s(i+1,j+1,k,icomp)
                - 7.0*(s(i-2,j-1,k,icomp) + s(i-2,j  ,k,icomp) + s(i-1,j-2,k,icomp) + s(i  ,j-2,k,icomp) +
                       s(i-1,j+1,k,icomp) + s(i  ,j+1,k,icomp) + s(i+1,j-1,k,icomp) + s(i+1,j  ,k,icomp))
               + 49.0*(s(i-1,j-1,k,icomp) + s(i  ,j-1,k,icomp) + s(i-1,j  ,k,icomp) + s(i  ,j  ,k,icomp)) ) / 144.0;
    });

    ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k){
        // compute initial estimates of slopes from unlimited corner points

        // local variables
        Real sumloc, redfac, redmax, div, kdp, sumdif, sgndif;

        Array1D<Real, 1, 4> diff;
        Array1D<Real, 1, 4> smin;
        Array1D<Real, 1, 4> smax;
        Array1D<Real, 1, 4> sc;

        Array1D<bool, 1, 4> allow_change;
        for (int mm=1; mm<=4; ++mm) {
            allow_change(mm) = true;
        }

        if ( i<=dlo.x && lo_x_physbc ) {
            allow_change(1) = false;
            allow_change(2) = false;
        }
        if ( i>=dhi.x+1 && hi_x_physbc ) {
            allow_change(3) = false;
            allow_change(4) = false;
        }
        if ( j<=dlo.y && lo_y_physbc ) {
            allow_change(1) = false;
            allow_change(3) = false;
        }
        if ( j>=dhi.y+1 && hi_y_physbc ) {
            allow_change(2) = false;
            allow_change(4) = false;
        }

        // compute initial estimates of slopes from unlimited corner points
        // sx
        slopes(i,j,k,0) = 0.5*(sint(i+1,j+1,k) + sint(i+1,j,k) - sint(i,j+1,k) - sint(i,j,k)) / hx;
        // sy
        slopes(i,j,k,1) = 0.5*(sint(i+1,j+1,k) - sint(i+1,j,k) + sint(i,j+1,k) - sint(i,j,k)) / hy;
        // sxy
        slopes(i,j,k,2) =     (sint(i+1,j+1,k) - sint(i+1,j,k) - sint(i,j+1,k) + sint(i,j,k)) / (hx*hy);

        if (limit_slopes) {

            // ++ / sint(i+1,j+1)
            sc(4) = s(i,j,k,icomp) + 0.5*(hx*slopes(i,j,k,0) + hy*slopes(i,j,k,1)) + 0.25*hx*hy*slopes(i,j,k,2);

            // +- / sint(i+1,j  )
            sc(3) = s(i,j,k,icomp) + 0.5*(hx*slopes(i,j,k,0) - hy*slopes(i,j,k,1)) - 0.25*hx*hy*slopes(i,j,k,2);

            // -+ / sint(i  ,j+1)
            sc(2) = s(i,j,k,icomp) - 0.5*(hx*slopes(i,j,k,0) - hy*slopes(i,j,k,1)) - 0.25*hx*hy*slopes(i,j,k,2);

            // -- / sint(i  ,j  )
            sc(1) = s(i,j,k,icomp) - 0.5*(hx*slopes(i,j,k,0) + hy*slopes(i,j,k,1)) + 0.25*hx*hy*slopes(i,j,k,2);

            // enforce max/min bounds
            smin(4) = amrex::min(s(i,j,k,icomp), s(i+1,j,k,icomp), s(i,j+1,k,icomp), s(i+1,j+1,k,icomp));
            smax(4) = amrex::max(s(i,j,k,icomp), s(i+1,j,k,icomp), s(i,j+1,k,icomp), s(i+1,j+1,k,icomp));

            smin(3) = amrex::min(s(i,j,k,icomp), s(i+1,j,k,icomp), s(i,j-1,k,icomp), s(i+1,j-1,k,icomp));
            smax(3) = amrex::max(s(i,j,k,icomp), s(i+1,j,k,icomp), s(i,j-1,k,icomp), s(i+1,j-1,k,icomp));

            smin(2) = amrex::min(s(i,j,k,icomp), s(i-1,j,k,icomp), s(i,j+1,k,icomp), s(i-1,j+1,k,icomp));
            smax(2) = amrex::max(s(i,j,k,icomp), s(i-1,j,k,icomp), s(i,j+1,k,icomp), s(i-1,j+1,k,icomp));

            smin(1) = amrex::min(s(i,j,k,icomp), s(i-1,j,k,icomp), s(i,j-1,k,icomp), s(i-1,j-1,k,icomp));
            smax(1) = amrex::max(s(i,j,k,icomp), s(i-1,j,k,icomp), s(i,j-1,k,icomp), s(i-1,j-1,k,icomp));

            for(int mm=1; mm<=4; ++mm){
               if (allow_change(mm)) {
                   sc(mm) = amrex::max(amrex::min(sc(mm), smax(mm)), smin(mm));
               }
            }

            // iterative loop
            for(int ll=1; ll<=3; ++ll){

               // compute the amount by which the average of the nodal values differs from cell-center value
               sumloc = 0.25*(sc(4) + sc(3) + sc(2) + sc(1));
               sumdif = (sumloc - s(i,j,k,icomp))*4.0;

               // sgndif = +(-)1 if the node average is too large(small)
               sgndif = std::copysign(1.0,sumdif);

               // compute how much each node is larger(smaller) than the cell-centered value
               for(int mm=1; mm<=4; ++mm){
                  diff(mm) = (sc(mm) - s(i,j,k,icomp))*sgndif;
               }

               kdp = 0;

               // count how many nodes are larger(smaller) than the cell-centered value
               for(int mm=1; mm<=4; ++mm){
                  if (diff(mm) > eps && allow_change(mm)) {
                     kdp = kdp+1;
                  }
               }

               // adjust node values
               for(int mm=1; mm<=4; ++mm){

                  // don't allow boundary nodes to change value
                  if (!allow_change(mm)) continue;

                  // how many node values are left to potentially adjust
                  if (kdp<1) {
                     div = 1.0;
                  } else {
                     div = kdp;
                  }

                  // if the node needs adjusting, figure out by how much the remaining sum is divy'ed up
                  if (diff(mm)>eps) {
                     redfac = sumdif*sgndif/div;
                     kdp = kdp-1;
                  } else {
                     redfac = 0.0;
                  }

                  // don't let the adjustment introduce any new extrema
                  if (sgndif > 0.0) {
                     redmax = sc(mm) - smin(mm);
                  } else {
                     redmax = smax(mm) - sc(mm);
                  }
                  redfac = amrex::min(redfac,redmax);

                  // adjust nodal value and decrement the excess
                  sumdif = sumdif - redfac*sgndif;
                  sc(mm) = sc(mm) - redfac*sgndif;
               }
            }

            // final slopes
            // sx
            slopes(i,j,k,0) = 0.5*( sc(4) + sc(3) -sc(1) - sc(2) )/hx;
            // sy
            slopes(i,j,k,1) = 0.5*( sc(4) + sc(2) -sc(1) - sc(3) )/hy;
            // sxy
            slopes(i,j,k,2) =     ( sc(1) + sc(4) -sc(2) - sc(3) )/(hx*hy);
        }
    });
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
 * \param [in]     bx          Current grid patch
 * \param [in]     geom        Level geometry.
 * \param [in]     icomp       Component of the Array4s.
 * \param [in]     s           Array4 of state.
 * \param [in,out] sedgex      Array4 containing x-edges.
 * \param [in,out] sedgey      Array4 containing y-edges.
 * \param [in]     slopes      Array4 containing slope information.
 * \param [in]     umac        Array4 for u-face velocity.
 * \param [in]     vmac        Array4 for v-face velocity.
 * \param [in]     force       Array4 for forces.
 * \param [in]     iconserv    Indicates conservative dimensions.
 * \param [in]     dt          Time step.
 * \param [in]     is_velocity Indicates a component is velocity so boundary conditions can
 *                             be properly addressed. The header hydro_constants.H
 *                             defines the component positon by [XY]VEL macro.
 *
 */

void
BDS::ComputeConc (Box const& bx,
                  const Geometry& geom,
                  int icomp,
                  Array4<Real const> const& s,
                  Array4<Real      > const& sedgex,
                  Array4<Real      > const& sedgey,
                  Array4<Real const> const& slopes,
                  Array4<Real const> const& umac,
                  Array4<Real const> const& vmac,
                  Array4<Real const> const& divu,
                  Array4<Real const> const& force,
                  int const* iconserv,
                  const Real dt, BCRec const* pbc,
                  const bool is_velocity)
{
    Box const& gbx = amrex::grow(bx,1);
    GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    FArrayBox ux_fab(gbx,1,The_Async_Arena());
    FArrayBox vy_fab(gbx,1,The_Async_Arena());
    auto const& ux    = ux_fab.array();
    auto const& vy    = vy_fab.array();

    Real hx = dx[0];
    Real hy = dx[1];

    Real dt2 = dt/2.0;
    Real dt3 = dt/3.0;

    Box const& domain = geom.Domain();
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    auto bc = pbc[icomp];
    bool lo_x_physbc = (bc.lo(0) == BCType::foextrap || bc.lo(0) == BCType::hoextrap || bc.lo(0) == BCType::ext_dir) ? true : false;
    bool hi_x_physbc = (bc.hi(0) == BCType::foextrap || bc.hi(0) == BCType::hoextrap || bc.hi(0) == BCType::ext_dir) ? true : false;
    bool lo_y_physbc = (bc.lo(1) == BCType::foextrap || bc.lo(1) == BCType::hoextrap || bc.lo(1) == BCType::ext_dir) ? true : false;
    bool hi_y_physbc = (bc.hi(1) == BCType::foextrap || bc.hi(1) == BCType::hoextrap || bc.hi(1) == BCType::ext_dir) ? true : false;

    // compute cell-centered ux, vy
    ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k){
        ux(i,j,k) = (umac(i+1,j,k) - umac(i,j,k)) / hx;
        vy(i,j,k) = (vmac(i,j+1,k) - vmac(i,j,k)) / hy;
    });

    // compute sedgex on x-faces
    Box const& xbx = amrex::surroundingNodes(bx,0);
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

        // set edge values equal to the ghost cell value since they store the physical condition on the boundary
        if ( i==dlo.x && lo_x_physbc ) {
            sedgex(i,j,k,icomp) = s(i-1,j,k,icomp);
            if (is_velocity && icomp == XVEL && (bc.lo(0) == BCType::foextrap ||  bc.lo(0) == BCType::hoextrap) ) {
                // make sure velocity is not blowing inward
                sedgex(i,j,k,icomp) = amrex::min(0._rt,sedgex(i,j,k,icomp));
            }
            // This will return from computing (i,j,k), not exit from ComputeConc.
            // Recall that ParallelFor is a macro that sets up a multidimensional loop
            // and then calls this lambda on the loop interior.
            return;
        }
        if ( i==dhi.x+1 && hi_x_physbc ) {
            sedgex(i,j,k,icomp) = s(i,j,k,icomp);
            if (is_velocity && icomp == XVEL && (bc.hi(0) == BCType::foextrap ||  bc.hi(0) == BCType::hoextrap) ) {
                // make sure velocity is not blowing inward
                sedgex(i,j,k,icomp) = amrex::max(0._rt,sedgex(i,j,k,icomp));
            }
            return;
        }

        // local variables
        Array1D<Real, 1, 2> del;
        Array1D<Real, 1, 2> p1;
        Array1D<Real, 1, 2> p2;
        Array1D<Real, 1, 2> p3;

        Array1D<Real, 1, 3> slope_tmp;

        int ioff, joff;
        Real isign,jsign;
        Real val1,val2,val3;
        Real u;
        Real gamma;

        // Since faces of xbx can overlap between tileboxes, need a temporary
        // to hold intermediate edgestate value
        Real xedge_tmp;

        ///////////////////////////////////////////////
        // compute sedgex without transverse corrections
        ///////////////////////////////////////////////

        if (umac(i,j,k) > 0) {
            isign = 1.;
            ioff = -1;
        } else {
            isign = -1.;
            ioff = 0;
        }

        for(int n=1; n<=3; ++n){
            slope_tmp(n) = slopes(i+ioff,j,k,n-1);
        }

        // centroid of rectangular volume
        del(1) = isign*0.5*hx - 0.5*umac(i,j,k)*dt;
        del(2) = 0.;
        xedge_tmp = eval(s(i+ioff,j,k,icomp),slope_tmp,del);

        // source term
        if (iconserv[icomp]) {
            xedge_tmp = xedge_tmp*(1. - dt2*ux(i+ioff,j,k));
        } else {
            xedge_tmp = xedge_tmp*(1. + dt2*vy(i+ioff,j,k));
        }
        if (force) {
            xedge_tmp += dt2*force(i+ioff,j,k,icomp);
        }

        ///////////////////////////////////////////////
        // compute \Gamma^{y+}
        ///////////////////////////////////////////////

        if (vmac(i+ioff,j+1,k) > 0.) {
            jsign = 1.;
            joff = 0;
        } else {
            jsign = -1.;
            joff = 1;
        }

        u = 0.;
        if (umac(i,j,k)*umac(i,j+joff,k) > 0.) {
            u = umac(i,j+joff,k);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;

        p2(1) = isign*0.5*hx - umac(i,j,k)*dt;
        p2(2) = jsign*0.5*hy;

        p3(1) = isign*0.5*hx - u*dt;
        p3(2) = jsign*0.5*hy - vmac(i+ioff,j+1,k)*dt;

        for(int n=1; n<=3; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k,n-1);
        }

        for (int ll=1; ll<=2; ++ll) {
            del(ll) = (p2(ll)+p3(ll))/2.;
        }
        val1 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        for (int ll=1; ll<=2; ++ll) {
            del(ll) = (p1(ll)+p3(ll))/2.;
        }
        val2 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        for (int ll=1; ll<=2; ++ll) {
            del(ll) = (p1(ll)+p2(ll))/2.;
        }
        val3 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        // average these centroid values to get the average value
        gamma = (val1+val2+val3)/3.;

        // source term
        if (iconserv[icomp]) {
            gamma = gamma*(1. - dt3*divu(i+ioff,j+joff,k));
        }

        ///////////////////////////////////////////////
        // correct sedgex with \Gamma^{y+}
        ///////////////////////////////////////////////

        gamma = gamma * vmac(i+ioff,j+1,k);

        xedge_tmp = xedge_tmp - dt*gamma/(2.*hy);

        ///////////////////////////////////////////////
        // compute \Gamma^{y-}
        ///////////////////////////////////////////////

        if (vmac(i+ioff,j,k) > 0.) {
            jsign = 1.;
            joff = -1;
        } else {
            jsign = -1.;
            joff = 0;
                }

        u = 0.;
        if (umac(i,j,k)*umac(i,j+joff,k) > 0.) {
            u = umac(i,j+joff,k);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;

        p2(1) = isign*0.5*hx - umac(i,j,k)*dt;
        p2(2) = jsign*0.5*hy;

        p3(1) = isign*0.5*hx - u*dt;
        p3(2) = jsign*0.5*hy - vmac(i+ioff,j,k)*dt;

        for(int n=1; n<=3; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k,n-1);
        }

        for (int ll=1; ll<=2; ++ll) {
            del(ll) = (p2(ll)+p3(ll))/2.;
        }
        val1 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        for (int ll=1; ll<=2; ++ll) {
            del(ll) = (p1(ll)+p3(ll))/2.;
        }
        val2 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        for (int ll=1; ll<=2; ++ll) {
            del(ll) = (p1(ll)+p2(ll))/2.;
        }
        val3 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        // average these centroid values to get the average value
        gamma = (val1+val2+val3)/3.;

        // source term
        if (iconserv[icomp]) {
            gamma = gamma*(1. - dt3*divu(i+ioff,j+joff,k));
        }

        ///////////////////////////////////////////////
        // correct sedgex with \Gamma^{y-}
        ///////////////////////////////////////////////

        gamma = gamma * vmac(i+ioff,j,k);
        sedgex(i,j,k,icomp) = xedge_tmp + dt*gamma/(2.*hy);
    });

    // compute sedgey on y-faces
    Box const& ybx = amrex::surroundingNodes(bx,1);
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

        // set edge values equal to the ghost cell value since they store the physical condition on the boundary
        if ( j==dlo.y && lo_y_physbc ) {
            sedgey(i,j,k,icomp) = s(i,j-1,k,icomp);
            if (is_velocity && icomp == YVEL && (bc.lo(1) == BCType::foextrap ||  bc.lo(1) == BCType::hoextrap) ) {
                // make sure velocity is not blowing inward
                sedgey(i,j,k,icomp) = amrex::min(0._rt,sedgey(i,j,k,icomp));
            }
            return;
        }
        if ( j==dhi.y+1 && hi_y_physbc ) {
            sedgey(i,j,k,icomp) = s(i,j,k,icomp);
            if (is_velocity && icomp == YVEL && (bc.hi(1) == BCType::foextrap ||  bc.hi(1) == BCType::hoextrap) ) {
                // make sure velocity is not blowing inward
                sedgey(i,j,k,icomp) = amrex::max(0._rt,sedgey(i,j,k,icomp));
            }
            return;
        }

        // local variables
        Array1D<Real, 1, 2> del;
        Array1D<Real, 1, 2> p1;
        Array1D<Real, 1, 2> p2;
        Array1D<Real, 1, 2> p3;

        Array1D<Real, 1, 3> slope_tmp;

        int ioff, joff;
        Real isign,jsign;
        Real val1,val2,val3;
        Real v;
        Real gamma;

        // To hold intermediate edgestate value
        Real yedge_tmp;

        ///////////////////////////////////////////////
        // compute sedgey without transverse corrections
        ///////////////////////////////////////////////

        // centroid of rectangular volume
        if (vmac(i,j,k) > 0.) {
            jsign = 1.;
            joff = -1;
        } else {
            jsign = -1.;
            joff = 0;
        }

        for(int n=1; n<=3; ++n){
            slope_tmp(n) = slopes(i,j+joff,k,n-1);
        }

        del(1) = 0.;
        del(2) = jsign*0.5*hy - 0.5*vmac(i,j,k)*dt;
        yedge_tmp = eval(s(i,j+joff,k,icomp),slope_tmp,del);

        // source term
        if (iconserv[icomp]) {
            yedge_tmp = yedge_tmp*(1. - dt2*vy(i,j+joff,k));
        } else {
            yedge_tmp = yedge_tmp*(1. + dt2*ux(i,j+joff,k));
        }
        if (force) {
            yedge_tmp += dt2*force(i,j+joff,k,icomp);
        }

        ///////////////////////////////////////////////
        // compute \Gamma^{x+}
        ///////////////////////////////////////////////

        if (umac(i+1,j+joff,k) > 0.) {
            isign = 1.;
            ioff = 0;
        } else {
            isign = -1.;
            ioff = 1;
        }

        v = 0.;
        if (vmac(i,j,k)*vmac(i+ioff,j,k) > 0.) {
            v = vmac(i+ioff,j,k);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;

        p2(1) = isign*0.5*hx;
        p2(2) = jsign*0.5*hy - vmac(i,j,k)*dt;

        p3(1) = isign*0.5*hx - umac(i+1,j+joff,k)*dt;
        p3(2) = jsign*0.5*hy - v*dt;

        for(int n=1; n<=3; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k,n-1);
        }

        for (int ll=1; ll<=2; ++ll) {
            del(ll) = (p2(ll)+p3(ll))/2.;
        }
        val1 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        for (int ll=1; ll<=2; ++ll) {
            del(ll) = (p1(ll)+p3(ll))/2.;
        }
        val2 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        for (int ll=1; ll<=2; ++ll) {
            del(ll) = (p1(ll)+p2(ll))/2.;
        }
        val3 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        // average these centroid values to get the average value
        gamma = (val1+val2+val3)/3.;

        // source term
        if (iconserv[icomp]) {
            gamma = gamma*(1. - dt3*divu(i+ioff,j+joff,k));
        }

        ///////////////////////////////////////////////
        // correct sedgey with \Gamma^{x+}
        ///////////////////////////////////////////////

        gamma = gamma * umac(i+1,j+joff,k);
        yedge_tmp = yedge_tmp - dt*gamma/(2.*hx);

        ///////////////////////////////////////////////
        // compute \Gamma^{x-}
        ///////////////////////////////////////////////

        if (umac(i,j+joff,k) > 0.) {
            isign = 1.;
            ioff = -1;
        } else {
            isign = -1.;
            ioff = 0;
        }

        v = 0.;
        if (vmac(i,j,k)*vmac(i+ioff,j,k) > 0.) {
            v = vmac(i+ioff,j,k);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;

        p2(1) = isign*0.5*hx;
        p2(2) = jsign*0.5*hy - vmac(i,j,k)*dt;

        p3(1) = isign*0.5*hx - umac(i,j+joff,k)*dt;
        p3(2) = jsign*0.5*hy - v*dt;

        for(int n=1; n<=3; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k,n-1);
        }

        for (int ll=1; ll<=2; ++ll) {
            del(ll) = (p2(ll)+p3(ll))/2.;
        }
        val1 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        for (int ll=1; ll<=2; ++ll) {
            del(ll) = (p1(ll)+p3(ll))/2.;
        }
        val2 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        for (int ll=1; ll<=2; ++ll) {
            del(ll) = (p1(ll)+p2(ll))/2.;
        }
        val3 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        // average these centroid values to get the average value
        gamma = (val1+val2+val3)/3.;

        // source term
        if (iconserv[icomp]) {
            gamma = gamma*(1. - dt3*divu(i+ioff,j+joff,k));
        }

        ///////////////////////////////////////////////
        // correct sedgey with \Gamma^{x-}
        ///////////////////////////////////////////////

        gamma = gamma * umac(i,j+joff,k);
        sedgey(i,j,k,icomp) = yedge_tmp + dt*gamma/(2.*hx);
    });
}

/** @} */
