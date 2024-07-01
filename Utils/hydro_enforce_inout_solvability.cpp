/** \addtogroup Utilities
 * @{
 */

#include <hydro_utils.H>

using namespace amrex;

namespace HydroUtils {

namespace {

void set_inout_masks(
    const int lev,
    const Vector<Array<MultiFab*, AMREX_SPACEDIM>>& vels_vec,
    Array<iMultiFab, AMREX_SPACEDIM>& inflow_masks,
    Array<iMultiFab, AMREX_SPACEDIM>& outflow_masks,
    const BCRec* bc_type,
    const Box& domain,
    const bool corners)
{
    // loop over the six orientations
    for (OrientationIter oit; oit != nullptr; ++oit) {
        const auto ori = oit();
        const auto side = ori.faceDir();
        const int dir = ori.coordDir();
        const auto islow = ori.isLow();
        const auto ishigh = ori.isHigh();

        // Multifab for normal mac velocity
        auto& vel_mf = vels_vec[lev][dir];
//Print() << vel_mf->boxArray() << std::endl;
        // mask iMFs for the respective velocity direction
        auto& inflow_mask = inflow_masks[dir];
        auto& outflow_mask = outflow_masks[dir];

        // domain extent indices for the velocities
        IndexType::CellIndex dir_index_type = (vel_mf->ixType()).ixType(dir);
        int dlo;
        if (dir_index_type == IndexType::CellIndex::CELL) {
            // lower boundary is at -1 for cell-centered velocity
            dlo = domain.smallEnd(dir) - 1;
        } else {
            // lower boundary is at  0 for face-centered velocity
            dlo = domain.smallEnd(dir);
        }
        int dhi = domain.bigEnd(dir) + 1;

        // get BCs for the normal velocity and set the boundary index
        // based on low or high side
        const BCRec ibcrec = bc_type[dir];
        int bc, bndry;
        if (side == Orientation::low) {
            bc = ibcrec.lo(dir);
            bndry = dlo;
        } else {
            bc = ibcrec.hi(dir);
            bndry = dhi;
        }

        // limit influx/outflux calculations to the in-out boundaries only
        // needs to change later?
        if (bc == BCType::direction_dependent) {
            for (MFIter mfi(*vel_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

                Box box = mfi.validbox();
//Print() << "validbox = " << box << std::endl;

                // include ghost cells for cell-centered
                // not for face-centered as boundary lies in valid region
                if (dir_index_type == IndexType::CellIndex::CELL) {
                    box.grow(dir, 1);
                }

                // include boundary corners if specified
                // this is relevant for cell-centered vels only
                // make this automatic based on cell-centered check ????
                if (corners) {
                    box.grow((dir+1)%AMREX_SPACEDIM, 1);
                    box.grow((dir+2)%AMREX_SPACEDIM, 1);
                }

                // Enter further only if the box bndry is at the domain bndry
                if ((islow && (box.smallEnd(dir) == dlo))
                 || (ishigh && (box.bigEnd(dir) == dhi))) {

                    // create a 2D box normal to dir at the low/high bndry
                    Box box2d(box); box2d.setRange(dir, bndry);

                    auto mac_vel = vel_mf->array(mfi);
                    auto in_mask = inflow_mask.array(mfi);
                    auto out_mask = outflow_mask.array(mfi);
Print() << "looping over 2d box: " << box2d << std::endl;

                    // tag cells as inflow or outflow by checking vel direction
                    ParallelFor(box2d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        if ((side == Orientation::low && mac_vel(i,j,k) >= 0)
                         || (side == Orientation::high && mac_vel(i,j,k) <= 0)) {
//Print() << "inflow at: " << i << " " << j << " " << k
//        << "  mac_vel = " << mac_vel(i,j,k) << std::endl;
                            in_mask(i,j,k) = 1;
                        } else {
//Print() << "outflow at: " << i << " " << j << " " << k
//        << "  mac_vel = " << mac_vel(i,j,k) << std::endl;
                            out_mask(i,j,k) = 1;
                        }
                    });
                }
            }
        }
    }
}

void compute_influx_outflux(
    const int lev,
    const Vector<Array<MultiFab*, AMREX_SPACEDIM>>& vels_vec,
    const Array<iMultiFab, AMREX_SPACEDIM>& inflow_masks,
    const Array<iMultiFab, AMREX_SPACEDIM>& outflow_masks,
    const Real* a_dx,
    Real& influx,
    Real& outflux,
    const bool corners)
{
    influx = 0.0, outflux = 0.0;

    // loop over the three dimensions
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {

        // normal face area
        const Real ds =
            a_dx[(idim+1) % AMREX_SPACEDIM] * a_dx[(idim+2) % AMREX_SPACEDIM];
//Print() << "ds is " << ds << std::endl;
        // Multifab for normal mac velocity
        auto& vel_mf = vels_vec[lev][idim];
//Print() << vel_mf->boxArray() << std::endl;

        // grow in the respective direction if vel is cell-centered
        IndexType index_type = vel_mf->ixType();
        index_type.flip(idim); IntVect ngrow = index_type.ixType();

        // grow in the transverse direction to include boundary corners
        // make this automatic based on cell-centered check ????
        if (corners) {
            ngrow[(idim+1)%AMREX_SPACEDIM] = 1;
            ngrow[(idim+2)%AMREX_SPACEDIM] = 1;
        }

        // mask iMFs for the respective velocity direction
        auto& inflow_mask = inflow_masks[idim];
        auto& outflow_mask = outflow_masks[idim];

        auto const& mac_vel_ma = vel_mf->const_arrays();
        auto const& inflow_mask_ma = inflow_mask.const_arrays();
        auto const& outflow_mask_ma = outflow_mask.const_arrays();

        influx += ds *
            ParReduce(TypeList<ReduceOpSum>{},
                     TypeList<Real>{},
                     *vel_mf, ngrow,
           [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
               noexcept -> GpuTuple<Real>
           {
               if (inflow_mask_ma[box_no](i,j,k)) {
//Print() << "counting inflow at: "<< i << " " << j << " " << k
//        << "  mac_vel = " << mac_vel_ma[box_no](i,j,k) << std::endl;
                   return { std::abs(mac_vel_ma[box_no](i,j,k)) };
               } else {
                   return { 0. };
               }
           });

        outflux += ds *
            ParReduce(TypeList<ReduceOpSum>{},
                     TypeList<Real>{},
                     *vel_mf, ngrow,
           [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
               noexcept -> GpuTuple<Real>
           {
               if (outflow_mask_ma[box_no](i,j,k)) {
//Print() << "counting outflow at: "<< i << " " << j << " " << k
//        << "  mac_vel = " << mac_vel_ma[box_no](i,j,k) << std::endl;
                   return { std::abs(mac_vel_ma[box_no](i,j,k)) };
               } else {
                   return { 0. };
               }
           });
    }

    ParallelDescriptor::ReduceRealSum(influx);
    ParallelDescriptor::ReduceRealSum(outflux);
Print() << "##### total influx is " << influx << std::endl;
Print() << "##### total outflux is " << outflux << std::endl;
}

void correct_outflow(
    const int lev,
    const Vector<Array<MultiFab*, AMREX_SPACEDIM>>& vels_vec,
    const Array<iMultiFab, AMREX_SPACEDIM>& outflow_masks,
    const BCRec* bc_type,
    const Box& domain,
    const Real alpha,
    const bool corners)
{
    // loop over the six orientations
    for (OrientationIter oit; oit != nullptr; ++oit) {
        const auto ori = oit();
        const auto side = ori.faceDir();
        const int dir = ori.coordDir();
        const auto islow = ori.isLow();
        const auto ishigh = ori.isHigh();

        // Multifab for normal mac velocity
        auto& vel_mf = vels_vec[lev][dir];
//Print() << vel_mf->boxArray() << std::endl;
        // mask iMFs for the respective velocity direction
        auto& outflow_mask = outflow_masks[dir];

        // domain extent indices for the velocities
        IndexType::CellIndex dir_index_type = (vel_mf->ixType()).ixType(dir);
        int dlo;
        if (dir_index_type == IndexType::CellIndex::CELL) {
            dlo = domain.smallEnd(dir) - 1; // cell-centered boundary
        } else {
            dlo = domain.smallEnd(dir);     // face-centered boundary
        }
        int dhi = domain.bigEnd(dir) + 1;

        // get BCs for the normal velocity and set the boundary index
        const BCRec ibcrec = bc_type[dir];
        int bc, bndry;
        if (side == Orientation::low) {
            bc = ibcrec.lo(dir);
            bndry = dlo;
        } else {
            bc = ibcrec.hi(dir);
            bndry = dhi;
        }

        if (bc == BCType::direction_dependent) {
            for (MFIter mfi(*vel_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

                Box box = mfi.validbox();
//Print() << "validbox = " << box << std::endl;
                if (dir_index_type == IndexType::CellIndex::CELL) {
                    box.grow(dir, 1);
                }
                if (corners) {
                    box.grow((dir+1)%AMREX_SPACEDIM, 1);
                    box.grow((dir+2)%AMREX_SPACEDIM, 1);
                }

                // Enter further only if the box boundary is at the domain boundary
                if ((islow && (box.smallEnd(dir) == dlo))
                 || (ishigh && (box.bigEnd(dir) == dhi))) {

                    // create a 2D box normal to dir at the low/high boundary
                    Box box2d(box); box2d.setRange(dir, bndry);

                    auto mac_vel = vel_mf->array(mfi);
                    auto out_mask = outflow_mask.array(mfi);
//Print() << "looping over 2d box: " << box2d << std::endl;
                    ParallelFor(box2d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        if (out_mask(i,j,k)) {
                            mac_vel(i,j,k) *= alpha;
                        }
                    });
                }
            }
        }
    }
}

} // file-local namespace

// !!!!!!! need to change mac-specific variable names
void enforceInOutSolvability (
    const Vector<Array<MultiFab*, AMREX_SPACEDIM>>& vels_vec,
    const BCRec* bc_type,
    const Vector<Geometry>& geom,
    const bool include_bndry_corners
)
{
    // get the level zero domain
    const Box domain = geom[0].Domain();

    const auto nlevs = int(vels_vec.size());
    for (int lev = 0; lev < nlevs; ++lev) {

        // masks to tag in/out flow at in-out boundaries
        // separate iMultifab for each velocity direction
        Array<iMultiFab, AMREX_SPACEDIM> inflow_masks;
        Array<iMultiFab, AMREX_SPACEDIM> outflow_masks;

        for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
        {
            auto& vel_mf = vels_vec[lev][idim];    // normal velocity multifab

            // grow in the respective direction if vel is cell-centered
            // to include the boundary cells
            IndexType index_type = vel_mf->ixType();
            index_type.flip(idim); IntVect ngrow = index_type.ixType();

            // grow in the transverse direction to include boundary corners
            if (include_bndry_corners) {
                ngrow[(idim+1)%AMREX_SPACEDIM] = 1;
                ngrow[(idim+2)%AMREX_SPACEDIM] = 1;
            }

            inflow_masks[idim].define(vel_mf->boxArray(), vel_mf->DistributionMap(), 1, ngrow);
            inflow_masks[idim].setVal(0);
            outflow_masks[idim].define(vel_mf->boxArray(), vel_mf->DistributionMap(), 1, ngrow);
            outflow_masks[idim].setVal(0);
        }
        set_inout_masks(lev, vels_vec, inflow_masks, outflow_masks, bc_type, domain, include_bndry_corners);

        const Real* a_dx = geom[lev].CellSize();
        Real influx = 0.0, outflux = 0.0;
        // now calculate the influx and outflux separately
        compute_influx_outflux(lev, vels_vec, inflow_masks, outflow_masks, a_dx, influx, outflux, include_bndry_corners);

        // apply correction factor to outflow
Print() << "##### Correcting outflow to match with inflow" << std::endl;
        const Real alpha = influx/outflux;
        correct_outflow(lev, vels_vec, outflow_masks, bc_type, domain, alpha, include_bndry_corners);

        // verify flux balance
        compute_influx_outflux(lev, vels_vec, inflow_masks, outflow_masks, a_dx, influx, outflux, include_bndry_corners);

    }   // levels loop
}

}
