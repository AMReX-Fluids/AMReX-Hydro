// The wall ghost cells: FillWallGhosts

#include "test_util.H"

#include <AMReX_Print.H>

namespace plic_test {

namespace {

constexpr int n = 32;
const std::string reconstruction = "myc_height";

// Walls at the low and high y faces; the other faces keep the ghost values they are given (ext_dir)
BCRec y_walls ()
{
    BCRec bc = uniform_bc(BCType::ext_dir);
    bc.setLo(1, BCType::foextrap);
    bc.setHi(1, BCType::foextrap);
    return bc;
}

// The ghost cells beyond the y walls, over the valid box of each fab
Vector<Box> beyond_y_walls (Box const& valid, Box const& domain, int ng)
{
    Box b = valid;
    b.grow(1, ng);
    return beyond_faces(b, domain, 1);
}

// Inflow and symmetry faces keep the ghost values of the solver: ext_dir on the low x side,
// reflect_even on the high x side, and walls elsewhere
void check_kept_faces ()
{
    const Geometry geom = unit_geometry(n, IntVect(0));
    const Box domain = geom.Domain();
    BoxArray ba(domain);
    ba.maxSize(8);
    DistributionMapping dm(ba);
    const int ng = PLIC::nGrowF(reconstruction);
    const Plane p = centered_plane({AMREX_D_DECL(1.0_rt, 0.4_rt, 0.1_rt)});

    MultiFab F(ba, dm, 1, ng);
    fill_plane(F, geom, p, IntVect(1));
    MultiFab Fkeep(ba, dm, 1, ng);
    MultiFab::Copy(Fkeep, F, 0, 0, 1, ng);
    BCRec bc = uniform_bc(BCType::foextrap);
    bc.setLo(0, BCType::ext_dir);
    bc.setHi(0, BCType::reflect_even);
    PLIC::FillWallGhosts(Fkeep, geom, bc, reconstruction);

    MultiFab::Subtract(Fkeep, F, 0, 0, 1, ng);
    const Real err = max_abs(Fkeep, [&] (Box const&, Box const& fab) {
        return beyond_faces(fab, domain, 0);
    });
    Print() << "  largest change beyond ext_dir and reflect_even faces " << err << "\n";
    check(err == 0.0_rt, "ext_dir and reflect_even ghost cells are kept");
}

// A plane meeting walls on every face continues as the same plane beyond them, away from the other
// walls: the stencils near an edge read the ghost cells of the other wall, and the search along the
// wall stops at it. Also with more ghost layers in F than the search reads.
void check_plane_through_walls ()
{
    const Geometry geom = unit_geometry(n, IntVect(0));
    const Box domain = geom.Domain();
    BoxArray ba(domain);
    ba.maxSize(8);
    DistributionMapping dm(ba);
    const int ng_min = PLIC::nGrowF(reconstruction);
    const Plane p = centered_plane({AMREX_D_DECL(1.0_rt, 0.4_rt, 0.1_rt)});

    for (const int ng : {ng_min, ng_min + PLIC::wall_nsearch + 1}) {
        MultiFab F(ba, dm, 1, ng);
        fill_plane(F, geom, p, IntVect(1));
        PLIC::FillWallGhosts(F, geom, uniform_bc(BCType::foextrap), reconstruction);

        const Real err = max_plane_error(F, geom, p, [&] (Box const&, Box const& fab) {
            Vector<Box> boxes;
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                const Box inner = away_from_other_faces(domain, dir, ng_min + PLIC::wall_nsearch, ng);
                for (Box const& b : beyond_faces(fab & inner, domain, dir)) { boxes.push_back(b); }
            }
            return boxes;
        });
        Print() << "  largest wall ghost error for a plane, " << ng << " ghost layers " << err << "\n";
        check(err < 100.0_rt * roundoff_tol, "wall ghost cells continue a planar interface");
    }
}

// A fine level covering part of the domain along the y walls: past its coarse/fine boundaries, F
// holds only nGrowF ghost layers of data (here the plane, as coarse interpolation would give), and
// the wall fill must not take the layers beyond them for data. The covered range in x is shifted so
// that the contact line on the low wall passes the coarse/fine boundary. Not periodic, as the plane
// is not.
void check_fine_level ()
{
    const Geometry geom = unit_geometry(n, IntVect(0));
    const Box domain = geom.Domain();
    const int ng = PLIC::nGrowF(reconstruction);
    // Liquid towards high x, so that layers read as zeros would be a gas region
    const Plane p = centered_plane({AMREX_D_DECL(-1.0_rt, -0.4_rt, -0.5_rt)});

    Real err = 0.0_rt;
    for (int xhi = 15; xhi <= 25; ++xhi) {
        Box covered = domain;
        covered.setSmall(0, xhi - 15);
        covered.setBig(0, xhi);
        BoxArray ba(covered);
        ba.maxSize(8);
        DistributionMapping dm(ba);
        MultiFab F(ba, dm, 1, ng);
        fill_plane(F, geom, p, IntVect::TheDimensionVector(1));
        PLIC::FillWallGhosts(F, geom, y_walls(), reconstruction);
        err = max(err, max_plane_error(F, geom, p, [&] (Box const& valid, Box const&) {
            return beyond_y_walls(valid, domain, ng);
        }));
    }
    Print() << "  largest wall ghost error for a plane, fine level " << err << "\n";
    check(err < 100.0_rt * roundoff_tol, "wall ghost cells of a fine level do not read past its ghost cells");
}

// Periodic in x, with the boundary condition of x left unset (bogus), as a solver may do: the fill
// must not read it, and must still continue a plane through the y walls
void check_periodic_unset_bc ()
{
    const Geometry geom = unit_geometry(n, IntVect::TheDimensionVector(0));
    const Box domain = geom.Domain();
    const int ng = PLIC::nGrowF(reconstruction);
    // Independent of x, so periodic in x
    const Plane p = centered_plane({AMREX_D_DECL(0.0_rt, 0.4_rt, 1.0_rt)});
    BCRec bc = y_walls();
    bc.setLo(0, BCType::bogus);
    bc.setHi(0, BCType::bogus);

    BoxArray ba(domain);
    ba.maxSize(8);
    DistributionMapping dm(ba);
    MultiFab F(ba, dm, 1, ng);
    fill_plane(F, geom, p, IntVect::TheDimensionVector(1));
    PLIC::FillWallGhosts(F, geom, bc, reconstruction);
    const Real err = max_plane_error(F, geom, p, [&] (Box const& valid, Box const&) {
        return beyond_y_walls(valid, domain, ng);
    });
    Print() << "  largest wall ghost error for a plane, x periodic with an unset condition " << err << "\n";
    check(err < 100.0_rt * roundoff_tol, "the condition of a periodic direction is not read");
}

// Each reconstruction continues a plane through a wall where its plane is exact, at each contact
// angle range of wall_contact_plane: the plane of the reconstruction (45 and 68 degrees), the MYC
// plane of a cell below the wall layer (11), and the height function inside the wall layer (3). The
// planes meet the low y wall at its center, and the cells checked lie within the search along the
// wall of that contact line.
void check_contact_angles ()
{
    const Geometry geom = unit_geometry(n, IntVect(0));
    const Box domain = geom.Domain();
    BoxArray ba(domain);
    ba.maxSize(8);
    DistributionMapping dm(ba);
    struct Case {
        const char* reconstruction;
        GpuArray<Real, AMREX_SPACEDIM> n;
        const char* angle;
    };
    // Angles between the interface and the wall in the xy plane. youngs is exact only at 45 degrees,
    // where the plane is symmetric about its normal.
    const Case cases[] = {
        {"youngs",     {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)}, "45 degrees"},
        {"myc",        {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)}, "45 degrees"},
        {"myc_height", {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)}, "45 degrees"},
        {"myc",        {AMREX_D_DECL(1.0_rt, 0.4_rt, 0.1_rt)}, "68 degrees"},
        {"myc_height", {AMREX_D_DECL(1.0_rt, 0.4_rt, 0.1_rt)}, "68 degrees"},
        {"myc",        {AMREX_D_DECL(0.2_rt, 1.0_rt, 0.1_rt)}, "11 degrees"},
        {"myc_height", {AMREX_D_DECL(0.2_rt, 1.0_rt, 0.1_rt)}, "11 degrees"},
        {"myc",        {AMREX_D_DECL(0.06_rt, 1.0_rt, 0.04_rt)}, "3 degrees"},
        {"myc_height", {AMREX_D_DECL(0.06_rt, 1.0_rt, 0.04_rt)}, "3 degrees"},
    };
    for (Case const& c : cases) {
        const std::string rec(c.reconstruction);
        const int ng = PLIC::nGrowF(rec);
        const Plane p{c.n, AMREX_D_TERM(0.5_rt * c.n[0], + 0.0_rt, + 0.5_rt * c.n[2])};
        MultiFab F(ba, dm, 1, ng);
        fill_plane(F, geom, p, IntVect::TheDimensionVector(1));
        PLIC::FillWallGhosts(F, geom, y_walls(), rec);
        const Real err = max_plane_error(F, geom, p, [&] (Box const& valid, Box const&) {
            const Box near_contact = away_from_other_faces(domain, 1, ng + PLIC::wall_nsearch, ng);
            return beyond_y_walls(valid & near_contact, domain, ng);
        });
        Print() << "  largest wall ghost error, " << rec << " at " << c.angle << " " << err << "\n";
        check(err < 100.0_rt * roundoff_tol, rec + " continues a plane through a wall at " + c.angle);
    }
}

// The fill does not depend on the box layout: the same curved interface, a sphere cutting the low y
// wall, with walls on every face, gives the same ghost values with one box, with boxes of 16 cells
// and of 8 cells. Each result is gathered on one box covering the domain and its ghost cells. Also
// periodic in x, with the sphere across the periodic boundary, where the search along the wall
// reads the periodic images of the boxes.
void check_layout_independence ()
{
    for (const int periodic_x : {0, 1}) {
        const Geometry geom = unit_geometry(n, IntVect(AMREX_D_DECL(periodic_x, 0, 0)));
        const Box domain = geom.Domain();
        const auto geomdata = geom.data();
        const Real cx = (periodic_x != 0) ? 0.05_rt : 0.4_rt;
        const std::string where = (periodic_x != 0) ? ", x periodic" : "";
        for (const char* name : {"youngs", "myc", "myc_height"}) {
            const std::string rec(name);
            const int ng = PLIC::nGrowF(rec);
            const BoxArray gathered_ba(grow(domain, ng));
            const DistributionMapping gathered_dm(gathered_ba);
            Vector<MultiFab> gathered;
            for (const int max_size : {n, 16, 8}) {
                BoxArray ba(domain);
                ba.maxSize(max_size);
                DistributionMapping dm(ba);
                MultiFab F(ba, dm, 1, ng);
                // A ramp of two cells across the sphere, extrapolated beyond the walls, and periodic
                // in x if x is
                auto const& fm = F.arrays();
                ParallelFor(F, IntVect(ng), [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept {
                    amrex::ignore_unused(k);
                    const int ii = (periodic_x != 0) ? i : Clamp(i, 0, n - 1);
                    const IntVect iv(AMREX_D_DECL(ii, Clamp(j, 0, n - 1), Clamp(k, 0, n - 1)));
                    const Real center[3] = {cx, 0.1_rt, 0.45_rt};
                    Real r2 = 0.0_rt;
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                        Real x = geomdata.ProbLo(d) + (static_cast<Real>(iv[d]) + 0.5_rt) * geomdata.CellSize(d)
                                 - center[d];
                        if (d == 0 && periodic_x != 0) {
                            x -= std::round(x);
                        }
                        r2 += x * x;
                    }
                    fm[b](i, j, k) = Clamp(0.5_rt + (0.3_rt - std::sqrt(r2)) / (2.0_rt * geomdata.CellSize(0)),
                                           0.0_rt, 1.0_rt);
                });
                Gpu::streamSynchronize();
                PLIC::FillWallGhosts(F, geom, uniform_bc(BCType::foextrap), rec);
                gathered.emplace_back(gathered_ba, gathered_dm, 1, 0);
                gathered.back().ParallelCopy(F, 0, 0, 1, ng, 0);
            }
            Real err = 0.0_rt;
            for (int l = 1; l < gathered.size(); ++l) {
                MultiFab::Subtract(gathered[l], gathered[0], 0, 0, 1, 0);
                err = max(err, max_abs(gathered[l], valid_only));
            }
            Print() << "  largest difference between box layouts, " << rec << where << " " << err << "\n";
            check(err == 0.0_rt, "wall ghost cells do not depend on the box layout, " + rec + where);
        }
    }
}

} // namespace

void check_walls ()
{
    Print() << "wall ghost cells\n";
    check_kept_faces();
    check_plane_through_walls();
    check_fine_level();
    check_periodic_unset_bc();
    check_contact_angles();
    check_layout_independence();
}

} // namespace plic_test
