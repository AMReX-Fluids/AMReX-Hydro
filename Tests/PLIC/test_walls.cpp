// The wall ghost cells: FillWallGhosts

#include "test_util.H"

#include <AMReX_Print.H>

#include <random>

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

// Largest |F - G| over every cell of every fab of F, ghost cells included, with G one box that holds
// them all: a ghost cell that several fabs hold is compared in each, where a ParallelCopy from F
// would keep only one of them
Real max_fab_difference (MultiFab const& F, MultiFab const& G)
{
    MultiFab R(F.boxArray(), F.DistributionMap(), 1, F.nGrowVect());
    R.ParallelCopy(G, 0, 0, 1, IntVect(0), F.nGrowVect());
    MultiFab::Subtract(R, F, 0, 0, 1, F.nGrowVect());
    return max_abs(R, [] (Box const& /*valid*/, Box const& fab) { return Vector<Box>{fab}; });
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

// Planes meeting walls on every face continue as the same plane in the nGrowF ghost layers, edges
// and corners included. Also with more ghost layers in F than the search reads (deeper, a plane
// meets ghost cells further along the wall from its contact line than the search reaches).
void check_plane_through_walls ()
{
    const Geometry geom = unit_geometry(n, IntVect(0));
    const Box domain = geom.Domain();
    BoxArray ba(domain);
    ba.maxSize(8);
    DistributionMapping dm(ba);
    const int ng_min = PLIC::nGrowF(reconstruction);
    for (const int ng : {ng_min, ng_min + PLIC::wall_nsearch + 1}) {
        Real err = 0.0_rt;
        for (auto const& nrm : {GpuArray<Real, AMREX_SPACEDIM>{AMREX_D_DECL(1.0_rt, 0.4_rt, 0.1_rt)},
                                GpuArray<Real, AMREX_SPACEDIM>{AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)},
                                GpuArray<Real, AMREX_SPACEDIM>{AMREX_D_DECL(0.3_rt, 1.0_rt, -0.6_rt)}}) {
            const Plane p = centered_plane(nrm);
            MultiFab F(ba, dm, 1, ng);
            fill_plane(F, geom, p, IntVect(1));
            PLIC::FillWallGhosts(F, geom, uniform_bc(BCType::foextrap), reconstruction);
            err = max(err, max_plane_error(F, geom, p, [&] (Box const&, Box const& fab) {
                return boxDiff(fab & grow(domain, ng_min), domain).data();
            }));
        }
        Print() << "  largest wall ghost error for planes, edges and corners included, " << ng
                << " ghost layers " << err << "\n";
        check(err < 100.0_rt * roundoff_tol, "wall ghost cells continue a planar interface, edges and corners included");
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
// angle range of wall_contact_planes: the plane of the reconstruction (45 and 68 degrees), the MYC
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

// myc_height continues random planes through a wall, each meeting the low y wall at its center:
// at any contact angle one of the planes of wall_contact_planes is exact. The cells checked lie
// beyond the low wall, within the search along it of the contact line.
void check_random_contact_planes ()
{
    const Geometry geom = unit_geometry(n, IntVect(0));
    const Box domain = geom.Domain();
    BoxArray ba(domain);
    ba.maxSize(8);
    DistributionMapping dm(ba);
    const int ng = PLIC::nGrowF(reconstruction);
    std::mt19937 gen(11);
    std::normal_distribution<Real> normal_dist;
    Real err = 0.0_rt;
    GpuArray<Real, AMREX_SPACEDIM> worst{};
    // Half in any orientation, half with slopes relative to the wall near 2a + b = 1, where the MYC
    // plane stops being exact and the search must reach a plane of myc_height (see wall_nsearch)
    std::uniform_real_distribution<Real> slope_a(0.45_rt, 0.6_rt);
    std::uniform_real_distribution<Real> slope_b(0.0_rt, 0.1_rt);
    for (int p = 0; p < 100; ++p) {
        GpuArray<Real, AMREX_SPACEDIM> nrm{AMREX_D_DECL(normal_dist(gen), normal_dist(gen), normal_dist(gen))};
        if (p % 2 == 1) {
            const Real sa = (p & 2) ? -1.0_rt : 1.0_rt;
            const Real a = slope_a(gen);
            const Real b = slope_b(gen);
            amrex::ignore_unused(b);  // 2D
            nrm = {AMREX_D_DECL(sa * a, (p & 4) ? -1.0_rt : 1.0_rt, b)};
        }
        const Plane pl{nrm, AMREX_D_TERM(0.5_rt * nrm[0], + 0.0_rt, + 0.5_rt * nrm[2])};
        MultiFab F(ba, dm, 1, ng);
        fill_plane(F, geom, pl, IntVect::TheDimensionVector(1));
        PLIC::FillWallGhosts(F, geom, y_walls(), reconstruction);
        const Real e = max_plane_error(F, geom, pl, [&] (Box const& valid, Box const&) {
            const Box near_contact = away_from_other_faces(domain, 1, ng + PLIC::wall_nsearch, ng);
            Box b = valid & near_contact;
            b.grow(1, ng);
            return Vector<Box>{beyond_faces(b, domain, 1)[0]};
        });
        if (e > err) {
            err = e;
            worst = nrm;
        }
    }
    Print() << "  largest wall ghost error over random contact planes " << err << ", normal";
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { Print() << " " << worst[d]; }
    Print() << "\n";
    check(err < 100.0_rt * roundoff_tol, "myc_height continues random planes through a wall");
}

// F on the domain of n^dim unit-length cells and its ghost cells: a ramp of two cells across the
// sphere of radius 0.3 about center, extrapolated beyond the walls, and periodic in x if periodic_x
void fill_sphere (MultiFab& F, Geometry const& geom, GpuArray<Real, 3> const& center, bool periodic_x)
{
    const auto geomdata = geom.data();
    auto const& fm = F.arrays();
    ParallelFor(F, F.nGrowVect(), [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept {
        amrex::ignore_unused(k);
        const int ii = periodic_x ? i : Clamp(i, 0, n - 1);
        const IntVect iv(AMREX_D_DECL(ii, Clamp(j, 0, n - 1), Clamp(k, 0, n - 1)));
        Real r2 = 0.0_rt;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            Real x = geomdata.ProbLo(d) + (static_cast<Real>(iv[d]) + 0.5_rt) * geomdata.CellSize(d) - center[d];
            if (d == 0 && periodic_x) {
                x -= std::round(x);
            }
            r2 += x * x;
        }
        fm[b](i, j, k) = Clamp(0.5_rt + (0.3_rt - std::sqrt(r2)) / (2.0_rt * geomdata.CellSize(0)),
                               0.0_rt, 1.0_rt);
    });
    Gpu::streamSynchronize();
}

// A plane that meets the low y wall only within the stencil of myc_height of the low x wall, nearly
// normal to x, continues through the y wall from the MYC planes of the cells there
void check_plane_near_edge ()
{
    const Geometry geom = unit_geometry(n, IntVect(0));
    const Box domain = geom.Domain();
    BoxArray ba(domain);
    ba.maxSize(8);
    DistributionMapping dm(ba);
    const int ng = PLIC::nGrowF(reconstruction);
    const Plane p{{AMREX_D_DECL(1.0_rt, -0.1_rt, 0.05_rt)}, 2.3_rt / static_cast<Real>(n)};
    MultiFab F(ba, dm, 1, ng);
    fill_plane(F, geom, p, IntVect(1));
    PLIC::FillWallGhosts(F, geom, uniform_bc(BCType::foextrap), reconstruction);
    Box below = domain;
    below.growLo(1, ng);
    const Real err = max_plane_error(F, geom, p, [&] (Box const&, Box const& fab) {
        return beyond_faces(fab & below, domain, 1);
    });
    Print() << "  largest wall ghost error for a plane near an edge " << err << "\n";
    check(err < 100.0_rt * roundoff_tol, "wall ghost cells continue a plane that meets the wall near an edge");
}

// The fill does not depend on the box layout: the same curved interface, a sphere cutting the low y
// wall, with walls on every face, gives the same ghost values with one box, with boxes of 16 cells
// and of 8 cells. Each result is gathered on one box covering the domain and its ghost cells. Also
// periodic in x, with the sphere across the periodic boundary, where the search along the wall
// reads the periodic images of the boxes, and with an ext_dir low x face, the sphere across its
// edge with the low y wall, whose ghost cells hold noise: the search reads no cell beyond it.
void check_layout_independence ()
{
    struct Case { bool periodic_x; bool ext_dir_x; Real cx; const char* where; };
    for (Case const& c : {Case{false, false, 0.4_rt, ""}, Case{true, false, 0.05_rt, ", x periodic"},
                          Case{false, true, 0.15_rt, ", ext_dir low x"}}) {
        const Geometry geom = unit_geometry(n, IntVect(AMREX_D_DECL(c.periodic_x ? 1 : 0, 0, 0)));
        const Box domain = geom.Domain();
        BCRec bc = uniform_bc(BCType::foextrap);
        if (c.ext_dir_x) {
            bc.setLo(0, BCType::ext_dir);
        }
        const std::string where(c.where);
        for (const char* name : {"youngs", "myc", "myc_height"}) {
            const std::string rec(name);
            const int ng = PLIC::nGrowF(rec);
            // The one-box layout, gathered on one box: each of its cells once
            const BoxArray gathered_ba(grow(domain, ng));
            MultiFab reference(gathered_ba, DistributionMapping(gathered_ba), 1, 0);
            Real err = 0.0_rt;
            for (const int max_size : {n, 16, 8}) {
                BoxArray ba(domain);
                ba.maxSize(max_size);
                DistributionMapping dm(ba);
                MultiFab F(ba, dm, 1, ng);
                fill_sphere(F, geom, {c.cx, 0.1_rt, 0.45_rt}, c.periodic_x);
                if (c.ext_dir_x) {
                    auto const& fm = F.arrays();
                    ParallelFor(F, F.nGrowVect(), [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept {
                        if (i < 0) {
                            // Noise, the same in every layout
                            const Real h = std::sin(12.9898_rt * static_cast<Real>(i) + 78.233_rt * static_cast<Real>(j)
                                                    + 37.719_rt * static_cast<Real>(k)) * 43758.5453_rt;
                            fm[b](i, j, k) = h - std::floor(h);
                        }
                    });
                    Gpu::streamSynchronize();
                }
                PLIC::FillWallGhosts(F, geom, bc, rec);
                if (max_size == n) {
                    reference.ParallelCopy(F, 0, 0, 1, ng, 0);
                } else {
                    err = max(err, max_fab_difference(F, reference));
                }
            }
            Print() << "  largest difference between box layouts, " << rec << where << " " << err << "\n";
            check(err == 0.0_rt, "wall ghost cells do not depend on the box layout, " + rec + where);
        }
    }
}

// The fill treats the walls and the directions along them alike: a curved interface with a symmetry
// (a sphere near the low x wall, centered on a plane of symmetry, with walls on every face) gives
// ghost values with the same symmetry, at edges and corners too. The symmetries: swapping x and y,
// x and z, and mirroring in x.
void check_wall_symmetry ()
{
    const Geometry geom = unit_geometry(n, IntVect(0));
    const Box domain = geom.Domain();
    BoxArray ba(domain);
    ba.maxSize(8);
    DistributionMapping dm(ba);
    Real err = 0.0_rt;
    Real fab_err = 0.0_rt;
    for (Symmetry const& sym : symmetries) {
        // On the plane of symmetry: x = y, x = z, or x = 1/2
        GpuArray<Real, 3> center{0.15_rt, 0.15_rt, 0.4_rt};
        if (sym.b == 0) { center[0] = 0.5_rt; }
        if (sym.b == 2) { amrex::Swap(center[1], center[2]); }
        for (const char* name : {"youngs", "myc", "myc_height"}) {
            const std::string rec(name);
            const int ng = PLIC::nGrowF(rec);
            MultiFab F(ba, dm, 1, ng);
            fill_sphere(F, geom, center, false);
            PLIC::FillWallGhosts(F, geom, uniform_bc(BCType::foextrap), rec);
            // Gathered on one box, then compared with its image
            const BoxArray gathered_ba(grow(domain, ng));
            MultiFab G(gathered_ba, DistributionMapping(gathered_ba), 1, 0);
            G.ParallelCopy(F, 0, 0, 1, ng, 0);
            // A ghost cell several boxes hold gets the same value in each, so G holds them all
            fab_err = max(fab_err, max_fab_difference(F, G));
            const Real e = max_over(G, valid_only,
                                    [=] AMREX_GPU_DEVICE (Array4<Real const> const& a, int i, int j, int k) noexcept -> Real {
                amrex::ignore_unused(k);
                const IntVect iv(AMREX_D_DECL(i, j, k));
                return std::abs(a(sym.cell(iv, n)) - a(iv));
            });
            Print() << "  largest difference from the image, " << sym.name << ", " << rec << " " << e << "\n";
            err = max(err, e);
        }
    }
    check(fab_err == 0.0_rt, "boxes sharing a wall ghost cell give it the same value");
    // To round-off: the normals do not treat the directions alike bit for bit
    check(err < roundoff_tol, "wall ghost cells keep the symmetries of the interface");
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
    check_random_contact_planes();
    check_plane_near_edge();
    check_layout_independence();
    check_wall_symmetry();
}

} // namespace plic_test
