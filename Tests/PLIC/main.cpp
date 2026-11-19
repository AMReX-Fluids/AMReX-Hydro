// Verification of the PLIC library, as a small solver would use it. The run fails (non-zero exit)
// when a check does not hold.
//
//   test.kind = kernels   the plane geometry, the normals, the plane fit (test_geometry.cpp), the
//                         wall ghost cells (test_walls.cpp) and the face density (test_density.cpp)
//   test.kind = momentum  the mass and momentum fluxes, with a density ratio of 1000
//                         (test_momentum.cpp)
//   test.kind = advect    a prescribed velocity on an AMR hierarchy (see inputs.*), advanced with
//                         PLIC::ComputeEdgeStateAndFluxes in the solver's own MFIter loop, or with
//                         HydroUtils::ComputeFluxesOnBoxFromState and advection_type "PLIC" if
//                         test.use_hydroutils = 1 (test_advect.cpp)

#include "test_util.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        amrex::ParmParse pp("test");
        std::string kind = "kernels";
        pp.query("kind", kind);
        if (kind == "kernels") {
            plic_test::check_geometry();
            plic_test::check_walls();
            plic_test::check_face_density();
            plic_test::check_fluxes();
        } else if (kind == "momentum") {
            plic_test::check_momentum();
        } else if (kind == "advect") {
            plic_test::check_advect();
        } else {
            amrex::Abort("test.kind must be kernels, momentum or advect");
        }
        amrex::Print() << ((plic_test::failures() == 0) ? "PASSED" : "FAILED") << "\n";
    }
    const int failures = plic_test::failures();
    amrex::Finalize();
    return (failures == 0) ? 0 : 1;
}
