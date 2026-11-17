// Verification of the PLIC library. The run fails (non-zero exit) when a check does not hold.
//
//   test.kind = kernels   the plane geometry, the normals and the plane fit (test_geometry.cpp)

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
        } else {
            amrex::Abort("test.kind must be kernels");
        }
        amrex::Print() << ((plic_test::failures() == 0) ? "PASSED" : "FAILED") << "\n";
    }
    const int failures = plic_test::failures();
    amrex::Finalize();
    return (failures == 0) ? 0 : 1;
}
