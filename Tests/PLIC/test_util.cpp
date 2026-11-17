#include "test_util.H"

#include <AMReX_Print.H>

namespace plic_test {

namespace {
int n_failures = 0;
} // namespace

void check (bool ok, std::string const& what)
{
    Print() << "  " << (ok ? "ok    " : "FAILED") << " " << what << "\n";
    if (!ok) { ++n_failures; }
}

int failures ()
{
    return n_failures;
}

} // namespace plic_test
