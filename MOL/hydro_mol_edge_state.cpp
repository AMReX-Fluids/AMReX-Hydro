/** \addtogroup MOL
 *  @{
 */

#include <hydro_mol.H>
#include <hydro_mol_edge_state_K.H>

using namespace amrex;

namespace {
    std::pair<bool,bool> has_extdir_or_ho (BCRec const* bcrec, int ncomp, int dir)
    {
        std::pair<bool,bool> r{false,false};
        for (int n = 0; n < ncomp; ++n)
        {
            r.first = r.first || bcrec[n].lo(dir) == BCType::ext_dir
                              || bcrec[n].lo(dir) == BCType::hoextrap;
            r.second = r.second || bcrec[n].hi(dir) == BCType::ext_dir
                                || bcrec[n].hi(dir) == BCType::hoextrap;
        }
        return r;
    }
}



//
// Compute edge state on REGULAR box
//
void
MOL::ComputeEdgeState (const Box& bx,
                       AMREX_D_DECL( Array4<Real> const& xedge,
                                     Array4<Real> const& yedge,
                                     Array4<Real> const& zedge),
                       Array4<Real const> const& q,
                       const int ncomp,
                       AMREX_D_DECL( Array4<Real const> const& umac,
                                     Array4<Real const> const& vmac,
                                     Array4<Real const> const& wmac),
                       const Box&       domain,
                       const Vector<BCRec>& bcs,
                       const        BCRec * d_bcrec_ptr,
                       bool         is_velocity)
{
    const int domain_ilo = domain.smallEnd(0);
    const int domain_ihi = domain.bigEnd(0);
    const int domain_jlo = domain.smallEnd(1);
    const int domain_jhi = domain.bigEnd(1);
#if (AMREX_SPACEDIM==3)
    const int domain_klo = domain.smallEnd(2);
    const int domain_khi = domain.bigEnd(2);
#endif

    AMREX_D_TERM( const Box& ubx = amrex::surroundingNodes(bx,0);,
                  const Box& vbx = amrex::surroundingNodes(bx,1);,
                  const Box& wbx = amrex::surroundingNodes(bx,2););

    // At an ext_dir boundary, the boundary value is on the face, not cell center.
    auto extdir_lohi = has_extdir_or_ho(bcs.dataPtr(), ncomp, 0);
    bool has_extdir_or_ho_lo = extdir_lohi.first;
    bool has_extdir_or_ho_hi = extdir_lohi.second;

    if ((has_extdir_or_ho_lo && domain_ilo >= ubx.smallEnd(0)-1) ||
        (has_extdir_or_ho_hi && domain_ihi <= ubx.bigEnd(0)))
    {
        amrex::ParallelFor(ubx, ncomp, [d_bcrec_ptr,q,domain_ilo,domain_ihi,umac,xedge,is_velocity]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            xedge(i,j,k,n) = MOL::hydro_mol_xedge_state_extdir( i, j, k, n, q, umac,
                                                                d_bcrec_ptr,
                                                                domain_ilo, domain_ihi,
                                                                is_velocity);
        });
    }
    else
    {
        amrex::ParallelFor(ubx, ncomp, [d_bcrec_ptr,q,domain_ilo,domain_ihi,umac,xedge,is_velocity]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            xedge(i,j,k,n) = MOL::hydro_mol_xedge_state( i, j, k, n, q, umac,
                                                         d_bcrec_ptr,
                                                         domain_ilo, domain_ihi,
                                                         is_velocity);
        });
    }

    extdir_lohi = has_extdir_or_ho(bcs.dataPtr(), ncomp, 1);
    has_extdir_or_ho_lo = extdir_lohi.first;
    has_extdir_or_ho_hi = extdir_lohi.second;
    if ((has_extdir_or_ho_lo && domain_jlo >= vbx.smallEnd(1)-1) ||
        (has_extdir_or_ho_hi && domain_jhi <= vbx.bigEnd(1)))
    {
        amrex::ParallelFor(vbx, ncomp, [d_bcrec_ptr,q,domain_jlo,domain_jhi,vmac,yedge,is_velocity]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            yedge(i,j,k,n) = MOL::hydro_mol_yedge_state_extdir( i, j, k, n, q, vmac,
                                                                d_bcrec_ptr,
                                                                domain_jlo, domain_jhi,
                                                                is_velocity);
        });
    }
    else
    {
        amrex::ParallelFor(vbx, ncomp, [d_bcrec_ptr,q,domain_jlo,domain_jhi,vmac,yedge,is_velocity]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            yedge(i,j,k,n) = MOL::hydro_mol_yedge_state( i, j, k, n, q, vmac,
                                                         d_bcrec_ptr,
                                                         domain_jlo, domain_jhi,
                                                         is_velocity);
        });
    }


#if ( AMREX_SPACEDIM ==3 )

    extdir_lohi = has_extdir_or_ho(bcs.dataPtr(), ncomp, 2);
    has_extdir_or_ho_lo = extdir_lohi.first;
    has_extdir_or_ho_hi = extdir_lohi.second;
    if ((has_extdir_or_ho_lo && domain_klo >= wbx.smallEnd(2)-1) ||
        (has_extdir_or_ho_hi && domain_khi <= wbx.bigEnd(2)))
    {
        amrex::ParallelFor(wbx, ncomp, [d_bcrec_ptr,q,domain_klo,domain_khi,wmac,zedge,is_velocity]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            zedge(i,j,k,n) = MOL::hydro_mol_zedge_state_extdir( i, j, k, n, q, wmac,
                                                                d_bcrec_ptr,
                                                                domain_klo, domain_khi,
                                                                is_velocity);
        });
    }
    else
    {
        amrex::ParallelFor(wbx, ncomp, [d_bcrec_ptr,q,domain_klo,domain_khi,wmac,zedge,is_velocity]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            zedge(i,j,k,n) = MOL::hydro_mol_zedge_state( i, j, k, n, q, wmac,
                                                         d_bcrec_ptr,
                                                         domain_klo, domain_khi,
                                                         is_velocity);
        });
    }

#endif
}
/** @}*/
