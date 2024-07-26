/*
 * \file hydro_godunov_ppm.H
 *
 * \addtogroup Godunov
 *  @{
 */

#ifndef HYDRO_PPM_GODUNOV_H
#define HYDRO_PPM_GODUNOV_H

#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <hydro_constants.H>

namespace PPM {

enum limiters {VanLeer, WENOZ, NoLimiter};

static constexpr int default_limiter = VanLeer;

struct nolimiter {

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real
    sedge1(const amrex::Real sm2,
           const amrex::Real sm1,
           const amrex::Real s0,
           const amrex::Real sp1,
           const amrex::Real /*sp2*/)
    {

        amrex::Real d1 = amrex::Real(0.5) * (sp1 - sm1);
        amrex::Real d2 = amrex::Real(0.5) * (s0 - sm2);

        amrex::Real sedge = amrex::Real(0.5)*(s0 + sm1) - amrex::Real(1./6.)*(d1 - d2);
        return sedge;

    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real
    sedge2(const amrex::Real /*sm2*/,
           const amrex::Real sm1,
           const amrex::Real s0,
           const amrex::Real sp1,
           const amrex::Real sp2)
    {
        amrex::Real d1 = amrex::Real(0.5) * (sp2 - s0);
        amrex::Real d2 = amrex::Real(0.5) * (sp1 - sm1);

        amrex::Real sedge = amrex::Real(0.5)*(sp1 + s0) - amrex::Real(1./6.)*(d1 - d2);
        return sedge;
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::GpuTuple<amrex::Real,amrex::Real>
    sm_sp(const amrex::Real /*s0*/,
         const amrex::Real sedge1,
         const amrex::Real sedge2)
    {
        return amrex::makeTuple(sedge1, sedge2);
    }
};

struct vanleer {

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real vanLeer(const amrex::Real a, const amrex::Real b, const amrex::Real c)
    {
        constexpr auto small_qty_sq = amrex::Real(1.e-10) * amrex::Real(1.e-10);

        amrex::Real dsc = amrex::Real(0.5)*(b - c);
        amrex::Real dsl = amrex::Real(2.0)*(a - c);
        amrex::Real dsr = amrex::Real(2.0)*(b - a);
        return (dsl*dsr > small_qty_sq) ?
            std::copysign(amrex::Real(1.0), dsc)*amrex::min(amrex::Math::abs(dsc),amrex::min(amrex::Math::abs(dsl),
                                                                                             amrex::Math::abs(dsr))) : amrex::Real(0.0);
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real
    sedge1(const amrex::Real sm2,
           const amrex::Real sm1,
           const amrex::Real s0,
           const amrex::Real sp1,
           const amrex::Real /*sp2*/)
    {

        amrex::Real d1 = vanLeer(s0,sp1,sm1);
        amrex::Real d2 = vanLeer(sm1,s0,sm2);

        amrex::Real sedge = amrex::Real(0.5)*(s0 + sm1) - amrex::Real(1./6.)*(d1 - d2);
        return amrex::min(amrex::max(sedge, amrex::min(s0, sm1)),amrex::max(s0,sm1));

    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real
    sedge2(const amrex::Real /*sm2*/,
           const amrex::Real sm1,
           const amrex::Real s0,
           const amrex::Real sp1,
           const amrex::Real sp2)
    {

        amrex::Real d1 = vanLeer(sp1,sp2,s0);
        amrex::Real d2 = vanLeer(s0,sp1,sm1);

        amrex::Real sedge = amrex::Real(0.5)*(sp1 + s0) - amrex::Real(1./6.)*(d1 - d2);
        return amrex::min(amrex::max(sedge, amrex::min(s0, sp1)),amrex::max(s0,sp1));


    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::GpuTuple<amrex::Real,amrex::Real>
    sm_sp(const amrex::Real s0,
         const amrex::Real sedge1,
         const amrex::Real sedge2)
    {
        amrex::Real sm_ = sedge1;
        amrex::Real sp_ = sedge2;

        if ((sedge2-s0)*(s0-sedge1) < 0.0) {
            sm_ = s0;
            sp_ = s0;
        } else if (amrex::Math::abs(sedge2-s0) >= amrex::Real(2.0)*amrex::Math::abs(sedge1-s0)) {
            sp_ = amrex::Real(3.0)*s0 - amrex::Real(2.0)*sedge1;
        } else if (amrex::Math::abs(sedge1-s0) >=  amrex::Real(2.0)*amrex::Math::abs(sedge2-s0)) {
            sm_ = amrex::Real(3.0)*s0 - amrex::Real(2.0)*sedge2;
        }
        return amrex::makeTuple(sm_,sp_);

        // this might be the correct way to do this in case the else if are both == 0.0
        // could also change the >= to just > for each
//        if ((sedge2-s0)*(s0-sedge1) < amrex::Real(0.0)) {
//            return amrex::makeTuple(s0,s0);
//        }
//        amrex::Real sm_ = sedge1;
//        amrex::Real sp_ = sedge2;
//        if (amrex::Math::abs(sedge2-s0) >= amrex::Real(2.0)*amrex::Math::abs(sedge1-s0)) {
//            sp_ = amrex::Real(3.0)*s0 - amrex::Real(2.0)*sedge1;
//        }
//        if (amrex::Math::abs(sedge1-s0) >=  amrex::Real(2.0)*amrex::Math::abs(sedge2-s0)) {
//            sm_ = amrex::Real(3.0)*s0 - amrex::Real(2.0)*sedge2;
//        }
//        return amrex::makeTuple(sm_,sp_);

    }
};

struct wenoz {

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real
    sedge2(const amrex::Real sm2,
           const amrex::Real sm1,
           const amrex::Real s0,
           const amrex::Real sp1,
           const amrex::Real sp2)
    {
        constexpr auto eps = amrex::Real(1.0e-6);

        const amrex::Real beta1 =
            amrex::Real(13.0) / amrex::Real(12.0) * (sm2 - amrex::Real(2.0) * sm1 + s0) * (sm2 - amrex::Real(2.0) * sm1 + s0) +
            amrex::Real(0.25) * (sm2 - amrex::Real(4.0) * sm1 + amrex::Real(3.0) * s0) * (sm2 - amrex::Real(4.0) * sm1 + amrex::Real(3.0) * s0);
        const amrex::Real beta2 =
            amrex::Real(13.0) / amrex::Real(12.0) * (sm1 - amrex::Real(2.0) * s0 + sp1) * (sm1 - amrex::Real(2.0) * s0 + sp1) +
            amrex::Real(0.25) * (sm1 - sp1) * (sm1 - sp1);
        const amrex::Real beta3 =
            amrex::Real(13.0) / amrex::Real(12.0) * (s0 - amrex::Real(2.0) * sp1 + sp2) * (s0 - amrex::Real(2.0) * sp1 + sp2) +
            amrex::Real(0.25) * (amrex::Real(3.0) * s0 - amrex::Real(4.0) * sp1 + sp2) * (amrex::Real(3.0) * s0 - amrex::Real(4.0) * sp1 + sp2);

        const amrex::Real t5 = amrex::Math::abs(beta3 - beta1);
        const amrex::Real omega1 = amrex::Real(0.1) * (amrex::Real(1.0) + t5 / (eps + beta1));
        const amrex::Real omega2 = amrex::Real(0.6) * (amrex::Real(1.0) + t5 / (eps + beta2));
        const amrex::Real omega3 = amrex::Real(0.3) * (amrex::Real(1.0) + t5 / (eps + beta3));

        const amrex::Real omega = omega1 + omega2 + omega3;

        const amrex::Real v_1 = amrex::Real(2.0) * sm2 - amrex::Real(7.0) * sm1 + amrex::Real(11.0) * s0;
        const amrex::Real v_2 = -sm1 + amrex::Real(5.0) * s0 + amrex::Real(2.0) * sp1;
        const amrex::Real v_3 = amrex::Real(2.0) * s0 + amrex::Real(5.0) * sp1 - sp2;

        return (omega1 * v_1 + omega2 * v_2 + omega3 * v_3) / (amrex::Real(6.0) * omega);
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real
    sedge1(const amrex::Real sm2,
           const amrex::Real sm1,
           const amrex::Real s0,
           const amrex::Real sp1,
           const amrex::Real sp2)
    {
        return sedge2(sp2,sp1,s0,sm1,sm2); // NOLINT(readability-suspicious-call-argument)
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::GpuTuple<amrex::Real,amrex::Real>
    sm_sp(const amrex::Real /*s0*/,
          const amrex::Real sedge1,
          const amrex::Real sedge2)
    {
        return amrex::makeTuple(sedge1, sedge2);
    }
};


template <typename Limiter>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SetXBCs ( const int i, const int j, const int k, const int n,
               amrex::Real &sm, amrex::Real &sp,
               amrex::Real &sedge1, amrex::Real &sedge2,
               const amrex::Array4<const amrex::Real> &s,
               const int bclo, const int bchi,
               const int domlo, const int domhi)
{
    using namespace amrex;

    if (bclo == BCType::ext_dir || bclo == BCType::hoextrap)
    {
        if (i == domlo)
        {
            sedge2 = -amrex::Real(0.2)*s(domlo-1,j,k,n) + amrex::Real(0.75)*s(domlo,j, k, n)
                     +amrex::Real(0.5)*s(domlo+1,j,k,n) - amrex::Real(0.05)*s(domlo+2,j,k,n);
            if constexpr (! std::is_same_v<nolimiter,Limiter>) {
                sedge2 = amrex::max(sedge2, amrex::min(s(domlo+1,j,k,n), s(domlo,j,k,n)));
                sedge2 = amrex::min(sedge2, amrex::max(s(domlo+1,j,k,n), s(domlo,j,k,n)));
            }

            sm = s(domlo-1,j,k,n);
            sp = sedge2;

        } else if (i == domlo+1) {

            sedge1 = -amrex::Real(0.2)*s(domlo-1,j,k,n) + amrex::Real(0.75)*s(domlo  ,j,k,n)
                   +  amrex::Real(0.5)*s(domlo+1,j,k,n) - amrex::Real(0.05)*s(domlo+2,j,k,n);
            if constexpr (! std::is_same_v<nolimiter,Limiter>) {
                sedge1 = amrex::max(sedge1, amrex::min(s(domlo+1,j,k,n), s(domlo,j,k,n)));
                sedge1 = amrex::min(sedge1, amrex::max(s(domlo+1,j,k,n), s(domlo,j,k,n)));
            }

            sp = sedge2;
            sm = sedge1;

            if constexpr (! std::is_same_v<nolimiter,Limiter>) {
                if ( (sp - s(domlo+1,j,k,n))*(s(domlo+1,j,k,n) - sm) <= amrex::Real(0.0))
                {
                    sp = s(domlo+1,j,k,n);
                    sm = s(domlo+1,j,k,n);
                } else if(amrex::Math::abs(sp - s(domlo+1,j,k,n)) >= amrex::Real(2.0)*amrex::Math::abs(sm - s(domlo+1,j,k,n))) {
                    sp = amrex::Real(3.0)*s(domlo+1,j,k,n) - amrex::Real(2.0)*sm;
                } else if(amrex::Math::abs(sm - s(domlo+1,j,k,n)) >= amrex::Real(2.0)*amrex::Math::abs(sp - s(domlo+1,j,k,n))) {
                    sm = amrex::Real(3.0)*s(domlo+1,j,k,n) - amrex::Real(2.0)*sp;
                }
            }
        }
    }

    if (bchi == BCType::ext_dir || bchi == BCType::hoextrap)
    {
        if (i == domhi)
        {
            sedge1 = -amrex::Real(0.2)*s(domhi+1,j,k,n) + amrex::Real(0.75)*s(domhi,j,k, n)
                     +amrex::Real(0.5)*s(domhi-1,j,k,n) - amrex::Real(0.05)*s(domhi-2,j,k,n);
            if constexpr (! std::is_same_v<nolimiter,Limiter>) {
                sedge1 = amrex::max(sedge1, amrex::min(s(domhi-1,j,k,n), s(domhi,j,k,n)));
                sedge1 = amrex::min(sedge1, amrex::max(s(domhi-1,j,k,n), s(domhi,j,k,n)));
            }

            sp = s(domhi+1,j,k,n);
            sm = sedge1;

        } else if (i == domhi-1) {

            sedge2 = -amrex::Real(0.2)*s(domhi+1,j,k,n) + amrex::Real(0.75)*s(domhi  ,j,k,n)
                     +amrex::Real(0.5)*s(domhi-1,j,k,n) - amrex::Real(0.05)*s(domhi-2,j,k,n);
            if constexpr (! std::is_same_v<nolimiter,Limiter>) {
                sedge2 = amrex::max(sedge2, amrex::min(s(domhi-1,j,k,n), s(domhi,j,k,n)));
                sedge2 = amrex::min(sedge2, amrex::max(s(domhi-1,j,k,n), s(domhi,j,k,n)));
            }

            sp = sedge2;
            sm = sedge1;

            if constexpr (! std::is_same_v<nolimiter,Limiter>) {
                if( (sp - s(domhi-1,j,k,n))*(s(domhi-1,j,k,n) - sm) <= amrex::Real(0.0))
                {
                    sp = s(domhi-1,j,k,n);
                    sm = s(domhi-1,j,k,n);
                } else if(amrex::Math::abs(sp - s(domhi-1,j,k,n)) >= 2.*amrex::Math::abs(sm - s(domhi-1,j,k,n))) {
                    sp = amrex::Real(3.0)*s(domhi-1,j,k,n) - amrex::Real(2.0)*sm;
                } else if(amrex::Math::abs(sm - s(domhi-1,j,k,n)) >= 2.*amrex::Math::abs(sp - s(domhi-1,j,k,n))) {
                    sm = amrex::Real(3.0)*s(domhi-1,j,k,n) - amrex::Real(2.0)*sp;
                }
            }
        }
    }
}

template <typename Limiter>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SetYBCs ( const int i, const int j, const int k, const int n,
               amrex::Real &sm, amrex::Real &sp,
               amrex::Real &sedge1, amrex::Real &sedge2,
               const amrex::Array4<const amrex::Real> &s,
               const int bclo, const int bchi,
               const int domlo, const int domhi)
{
    using namespace amrex;

    if (bclo == BCType::ext_dir || bclo == BCType::hoextrap)
    {
        if (j == domlo)
        {
            sedge2 = -amrex::Real(0.2)*s(i,domlo-1,k,n) + amrex::Real(0.75)*s(i,domlo  ,k,n)
                     +amrex::Real(0.5)*s(i,domlo+1,k,n) - amrex::Real(0.05)*s(i,domlo+2,k,n);
            if constexpr (! std::is_same_v<nolimiter,Limiter>) {
                sedge2 = amrex::max(sedge2, amrex::min(s(i,domlo+1,k,n), s(i,domlo,k,n)));
                sedge2 = amrex::min(sedge2, amrex::max(s(i,domlo+1,k,n), s(i,domlo,k,n)));
            }

            sm = s(i,domlo-1,k,n);
            sp = sedge2;

        } else if (j == domlo+1) {

            sedge1 = -amrex::Real(0.2)*s(i,domlo-1,k,n) + amrex::Real(0.75)*s(i,domlo  ,k,n)
                     +amrex::Real(0.5)*s(i,domlo+1,k,n) - amrex::Real(0.05)*s(i,domlo+2,k,n);
            if constexpr (! std::is_same_v<nolimiter,Limiter>) {
                sedge1 = amrex::max(sedge1, amrex::min(s(i,domlo+1,k,n), s(i,domlo,k,n)));
                sedge1 = amrex::min(sedge1, amrex::max(s(i,domlo+1,k,n), s(i,domlo,k,n)));
            }

            sp = sedge2;
            sm = sedge1;

            if constexpr (! std::is_same_v<nolimiter,Limiter>) {
                if ( (sp - s(i,domlo+1,k,n))*(s(i,domlo+1,k,n) - sm) <= amrex::Real(0.0))
                {
                    sp = s(i,domlo+1,k,n);
                    sm = s(i,domlo+1,k,n);
                } else if(amrex::Math::abs(sp - s(i,domlo+1,k,n)) >= amrex::Real(2.0)*amrex::Math::abs(sm - s(i,domlo+1,k,n))) {
                    sp = amrex::Real(3.0)*s(i,domlo+1,k,n) - amrex::Real(2.0)*sm;
                } else if(amrex::Math::abs(sm - s(i,domlo+1,k,n)) >= amrex::Real(2.0)*amrex::Math::abs(sp - s(i,domlo+1,k,n))) {
                    sm = amrex::Real(3.0)*s(i,domlo+1,k,n) - amrex::Real(2.0)*sp;
                }
            }
        }
    }

    if (bchi == BCType::ext_dir || bchi == BCType::hoextrap)
    {
        if (j == domhi)
        {
            sedge1 = -amrex::Real(0.2)*s(i,domhi+1,k,n) + amrex::Real(0.75)*s(i,domhi  ,k,n)
                     +amrex::Real(0.5)*s(i,domhi-1,k,n) - amrex::Real(0.05)*s(i,domhi-2,k,n);
            if constexpr (! std::is_same_v<nolimiter,Limiter>) {
                sedge1 = amrex::max(sedge1, amrex::min(s(i,domhi-1,k,n), s(i,domhi,k,n)));
                sedge1 = amrex::min(sedge1, amrex::max(s(i,domhi-1,k,n), s(i,domhi,k,n)));
            }

            sp = s(i,domhi+1,k, n);
            sm = sedge1;

        } else if (j == domhi-1) {

            sedge2 = -amrex::Real(0.2)*s(i,domhi+1,k,n) + amrex::Real(0.75)*s(i,domhi  ,k,n)
                     +amrex::Real(0.5)*s(i,domhi-1,k,n) - amrex::Real(0.05)*s(i,domhi-2,k,n);
            if constexpr (! std::is_same_v<nolimiter,Limiter>) {
                sedge2 = amrex::max(sedge2, amrex::min(s(i,domhi-1,k,n), s(i,domhi,k,n)));
                sedge2 = amrex::min(sedge2, amrex::max(s(i,domhi-1,k,n), s(i,domhi,k,n)));
            }

            sp = sedge2;
            sm = sedge1;

            if constexpr (! std::is_same_v<nolimiter,Limiter>) {
                if( (sp - s(i,domhi-1,k,n))*(s(i,domhi-1,k,n) - sm) <= amrex::Real(0.0)){
                    sp = s(i,domhi-1,k,n);
                    sm = s(i,domhi-1,k,n);
                } else if(amrex::Math::abs(sp - s(i,domhi-1,k,n)) >= 2.*amrex::Math::abs(sm - s(i,domhi-1,k,n))) {
                    sp = amrex::Real(3.0)*s(i,domhi-1,k,n) - amrex::Real(2.0)*sm;
                } else if(amrex::Math::abs(sm - s(i,domhi-1,k,n)) >= 2.*amrex::Math::abs(sp - s(i,domhi-1,k,n))) {
                    sm = amrex::Real(3.0)*s(i,domhi-1,k,n) - amrex::Real(2.0)*sp;
                }
            }
        }
    }
}

#if (AMREX_SPACEDIM==3)
template <typename Limiter>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SetZBCs ( const int i, const int j, const int k, const int n,
               amrex::Real &sm, amrex::Real &sp,
               amrex::Real &sedge1, amrex::Real &sedge2,
               const amrex::Array4<const amrex::Real> &s,
               const int bclo, const int bchi,
               const int domlo, const int domhi)
{
    using namespace amrex;

    if (bclo == BCType::ext_dir || bclo == BCType::hoextrap)
    {
        if (k == domlo)
        {
            sedge2 = -amrex::Real(0.2)*s(i,j,domlo-1,n) + amrex::Real(0.75)*s(i,j,domlo  ,n)
                     +amrex::Real(0.5)*s(i,j,domlo+1,n) - amrex::Real(0.05)*s(i,j,domlo+2,n);
            if constexpr (! std::is_same_v<nolimiter,Limiter>) {
                sedge2 = amrex::max(sedge2, amrex::min(s(i,j,domlo+1,n), s(i,j,domlo,n)));
                sedge2 = amrex::min(sedge2, amrex::max(s(i,j,domlo+1,n), s(i,j,domlo,n)));
            }

            sm = s(i,j,domlo-1,n);
            sp = sedge2;

        } else if (k == domlo+1) {

            sedge1 = -amrex::Real(0.2)*s(i,j,domlo-1,n) + amrex::Real(0.75)*s(i,j,domlo  ,n)
                     +amrex::Real(0.5)*s(i,j,domlo+1,n) - amrex::Real(0.05)*s(i,j,domlo+2,n);
            if constexpr (! std::is_same_v<nolimiter,Limiter>) {
                sedge1 = amrex::max(sedge1, amrex::min(s(i,j,domlo+1,n), s(i,j,domlo,n)));
                sedge1 = amrex::min(sedge1, amrex::max(s(i,j,domlo+1,n), s(i,j,domlo,n)));
            }

            sp = sedge2;
            sm = sedge1;

            if constexpr (! std::is_same_v<nolimiter,Limiter>) {
                if ( (sp - s(i,j,domlo+1,n))*(s(i,j,domlo+1,n) - sm) <= 0. )
                {
                    sp = s(i,j,domlo+1,n);
                    sm = s(i,j,domlo+1,n);
                } else if(amrex::Math::abs(sp - s(i,j,domlo+1,n)) >= 2.*amrex::Math::abs(sm - s(i,j,domlo+1,n))) {
                    sp = amrex::Real(3.0)*s(i,j,domlo+1,n) - amrex::Real(2.0)*sm;
                } else if(amrex::Math::abs(sm - s(i,j,domlo+1,n)) >= 2.*amrex::Math::abs(sp - s(i,j,domlo+1,n))) {
                    sm = amrex::Real(3.0)*s(i,j,domlo+1,n) - amrex::Real(2.0)*sp;
                }
            }
        }
    }

    if (bchi == BCType::ext_dir || bchi == BCType::hoextrap)
    {
        if (k == domhi)
        {
            sedge1 =  -amrex::Real(0.2)*s(i,j,domhi+1,n) + amrex::Real(0.75)*s(i,j,domhi  ,n)
                      +amrex::Real(0.5)*s(i,j,domhi-1,n) - amrex::Real(0.05)*s(i,j,domhi-2,n);
            if constexpr (! std::is_same_v<nolimiter,Limiter>) {
                sedge1 = amrex::max(sedge1, amrex::min(s(i,j,domhi-1,n), s(i,j,domhi,n)));
                sedge1 = amrex::min(sedge1, amrex::max(s(i,j,domhi-1,n), s(i,j,domhi,n)));
            }

            sp = s(i,j,domhi+1,n);
            sm = sedge1;

        } else if (k == domhi-1) {

            sedge2 = -amrex::Real(0.2)*s(i,j,domhi+1,n) + amrex::Real(0.75)*s(i,j,domhi  ,n)
                     +amrex::Real(0.5)*s(i,j,domhi-1,n) - amrex::Real(0.05)*s(i,j,domhi-2,n);
            if constexpr (! std::is_same_v<nolimiter,Limiter>) {
                sedge2 = amrex::max(sedge2, amrex::min(s(i,j,domhi-1,n), s(i,j,domhi,n)));
                sedge2 = amrex::min(sedge2, amrex::max(s(i,j,domhi-1,n), s(i,j,domhi,n)));
            }

            sp = sedge2;
            sm = sedge1;

            if constexpr (! std::is_same_v<nolimiter,Limiter>) {
                if ( (sp - s(i,j,domhi-1,n))*(s(i,j,domhi-1,n) - sm) <= 0. )
                {
                    sp = s(i,j,domhi-1,n);
                    sm = s(i,j,domhi-1,n);
                } else if(amrex::Math::abs(sp - s(i,j,domhi-1,n)) >= 2.*amrex::Math::abs(sm - s(i,j,domhi-1,n))) {
                    sp = amrex::Real(3.0)*s(i,j,domhi-1,n) - amrex::Real(2.0)*sm;
                } else if(amrex::Math::abs(sm - s(i,j,domhi-1,n)) >= 2.*amrex::Math::abs(sp - s(i,j,domhi-1,n))) {
                    sm = amrex::Real(3.0)*s(i,j,domhi-1,n) - amrex::Real(2.0)*sp;
                }
            }

        }
    }
}
#endif


// Right now only ppm type 1 is supported on GPU
// This version is called before the MAC projection, when we use the cell-centered velocity
//      for upwinding
template <typename Limiter>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void PredictVelOnXFace ( const int i, const int j, const int k, const int n,
                         const amrex::Real dtdx, const amrex::Real v_ad,
                         const amrex::Array4<const amrex::Real> &S,
                         const amrex::Array4<amrex::Real> &Im,
                         const amrex::Array4<amrex::Real> &Ip,
                         const amrex::BCRec bc, const int domlo, const int domhi,
                         const Limiter& /*limiter*/)
{
    constexpr auto half{amrex::Real(0.5)};
    constexpr auto two3rds{amrex::Real(2.0/3.0)};

    amrex::Real sm2 = S(i-2,j,k,n);
    amrex::Real sm1 = S(i-1,j,k,n);
    amrex::Real s0  = S(i  ,j,k,n);
    amrex::Real sp1 = S(i+1,j,k,n);
    amrex::Real sp2 = S(i+2,j,k,n);

    amrex::Real sedge1 = Limiter::sedge1(sm2,sm1,s0,sp1,sp2);
    amrex::Real sedge2 = Limiter::sedge2(sm2,sm1,s0,sp1,sp2);

    amrex::Real sm, sp;
    amrex::Tie(sm,sp) = Limiter::sm_sp(s0,sedge1,sedge2);

    SetXBCs<Limiter>(i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(0), bc.hi(0), domlo, domhi);

    amrex::Real s6 = amrex::Real(6.0)*s0 - amrex::Real(3.0)*(sm + sp);

    amrex::Real sigma = amrex::Math::abs(v_ad)*dtdx;

    if (v_ad > small_vel)
    {
        Ip(i,j,k,n) = sp - (half*sigma)*((sp-sm) - (amrex::Real(1.0) - two3rds*sigma)*s6);
        Im(i,j,k,n) = S(i,j,k,n);
    }
    else if (v_ad < -small_vel)
    {
        Ip(i,j,k,n) = S(i,j,k,n);
        Im(i,j,k,n) = sm + (half*sigma)*((sp-sm) + (amrex::Real(1.0) - two3rds*sigma)*s6);
    } else
    {
        Ip(i,j,k,n) = S(i,j,k,n);
        Im(i,j,k,n) = S(i,j,k,n);
    }
}

template <typename Limiter>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void PredictVelOnYFace ( const int i, const int j, const int k, const int n,
                         const amrex::Real dtdy, const amrex::Real v_ad,
                         const amrex::Array4<const amrex::Real> &S,
                         const amrex::Array4<amrex::Real> &Im,
                         const amrex::Array4<amrex::Real> &Ip,
                         const amrex::BCRec bc, const int domlo, const int domhi,
                         const Limiter& /*limiter*/)
{
    constexpr auto half{amrex::Real(0.5)};
    constexpr auto two3rds{amrex::Real(2.0/3.0)};

    amrex::Real sm2 = S(i,j-2,k,n);
    amrex::Real sm1 = S(i,j-1,k,n);
    amrex::Real s0  = S(i,j  ,k,n);
    amrex::Real sp1 = S(i,j+1,k,n);
    amrex::Real sp2 = S(i,j+2,k,n);

    amrex::Real sedge1 = Limiter::sedge1(sm2,sm1,s0,sp1,sp2);
    amrex::Real sedge2 = Limiter::sedge2(sm2,sm1,s0,sp1,sp2);

    amrex::Real sm, sp;
    amrex::Tie(sm,sp) = Limiter::sm_sp(s0,sedge1,sedge2);

    SetYBCs<Limiter>(i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(1), bc.hi(1), domlo, domhi);

    amrex::Real s6 = amrex::Real(6.0)*s0 - amrex::Real(3.0)*(sm + sp);

    amrex::Real sigma = amrex::Math::abs(v_ad)*dtdy;

    if (v_ad > small_vel)
    {
        Ip(i,j,k,n) = sp - (half*sigma)*((sp-sm) - (amrex::Real(1.0) - two3rds*sigma)*s6);
        Im(i,j,k,n) = S(i,j,k,n);
    }
    else if (v_ad < -small_vel)
    {
        Ip(i,j,k,n) = S(i,j,k,n);
        Im(i,j,k,n) = sm + (half*sigma)*((sp-sm) + (amrex::Real(1.0) - two3rds*sigma)*s6);
    } else
    {
        Ip(i,j,k,n) = S(i,j,k,n);
        Im(i,j,k,n) = S(i,j,k,n);
    }
}

#if (AMREX_SPACEDIM==3)
template <typename Limiter>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void PredictVelOnZFace ( const int i, const int j, const int k, const int n,
                         const amrex::Real dtdz, const amrex::Real v_ad,
                         const amrex::Array4<const amrex::Real> &S,
                         const amrex::Array4<amrex::Real> &Im,
                         const amrex::Array4<amrex::Real> &Ip,
                         const amrex::BCRec bc, const int domlo, const int domhi,
                         const Limiter& /*limiter*/)
{
    constexpr auto half{amrex::Real(0.5)};
    constexpr auto two3rds{amrex::Real(2.0/3.0)};

    amrex::Real sm2 = S(i,j,k-2,n);
    amrex::Real sm1 = S(i,j,k-1,n);
    amrex::Real s0  = S(i,j,k  ,n);
    amrex::Real sp1 = S(i,j,k+1,n);
    amrex::Real sp2 = S(i,j,k+2,n);

    amrex::Real sedge1 = Limiter::sedge1(sm2,sm1,s0,sp1,sp2);
    amrex::Real sedge2 = Limiter::sedge2(sm2,sm1,s0,sp1,sp2);

    amrex::Real sm, sp;
    amrex::Tie(sm,sp) = Limiter::sm_sp(s0,sedge1,sedge2);

    SetZBCs<Limiter>(i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(2), bc.hi(2), domlo, domhi);

    amrex::Real s6 = amrex::Real(6.0)*s0 - amrex::Real(3.0)*(sm + sp);

    amrex::Real sigma = amrex::Math::abs(v_ad)*dtdz;

    if (v_ad > small_vel)
    {
        Ip(i,j,k,n) = sp - (half*sigma)*((sp-sm) - (amrex::Real(1.0) - two3rds*sigma)*s6);
        Im(i,j,k,n) = S(i,j,k,n);
    }
    else if (v_ad < -small_vel)
    {
        Ip(i,j,k,n) = S(i,j,k,n);
        Im(i,j,k,n) = sm + (half*sigma)*((sp-sm) + (amrex::Real(1.0) - two3rds*sigma)*s6);
    } else
    {
        Ip(i,j,k,n) = S(i,j,k,n);
        Im(i,j,k,n) = S(i,j,k,n);
    }
}

#endif

// Right now only ppm type 1 is supported on GPU
// This version is called after the MAC projection, when we use the MAC-projected velocity
//      for upwinding
template <typename Limiter>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void PredictStateOnXFace ( const int i, const int j, const int k, const int n,
                           const amrex::Real dt, const amrex::Real dx,
                           amrex::Real& Im, amrex::Real& Ip,
                           const amrex::Array4<const amrex::Real> &S,
                           const amrex::Array4<const amrex::Real> &vel_edge,
                           const amrex::BCRec bc,
                           const int domlo, const int domhi,
                           const Limiter& /*limiter*/)
{
    constexpr auto half{amrex::Real(0.5)};
    constexpr auto two3rds{amrex::Real(2.0/3.0)};

    amrex::Real sm2 = S(i-2,j,k,n);
    amrex::Real sm1 = S(i-1,j,k,n);
    amrex::Real s0  = S(i  ,j,k,n);
    amrex::Real sp1 = S(i+1,j,k,n);
    amrex::Real sp2 = S(i+2,j,k,n);

    amrex::Real sedge1 = Limiter::sedge1(sm2,sm1,s0,sp1,sp2);
    amrex::Real sedge2 = Limiter::sedge2(sm2,sm1,s0,sp1,sp2);

    amrex::Real sm, sp;
    amrex::Tie(sm,sp) = Limiter::sm_sp(s0,sedge1,sedge2);

    SetXBCs<Limiter>(i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(0), bc.hi(0), domlo, domhi);

    amrex::Real s6 = 6.0*s0 - 3.0*(sm + sp);

    amrex::Real sigmap = amrex::Math::abs(vel_edge(i+1,j,k))*dt/dx;
    amrex::Real sigmam = amrex::Math::abs(vel_edge(i  ,j,k))*dt/dx;

    if (vel_edge(i+1,j,k) > small_vel)
        Ip = sp - (half*sigmap)*((sp - sm) - (amrex::Real(1.0) -two3rds*sigmap)*s6);
    else
        Ip = S(i,j,k,n);

    if(vel_edge(i,j,k) < -small_vel)
        Im = sm + (half*sigmam)*((sp-sm) + (amrex::Real(1.0) - two3rds*sigmam)*s6);
    else
        Im = S(i,j,k,n);
}

template <typename Limiter>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void PredictStateOnYFace ( const int i, const int j, const int k, const int n,
                           const amrex::Real dt, const amrex::Real dx,
                           amrex::Real& Im, amrex::Real& Ip,
                           const amrex::Array4<const amrex::Real> &S,
                           const amrex::Array4<const amrex::Real> &vel_edge,
                           const amrex::BCRec bc,
                           const int domlo, const int domhi,
                           const Limiter& /*limiter*/)
{
    constexpr auto half{amrex::Real(0.5)};
    constexpr auto two3rds{amrex::Real(2.0/3.0)};

    amrex::Real sm2 = S(i,j-2,k,n);
    amrex::Real sm1 = S(i,j-1,k,n);
    amrex::Real s0  = S(i,j  ,k,n);
    amrex::Real sp1 = S(i,j+1,k,n);
    amrex::Real sp2 = S(i,j+2,k,n);

    amrex::Real sedge1 = Limiter::sedge1(sm2,sm1,s0,sp1,sp2);
    amrex::Real sedge2 = Limiter::sedge2(sm2,sm1,s0,sp1,sp2);

    amrex::Real sm, sp;
    amrex::Tie(sm,sp) = Limiter::sm_sp(s0,sedge1,sedge2);

    SetYBCs<Limiter>(i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(1), bc.hi(1), domlo, domhi);

    amrex::Real s6 = 6.0*s0- 3.0*(sm + sp);

    amrex::Real sigmap = amrex::Math::abs(vel_edge(i,j+1,k))*dt/dx;
    amrex::Real sigmam = amrex::Math::abs(vel_edge(i,j  ,k))*dt/dx;

    if (vel_edge(i,j+1,k) > small_vel)
        Ip = sp - (half*sigmap)*((sp - sm) - (amrex::Real(1.0) -two3rds*sigmap)*s6);
    else
        Ip = S(i,j,k,n);

    if (vel_edge(i,j,k) < -small_vel)
        Im = sm + (half*sigmam)*((sp-sm) + (amrex::Real(1.0) - two3rds*sigmam)*s6);
    else
        Im = S(i,j,k,n);
}


#if (AMREX_SPACEDIM==3)
template <typename Limiter>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void PredictStateOnZFace ( const int i, const int j, const int k, const int n,
                           const amrex::Real dt, const amrex::Real dx,
                           amrex::Real& Im, amrex::Real& Ip,
                           const amrex::Array4<const amrex::Real> &S,
                           const amrex::Array4<const amrex::Real> &vel_edge,
                           const amrex::BCRec bc,
                           const int domlo, const int domhi,
                           const Limiter& /*limiter*/)
{

    constexpr auto half{amrex::Real(0.5)};
    constexpr auto two3rds{amrex::Real(2.0/3.0)};

    amrex::Real sm2 = S(i,j,k-2,n);
    amrex::Real sm1 = S(i,j,k-1,n);
    amrex::Real s0  = S(i,j,k  ,n);
    amrex::Real sp1 = S(i,j,k+1,n);
    amrex::Real sp2 = S(i,j,k+2,n);

    amrex::Real sedge1 = Limiter::sedge1(sm2,sm1,s0,sp1,sp2);
    amrex::Real sedge2 = Limiter::sedge2(sm2,sm1,s0,sp1,sp2);

    amrex::Real sm, sp;
    amrex::Tie(sm,sp) = Limiter::sm_sp(s0,sedge1,sedge2);

    SetZBCs<Limiter>(i, j, k, n, sm, sp, sedge1, sedge2, S, bc.lo(2), bc.hi(2), domlo, domhi);

    amrex::Real s6 = 6.0*s0- 3.0*(sm + sp);
    amrex::Real sigmap = amrex::Math::abs(vel_edge(i,j,k+1))*dt/dx;
    amrex::Real sigmam = amrex::Math::abs(vel_edge(i,j,k  ))*dt/dx;

    if(vel_edge(i,j,k+1) > small_vel)
        Ip = sp - (half*sigmap)*((sp-sm) - (amrex::Real(1.0) -two3rds*sigmap)*s6);
    else
        Ip = S(i,j,k,n);

    if(vel_edge(i,j,k) < -small_vel)
        Im = sm + (half*sigmam)*((sp-sm) + (amrex::Real(1.0) - two3rds*sigmam)*s6);
    else
        Im = S(i,j,k,n);
}
#endif

template <typename Limiter>
void PredictVelOnFaces (amrex::Box const& bx,
                        AMREX_D_DECL(amrex::Array4<amrex::Real> const& Imx,
                                     amrex::Array4<amrex::Real> const& Imy,
                                     amrex::Array4<amrex::Real> const& Imz),
                        AMREX_D_DECL(amrex::Array4<amrex::Real> const& Ipx,
                                     amrex::Array4<amrex::Real> const& Ipy,
                                     amrex::Array4<amrex::Real> const& Ipz),
                        amrex::Array4<amrex::Real const> const& q,
                        amrex::Array4<amrex::Real const> const& vel,
                        amrex::Geometry geom,
                        amrex::Real dt,
                        amrex::BCRec const* pbc,
                        const Limiter& limiter)
{
    const amrex::Box& domain = geom.Domain();
    const amrex::Dim3 dlo = amrex::lbound(domain);
    const amrex::Dim3 dhi = amrex::ubound(domain);

    const auto dx = geom.CellSizeArray();
    AMREX_D_TERM( amrex::Real l_dtdx = dt / dx[0];,
                  amrex::Real l_dtdy = dt / dx[1];,
                  amrex::Real l_dtdz = dt / dx[2];);

    amrex::ParallelFor(bx, AMREX_SPACEDIM,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        PredictVelOnXFace(i,j,k,n,l_dtdx,vel(i,j,k,0),q,Imx,Ipx,pbc[n],dlo.x,dhi.x,limiter);
        PredictVelOnYFace(i,j,k,n,l_dtdy,vel(i,j,k,1),q,Imy,Ipy,pbc[n],dlo.y,dhi.y,limiter);
#if (AMREX_SPACEDIM==3)
        PredictVelOnZFace(i,j,k,n,l_dtdz,vel(i,j,k,2),q,Imz,Ipz,pbc[n],dlo.z,dhi.z,limiter);
#endif
    });
}

template <typename Limiter>
void PredictStateOnFaces (amrex::Box const& bx,
                        AMREX_D_DECL(amrex::Array4<amrex::Real> const& Imx,
                                     amrex::Array4<amrex::Real> const& Imy,
                                     amrex::Array4<amrex::Real> const& Imz),
                        AMREX_D_DECL(amrex::Array4<amrex::Real> const& Ipx,
                                     amrex::Array4<amrex::Real> const& Ipy,
                                     amrex::Array4<amrex::Real> const& Ipz),
                        AMREX_D_DECL(amrex::Array4<amrex::Real const> const& umac,
                                     amrex::Array4<amrex::Real const> const& vmac,
                                     amrex::Array4<amrex::Real const> const& wmac),
                        amrex::Array4<amrex::Real const> const& q,
                        amrex::Geometry geom,
                        amrex::Real l_dt,
                        amrex::BCRec const* pbc,
                        const int ncomp,
                        const Limiter& limiter)
{
    const amrex::Box& domain = geom.Domain();
    const amrex::Dim3 dlo = amrex::lbound(domain);
    const amrex::Dim3 dhi = amrex::ubound(domain);

    AMREX_D_TERM (const auto dx = geom.CellSize(0);,
                  const auto dy = geom.CellSize(1);,
                  const auto dz = geom.CellSize(2););

    amrex::ParallelFor(bx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        PPM::PredictStateOnXFace(i, j, k, n, l_dt, dx, Imx(i,j,k,n), Ipx(i,j,k,n),
                               q, umac, pbc[n], dlo.x, dhi.x,limiter);
        PPM::PredictStateOnYFace(i, j, k, n, l_dt, dy, Imy(i,j,k,n), Ipy(i,j,k,n),
                               q, vmac, pbc[n], dlo.y, dhi.y,limiter);
#if (AMREX_SPACEDIM==3)
        PPM::PredictStateOnZFace(i, j, k, n, l_dt, dz, Imz(i,j,k,n), Ipz(i,j,k,n),
                               q, wmac, pbc[n], dlo.z, dhi.z,limiter);
#endif
    });
}


}

#endif
/** @} */
