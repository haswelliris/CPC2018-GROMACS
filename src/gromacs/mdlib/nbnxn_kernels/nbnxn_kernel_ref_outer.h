/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#define UNROLLI    NBNXN_CPU_CLUSTER_I_SIZE
#define UNROLLJ    NBNXN_CPU_CLUSTER_I_SIZE

/* We could use para_nbat->xstride and para_nbat->fstride, but macros might be faster */
#define X_STRIDE   3
#define F_STRIDE   3
/* Local i-atom buffer strides */
#define XI_STRIDE  3
#define FI_STRIDE  3


/* All functionality defines are set here, except for:
 * CALC_ENERGIES, ENERGY_GROUPS which are defined before.
 * CHECK_EXCLS, which is set just before including the inner loop contents.
 */

/* We always calculate shift forces, because it's cheap anyhow */
#define CALC_SHIFTFORCES

#ifdef CALC_COUL_RF
#define NBK_FUNC_NAME2(ljt, feg) nbnxn_kernel ## _ElecRF ## ljt ## feg ## _ref
#endif
#ifdef CALC_COUL_TAB
#ifndef VDW_CUTOFF_CHECK
#define NBK_FUNC_NAME2(ljt, feg) nbnxn_kernel ## _ElecQSTab ## ljt ## feg ## _ref
#else
#define NBK_FUNC_NAME2(ljt, feg) nbnxn_kernel ## _ElecQSTabTwinCut ## ljt ## feg ## _ref
#endif
#endif

#if defined LJ_CUT && !defined LJ_EWALD
#define NBK_FUNC_NAME(feg) NBK_FUNC_NAME2(_VdwLJ, feg)
#elif defined LJ_FORCE_SWITCH
#define NBK_FUNC_NAME(feg) NBK_FUNC_NAME2(_VdwLJFsw, feg)
#elif defined LJ_POT_SWITCH
#define NBK_FUNC_NAME(feg) NBK_FUNC_NAME2(_VdwLJPsw, feg)
#elif defined LJ_EWALD
#ifdef LJ_EWALD_COMB_GEOM
#define NBK_FUNC_NAME(feg) NBK_FUNC_NAME2(_VdwLJEwCombGeom, feg)
#else
#define NBK_FUNC_NAME(feg) NBK_FUNC_NAME2(_VdwLJEwCombLB, feg)
#endif
#else
#error "No VdW type defined"
#endif

static void
#ifndef CALC_ENERGIES
NBK_FUNC_NAME(_F)
#else
#ifndef ENERGY_GROUPS
NBK_FUNC_NAME(_VF)
#else
NBK_FUNC_NAME(_VgrpF)
#endif
#endif
#undef NBK_FUNC_NAME
#undef NBK_FUNC_NAME2
(const nbnxn_pairlist_t     *para_nbl,
 const nbnxn_atomdata_t     *para_nbat,
 const interaction_const_t  *para_ic,
 rvec                       *para_shift_vec,
 real                       *para_f
#ifdef CALC_SHIFTFORCES
 ,
 real                       *para_fshift
#endif
#ifdef CALC_ENERGIES
 ,
 real                       *para_Vvdw,
 real                       *para_Vc
#endif
)
{
    // =========== DEF DATA =============
    const nbnxn_ci_t   *nbln;
    const nbnxn_cj_t   *l_cj;
    const int          *type;
    const real         *q;
    const real         *para_shiftvec;
    const real         *x;
    const real         *nbfp;
    real                rcut2;
#ifdef VDW_CUTOFF_CHECK
    real                rvdw2;
#endif
    int                 ntype2;
    real                facel;
    real               *nbfp_i;
    int                 n, ci, ci_sh;
    int                 ish, ishf;
    gmx_bool            do_LJ, half_LJ, do_coul, do_self;
    int                 cjind0, cjind1, cjind;
    int                 ip, jp;

    real                xi[UNROLLI*XI_STRIDE];
    real                fi[UNROLLI*FI_STRIDE];
    real                qi[UNROLLI];

#ifdef CALC_ENERGIES
#ifndef ENERGY_GROUPS

    real       Vvdw_ci, Vc_ci;
#else
    int        egp_mask;
    int        egp_sh_i[UNROLLI];
#endif
#endif
#ifdef LJ_POT_SWITCH
    real       swV3, swV4, swV5;
    real       swF2, swF3, swF4;
#endif
#ifdef LJ_EWALD
    real        lje_coeff2, lje_coeff6_6, lje_vc;
    const real *ljc;
#endif

#ifdef CALC_COUL_RF
    real       k_rf2;
#ifdef CALC_ENERGIES
    real       k_rf, c_rf;
#endif
#endif
#ifdef CALC_COUL_TAB
    real       tabscale;
#ifdef CALC_ENERGIES
    real       halfsp;
#endif
#ifndef GMX_DOUBLE
    const real *tab_coul_FDV0;
#else
    const real *tab_coul_F;
    const real *tab_coul_V;
#endif
#endif
    // =========== DEF DATA =============

    // =========== INIT DATA =============
#ifdef LJ_POT_SWITCH
    swV3 = para_ic->vdw_switch.c3;
    swV4 = para_ic->vdw_switch.c4;
    swV5 = para_ic->vdw_switch.c5;
    swF2 = 3*para_ic->vdw_switch.c3;
    swF3 = 4*para_ic->vdw_switch.c4;
    swF4 = 5*para_ic->vdw_switch.c5;
#endif

#ifdef LJ_EWALD
    lje_coeff2   = para_ic->ewaldcoeff_lj*para_ic->ewaldcoeff_lj;
    lje_coeff6_6 = lje_coeff2*lje_coeff2*lje_coeff2/6.0;
    lje_vc       = para_ic->sh_lj_ewald;

    ljc          = para_nbat->nbfp_comb;
#endif

#ifdef CALC_COUL_RF
    k_rf2 = 2*para_ic->k_rf;
#ifdef CALC_ENERGIES
    k_rf = para_ic->k_rf;
    c_rf = para_ic->c_rf;
#endif
#endif
#ifdef CALC_COUL_TAB
    tabscale = para_ic->tabq_scale;
#ifdef CALC_ENERGIES
    halfsp = 0.5/para_ic->tabq_scale;
#endif

#ifndef GMX_DOUBLE
    tab_coul_FDV0 = para_ic->tabq_coul_FDV0;
#else
    tab_coul_F    = para_ic->tabq_coul_F;
    tab_coul_V    = para_ic->tabq_coul_V;
#endif
#endif

#ifdef ENERGY_GROUPS
    egp_mask = (1<<para_nbat->neg_2log) - 1;
#endif


    rcut2               = para_ic->rcoulomb*para_ic->rcoulomb;
#ifdef VDW_CUTOFF_CHECK
    rvdw2               = para_ic->rvdw*para_ic->rvdw;
#endif

    ntype2              = para_nbat->ntype*2;
    nbfp                = para_nbat->nbfp;
    q                   = para_nbat->q;
    type                = para_nbat->type;
    facel               = para_ic->epsfac;
    para_shiftvec            = para_shift_vec[0];
    x                   = para_nbat->x;

    l_cj = para_nbl->cj;
    // =========== INIT DATA =============

    /* =========== DATA'S STAT =========== //
    int ci_SZ = para_nbl->nci;
    int cj_SZ = para_nbl->ncj;
    int i, j, k;

    for(k = 0; k < 64; ++k)
    {
        if(k == host_param.host_rank)
        {
            printf("\n\n=====RANK %d=======\n", host_param.host_rank);
            printf("ci_SZ =%d, cj_SZ =%d\n\n", ci_SZ, cj_SZ);
            for(i = 0; i < ci_SZ; ++i)
            {
                nbln = &para_nbl->ci[i];
                printf("ci[%d]'s info:\n", i);
                printf("ci =%d, cjind0 =%d, cjind1 =%d\n", nbln->ci, nbln->cj_ind_start, nbln->cj_ind_end);
                printf("cj list = [");
                for(j = nbln->cj_ind_start; j < nbln->cj_ind_end; ++j)
                {
                    printf("%d, ", para_nbl->cj[j].cj);
                }
                printf("]\n\n");
            }
            printf("\n==================\n\n");
        }
        else
        {
            sleep(1);
        }
    }


    // =========== DATA'S STAT =========== */

    for (n = 0; n < para_nbl->nci; n++)
    {
        int i, d;

        nbln = &para_nbl->ci[n];

        ish              = (nbln->shift & NBNXN_CI_SHIFT);
        /* x, para_f and para_fshift are assumed to be stored with stride 3 */
        ishf             = ish*DIM;
        cjind0           = nbln->cj_ind_start;
        cjind1           = nbln->cj_ind_end;
        /* Currently only works super-cells equal to sub-cells */
        ci               = nbln->ci;
        ci_sh            = (ish == CENTRAL ? ci : -1);

        /* We have 5 LJ/C combinations, but use only three inner loops,
         * as the other combinations are unlikely and/or not much faster:
         * inner half-LJ + C for half-LJ + C / no-LJ + C
         * inner LJ + C      for full-LJ + C
         * inner LJ          for full-LJ + no-C / half-LJ + no-C
         */
        do_LJ   = (nbln->shift & NBNXN_CI_DO_LJ(0));
        do_coul = (nbln->shift & NBNXN_CI_DO_COUL(0));
        half_LJ = ((nbln->shift & NBNXN_CI_HALF_LJ(0)) || !do_LJ) && do_coul;
#ifdef LJ_EWALD
        do_self = TRUE;
#else
        do_self = do_coul;
#endif

#ifdef CALC_ENERGIES
#ifndef ENERGY_GROUPS
        Vvdw_ci = 0;
        Vc_ci   = 0;
#else
        for (i = 0; i < UNROLLI; i++)
        {
            egp_sh_i[i] = ((para_nbat->energrp[ci]>>(i*para_nbat->neg_2log)) & egp_mask)*para_nbat->nenergrp;
        }
#endif
#endif

        for (i = 0; i < UNROLLI; i++)
        {
            for (d = 0; d < DIM; d++)
            {
                xi[i*XI_STRIDE+d] = x[(ci*UNROLLI+i)*X_STRIDE+d] + para_shiftvec[ishf+d];
                fi[i*FI_STRIDE+d] = 0;
            }

            qi[i] = facel*q[ci*UNROLLI+i];
        }

#ifdef CALC_ENERGIES
        if (do_self)
        {
            real Vc_sub_self;

#ifdef CALC_COUL_RF
            Vc_sub_self = 0.5*c_rf;
#endif
#ifdef CALC_COUL_TAB
#ifdef GMX_DOUBLE
            Vc_sub_self = 0.5*tab_coul_V[0];
#else
            Vc_sub_self = 0.5*tab_coul_FDV0[2];
#endif
#endif

            if (l_cj[nbln->cj_ind_start].cj == ci_sh)
            {
                for (i = 0; i < UNROLLI; i++)
                {
                    int egp_ind;
#ifdef ENERGY_GROUPS
                    egp_ind = egp_sh_i[i] + ((para_nbat->energrp[ci]>>(i*para_nbat->neg_2log)) & egp_mask);
#else
                    egp_ind = 0;
#endif
                    /* Coulomb self interaction */
                    para_Vc[egp_ind]   -= qi[i]*q[ci*UNROLLI+i]*Vc_sub_self;

#ifdef LJ_EWALD
                    /* LJ Ewald self interaction */
                    para_Vvdw[egp_ind] += 0.5*para_nbat->nbfp[para_nbat->type[ci*UNROLLI+i]*(para_nbat->ntype + 1)*2]/6*lje_coeff6_6;
#endif
                }
            }
        }
#endif  /* CALC_ENERGIES */

        cjind = cjind0;
        while (cjind < cjind1 && para_nbl->cj[cjind].excl != 0xffff)
        {
#define CHECK_EXCLS
            if (half_LJ)
            {
#define CALC_COULOMB
#define HALF_LJ
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_inner.h"
#undef HALF_LJ
#undef CALC_COULOMB
            }
            else if (do_coul)
            {
#define CALC_COULOMB
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_inner.h"
#undef CALC_COULOMB
            }
            else
            {
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_inner.h"
            }
#undef CHECK_EXCLS
            cjind++;
        }

        for (; (cjind < cjind1); cjind++)
        {
            if (half_LJ)
            {
#define CALC_COULOMB
#define HALF_LJ
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_inner.h"
#undef HALF_LJ
#undef CALC_COULOMB
            }
            else if (do_coul)
            {
#define CALC_COULOMB
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_inner.h"
#undef CALC_COULOMB
            }
            else
            {
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_inner.h"
            }
        }

        /* Add accumulated i-forces to the force array */
        for (i = 0; i < UNROLLI; i++)
        {
            for (d = 0; d < DIM; d++)
            {
                para_f[(ci*UNROLLI+i)*F_STRIDE+d] += fi[i*FI_STRIDE+d];
            }
        }
#ifdef CALC_SHIFTFORCES
        if (para_fshift != NULL)
        {
            /* Add i forces to shifted force list */
            for (i = 0; i < UNROLLI; i++)
            {
                for (d = 0; d < DIM; d++)
                {
                    para_fshift[ishf+d] += fi[i*FI_STRIDE+d];
                }
            }
        }
#endif

#ifdef CALC_ENERGIES
#ifndef ENERGY_GROUPS
        *para_Vvdw += Vvdw_ci;
        *para_Vc   += Vc_ci;
#endif
#endif
    }
}

#undef CALC_SHIFTFORCES

#undef X_STRIDE
#undef F_STRIDE
#undef XI_STRIDE
#undef FI_STRIDE

#undef UNROLLI
#undef UNROLLJ
