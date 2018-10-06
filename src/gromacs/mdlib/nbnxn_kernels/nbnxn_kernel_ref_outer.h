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

/* We could use host_func_para.nbat->xstride and host_func_para.nbat->fstride, but macros might be faster */
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
(int device_core_id)
{
    // =========== DEF DATA =============
    const nbnxn_ci_t   *nbln;
    const nbnxn_cj_t   *l_cj;
    const int          *type;
    const real         *q;
    const real         *func_para_shiftvec;
    const real         *x;
    const real         *nbfp;
    real                rcut2;

    int                 ntype2;
    real                facel;
    real               *nbfp_i; // UNUSED
    int                 n, ci, ci_sh;
    int                 ish, ishf;
    gmx_bool            do_LJ, half_LJ, do_coul, do_self;
    int                 cjind0, cjind1, cjind;
    int                 ip, jp;

    real                xi[UNROLLI*XI_STRIDE];
    real                fi[UNROLLI*FI_STRIDE];
    real                qi[UNROLLI];

#ifdef CALC_ENERGIES
    real       Vvdw_ci, Vc_ci;
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
#ifdef CALC_COUL_RF
    k_rf2 = 2*host_func_para.ic->k_rf;
#ifdef CALC_ENERGIES
    k_rf = host_func_para.ic->k_rf;
    c_rf = host_func_para.ic->c_rf;
#endif
#endif
#ifdef CALC_COUL_TAB
    tabscale = host_func_para.ic->tabq_scale;
#ifdef CALC_ENERGIES
    halfsp = 0.5/host_func_para.ic->tabq_scale;
#endif

#ifndef GMX_DOUBLE
    tab_coul_FDV0 = host_func_para.ic->tabq_coul_FDV0;
#else
    tab_coul_F    = host_func_para.ic->tabq_coul_F;
    tab_coul_V    = host_func_para.ic->tabq_coul_V;
#endif
#endif


    rcut2               = host_func_para.ic->rcoulomb*host_func_para.ic->rcoulomb;

    ntype2              = host_func_para.nbat->ntype*2;
    nbfp                = host_func_para.nbat->nbfp;
    q                   = host_func_para.nbat->q;
    type                = host_func_para.nbat->type;
    facel               = host_func_para.ic->epsfac;
    func_para_shiftvec            = host_func_para.shift_vec[0];
    x                   = host_func_para.nbat->x;

    l_cj = host_func_para.nbl->cj;
    // =========== INIT DATA =============

    /* =========== DATA'S STAT =========== //
    int ci_SZ = host_func_para.nbl->nci;
    int cj_SZ = host_func_para.nbl->ncj;
    int i, j, k;

    for(k = 0; k < 64; ++k)
    {
        if(k == host_param.host_rank)
        {
            printf("\n\n=====RANK %d=======\n", host_param.host_rank);
            printf("ci_SZ =%d, cj_SZ =%d\n\n", ci_SZ, cj_SZ);
            for(i = 0; i < ci_SZ; ++i)
            {
                nbln = &host_func_para.nbl->ci[i];
                printf("ci[%d]'s info:\n", i);
                printf("ci =%d, cjind0 =%d, cjind1 =%d\n", nbln->ci, nbln->cj_ind_start, nbln->cj_ind_end);
                printf("cj list = [");
                for(j = nbln->cj_ind_start; j < nbln->cj_ind_end; ++j)
                {
                    printf("%d, ", host_func_para.nbl->cj[j].cj);
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

#ifdef SW_NEW_ALG

    int natoms = host_func_para.nbat->natoms;
    int fstride = host_func_para.nbat->fstride;
    int sizeof_f = natoms*fstride;
    int sizeof_fshift = SHIFTS*DIM; // 135 * sizeof(float)

    /* =========== DATA'S STAT =========== //
    if(sizeof_f % 12 != 0)
    {
        OLOG("f cannot mod 12! sizeof_f =%d.", sizeof_f);
    }
    // =========== DATA'S STAT =========== */

    int start_f_div_12 = BLOCK_HEAD(device_core_id, 64, sizeof_f/12);
    int start_f = start_f_div_12*12;
    int sz_f_div_12 = BLOCK_SIZE(device_core_id, 64, sizeof_f/12);
    int sz_f = sz_f_div_12*12;
    int end_f_div_12 = start_f_div_12 + sz_f_div_12;
    int end_f = end_f_div_12*12;

#ifdef CALC_SHIFTFORCES
    real  ldm_fshift[SHIFTS*DIM];
#endif
    real* ldm_f;
    real  ldm_Vvdw = 0;
    real  ldm_Vc = 0;

    ldm_f = (real*)malloc(sz_f*sizeof(real));

#ifdef CALC_SHIFTFORCES
    memset(ldm_fshift, 0, SHIFTS*DIM*sizeof(real));
#endif
    memcpy(ldm_f, host_func_para.f + start_f, sz_f*sizeof(real));

#define IN_F_BLOCK(idx) ((idx) >= start_f_div_12 && (idx) < end_f_div_12)

    for (n = 0; n < host_func_para.nbl->nci; n++)
    {
        int i, d, write_ci;

        nbln = &host_func_para.nbl->ci[n];

        ish              = (nbln->shift & NBNXN_CI_SHIFT);
        /* x, host_func_para.f and host_func_para.fshift are assumed to be stored with stride 3 */
        ishf             = ish*DIM;
        cjind0           = nbln->cj_ind_start;
        cjind1           = nbln->cj_ind_end;
        /* Currently only works super-cells equal to sub-cells */
        ci               = nbln->ci;
        ci_sh            = (ish == CENTRAL ? ci : -1);
        write_ci         = IN_F_BLOCK(ci);

        /* We have 5 LJ/C combinations, but use only three inner loops,
         * as the other combinations are unlikely and/or not much faster:
         * inner half-LJ + C for half-LJ + C / no-LJ + C
         * inner LJ + C      for full-LJ + C
         * inner LJ          for full-LJ + no-C / half-LJ + no-C
         */
        do_LJ   = (nbln->shift & NBNXN_CI_DO_LJ(0));
        do_coul = (nbln->shift & NBNXN_CI_DO_COUL(0));
        half_LJ = ((nbln->shift & NBNXN_CI_HALF_LJ(0)) || !do_LJ) && do_coul;
    
        do_self = do_coul;

#ifdef CALC_ENERGIES
        Vvdw_ci = 0;
        Vc_ci   = 0;
#endif

        //TODO: ldm load: x, q， func_para_shiftvec
        for (i = 0; i < UNROLLI; i++)
        {
            for (d = 0; d < DIM; d++)
            {
                xi[i*XI_STRIDE+d] = x[(ci*UNROLLI+i)*X_STRIDE+d] + func_para_shiftvec[ishf+d];
                fi[i*FI_STRIDE+d] = 0;
            }

            qi[i] = facel*q[ci*UNROLLI+i];
        }
#ifdef DENUG_F
        TLOG("Miao 1.\n");
#endif
#ifdef CALC_ENERGIES
        if (do_self && write_ci)
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
                    //TODO: REDUCE SUM
                    /* Coulomb self interaction */
                    ldm_Vc   -= qi[i]*q[ci*UNROLLI+i]*Vc_sub_self;
                }
            }
        }
#endif  /* CALC_ENERGIES */
#ifdef DENUG_F
        TLOG("Miao 2.\n");
#endif
        cjind = cjind0;
        while (cjind < cjind1 && host_func_para.nbl->cj[cjind].excl != 0xffff)
        {
#define CHECK_EXCLS
            if (half_LJ)
            {
#define CALC_COULOMB
#define HALF_LJ
#include "nbnxn_kernel_ref_inner.h"
#undef HALF_LJ
#undef CALC_COULOMB
            }
            else if (do_coul)
            {
#define CALC_COULOMB
#include "nbnxn_kernel_ref_inner.h"
#undef CALC_COULOMB
            }
            else
            {
#include "nbnxn_kernel_ref_inner.h"
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
#include "nbnxn_kernel_ref_inner.h"
#undef HALF_LJ
#undef CALC_COULOMB
            }
            else if (do_coul)
            {
#define CALC_COULOMB
#include "nbnxn_kernel_ref_inner.h"
#undef CALC_COULOMB
            }
            else
            {
#include "nbnxn_kernel_ref_inner.h"
            }
        }
#ifdef DENUG_F
        TLOG("Miao 3.\n");
#endif
        if(write_ci)
        {
            /* Add accumulated i-forces to the force array */
            for (i = 0; i < UNROLLI; i++)
            {
                for (d = 0; d < DIM; d++)
                {
                    //TODO: REDUCE SUM
                    ldm_f[(ci*UNROLLI+i)*F_STRIDE+d-start_f] += fi[i*FI_STRIDE+d];
                }
            }
#ifdef CALC_SHIFTFORCES
            if (host_func_para.expand_fshift != NULL)
            {
                /* Add i forces to shifted force list */
                for (i = 0; i < UNROLLI; i++)
                {
                    for (d = 0; d < DIM; d++)
                    {
                        //TODO: REDUCE SUM
                        ldm_fshift[ishf+d] += fi[i*FI_STRIDE+d];
                    }
                }
            }
#endif
#ifdef CALC_ENERGIES
            //TODO: REDUCE SUM
            ldm_Vvdw += Vvdw_ci;
            ldm_Vc   += Vc_ci;
#endif
        }
#ifdef DENUG_F
        TLOG("Miao 4.\n");
#endif
    }

    memcpy(host_func_para.f + start_f, ldm_f, sz_f*sizeof(real));
    free(ldm_f);
#ifdef CALC_SHIFTFORCES
    if (host_func_para.expand_fshift != NULL)
    {
        memcpy(host_func_para.expand_fshift+device_core_id*SHIFTS*DIM, ldm_fshift, SHIFTS*DIM*sizeof(real));
    }
#endif
#ifdef CALC_ENERGIES
    host_func_para.expand_Vvdw[device_core_id] = ldm_Vvdw;
    host_func_para.expand_Vc  [device_core_id] = ldm_Vc;
#endif
#ifdef DENUG_F
    TLOG("Miao 5.\n");
#endif

#undef IN_F_BLOCK

#else  /* SW_NEW_ALG */

    int start_nci = BLOCK_HEAD(device_core_id, 64, host_func_para.nbl->nci);
    int end_nci = start_nci + BLOCK_SIZE(device_core_id, 64, host_func_para.nbl->nci);

    for (n = start_nci; n < end_nci; n++)
    {
        int i, d;

        nbln = &host_func_para.nbl->ci[n];

        ish              = (nbln->shift & NBNXN_CI_SHIFT);
        /* x, host_func_para.f and host_func_para.fshift are assumed to be stored with stride 3 */
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
    
        do_self = do_coul;

#ifdef CALC_ENERGIES
        Vvdw_ci = 0;
        Vc_ci   = 0;
#endif

        //TODO: ldm load: x, q， func_para_shiftvec
        for (i = 0; i < UNROLLI; i++)
        {
            for (d = 0; d < DIM; d++)
            {
                xi[i*XI_STRIDE+d] = x[(ci*UNROLLI+i)*X_STRIDE+d] + func_para_shiftvec[ishf+d];
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
                    //TODO: REDUCE SUM
                    /* Coulomb self interaction */
                    host_func_para.Vc[0]   -= qi[i]*q[ci*UNROLLI+i]*Vc_sub_self;
                }
            }
        }
#endif  /* CALC_ENERGIES */
// ==============================================
        cjind = cjind0;
        while (cjind < cjind1 && host_func_para.nbl->cj[cjind].excl != 0xffff)
        {
#define CHECK_EXCLS
            if (half_LJ)
            {
#define CALC_COULOMB
#define HALF_LJ
#include "nbnxn_kernel_ref_inner.h"
#undef HALF_LJ
#undef CALC_COULOMB
            }
            else if (do_coul)
            {
#define CALC_COULOMB
#include "nbnxn_kernel_ref_inner.h"
#undef CALC_COULOMB
            }
            else
            {
#include "nbnxn_kernel_ref_inner.h"
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
#include "nbnxn_kernel_ref_inner.h"
#undef HALF_LJ
#undef CALC_COULOMB
            }
            else if (do_coul)
            {
#define CALC_COULOMB
#include "nbnxn_kernel_ref_inner.h"
#undef CALC_COULOMB
            }
            else
            {
#include "nbnxn_kernel_ref_inner.h"
            }
        }
//===========================================
        /* Add accumulated i-forces to the force array */
        for (i = 0; i < UNROLLI; i++)
        {
            for (d = 0; d < DIM; d++)
            {
                //TODO: REDUCE SUM
                host_func_para.f[(ci*UNROLLI+i)*F_STRIDE+d] += fi[i*FI_STRIDE+d];
            }
        }
#ifdef CALC_SHIFTFORCES
        if (host_func_para.fshift != NULL)
        {
            /* Add i forces to shifted force list */
            for (i = 0; i < UNROLLI; i++)
            {
                for (d = 0; d < DIM; d++)
                {
                    //TODO: REDUCE SUM
                    host_func_para.fshift[ishf+d] += fi[i*FI_STRIDE+d];
                }
            }
        }
#endif
// ============================================
#ifdef CALC_ENERGIES
        //TODO: REDUCE SUM
        *host_func_para.Vvdw += Vvdw_ci;
        *host_func_para.Vc   += Vc_ci;
#endif
    }

#endif  /* SW_NEW_ALG */
}

#undef CALC_SHIFTFORCES

#undef X_STRIDE
#undef F_STRIDE
#undef XI_STRIDE
#undef FI_STRIDE

#undef UNROLLI
#undef UNROLLJ

