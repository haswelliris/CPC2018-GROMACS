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
#include "gmxpre.h"

#include "nbnxn_kernel_ref.h"

#include "config.h"

#include <assert.h>
#include <math.h>
#include <string.h>

#include "gromacs/legacyheaders/force.h"
#include "gromacs/legacyheaders/gmx_omp_nthreads.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_common.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/utility/smalloc.h"

#ifdef __cplusplus
extern "C" {
#endif
#include "gromacs/mdlib/nbnxn_kernels/sw_subcore/sw/SwHost.h"
#ifdef __cplusplus
}
#endif

/*! \brief Typedefs for declaring lookup tables of kernel functions.
 */

typedef void (*p_nbk_func_noener)(int device_core_id);

typedef void (*p_nbk_func_ener)(int device_core_id);

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    nbnxn_pairlist_t     *nbl;
    nbnxn_atomdata_t     *nbat;
    interaction_const_t  *ic;
    rvec                       *shift_vec;
    real                       *f;
    real                       *expand_Vvdw;
    real                       *expand_Vc;
    real                       *expand_fshift;

    real *tabq_coul_F;
    real *tabq_coul_V;
    real *tabq_coul_FDV0;

    int *f_start;
    int *f_end;
} func_para_t;

func_para_t host_func_para;

#ifdef __cplusplus
}
#endif

/* Analytical reaction-field kernels */
//++++++++++++++++++++++++++++++++++++
#define CALC_COUL_RF
    //===================================
    #define LJ_CUT
    #include "nbnxn_kernel_ref_includes.h"
    #undef LJ_CUT

    #ifndef SW_TEST_FUNC /* in SwConfig.h */
        #define LJ_FORCE_SWITCH
        #include "nbnxn_kernel_ref_includes.h"
        #undef LJ_FORCE_SWITCH

        #define LJ_POT_SWITCH
        #include "nbnxn_kernel_ref_includes.h"
        #undef LJ_POT_SWITCH

        #define LJ_EWALD
            #define LJ_CUT
                #define LJ_EWALD_COMB_GEOM
                #include "nbnxn_kernel_ref_includes.h"
                #undef LJ_EWALD_COMB_GEOM

                #define LJ_EWALD_COMB_LB
                #include "nbnxn_kernel_ref_includes.h"
                #undef LJ_EWALD_COMB_LB
            #undef LJ_CUT
        #undef LJ_EWALD
    #endif
    //====================================
#undef CALC_COUL_RF
//++++++++++++++++++++++++++++++++++++


/* Tabulated exclusion interaction electrostatics kernels */
//++++++++++++++++++++++++++++++++++++
#define CALC_COUL_TAB
    //====================================
    #define LJ_CUT
    #include "nbnxn_kernel_ref_includes.h"
    #undef LJ_CUT

    #ifndef SW_TEST_FUNC /* in SwConfig.h */
        #define LJ_FORCE_SWITCH
        #include "nbnxn_kernel_ref_includes.h"
        #undef LJ_FORCE_SWITCH
    
        #define LJ_POT_SWITCH
        #include "nbnxn_kernel_ref_includes.h"
        #undef LJ_POT_SWITCH

        #define LJ_EWALD
            #define LJ_CUT
                #define LJ_EWALD_COMB_GEOM
                #include "nbnxn_kernel_ref_includes.h"
                #undef LJ_EWALD_COMB_GEOM

                #define LJ_EWALD_COMB_LB
                #include "nbnxn_kernel_ref_includes.h"
                #undef LJ_EWALD_COMB_LB
            #undef LJ_CUT
        #undef LJ_EWALD
    #endif
    //====================================



    /* Twin-range cut-off kernels */
    //====================================
    #ifndef SW_TEST_FUNC /* in SwConfig.h */

    #define VDW_CUTOFF_CHECK
        #define LJ_CUT
        #include "nbnxn_kernel_ref_includes.h"
        #undef LJ_CUT

        #define LJ_FORCE_SWITCH
        #include "nbnxn_kernel_ref_includes.h"
        #undef LJ_FORCE_SWITCH

        #define LJ_POT_SWITCH
        #include "nbnxn_kernel_ref_includes.h"
        #undef LJ_POT_SWITCH

        #define LJ_EWALD
            #define LJ_CUT
                #define LJ_EWALD_COMB_GEOM
                #include "nbnxn_kernel_ref_includes.h"
                #undef LJ_EWALD_COMB_GEOM

                #define LJ_EWALD_COMB_LB
                #include "nbnxn_kernel_ref_includes.h"
                #undef LJ_EWALD_COMB_LB
            #undef LJ_CUT
        #undef LJ_EWALD
    #undef VDW_CUTOFF_CHECK

    #endif
    //====================================
#undef CALC_COUL_TAB
//++++++++++++++++++++++++++++++++++++

enum {
    coultRF, coultTAB, coultTAB_TWIN, coultNR
};

enum {
    vdwtCUT, vdwtFSWITCH, vdwtPSWITCH, vdwtEWALDGEOM, vdwtEWALDLB, vdwtNR
};

p_nbk_func_noener p_nbk_c_noener[coultNR][vdwtNR] =
{
#ifdef SW_TEST_FUNC /* in SwConfig.h */
    { nbnxn_kernel_ElecRF_VdwLJ_F_ref,           NULL, NULL, NULL, NULL           },
    { nbnxn_kernel_ElecQSTab_VdwLJ_F_ref,        NULL, NULL, NULL, NULL           },
    { NULL, NULL, NULL, NULL, NULL  }
#else
    { nbnxn_kernel_ElecRF_VdwLJ_F_ref,           nbnxn_kernel_ElecRF_VdwLJFsw_F_ref,           nbnxn_kernel_ElecRF_VdwLJPsw_F_ref,           nbnxn_kernel_ElecRF_VdwLJEwCombGeom_F_ref,           nbnxn_kernel_ElecRF_VdwLJEwCombLB_F_ref           },
    { nbnxn_kernel_ElecQSTab_VdwLJ_F_ref,        nbnxn_kernel_ElecQSTab_VdwLJFsw_F_ref,        nbnxn_kernel_ElecQSTab_VdwLJPsw_F_ref,        nbnxn_kernel_ElecQSTab_VdwLJEwCombGeom_F_ref,        nbnxn_kernel_ElecQSTab_VdwLJEwCombLB_F_ref        },
    { nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_F_ref, nbnxn_kernel_ElecQSTabTwinCut_VdwLJFsw_F_ref, nbnxn_kernel_ElecQSTabTwinCut_VdwLJPsw_F_ref, nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombGeom_F_ref, nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombLB_F_ref }
#endif
};

p_nbk_func_ener p_nbk_c_ener[coultNR][vdwtNR] =
{
#ifdef SW_TEST_FUNC /* in SwConfig.h */
    { nbnxn_kernel_ElecRF_VdwLJ_VF_ref,           NULL, NULL, NULL, NULL           },
    { nbnxn_kernel_ElecQSTab_VdwLJ_VF_ref,        NULL, NULL, NULL, NULL           },
    { NULL, NULL, NULL, NULL, NULL  }
#else
    { nbnxn_kernel_ElecRF_VdwLJ_VF_ref,           nbnxn_kernel_ElecRF_VdwLJFsw_VF_ref,           nbnxn_kernel_ElecRF_VdwLJPsw_VF_ref,           nbnxn_kernel_ElecRF_VdwLJEwCombGeom_VF_ref,           nbnxn_kernel_ElecRF_VdwLJEwCombLB_VF_ref            },
    { nbnxn_kernel_ElecQSTab_VdwLJ_VF_ref,        nbnxn_kernel_ElecQSTab_VdwLJFsw_VF_ref,        nbnxn_kernel_ElecQSTab_VdwLJPsw_VF_ref,        nbnxn_kernel_ElecQSTab_VdwLJEwCombGeom_VF_ref,        nbnxn_kernel_ElecQSTab_VdwLJEwCombLB_VF_ref         },
    { nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VF_ref, nbnxn_kernel_ElecQSTabTwinCut_VdwLJFsw_VF_ref, nbnxn_kernel_ElecQSTabTwinCut_VdwLJPsw_VF_ref, nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombGeom_VF_ref, nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombLB_VF_ref  }
#endif
};
#ifdef SW_ENERGRP /* in SwConfig */

p_nbk_func_ener p_nbk_c_energrp[coultNR][vdwtNR] =
{
#ifdef SW_TEST_FUNC /* in SwConfig.h */
    { nbnxn_kernel_ElecRF_VdwLJ_VgrpF_ref,           NULL, NULL, NULL, NULL           },
    { nbnxn_kernel_ElecQSTab_VdwLJ_VgrpF_ref,        NULL, NULL, NULL, NULL           },
    { NULL, NULL, NULL, NULL, NULL  }
#else
    { nbnxn_kernel_ElecRF_VdwLJ_VgrpF_ref,           nbnxn_kernel_ElecRF_VdwLJFsw_VgrpF_ref,           nbnxn_kernel_ElecRF_VdwLJPsw_VgrpF_ref,           nbnxn_kernel_ElecRF_VdwLJEwCombGeom_VgrpF_ref,           nbnxn_kernel_ElecRF_VdwLJEwCombLB_VgrpF_ref           },
    { nbnxn_kernel_ElecQSTab_VdwLJ_VgrpF_ref,        nbnxn_kernel_ElecQSTab_VdwLJFsw_VgrpF_ref,        nbnxn_kernel_ElecQSTab_VdwLJPsw_VgrpF_ref,        nbnxn_kernel_ElecQSTab_VdwLJEwCombGeom_VgrpF_ref,        nbnxn_kernel_ElecQSTab_VdwLJEwCombLB_VgrpF_ref        },
    { nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VgrpF_ref, nbnxn_kernel_ElecQSTabTwinCut_VdwLJFsw_VgrpF_ref, nbnxn_kernel_ElecQSTabTwinCut_VdwLJPsw_VgrpF_ref, nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombGeom_VgrpF_ref, nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombLB_VgrpF_ref }
#endif
};

#endif

void fake_device_run()
{
    int func_type = host_out_param[FUNC_TYPE];
    int func_i    = host_out_param[FUNC_I];
    int func_j    = host_out_param[FUNC_J];
    int k;

    if(func_type == FUNC_NO_ENER)
    {
        for(k = 0; k < 64; ++k)
        {
            p_nbk_c_noener[func_i][func_j](k);
        }
    }
    else if(func_type == FUNC_ENER)
    {
        for(k = 0; k < 64; ++k)
        {
            p_nbk_c_ener[func_i][func_j](k);
        }
    }
    else
    {
        gmx_incons("Unknown FUNCTYPE");
    }
}

void deep_copy_nbl(nbnxn_pairlist_t *dst, nbnxn_pairlist_t *src, int new_dst, int del_src_and_no_copy)
{
    if(new_dst)
    {
        dst->ci = (nbnxn_ci_t*)malloc(src->nci*sizeof(nbnxn_ci_t));
        dst->cj = (nbnxn_cj_t*)malloc(src->ncj*sizeof(nbnxn_cj_t));

        dst->nci = src->nci;
        dst->ncj = src->ncj;
        memcpy(dst->ci, src->ci, src->nci*sizeof(nbnxn_ci_t));
        memcpy(dst->cj, src->cj, src->ncj*sizeof(nbnxn_cj_t));
    }
    if(del_src_and_no_copy)
    {
        free(src->ci);
        free(src->cj);
    }
}

void deep_copy_nbat(nbnxn_atomdata_t *dst, nbnxn_atomdata_t *src, int new_dst, int del_src_and_no_copy)
{
    if(new_dst)
    {
        dst->x = (real*)malloc(src->natoms*src->xstride*sizeof(real));
        dst->q = (real*)malloc(src->natoms*sizeof(real));
        dst->nbfp = (real*)malloc(src->ntype*src->ntype*2*sizeof(real));
        dst->type = (int*)malloc(src->natoms*sizeof(int));
        // if the nbfp can load to LDM?
        // TLOG("nbfpSZ =%d Byte\n", src->ntype*src->ntype*2*sizeof(real));
        // watch VAR of type
        // int i;
        // TLOG("type = [ ");
        // for(i = 0; i < src->natoms; ++i)
        // {
        //     if(host_param.host_rank == 0)
        //     {
        //         printf(" %d,", src->type[i]);
        //     }
        // }
        // if(host_param.host_rank == 0)
        // {
        //         printf(" ]\n");
        // }
        dst->natoms = src->natoms;
        dst->ntype = src->ntype;
        dst->xstride = src->xstride;
        dst->fstride = src->fstride;
        memcpy(dst->x, src->x, src->natoms*src->xstride*sizeof(real));
        memcpy(dst->q, src->q, src->natoms*sizeof(real));
        memcpy(dst->nbfp, src->nbfp, src->ntype*src->ntype*2*sizeof(real));
        memcpy(dst->type, src->type, src->natoms*sizeof(int));
    }
    if(del_src_and_no_copy)
    {
        free(src->x);
        free(src->q);
        free(src->nbfp);
        free(src->type);
    }
}

void subcore_loadbalance(nbnxn_atomdata_t *nbat, nbnxn_pairlist_t *nbl, int *f_start, int *f_end)
{
    int n, i, device_core_id;
    int f_count_len = nbat->natoms/4;

	int *f_count = (int*)malloc(f_count_len*sizeof(int));
	memset(f_count, 0, f_count_len*sizeof(int));

	int tot_work_count = 0, avg_work_load = 0;

	int *workload = (int*)malloc(64*sizeof(int));
	memset(workload, 0, 64*sizeof(int));

	for (n = 0; n < nbl->nci; n++) {
	    f_count[nbl->ci[n].ci] += nbl->ci[n].cj_ind_end - nbl->ci[n].cj_ind_start;
	    tot_work_count += nbl->ci[n].cj_ind_end - nbl->ci[n].cj_ind_start;
	    for (i = nbl->ci[n].cj_ind_start; i < nbl->ci[n].cj_ind_end; i++) {
	    	f_count[nbl->cj[i].cj]++;
	    	tot_work_count++;
	    }
	}
	avg_work_load = tot_work_count/64+1;
	int p = 0;
	while (f_count[f_count_len-1] == 0)
	    f_count_len--;
    for (device_core_id = 0; device_core_id < 64; device_core_id++) {
	    while (f_count[p] == 0) p++;
    	f_start[device_core_id] = p;
		while (p < f_count_len && 
    			(workload[device_core_id] < avg_work_load || (f_count_len - p) > (63 - device_core_id)*MAX_F_LDM_SIZE ) && 
    			(p - f_start[device_core_id]) <= MAX_F_LDM_SIZE)
    	{
    		workload[device_core_id]+=f_count[p];
    		p++;
    	}
    	f_end[device_core_id] = p;
	}
    // for (device_core_id = 0; device_core_id < 64; device_core_id++) {
    //     f_start[device_core_id] = BLOCK_HEAD(device_core_id, 64, f_count_len);
    //     f_end  [device_core_id] = f_start[device_core_id] + BLOCK_SIZE(device_core_id, 64, f_count_len);
    // }
    free(f_count);
    free(workload);
}

void
nbnxn_kernel_ref(const nbnxn_pairlist_set_t *nbl_list,
                 const nbnxn_atomdata_t     *nbat,
                 const interaction_const_t  *ic,
                 rvec                       *shift_vec,
                 int                         force_flags,
                 int                         clearF,
                 real                       *fshift,
                 real                       *Vc,
                 real                       *Vvdw,
                 gmx_wallcycle_t             wcycle)
{
    int                nnbl;
    nbnxn_pairlist_t **nbl;
    int                coult;
    int                vdwt;
    int                nb;
    int                nthreads gmx_unused;

    nnbl = nbl_list->nnbl;
    nbl  = nbl_list->nbl;

    if (EEL_RF(ic->eeltype) || ic->eeltype == eelCUT)
    {
        coult = coultRF;
    }
    else
    {
        if (ic->rcoulomb == ic->rvdw)
        {
            coult = coultTAB;
        }
        else
        {
            coult = coultTAB_TWIN;
        }
    }

    if (ic->vdwtype == evdwCUT)
    {
        switch (ic->vdw_modifier)
        {
            case eintmodPOTSHIFT:
            case eintmodNONE:
                vdwt = vdwtCUT;
                break;
            case eintmodFORCESWITCH:
                vdwt = vdwtFSWITCH;
                break;
            case eintmodPOTSWITCH:
                vdwt = vdwtPSWITCH;
                break;
            default:
                gmx_incons("Unsupported VdW modifier");
                break;
        }
    }
    else if (ic->vdwtype == evdwPME)
    {
        if (ic->ljpme_comb_rule == ljcrGEOM)
        {
            assert(nbat->comb_rule == ljcrGEOM);
            vdwt = vdwtEWALDGEOM;
        }
        else
        {
            assert(nbat->comb_rule == ljcrLB);
            vdwt = vdwtEWALDLB;
        }
    }
    else
    {
        gmx_incons("Unsupported vdwtype in nbnxn reference kernel");
    }

    nthreads = gmx_omp_nthreads_get(emntNonbonded);
#pragma omp parallel for schedule(static) num_threads(nthreads)
    for (nb = 0; nb < nnbl; nb++)
    {
        int *f_start = (int*)malloc(64*sizeof(int));
    	memset(f_start, 0, 64*sizeof(int));
	    int *f_end = (int*)malloc(64*sizeof(int));
	    memset(f_end, 0, 64*sizeof(int));
        //wallcycle_sub_start(wcycle, ewcsMEMCPY);
        subcore_loadbalance(nbat, nbl[nb], f_start, f_end);
        //wallcycle_sub_stop(wcycle, ewcsMEMCPY);
        nbnxn_atomdata_output_t *out;
        real                    *fshift_p;

        out = &nbat->out[nb];

        if (clearF == enbvClearFYes)
        {
            clear_f(nbat, nb, out->f);
        }

        if ((force_flags & GMX_FORCE_VIRIAL) && nnbl == 1)
        {
            fshift_p = fshift;
        }
        else
        {
            fshift_p = out->fshift;

            if (clearF == enbvClearFYes)
            {
                clear_fshift(fshift_p);
            }
        }
        //wallcycle_sub_start(wcycle, ewcsMEMCPY);
        real *expand_fshift = (real*)malloc(SHIFTS*DIM*64*sizeof(real));

        real *other_f       = (real*)malloc(nbat->natoms*nbat->fstride*sizeof(real));
        memcpy(other_f, out->f, nbat->natoms*nbat->fstride*sizeof(real));

        nbnxn_pairlist_t other_nbl;
        deep_copy_nbl(&other_nbl, nbl[nb], 1, 0);

        nbnxn_atomdata_t other_nbat;
        deep_copy_nbat(&other_nbat, nbat, 1, 0);

        rvec *other_shift_vec = (rvec*)malloc(SHIFTS*DIM*sizeof(real));
        memcpy(other_shift_vec, shift_vec, SHIFTS*DIM*sizeof(real));

        real *other_tabq_coul_F = NULL;
        real *other_tabq_coul_V = NULL;
        real *other_tabq_coul_FDV0 = NULL;

#ifndef GMX_DOUBLE
        other_tabq_coul_FDV0 = (real*)malloc(ic->tabq_size*4*sizeof(real));
        memcpy(other_tabq_coul_FDV0, ic->tabq_coul_FDV0, ic->tabq_size*4*sizeof(real));
#else
        other_tabq_coul_F    = (real*)malloc(ic->tabq_size*sizeof(real));
        other_tabq_coul_V    = (real*)malloc(ic->tabq_size*sizeof(real));
        memcpy(other_tabq_coul_F, ic->tabq_coul_F, ic->tabq_size*sizeof(real));
        memcpy(other_tabq_coul_V, ic->tabq_coul_V, ic->tabq_size*sizeof(real));
#endif
        //wallcycle_sub_stop(wcycle, ewcsMEMCPY);
        // if the tabq_coul_FDV0 can load to LDM?
        // TLOG("tabq_coul_FDV0_SZ =%d Byte\n", ic->tabq_size*4*sizeof(real));
        // TLOG("tabq_size =%d, ntype =%d, natoms =%d\n", ic->tabq_size, nbat->ntype, nbat->natoms);
#ifdef DEBUG_SDLB
        TLOG("kaCHI sizeof(nbnxn_pairlist_t) =%d\n", sizeof(nbnxn_pairlist_t));
#endif
#ifdef HAHAHAHAHAHAHAHA
        if(coult == 1)
        {
            usleep(40000);
        }
        else
        {
            usleep(400000);
        }
#endif
        if (!(force_flags & GMX_FORCE_ENERGY))
        {
            host_out_param[PARAM_DEVICE_ACTION] = DEVICE_ACTION_RUN;
            host_out_param[FUNC_TYPE] = FUNC_NO_ENER;
            host_out_param[FUNC_I] = coult;
            host_out_param[FUNC_J] = vdwt;
            host_out_param[FUNC_PARAM_PTR] = (long)&host_func_para;
            host_func_para.nbl = &other_nbl; // read only
            host_func_para.nbat = &other_nbat; // read only
            host_func_para.ic = ic;      // read only
            host_func_para.shift_vec = other_shift_vec; // read only
            host_func_para.f = other_f; // write only, reduce FIN
            host_func_para.expand_Vvdw = NULL;// write only, reduce FIN
            host_func_para.expand_Vc = NULL;// write only, reduce FIN
            host_func_para.expand_fshift = expand_fshift;// write only, reduce FIN
            host_func_para.tabq_coul_F = other_tabq_coul_F;// read only
            host_func_para.tabq_coul_V = other_tabq_coul_V;// read only
            host_func_para.tabq_coul_FDV0 = other_tabq_coul_FDV0;// read only
            host_func_para.f_start = f_start;
            host_func_para.f_end   = f_end;
#ifdef SW_HOST_LOG /* in SwConfig */
            if((host_param.host_rank + host_notice_counter) % 64 == 0)
            {
                OLOG("FuncType =%d, I =%d, J =%d\n", FUNC_NO_ENER, coult, vdwt);
            }
#endif
            notice_device();
            wait_device();

            /* Don't calculate energies */
            //fake_device_run();
            //p_nbk_c_noener[coult][vdwt]();

        }
        else if (out->nV == 1)
        {
            /* No energy groups */
            out->Vvdw[0] = 0;
            out->Vc[0]   = 0;

            real *expand_Vvdw   = (real*)malloc(64*sizeof(real));
            real *expand_Vc     = (real*)malloc(64*sizeof(real));

            host_out_param[PARAM_DEVICE_ACTION] = DEVICE_ACTION_RUN;
            host_out_param[FUNC_TYPE] = FUNC_ENER;
            host_out_param[FUNC_I] = coult;
            host_out_param[FUNC_J] = vdwt;
            host_out_param[FUNC_PARAM_PTR] = (long)&host_func_para;
            host_func_para.nbl = &other_nbl;
            host_func_para.nbat = &other_nbat;
            host_func_para.ic = ic;
            host_func_para.shift_vec = other_shift_vec;
            host_func_para.f = other_f;
            host_func_para.expand_Vvdw = expand_Vvdw;
            host_func_para.expand_Vc = expand_Vc;
            host_func_para.expand_fshift = expand_fshift;
            host_func_para.tabq_coul_F = other_tabq_coul_F;
            host_func_para.tabq_coul_V = other_tabq_coul_V;
            host_func_para.tabq_coul_FDV0 = other_tabq_coul_FDV0;
            host_func_para.f_start = f_start;
            host_func_para.f_end   = f_end;
#ifdef SW_HOST_LOG /* in SwConfig */
            if((host_param.host_rank + host_notice_counter) % 64 == 0)
            {
                OLOG("FuncType =%d, I =%d, J =%d\n", FUNC_ENER, coult, vdwt);
            }
#endif
            notice_device();
            wait_device();

            //fake_device_run();
            //p_nbk_c_ener[coult][vdwt]();

#ifdef DEBUG_FPEX
            TLOG("MOee 0.\n");
            //wait_device();
#endif
            // reduce Vvdw
            int i;
            for(i = 0; i < 64 - 1; ++i)
            {
                expand_Vvdw[i+1] += expand_Vvdw[i];
            }
            out->Vvdw[0] = expand_Vvdw[63];
            free(expand_Vvdw);
#ifdef DEBUG_FPEX
            {
                TLOG("MOee 0.2\n");
                if(host_param.host_rank == 0)
                {
                    printf("Vc= [");
                    for(i = 0; i < 64; ++i)
                    {
                        printf(" %f", expand_Vc[i]);
                    }
                    printf("]\n");
                    wait_device();
                }
                else
                {
                    wait_device();
                }
            }
#endif
            // reduce Vc
            for(i = 0; i < 64 - 1; ++i)
            {
                expand_Vc[i+1] += expand_Vc[i];
            }
            out->Vc[0] = expand_Vc[63];
            free(expand_Vc);
#ifdef DEBUG_FPEX
            TLOG("MOee 1.\n");
            wait_device();
#endif
        }
        else
        {
            /* Calculate energy group contributions */
            int i;

            for (i = 0; i < out->nV; i++)
            {
                out->Vvdw[i] = 0;
            }
            for (i = 0; i < out->nV; i++)
            {
                out->Vc[i] = 0;
            }

            real *expand_Vvdw   = (real*)malloc(out->nV*64*sizeof(real));
            real *expand_Vc     = (real*)malloc(out->nV*64*sizeof(real));

            host_out_param[PARAM_DEVICE_ACTION] = DEVICE_ACTION_RUN;
            host_out_param[FUNC_TYPE] = FUNC_ENERGRP;
            host_out_param[FUNC_I] = coult;
            host_out_param[FUNC_J] = vdwt;
            host_out_param[FUNC_PARAM_PTR] = (long)&host_func_para;
            host_func_para.nbl = &other_nbl;
            host_func_para.nbat = &other_nbat;
            host_func_para.ic = ic;
            host_func_para.shift_vec = other_shift_vec;
            host_func_para.f = other_f;
            host_func_para.expand_Vvdw = expand_Vvdw;
            host_func_para.expand_Vc = expand_Vc;
            host_func_para.expand_fshift = expand_fshift;
            host_func_para.tabq_coul_F = other_tabq_coul_F;
            host_func_para.tabq_coul_V = other_tabq_coul_V;
            host_func_para.tabq_coul_FDV0 = other_tabq_coul_FDV0;
            host_func_para.f_start = f_start;
            host_func_para.f_end   = f_end;
#ifdef SW_HOST_LOG /* in SwConfig */
            if((host_param.host_rank + host_notice_counter) % 64 == 0)
            {
                OLOG("FuncType =%d, I =%d, J =%d\n", FUNC_ENERGRP, coult, vdwt);
            }
#endif
            notice_device();
            wait_device();

            //fake_device_run();
#ifdef SW_ENERGRP /* in SwConfig */
            //p_nbk_c_energrp[coult][vdwt]();
#else
            static int grp_call=0;
            if(grp_call == 0)
            {
                OLOG("No Energy Group Function.\n");
                grp_call = 1;
            }

#endif
            // reduce Vvdw
            int j;
            for(i = 0; i < 64 - 1; ++i)
            {
                for (j = 0; j < out->nV; j++)
                {
                    expand_Vvdw[(i+1)*out->nV+j] += expand_Vvdw[(i)*out->nV+j];
                }
            }
            for (j = 0; j < out->nV; j++)
            {
                out->Vvdw[j] = expand_Vvdw[(63)*out->nV+j];
            }
            free(expand_Vvdw);

            // reduce Vc
            for(i = 0; i < 64 - 1; ++i)
            {
                for (j = 0; j < out->nV; j++)
                {
                    expand_Vc[(i+1)*out->nV+j] += expand_Vc[(i)*out->nV+j];
                }
            }
            for (j = 0; j < out->nV; j++)
            {
                out->Vc[j] = expand_Vc[(63)*out->nV+j];
            }
            free(expand_Vc);
        }
#ifdef DEBUG_FPEX
        TLOG("MOee 2.\n");
        wait_device();
#endif
        // reduce fshift
        int i, j;
        for(i = 0; i < 64 - 1; ++i)
        {
            for(j = 0; j < SHIFTS*DIM; ++j)
            {
                expand_fshift[(i+1)*SHIFTS*DIM+j] += expand_fshift[(i)*SHIFTS*DIM+j];
            }
        }
#ifdef DEBUG_FPEX
        TLOG("MOee 3.\n");
        wait_device();
#endif
        for(j = 0; j < SHIFTS*DIM; ++j)
        {
            fshift_p[j] += expand_fshift[63*SHIFTS*DIM+j];
        }
        free(expand_fshift);
#ifdef DEBUG_FPEX
        TLOG("MOee 4.\n");
        wait_device();
#endif
        //wallcycle_sub_start(wcycle, ewcsMEMCPY);
        memcpy(out->f, other_f, nbat->natoms*nbat->fstride*sizeof(real));
        free(other_f);

        deep_copy_nbl(nbl[nb], &other_nbl, 0, 1);

        deep_copy_nbat(nbat, &other_nbat, 0, 1);
#ifdef DEBUG_FPEX
        TLOG("MOee 5.\n");
        wait_device();
#endif
        free(other_shift_vec);

#ifndef GMX_DOUBLE
        free(other_tabq_coul_FDV0);
#else
        free(other_tabq_coul_F);
        free(other_tabq_coul_V);
#endif
        //wallcycle_sub_stop(wcycle, ewcsMEMCPY);
        free(f_start);
        free(f_end);
#ifdef DEBUG_FPEX
        TLOG("MOee 6.\n");
        wait_device();
#endif
    }

    if (force_flags & GMX_FORCE_ENERGY)
    {
        reduce_energies_over_lists(nbat, nnbl, Vvdw, Vc);
    }
}
