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
	    for (i = nbl->ci[n].cj_ind_start; i < nbl->ci[n].cj_ind_end; i+=load_balance_step) {
	    	f_count[nbl->cj[i].cj]+=load_balance_step;
	    	tot_work_count+=load_balance_step;
	    }
	}
	avg_work_load = tot_work_count/64+1;
	int p = 0;
	while (f_count[f_count_len-1] == 0)
	    f_count_len--;
    for (device_core_id = 0; device_core_id < 64; device_core_id++)
    {
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
        if (device_core_id != 0) {
            tot_work_count -= workload[device_core_id];
            avg_work_load = tot_work_count/(63-device_core_id)+1;
        }
	}
    for (device_core_id = 63; workload[device_core_id] == 0 && device_core_id >= 0; device_core_id--)
    {
        int max_workload = workload[0], max_workload_index = 0;
        for (i = 1; i < device_core_id; i++)
            if (max_workload < workload[i]) {
                max_workload = workload[i];
                max_workload_index = i;
            }
        f_end[device_core_id] = f_end[max_workload_index];
        int mid = (f_start[max_workload_index]+f_end[max_workload_index])/2;
        f_end[max_workload_index] = mid;
        f_start[device_core_id] = mid;
    }
    // for (device_core_id = 0; device_core_id < 64; device_core_id++)
    //     if (f_end[device_core_id]-f_start[device_core_id] == 0)
    // printf("!");
    // for (device_core_id = 0; device_core_id < 64; device_core_id++) {
    //     f_start[device_core_id] = BLOCK_HEAD(device_core_id, 64, f_count_len);
    //     f_end  [device_core_id] = f_start[device_core_id] + BLOCK_SIZE(device_core_id, 64, f_count_len);
    // }
    free(f_count);
    free(workload);
}

struct 
{
    nbnxn_pairlist_set_t *nbl_list;
    nbnxn_atomdata_t     *nbat;
    interaction_const_t  *ic;
    rvec                 *shift_vec;
    int                   force_flags;
    int                   clearF;
    real                 *fshift;
    real                 *Vc;
    real                 *Vvdw;
    gmx_wallcycle_t      wcycle;

    // ===============================

    int              *f_start;
    int              *f_end;

    real             *expand_fshift;
    real             *expand_Vvdw;
    real             *expand_Vc;

    real             *other_f;
    nbnxn_pairlist_t  other_nbl;
    nbnxn_atomdata_t  other_nbat;
    rvec             *other_shift_vec;
    real             *other_tabq_coul_F;
    real             *other_tabq_coul_V;
    real             *other_tabq_coul_FDV0;

    real             *fshift_p;
} _kernel_para;
int _kernel_launched = 0;

void nbnxn_kernel_ref_launch()
{
    int                nnbl;
    nbnxn_pairlist_t **nbl;
    int                coult;
    int                vdwt;
    int                nb;
    int                nthreads gmx_unused;

    nnbl = _kernel_para.nbl_list->nnbl;
    nbl  = _kernel_para.nbl_list->nbl;

    if (EEL_RF(_kernel_para.ic->eeltype) || _kernel_para.ic->eeltype == eelCUT)
    {
        coult = coultRF;
    }
    else
    {
        if (_kernel_para.ic->rcoulomb == _kernel_para.ic->rvdw)
        {
            coult = coultTAB;
        }
        else
        {
            coult = coultTAB_TWIN;
        }
    }

    if (_kernel_para.ic->vdwtype == evdwCUT)
    {
        switch (_kernel_para.ic->vdw_modifier)
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
    else if (_kernel_para.ic->vdwtype == evdwPME)
    {
        if (_kernel_para.ic->ljpme_comb_rule == ljcrGEOM)
        {
            assert(_kernel_para.nbat->comb_rule == ljcrGEOM);
            vdwt = vdwtEWALDGEOM;
        }
        else
        {
            assert(_kernel_para.nbat->comb_rule == ljcrLB);
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
        _kernel_para.f_start = (int*)malloc(64*sizeof(int));
    	memset(_kernel_para.f_start, 0, 64*sizeof(int));
	    _kernel_para.f_end = (int*)malloc(64*sizeof(int));
	    memset(_kernel_para.f_end, 0, 64*sizeof(int));
        //wallcycle_sub_start(wcycle, ewcsMEMCPY);
        subcore_loadbalance(_kernel_para.nbat, nbl[nb], _kernel_para.f_start, _kernel_para.f_end);
        //wallcycle_sub_stop(wcycle, ewcsMEMCPY);
        nbnxn_atomdata_output_t *out;

        out = &_kernel_para.nbat->out[nb];

        if (_kernel_para.clearF == enbvClearFYes)
        {
            clear_f(_kernel_para.nbat, nb, out->f);
        }

        if ((_kernel_para.force_flags & GMX_FORCE_VIRIAL) && nnbl == 1)
        {
            _kernel_para.fshift_p = _kernel_para.fshift;
        }
        else
        {
            _kernel_para.fshift_p = out->fshift;

            if (_kernel_para.clearF == enbvClearFYes)
            {
                clear_fshift(_kernel_para.fshift_p);
            }
        }
        //wallcycle_sub_start(wcycle, ewcsMEMCPY);
        _kernel_para.expand_fshift = (real*)malloc(SHIFTS*DIM*64*sizeof(real));

        _kernel_para.other_f       = (real*)malloc(
            _kernel_para.nbat->natoms*_kernel_para.nbat->fstride*sizeof(real));
        memcpy(_kernel_para.other_f, out->f, 
               _kernel_para.nbat->natoms*_kernel_para.nbat->fstride*sizeof(real));

        deep_copy_nbl(&_kernel_para.other_nbl, nbl[nb], 1, 0);
        deep_copy_nbat(&_kernel_para.other_nbat, _kernel_para.nbat, 1, 0);

        _kernel_para.other_shift_vec = (rvec*)malloc(SHIFTS*DIM*sizeof(real));
        memcpy(_kernel_para.other_shift_vec, _kernel_para.shift_vec, SHIFTS*DIM*sizeof(real));

        _kernel_para.other_tabq_coul_F = NULL;
        _kernel_para.other_tabq_coul_V = NULL;
        _kernel_para.other_tabq_coul_FDV0 = NULL;

#ifndef GMX_DOUBLE
        _kernel_para.other_tabq_coul_FDV0 = (real*)malloc(_kernel_para.ic->tabq_size*4*sizeof(real));
        memcpy(_kernel_para.other_tabq_coul_FDV0, _kernel_para.ic->tabq_coul_FDV0, 
               _kernel_para.ic->tabq_size*4*sizeof(real));
#else
        _kernel_para.other_tabq_coul_F    = (real*)malloc(_kernel_para.ic->tabq_size*sizeof(real));
        _kernel_para.other_tabq_coul_V    = (real*)malloc(_kernel_para.ic->tabq_size*sizeof(real));
        memcpy(_kernel_para.other_tabq_coul_F, _kernel_para.ic->tabq_coul_F, 
               _kernel_para.ic->tabq_size*sizeof(real));
        memcpy(_kernel_para.other_tabq_coul_V, _kernel_para.ic->tabq_coul_V, 
               _kernel_para.ic->tabq_size*sizeof(real));
#endif
        //wallcycle_sub_stop(wcycle, ewcsMEMCPY);
        // if the tabq_coul_FDV0 can load to LDM?
        // TLOG("tabq_coul_FDV0_SZ =%d Byte\n", ic->tabq_size*4*sizeof(real));
        // TLOG("tabq_size =%d, ntype =%d, natoms =%d\n", ic->tabq_size, nbat->ntype, nbat->natoms);
        if (!(_kernel_para.force_flags & GMX_FORCE_ENERGY))
        {
            host_out_param[PARAM_DEVICE_ACTION] = DEVICE_ACTION_RUN;
            host_out_param[FUNC_TYPE] = FUNC_NO_ENER;
            host_out_param[FUNC_I] = coult;
            host_out_param[FUNC_J] = vdwt;
            host_out_param[FUNC_PARAM_PTR] = (long)&host_func_para;
            host_func_para.nbl = &_kernel_para.other_nbl; // read only
            host_func_para.nbat = &_kernel_para.other_nbat; // read only
            host_func_para.ic = _kernel_para.ic;      // read only
            host_func_para.shift_vec = _kernel_para.other_shift_vec; // read only
            host_func_para.f = _kernel_para.other_f; // write only, reduce FIN
            host_func_para.expand_Vvdw = NULL;// write only, reduce FIN
            host_func_para.expand_Vc = NULL;// write only, reduce FIN
            host_func_para.expand_fshift = _kernel_para.expand_fshift;// write only, reduce FIN
            host_func_para.tabq_coul_F = _kernel_para.other_tabq_coul_F;// read only
            host_func_para.tabq_coul_V = _kernel_para.other_tabq_coul_V;// read only
            host_func_para.tabq_coul_FDV0 = _kernel_para.other_tabq_coul_FDV0;// read only
            host_func_para.f_start = _kernel_para.f_start;
            host_func_para.f_end   = _kernel_para.f_end;
            notice_device();
            //wait_device();
            _kernel_launched = 1;

            /* Don't calculate energies */
            //fake_device_run();
            //p_nbk_c_noener[coult][vdwt]();

        }
        else if (out->nV == 1)
        {
            /* No energy groups */
            out->Vvdw[0] = 0;
            out->Vc[0]   = 0;

            _kernel_para.expand_Vvdw   = (real*)malloc(64*sizeof(real));
            _kernel_para.expand_Vc     = (real*)malloc(64*sizeof(real));

            host_out_param[PARAM_DEVICE_ACTION] = DEVICE_ACTION_RUN;
            host_out_param[FUNC_TYPE] = FUNC_ENER;
            host_out_param[FUNC_I] = coult;
            host_out_param[FUNC_J] = vdwt;
            host_out_param[FUNC_PARAM_PTR] = (long)&host_func_para;
            host_func_para.nbl = &_kernel_para.other_nbl;
            host_func_para.nbat = &_kernel_para.other_nbat;
            host_func_para.ic = _kernel_para.ic;
            host_func_para.shift_vec = _kernel_para.other_shift_vec;
            host_func_para.f = _kernel_para.other_f;
            host_func_para.expand_Vvdw = _kernel_para.expand_Vvdw;
            host_func_para.expand_Vc = _kernel_para.expand_Vc;
            host_func_para.expand_fshift = _kernel_para.expand_fshift;
            host_func_para.tabq_coul_F = _kernel_para.other_tabq_coul_F;
            host_func_para.tabq_coul_V = _kernel_para.other_tabq_coul_V;
            host_func_para.tabq_coul_FDV0 = _kernel_para.other_tabq_coul_FDV0;
            host_func_para.f_start = _kernel_para.f_start;
            host_func_para.f_end   = _kernel_para.f_end;
            notice_device();
            //wait_device();
            _kernel_launched = 1;

            //fake_device_run();
            //p_nbk_c_ener[coult][vdwt]();
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

            _kernel_para.expand_Vvdw   = (real*)malloc(out->nV*64*sizeof(real));
            _kernel_para.expand_Vc     = (real*)malloc(out->nV*64*sizeof(real));

            host_out_param[PARAM_DEVICE_ACTION] = DEVICE_ACTION_RUN;
            host_out_param[FUNC_TYPE] = FUNC_ENERGRP;
            host_out_param[FUNC_I] = coult;
            host_out_param[FUNC_J] = vdwt;
            host_out_param[FUNC_PARAM_PTR] = (long)&host_func_para;
            host_func_para.nbl = &_kernel_para.other_nbl;
            host_func_para.nbat = &_kernel_para.other_nbat;
            host_func_para.ic = _kernel_para.ic;
            host_func_para.shift_vec = _kernel_para.other_shift_vec;
            host_func_para.f = _kernel_para.other_f;
            host_func_para.expand_Vvdw = _kernel_para.expand_Vvdw;
            host_func_para.expand_Vc = _kernel_para.expand_Vc;
            host_func_para.expand_fshift = _kernel_para.expand_fshift;
            host_func_para.tabq_coul_F = _kernel_para.other_tabq_coul_F;
            host_func_para.tabq_coul_V = _kernel_para.other_tabq_coul_V;
            host_func_para.tabq_coul_FDV0 = _kernel_para.other_tabq_coul_FDV0;
            host_func_para.f_start = _kernel_para.f_start;
            host_func_para.f_end   = _kernel_para.f_end;

            notice_device();
            //wait_device();
            _kernel_launched = 1;

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
        }
    }
}

void nbnxn_kernel_ref_reduce()
{
    if(_kernel_launched != 1)
    {
        TLOG("FORM reducer: ERROR! NO KERNEL LAUNCHED!");
        return;
    }

    int                nnbl;
    nbnxn_pairlist_t **nbl;
    int                coult;
    int                vdwt;
    int                nb;
    int                nthreads gmx_unused;

    nnbl = _kernel_para.nbl_list->nnbl;
    nbl  = _kernel_para.nbl_list->nbl;

    if (EEL_RF(_kernel_para.ic->eeltype) || _kernel_para.ic->eeltype == eelCUT)
    {
        coult = coultRF;
    }
    else
    {
        if (_kernel_para.ic->rcoulomb == _kernel_para.ic->rvdw)
        {
            coult = coultTAB;
        }
        else
        {
            coult = coultTAB_TWIN;
        }
    }

    if (_kernel_para.ic->vdwtype == evdwCUT)
    {
        switch (_kernel_para.ic->vdw_modifier)
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
    else if (_kernel_para.ic->vdwtype == evdwPME)
    {
        if (_kernel_para.ic->ljpme_comb_rule == ljcrGEOM)
        {
            assert(_kernel_para.nbat->comb_rule == ljcrGEOM);
            vdwt = vdwtEWALDGEOM;
        }
        else
        {
            assert(_kernel_para.nbat->comb_rule == ljcrLB);
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
        nbnxn_atomdata_output_t *out;

        out = &_kernel_para.nbat->out[nb];
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
        if (!(_kernel_para.force_flags & GMX_FORCE_ENERGY))
        {
            wait_device();

            /* Don't calculate energies */

        }
        else if (out->nV == 1)
        {
            /* No energy groups */

            wait_device();

            // reduce Vvdw
            int i;
            for(i = 0; i < 64 - 1; ++i)
            {
                _kernel_para.expand_Vvdw[i+1] += _kernel_para.expand_Vvdw[i];
            }
            out->Vvdw[0] = _kernel_para.expand_Vvdw[63];
            free(_kernel_para.expand_Vvdw);
            // reduce Vc
            for(i = 0; i < 64 - 1; ++i)
            {
                _kernel_para.expand_Vc[i+1] += _kernel_para.expand_Vc[i];
            }
            out->Vc[0] = _kernel_para.expand_Vc[63];
            free(_kernel_para.expand_Vc);
        }
        else
        {
            /* Calculate energy group contributions */
            int i;

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
                    _kernel_para.expand_Vvdw[(i+1)*out->nV+j] += _kernel_para.expand_Vvdw[(i)*out->nV+j];
                }
            }
            for (j = 0; j < out->nV; j++)
            {
                out->Vvdw[j] = _kernel_para.expand_Vvdw[(63)*out->nV+j];
            }
            free(_kernel_para.expand_Vvdw);

            // reduce Vc
            for(i = 0; i < 64 - 1; ++i)
            {
                for (j = 0; j < out->nV; j++)
                {
                    _kernel_para.expand_Vc[(i+1)*out->nV+j] += _kernel_para.expand_Vc[(i)*out->nV+j];
                }
            }
            for (j = 0; j < out->nV; j++)
            {
                out->Vc[j] = _kernel_para.expand_Vc[(63)*out->nV+j];
            }
            free(_kernel_para.expand_Vc);
        }

        // reduce fshift
        int i, j;
        for(i = 0; i < 64 - 1; ++i)
        {
            for(j = 0; j < SHIFTS*DIM; ++j)
            {
                _kernel_para.expand_fshift[(i+1)*SHIFTS*DIM+j] += _kernel_para.expand_fshift[(i)*SHIFTS*DIM+j];
            }
        }
        for(j = 0; j < SHIFTS*DIM; ++j)
        {
            _kernel_para.fshift_p[j] += _kernel_para.expand_fshift[63*SHIFTS*DIM+j];
        }
        

        memcpy(out->f, _kernel_para.other_f, 
               _kernel_para.nbat->natoms*_kernel_para.nbat->fstride*sizeof(real));
        deep_copy_nbl(nbl[nb], &_kernel_para.other_nbl, 0, 1);
        deep_copy_nbat(_kernel_para.nbat, &_kernel_para.other_nbat, 0, 1);

        free(_kernel_para.expand_fshift);
        free(_kernel_para.other_f);
        free(_kernel_para.other_shift_vec);
#ifndef GMX_DOUBLE
        free(_kernel_para.other_tabq_coul_FDV0);
#else
        free(_kernel_para.other_tabq_coul_F);
        free(_kernel_para.other_tabq_coul_V);
#endif
        free(_kernel_para.f_start);
        free(_kernel_para.f_end);
    }

    if (_kernel_para.force_flags & GMX_FORCE_ENERGY)
    {
        reduce_energies_over_lists(_kernel_para.nbat, nnbl, _kernel_para.Vvdw, _kernel_para.Vc);
    }

    _kernel_launched = 0;
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
    if(_kernel_launched == 1)
    {
        TLOG("FORM launcher: ERROR! THERE IS A KERNEL LAUNCHED.\n");
        return;
    }

    _kernel_para.nbl_list    = nbl_list;
    _kernel_para.nbat        = nbat;
    _kernel_para.ic          = ic;
    _kernel_para.shift_vec   = shift_vec;
    _kernel_para.force_flags = force_flags;
    _kernel_para.clearF      = clearF;
    _kernel_para.fshift      = fshift;
    _kernel_para.Vc          = Vc;
    _kernel_para.Vvdw        = Vvdw;
    _kernel_para.wcycle      = wcycle;

    nbnxn_kernel_ref_launch();
}
