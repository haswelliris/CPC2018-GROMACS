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

/* We could use nbat->xstride and nbat->fstride, but macros might be faster */
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
(const nbnxn_pairlist_t     *nbl,
 const nbnxn_atomdata_t     *nbat,
 const interaction_const_t  *ic,
 rvec                       *shift_vec,
 real                       *f,
 real                       *fshift,
 real                       *Vvdw,
 real                       *Vc
)
{

    int macro_para = 0;
    enum {para_CALC_COUL_RF, para_CALC_COUL_TAB, para_CALC_ENERGIES, para_ENERGY_GROUPS, 
        para_LJ_CUT, para_LJ_EWALD, para_LJ_EWALD_COMB_GEOM, para_LJ_EWALD_COMB_LB, 
        para_LJ_FORCE_SWITCH, para_LJ_POT_SWITCH, para_VDW_CUTOFF_CHECK, 
        para_EXCL_FORCES, para_count
    };

    #ifdef CALC_COUL_RF
        macro_para |= 1 << para_CALC_COUL_RF;
    #endif
    #ifdef CALC_COUL_TAB
        macro_para |= 1 << para_CALC_COUL_TAB;
    #endif
    #ifdef CALC_ENERGIES
        macro_para |= 1 << para_CALC_ENERGIES;
    #endif
    #ifdef ENERGY_GROUPS
        macro_para |= 1 << para_ENERGY_GROUPS;
    #endif
    #ifdef LJ_CUT
        macro_para |= 1 << para_LJ_CUT;
    #endif
    #ifdef LJ_EWALD
        macro_para |= 1 << para_LJ_EWALD;
    #endif
    #ifdef LJ_EWALD_COMB_GEOM
        macro_para |= 1 << para_LJ_EWALD_COMB_GEOM;
    #endif
    #ifdef LJ_EWALD_COMB_LB
        macro_para |= 1 << para_LJ_EWALD_COMB_LB;
    #endif
    #ifdef LJ_FORCE_SWITCH
        macro_para |= 1 << para_LJ_FORCE_SWITCH;
    #endif
    #ifdef LJ_POT_SWITCH
        macro_para |= 1 << para_LJ_POT_SWITCH;
    #endif
    #ifdef VDW_CUTOFF_CHECK
        macro_para |= 1 << para_VDW_CUTOFF_CHECK;
    #endif

    #define macro_has(para_name) ((macro_para >> para_name) & 1)

    real fshift_host[64*SHIFTS*DIM];
    real Vvdw_host[64];
    real Vc_host[64];

    real *f_host = (real*)malloc(nbat->natoms*F_STRIDE*sizeof(real));
    memcpy(f_host, f, nbat->natoms*F_STRIDE*sizeof(real));
    
    workLoadPara_host.macro_para = macro_para;
    workLoadPara_host.nbl        = nbl;
    workLoadPara_host.nbat       = nbat;
    workLoadPara_host.ic         = ic;
    workLoadPara_host.shift_vec  = shift_vec;
    workLoadPara_host.f          = f_host;
    workLoadPara_host.fshift     = fshift;
    workLoadPara_host.fshift_host= fshift_host;
    workLoadPara_host.Vvdw       = Vvdw;
    workLoadPara_host.Vvdw_host  = Vvdw_host;
    workLoadPara_host.Vc         = Vc;
    workLoadPara_host.Vc_host    = Vc_host;

    int device_core_id;

    // fake subcore
    for (device_core_id = 0; device_core_id < 64; device_core_id++)
        subcore_func(&workLoadPara_host, device_core_id);

    // real subcore
    // host_param.host_to_device[WORKLOADPARA] = (long)&workLoadPara_host;
    // host_param.host_to_device[PARAM_DEVICE_ACTION] = DEVICE_ACTION_RUN;
    // notice_device();
    // wait_device();
    // printf("done!\n");
    // exit(0);

    memcpy(f, f_host, nbat->natoms*F_STRIDE*sizeof(real));

    int my_i;
    for (device_core_id = 0; device_core_id < 64-1; device_core_id++)
        for (my_i = 0; my_i < SHIFTS*DIM; my_i++)
            fshift_host[(device_core_id+1)*SHIFTS*DIM+my_i] += fshift_host[device_core_id*SHIFTS*DIM+my_i];
    for (my_i = 0; my_i < SHIFTS*DIM; my_i++)
        fshift[my_i] = fshift_host[63*SHIFTS*DIM+my_i];

    if (macro_has(para_CALC_ENERGIES)) {
        for (device_core_id = 0; device_core_id < 64; device_core_id++) {
            *Vvdw += Vvdw_host[device_core_id];
            *Vc += Vc_host[device_core_id];
        }
    }

#ifdef COUNT_PAIRS
    printf("atom pairs %d\n", npair);
#endif
}

#undef CALC_SHIFTFORCES

#undef X_STRIDE
#undef F_STRIDE
#undef XI_STRIDE
#undef FI_STRIDE

#undef UNROLLI
#undef UNROLLJ
