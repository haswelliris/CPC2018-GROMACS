#include "nbnxn_kernel_ref_outer_fun.h"
#include "SwConfig.h"

#ifdef HOST_RUN
	#include <string.h>
	#define aget_mem(A, B, C) memcpy(A, B, C)
	#define sget_mem(A, B, C) memcpy(A, B, C)
	#define put_mem(A, B, C) memcpy(A, B, C)
	#define alloc(A) malloc(A)
#else
	#include "SwDevice.h"
	#define aget_mem(A, B, C) async_get(A, B, C)
	#define sget_mem(A, B, C) sync_get(A, B, C)
	#define put_mem(A, B, C) sync_put(A, B, C)
	#define alloc(A) device_malloc(A)
	#define free(A) ldm_free(A)
#endif


// asign

int 						macro_para;
nbnxn_pairlist_t			nbl; // 会用到这个结构体里少量值
// nbl.nci = 2000
// nbl.ncj = 40000
nbnxn_atomdata_t			nbat; // 会用到这个结构体里比较多的东西
interaction_const_t			ic; // 用于大量初始化，inner中也有很多使用
rvec					 	*shift_vec; // 不需要传入，已提取到shiftvec
real					 	*f; // 会被更改的量，有规约！
#define F_LOCAL_SIZE		5000
real						f_local[F_LOCAL_SIZE];
real					 	*fshift; // 会被更改的量，有规约！
real						fshift_local[SHIFTS*DIM];
real					 	*Vvdw; // 会被更改的量，有规约！
real					 	*Vc; // 会被更改的量，有规约！

nbnxn_ci_t 					nbln; // 每次迭代值开始初始化，只在从核使用，使用其对象成员
nbnxn_cj_t 					l_cj; // 局部变量，初始化后在循环内不变，使用数组内容
const int					*type; // 初始化后在循环内不变，使用其数组内容
const real					*q; // 初始化后在循环内不变，使用其数组内容
const real					*shiftvec; // 初始化后在循环内不变，使用其数组内容
const real					*x; // 初始化后在循环内不变，使用其数组内容
const real					*nbfp; // 初始化后保持不变，使用其数组内容
real						rcut2; // 初始化后保持不变
real						rvdw2; // 初始化后保持不变
int							ntype2; // 初始化后保持不变
real						facel; // 初始化后保持不变
real						*nbfp_i; // 没有用到
int							n; // 循环变量
int 						ci; // 每次循环开始初始化
int 						ci_sh; // 每次循环开始初始化
int							ish; // 每次循环开始初始化
int 						ishf; // 每次循环开始初始化
gmx_bool					do_LJ; // 每次循环开始初始化
gmx_bool 					half_LJ; // 每次循环开始初始化
gmx_bool 					do_coul; // 每次循环开始初始化
gmx_bool 					do_self; // 每次循环开始初始化
int							cjind0; // 子循环用
int 						cjind1; // 子循环用
int 						cjind; // 子循环用
int							ip; // 好像没用
int 						jp; // 好像也没用

real						xi[UNROLLI*XI_STRIDE]; // 初始化后保持不变，貌似用于暂存
real						fi[UNROLLI*FI_STRIDE]; // 存放中间结果的变量
real						qi[UNROLLI]; // 初始化后保持不变

real	   					Vvdw_ci; // 每次循环开始时清零，最后用于加到规约变量Vvdw中
real						Vc_ci; // 每次循环开始时清零，最后用于加到规约变量Vc中

int							egp_mask; // 初始化后保持不变
int							egp_sh_i[UNROLLI]; // 每次循环开始初始化

real						swV3; // 初始化后保持不变
real 						swV4; // 初始化后保持不变
real 						swV5; // 初始化后保持不变
real						swF2; // 初始化后保持不变
real 						swF3; // 初始化后保持不变
real 						swF4; // 初始化后保持不变

real						lje_coeff2; // 初始化后保持不变
real 						lje_coeff6_6; // 初始化后保持不变
real 						lje_vc; // 初始化后保持不变
const 						real *ljc; // 初始化后保持不变，使用其数组内容 

real	   					k_rf2; // 初始化后保持不变

real	   					k_rf; // 初始化后保持不变
real						c_rf; // 初始化后保持不变

real						tabscale; // 没有用

real						halfsp; // 初始化后保持不变

#ifndef GMX_DOUBLE
	const real 				*tab_coul_FDV0; // 初始化后在循环内不变，使用其数组内容
#else
	const real 				*tab_coul_F; // 初始化后在循环内不变，使用其数组内容
	const real 				*tab_coul_V; // 初始化后在循环内不变，使用其数组内容
#endif

struct WorkLoadPara workLoadPara;

// 注意上从核这里要加extern
#ifdef HOST_RUN
void subcore_func(struct WorkLoadPara *workLoadPara_pass, int device_core_id)
#else
void subcore_func()
#endif
{
	sget_mem(&workLoadPara, workLoadPara_pass, sizeof(struct WorkLoadPara));

	macro_para = workLoadPara.macro_para;
	// nbl = workLoadPara.nbl;
	// nbat = workLoadPara.nbat;
	// ic = workLoadPara.ic;
	shift_vec = workLoadPara.shift_vec;
	f = workLoadPara.f;
	fshift = workLoadPara.fshift;
	Vvdw = workLoadPara.Vvdw;
	Vc = workLoadPara.Vc;

	int f_start = BLOCK_HEAD(device_core_id, 64, nbat.natoms/4)*12;
    int f_end = f_start + BLOCK_SIZE(device_core_id, 64, nbat.natoms/4)*12;
    if (f_end - f_start > F_LOCAL_SIZE)
    	printf("F_LOCAL_SIZE is not big enough!\n");

	// load obj

	aget_mem(&nbl, workLoadPara.nbl, sizeof(nbnxn_pairlist_t));
	aget_mem(&nbat, workLoadPara.nbat, sizeof(nbnxn_atomdata_t));
	aget_mem(&ic, workLoadPara.ic, sizeof(interaction_const_t));
	aget_mem(f_local, f+f_start, (f_end-f_start)*sizeof(real));

	memset(fshift_local, 0, sizeof(fshift_local));

	#ifndef HOST_RUN
		wait_all_async_get();
	#endif

	// end load obj

	// load sub obj

	

	// end load sub obj

	// init

	if (macro_has(para_LJ_POT_SWITCH)) {
	    swV3 = ic.vdw_switch.c3;
	    swV4 = ic.vdw_switch.c4;
	    swV5 = ic.vdw_switch.c5;
	    swF2 = 3*ic.vdw_switch.c3;
	    swF3 = 4*ic.vdw_switch.c4;
	    swF4 = 5*ic.vdw_switch.c5;
	}

	if (macro_has(para_LJ_EWALD)) {
	    lje_coeff2   = ic.ewaldcoeff_lj*ic.ewaldcoeff_lj;
	    lje_coeff6_6 = lje_coeff2*lje_coeff2*lje_coeff2/6.0;
	    lje_vc       = ic.sh_lj_ewald;

	    ljc          = nbat.nbfp_comb;
	}

	if (macro_has(para_CALC_COUL_RF)) {
	    k_rf2 = 2*ic.k_rf;
	    if (macro_has(para_CALC_ENERGIES)) {
		    k_rf = ic.k_rf;
		    c_rf = ic.c_rf;
		}
	}


	if (macro_has(para_CALC_COUL_TAB)) {
	    tabscale = ic.tabq_scale;
	    if (macro_has(para_CALC_ENERGIES)) {
	    	halfsp = 0.5/ic.tabq_scale;
		}

		#ifndef GMX_DOUBLE
	    	tab_coul_FDV0 = ic.tabq_coul_FDV0;
		#else
	    	tab_coul_F    = ic.tabq_coul_F;
	    	tab_coul_V    = ic.tabq_coul_V;
		#endif
	}

	if (macro_has(para_ENERGY_GROUPS)) {
    	egp_mask = (1<<nbat.neg_2log) - 1;
	}

    rcut2               = ic.rcoulomb*ic.rcoulomb;

    if (macro_has(para_VDW_CUTOFF_CHECK)) {
    	rvdw2           = ic.rvdw*ic.rvdw;
	}


    ntype2              = nbat.ntype*2;
    nbfp                = nbat.nbfp;
    q                   = nbat.q;
    type                = nbat.type;
    facel               = ic.epsfac;
    shiftvec            = shift_vec[0];
    x                   = nbat.x;

    // l_cj = nbl.cj;

	// end init

    // printf("%d\n", nbl.ncj);

	for (n = 0; n < nbl.nci; n++)
	// int task_num = BLOCK_SIZE(device_core_id, 64, nbl.nci);
	// for (n = BLOCK_HEAD(device_core_id, 64, nbl.nci); task_num; task_num--, n++)
	{
		int i, d;

		// nbln = &nbl.ci[n];
		sget_mem(&nbln, &nbl.ci[n], sizeof(nbnxn_ci_t));

		sget_mem(&l_cj, &nbl.cj[nbln.cj_ind_start], sizeof(nbnxn_cj_t));

		ish			  = (nbln.shift & NBNXN_CI_SHIFT);
		/* x, f and fshift are assumed to be stored with stride 3 */
		ishf			 = ish*DIM;
		cjind0		   = nbln.cj_ind_start;
		cjind1		   = nbln.cj_ind_end;
		/* Currently only works super-cells equal to sub-cells */
		ci			   = nbln.ci;
		ci_sh			= (ish == CENTRAL ? ci : -1);

		/* We have 5 LJ/C combinations, but use only three inner loops,
		 * as the other combinations are unlikely and/or not much faster:
		 * inner half-LJ + C for half-LJ + C / no-LJ + C
		 * inner LJ + C	  for full-LJ + C
		 * inner LJ		  for full-LJ + no-C / half-LJ + no-C
		 */
		do_LJ   = (nbln.shift & NBNXN_CI_DO_LJ(0));
		do_coul = (nbln.shift & NBNXN_CI_DO_COUL(0));
		half_LJ = ((nbln.shift & NBNXN_CI_HALF_LJ(0)) || !do_LJ) && do_coul;

		if (macro_has(para_LJ_EWALD))
			do_self = TRUE;
		else 
			do_self = do_coul;


		if (macro_has(para_CALC_ENERGIES)) {
			if (!macro_has(para_ENERGY_GROUPS)) {
				Vvdw_ci = 0;
				Vc_ci   = 0;
			}
			else {
				for (i = 0; i < UNROLLI; i++)
				{
					egp_sh_i[i] = ((nbat.energrp[ci]>>(i*nbat.neg_2log)) & egp_mask)*nbat.nenergrp;
				}
			}
		}

		for (i = 0; i < UNROLLI; i++)
		{
			for (d = 0; d < DIM; d++)
			{
				xi[i*XI_STRIDE+d] = x[(ci*UNROLLI+i)*X_STRIDE+d] + shiftvec[ishf+d];
				fi[i*FI_STRIDE+d] = 0;
			}

			qi[i] = facel*q[ci*UNROLLI+i];
		}

		if (macro_has(para_CALC_ENERGIES)) {
			if (do_self)
			{
				real Vc_sub_self;

				if (macro_has(para_CALC_COUL_RF))
					Vc_sub_self = 0.5*c_rf;

				if (macro_has(para_CALC_COUL_TAB)) {
					#ifdef GMX_DOUBLE
								Vc_sub_self = 0.5*tab_coul_V[0];
					#else
								Vc_sub_self = 0.5*tab_coul_FDV0[2];
					#endif
				}


				if (l_cj.cj == ci_sh)
				{
					for (i = 0; i < UNROLLI; i++)
					{
						int egp_ind;

						if (macro_has(para_ENERGY_GROUPS))
							egp_ind = egp_sh_i[i] + ((nbat.energrp[ci]>>(i*nbat.neg_2log)) & egp_mask);
						else
							egp_ind = 0;

						/* Coulomb self interaction */
						if (BLOCK_HINT(ci*UNROLLI*F_STRIDE, f_start, f_end))
							Vc[egp_ind]   -= qi[i]*q[ci*UNROLLI+i]*Vc_sub_self;

						if (macro_has(para_LJ_EWALD)) {
							/* LJ Ewald self interaction */
							Vvdw[egp_ind] += 0.5*nbat.nbfp[nbat.type[ci*UNROLLI+i]*(nbat.ntype + 1)*2]/6*lje_coeff6_6;
						}
					}
				}
			}
		}

		cjind = cjind0;
		while (cjind < cjind1 && nbl.cj[cjind].excl != 0xffff)
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
			sget_mem(&l_cj, &nbl.cj[cjind], sizeof(nbnxn_cj_t));
		}

		for (; (cjind < cjind1); cjind++)
		{
			sget_mem(&l_cj, &nbl.cj[cjind], sizeof(nbnxn_cj_t));
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
		// ninner += cjind1 - cjind0;
		if (BLOCK_HINT(ci*UNROLLI*F_STRIDE, f_start, f_end)) {
		/* Add accumulated i-forces to the force array */
		for (i = 0; i < UNROLLI; i++)
		{
			for (d = 0; d < DIM; d++)
			{
				f_local[(ci*UNROLLI+i)*F_STRIDE+d-f_start] += fi[i*FI_STRIDE+d];
			}
		}
		// #ifdef CALC_SHIFTFORCES
		if (fshift != NULL)
		{
			/* Add i forces to shifted force list */
			for (i = 0; i < UNROLLI; i++)
			{
				for (d = 0; d < DIM; d++)
				{
					fshift_local[ishf+d] += fi[i*FI_STRIDE+d];
				}
			}
		}
		// #endif

		if (macro_has(para_CALC_ENERGIES)) {
			if (!macro_has(para_ENERGY_GROUPS)) {
				*Vvdw += Vvdw_ci;
				*Vc   += Vc_ci;
			}
		}
		} // end of ci write test
	}
	put_mem(f+f_start, f_local, (f_end-f_start)*sizeof(real));
	put_mem(workLoadPara.fshift_host + device_core_id*SHIFTS*DIM, fshift_local, sizeof(fshift_local));
}
