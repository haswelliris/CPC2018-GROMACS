#include "nbnxn_kernel_ref_outer_fun.h"

void subcore_func(
		int 						macro_para,
		const nbnxn_pairlist_t		*nbl,
		const nbnxn_atomdata_t		*nbat,
		const interaction_const_t	*ic,
		rvec					 	*shift_vec,
		real					 	*f,
		real					 	*fshift,
		real					 	*Vvdw,
		real					 	*Vc,
		const nbnxn_ci_t 			*nbln,
		const nbnxn_cj_t 			*l_cj,
		const int					*type,
		const real					*q,
		const real					*shiftvec,
		const real					*x,
		const real					*nbfp,
		real						rcut2,
		real						rvdw2,
		int							ntype2,
		real						facel,
		real						*nbfp_i,
		int							n,
		int 						ci,
		int 						ci_sh,
		int							ish,
		int 						ishf,
		gmx_bool					do_LJ,
		gmx_bool 					half_LJ,
		gmx_bool 					do_coul,
		gmx_bool 					do_self,
		int							cjind0,
		int 						cjind1,
		int 						cjind,
		int							ip,
		int 						jp,

		real						*xi,
		real						*fi,
		real						*qi,

		real	   					Vvdw_ci,
		real						Vc_ci,

		int							egp_mask,
		int							*egp_sh_i,

		real						swV3,
		real 						swV4,
		real 						swV5,
		real						swF2,
		real 						swF3,
		real 						swF4,

		real						lje_coeff2,
		real 						lje_coeff6_6,
		real 						lje_vc,
		const 						real *ljc,

		real	   					k_rf2,

		real	   					k_rf,
		real						c_rf,

		real						tabscale,

		real						halfsp,

		#ifndef GMX_DOUBLE
			const real 				*tab_coul_FDV0,
		#else
			const real 				*tab_coul_F,
			const real 				*tab_coul_V,
		#endif
		int 						ninner

	) {
	for (n = nbl->nci -1; n >= 0; n--)
	{
		int i, d;

		nbln = &nbl->ci[n];

		ish			  = (nbln->shift & NBNXN_CI_SHIFT);
		/* x, f and fshift are assumed to be stored with stride 3 */
		ishf			 = ish*DIM;
		cjind0		   = nbln->cj_ind_start;
		cjind1		   = nbln->cj_ind_end;
		/* Currently only works super-cells equal to sub-cells */
		ci			   = nbln->ci;
		ci_sh			= (ish == CENTRAL ? ci : -1);

		/* We have 5 LJ/C combinations, but use only three inner loops,
		 * as the other combinations are unlikely and/or not much faster:
		 * inner half-LJ + C for half-LJ + C / no-LJ + C
		 * inner LJ + C	  for full-LJ + C
		 * inner LJ		  for full-LJ + no-C / half-LJ + no-C
		 */
		do_LJ   = (nbln->shift & NBNXN_CI_DO_LJ(0));
		do_coul = (nbln->shift & NBNXN_CI_DO_COUL(0));
		half_LJ = ((nbln->shift & NBNXN_CI_HALF_LJ(0)) || !do_LJ) && do_coul;

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
					egp_sh_i[i] = ((nbat->energrp[ci]>>(i*nbat->neg_2log)) & egp_mask)*nbat->nenergrp;
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


				if (l_cj[nbln->cj_ind_start].cj == ci_sh)
				{
					for (i = 0; i < UNROLLI; i++)
					{
						int egp_ind;

						if (macro_has(para_ENERGY_GROUPS))
							egp_ind = egp_sh_i[i] + ((nbat->energrp[ci]>>(i*nbat->neg_2log)) & egp_mask);
						else
							egp_ind = 0;

						/* Coulomb self interaction */
						Vc[egp_ind]   -= qi[i]*q[ci*UNROLLI+i]*Vc_sub_self;

						if (macro_has(para_LJ_EWALD)) {
							/* LJ Ewald self interaction */
							Vvdw[egp_ind] += 0.5*nbat->nbfp[nbat->type[ci*UNROLLI+i]*(nbat->ntype + 1)*2]/6*lje_coeff6_6;
						}
					}
				}
			}
		}

		cjind = cjind0;
		while (cjind < cjind1 && nbl->cj[cjind].excl != 0xffff)
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
		ninner += cjind1 - cjind0;

		/* Add accumulated i-forces to the force array */
		for (i = 0; i < UNROLLI; i++)
		{
			for (d = 0; d < DIM; d++)
			{
				f[(ci*UNROLLI+i)*F_STRIDE+d] += fi[i*FI_STRIDE+d];
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
					fshift[ishf+d] += fi[i*FI_STRIDE+d];
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
	}
}
