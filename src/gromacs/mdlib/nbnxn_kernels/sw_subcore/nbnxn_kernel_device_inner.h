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

/* When calculating RF or Ewald interactions we calculate the electrostatic
 * forces and energies on excluded atom pairs here in the non-bonded loops.
 */
#if defined CHECK_EXCLS && (defined CALC_COULOMB || defined LJ_EWALD)
#define EXCL_FORCES
#endif

{
    int cj;
    int i;
    int write_cj;

    //TODO: ldm load: l_cj
    //cj               = l_cj[cjind].cj;
    cj               = Ccj_p->cj;
    write_cj         = IN_F_BLOCK(cj);

    if(write_ci || write_cj)
    {

    Cxj_p = xj_C(cj);
    Cqj_p = qj_C(cj);
    Ctj_p = tj_C(cj);
    for (i = 0; i < UNROLLI; i++)
    {
        int ai;
        int type_i_off;
        int j;

        //ai = ci*UNROLLI + i;

        //TODO: ldm load: type
        //type_i_off = type[ai]*ntype2;
        type_i_off = Cti_p[i]*ntype2;

#ifdef SIMD_INNER // ============================================
        int aj[4];
        realv4 dxV, dyV, dzV;
        realv4 rsqV, rinvV;
        realv4 rinvsqV, rinvsixV;
        realv4 c6V, c12V;
        realv4 FrLJ6V, FrLJ12V, frLJV, VLJV;
        FrLJ6V.v = 0.0, FrLJ12V.v = 0.0, frLJV.v = 0.0, VLJV.v = 0.0;

#ifdef CALC_COULOMB
        realv4 qqV;
        realv4 fcoulV;
#ifdef CALC_COUL_TAB
        realv4 rsV, fracV;
        int  ri[4];
        realv4 fexclV;
#endif
#ifdef CALC_ENERGIES
        realv4 vcoulV;
#endif
#endif
        realv4 fscalV;
        realv4 fxV, fyV, fzV;

        realv4 skipmaskV;

#ifdef CHECK_EXCLS
        realv4 interactV;
        //DEVICE_CODE_FENCE();
        interactV.p[0] = (real)((Ccj_p->excl>>(i*UNROLLI + 0)) & 1);
        interactV.p[1] = (real)((Ccj_p->excl>>(i*UNROLLI + 1)) & 1);
        interactV.p[2] = (real)((Ccj_p->excl>>(i*UNROLLI + 2)) & 1);
        interactV.p[3] = (real)((Ccj_p->excl>>(i*UNROLLI + 3)) & 1);
#ifndef EXCL_FORCES
        skipmaskV.v = interactV.v;
#else   
        //DEVICE_CODE_FENCE();
        skipmaskV.p[0] = (real)!(cj == ci_sh && 0 <= i);
        skipmaskV.p[1] = (real)!(cj == ci_sh && 1 <= i);
        skipmaskV.p[2] = (real)!(cj == ci_sh && 2 <= i);
        skipmaskV.p[3] = (real)!(cj == ci_sh && 3 <= i);
#endif  // EXCL_FORCES
#else
#define interactV oneV
        skipmaskV.v = 1.0;
#endif  // CHECK_EXCLS
        realv4 xiXV;
        realv4 xiYV;
        realv4 xiZV;
        xiXV.v = xi[i*XI_STRIDE+XX];
        xiYV.v = xi[i*XI_STRIDE+YY];
        xiZV.v = xi[i*XI_STRIDE+ZZ];

        realv4 xjXV, xjYV, xjZV;
        //DEVICE_CODE_FENCE();
        xjXV.p[0] = Cxj_p[0*X_STRIDE+XX];
        xjYV.p[0] = Cxj_p[0*X_STRIDE+YY];
        xjZV.p[0] = Cxj_p[0*X_STRIDE+ZZ];
        xjXV.p[1] = Cxj_p[1*X_STRIDE+XX];
        xjYV.p[1] = Cxj_p[1*X_STRIDE+YY];
        xjZV.p[1] = Cxj_p[1*X_STRIDE+ZZ];
        xjXV.p[2] = Cxj_p[2*X_STRIDE+XX];
        xjYV.p[2] = Cxj_p[2*X_STRIDE+YY];
        xjZV.p[2] = Cxj_p[2*X_STRIDE+ZZ];
        xjXV.p[3] = Cxj_p[3*X_STRIDE+XX];
        xjYV.p[3] = Cxj_p[3*X_STRIDE+YY];
        xjZV.p[3] = Cxj_p[3*X_STRIDE+ZZ];

        aj[0] = cj*UNROLLJ + 0;
        aj[1] = cj*UNROLLJ + 1;
        aj[2] = cj*UNROLLJ + 2;
        aj[3] = cj*UNROLLJ + 3;

        dxV.v  = xiXV.v - xjXV.v;
        dyV.v  = xiYV.v - xjYV.v;
        dzV.v  = xiZV.v - xjZV.v;

        rsqV.v = dxV.v*dxV.v + dyV.v*dyV.v + dzV.v*dzV.v;

        skipmaskV.v = simd_vsellt(rsqV.v - rcut2, skipmaskV.v, zeroV.v);
        //skipmask = (rsq -  rcut2 < 0) ? skipmask : 0;
        //skipmask = (rsq >= rcut2 + 0) ? 0 : skipmask;

#ifdef CHECK_EXCLS
        rsqV.v = rsqV.v + (1.0-interactV.v)*NBNXN_AVOID_SING_R2_INC;
#endif

        rinvV.v = 1.0/simd_vsqrt(rsqV.v);
        rinvV.v = rinvV.v*skipmaskV.v;
        rinvsqV.v  = rinvV.v*rinvV.v;

        //DEVICE_CODE_FENCE();
#ifdef HALF_LJ
        if (i < UNROLLI/2)
#endif
        {
            // c6      = nbfp[type_i_off+Ctj_p[j]*2  ];
            // c12     = nbfp[type_i_off+Ctj_p[j]*2+1];
            c6V.p [0]   = nbfp[type_i_off+Ctj_p[0]*2  ];
            c12V.p[0]   = nbfp[type_i_off+Ctj_p[0]*2+1];
            c6V.p [1]   = nbfp[type_i_off+Ctj_p[1]*2  ];
            c12V.p[1]   = nbfp[type_i_off+Ctj_p[1]*2+1];
            c6V.p [2]   = nbfp[type_i_off+Ctj_p[2]*2  ];
            c12V.p[2]   = nbfp[type_i_off+Ctj_p[2]*2+1];
            c6V.p [3]   = nbfp[type_i_off+Ctj_p[3]*2  ];
            c12V.p[3]   = nbfp[type_i_off+Ctj_p[3]*2+1];

#if defined LJ_CUT
            rinvsixV.v = interactV.v*rinvsqV.v*rinvsqV.v*rinvsqV.v;
            FrLJ6V.v   = c6V.v*rinvsixV.v;
            FrLJ12V.v  = c12V.v*rinvsixV.v*rinvsixV.v;
            frLJV.v    = FrLJ12V.v - FrLJ6V.v;
#if defined CALC_ENERGIES
            VLJV.v     = (FrLJ12V.v + c12V.v*cpot)/12 - (FrLJ6V.v + c6V.v*cpot)/6;
#endif
#endif

#if defined CALC_ENERGIES
            VLJV.v     = VLJV.v * interactV.v;
            VLJV.v     = VLJV.v * skipmaskV.v;
#endif
#ifdef CALC_ENERGIES
            Vvdw_ci += (VLJV.p[0] + VLJV.p[1] + VLJV.p[2] + VLJV.p[3]);
#endif
        }

#ifdef CALC_COULOMB 
        //qq = skipmask * qi[i] * Cqj_p[j];
        realv4 qiV;
        qiV.v = qi[i];
        realv4 qjV;
        simd_load(qjV.v, &Cqj_p[0]);
        qqV.v = skipmaskV.v * qiV.v * qjV.v;
#ifdef CALC_COUL_RF
        fcoulV.v  = qqV.v*(interactV.v*rinvV.v*rinvsqV.v - k_rf2);
#ifdef CALC_ENERGIES
        vcoulV.v  = qqV.v*(interactV.v*rinvV.v + k_rf*rsqV.v - c_rf);
#endif // CALC_ENERGIES
#endif // CALC_COUL_RF

#ifdef CALC_COUL_TAB
        rsV.v    = rsqV.v*rinvV.v*ic.tabq_scale;

        //DEVICE_CODE_FENCE();
        ri[0]  = (int)rsV.p[0];
        ri[1]  = (int)rsV.p[1];
        ri[2]  = (int)rsV.p[2];
        ri[3]  = (int)rsV.p[3];
        
        //DEVICE_CODE_FENCE();
        fracV.p[0]   = rsV.p[0] - (float)ri[0];
        fracV.p[1]   = rsV.p[1] - (float)ri[1];
        fracV.p[2]   = rsV.p[2] - (float)ri[2];
        fracV.p[3]   = rsV.p[3] - (float)ri[3];
#ifndef GMX_DOUBLE
        realv4 tab_coul_FDV0_FV, tab_coul_FDV0_DV;
#ifdef CALC_ENERGIES
        realv4 tab_coul_FDV0_VV;
#endif
        //DEVICE_CODE_FENCE();
        tab_coul_FDV0_FV.p[0] = tab_coul_FDV0[ri[0]*4];
        tab_coul_FDV0_DV.p[0] = tab_coul_FDV0[ri[0]*4+1];
#ifdef CALC_ENERGIES
        tab_coul_FDV0_VV.p[0] = tab_coul_FDV0[ri[0]*4+2];
#endif
        tab_coul_FDV0_FV.p[1] = tab_coul_FDV0[ri[1]*4];
        tab_coul_FDV0_DV.p[1] = tab_coul_FDV0[ri[1]*4+1];
#ifdef CALC_ENERGIES
        tab_coul_FDV0_VV.p[1] = tab_coul_FDV0[ri[1]*4+2];
#endif
        tab_coul_FDV0_FV.p[2] = tab_coul_FDV0[ri[2]*4];
        tab_coul_FDV0_DV.p[2] = tab_coul_FDV0[ri[2]*4+1];
#ifdef CALC_ENERGIES
        tab_coul_FDV0_VV.p[2] = tab_coul_FDV0[ri[2]*4+2];
#endif
        tab_coul_FDV0_FV.p[3] = tab_coul_FDV0[ri[3]*4];
        tab_coul_FDV0_DV.p[3] = tab_coul_FDV0[ri[3]*4+1];
#ifdef CALC_ENERGIES
        tab_coul_FDV0_VV.p[3] = tab_coul_FDV0[ri[3]*4+2];
#endif

        fexclV.v  = tab_coul_FDV0_FV.v + fracV.v*tab_coul_FDV0_DV.v;
        //fexcl  = tab_coul_FDV0[ri*4] + frac*tab_coul_FDV0[ri*4+1];
        //fexcl  = (1 - frac)*tab_coul_FDV0[ri*2] + frac*tab_coul_FDV0[(ri+1)*2];
#else
        realv4 tab_coul_F0V, tab_coul_F1V;
        //DEVICE_CODE_FENCE();
        tab_coul_F0V.p[0] = tab_coul_F[ri[0]];
        tab_coul_F1V.p[0] = tab_coul_F[ri[0]+1];
        tab_coul_F0V.p[1] = tab_coul_F[ri[1]];
        tab_coul_F1V.p[1] = tab_coul_F[ri[1]+1];
        tab_coul_F0V.p[2] = tab_coul_F[ri[2]];
        tab_coul_F1V.p[2] = tab_coul_F[ri[2]+1];
        tab_coul_F0V.p[3] = tab_coul_F[ri[3]];
        tab_coul_F1V.p[3] = tab_coul_F[ri[3]+1];
        /* fexcl = (1-frac) * F_i + frac * F_(i+1) */
        //fexcl  = (1 - frac)*tab_coul_F[ri] + frac*tab_coul_F[ri+1];
        fexclV  = (1 - fracV)*tab_coul_F0V + fracV*tab_coul_F1V;
#endif // GMX_DOUBLE
        fcoulV.v  = interactV.v*rinvsqV.v - fexclV.v;
#ifdef CALC_ENERGIES
#ifndef GMX_DOUBLE
        vcoulV.v  = qqV.v*(interactV.v*(rinvV.v - ic.sh_ewald)
                         -(tab_coul_FDV0_VV.v
                           -halfsp*fracV.v*(tab_coul_FDV0_FV.v + fexclV.v)));
        // vcoul  = qq*(interact*(rinv - ic.sh_ewald)
        //              -(tab_coul_FDV0[ri*2+1]
        //                -halfsp*frac*(tab_coul_FDV0[ri*2] + fexcl)));
        /* 7 flops for float 1/r-table energy (8 with excls) */
#else
        realv4 tab_coul_VV;
        //DEVICE_CODE_FENCE();
        tab_coul_FVV.p[0] = tab_coul_V[ri[0]];
        tab_coul_FVV.p[1] = tab_coul_V[ri[1]];
        tab_coul_FVV.p[2] = tab_coul_V[ri[2]];
        tab_coul_FVV.p[3] = tab_coul_V[ri[3]];
        vcoulV.v  = qqV.v*(interactV.v*(rinvV.v - ic.sh_ewald)
                         -(tab_coul_VV.v
                           -halfsp*fracV.v*(tab_coul_F0V.v + fexclV.v)));
#endif // GMX_DOUBLE
#endif // CALC_ENERGIES
        fcoulV.v = fcoulV.v*qqV.v*rinvV.v;
#endif // CALC_COUL_TAB
#ifdef CALC_ENERGIES
        Vc_ci += (vcoulV.p[0] + vcoulV.p[1] + vcoulV.p[2] + vcoulV.p[3]);
        /* 1 flop for Coulomb energy addition */
#endif // CALC_ENERGIES
#endif // CALC_COULOMB

#ifdef CALC_COULOMB
#ifdef HALF_LJ
        if (i < UNROLLI/2)
#endif
        {
            fscalV.v = frLJV.v*rinvsqV.v + fcoulV.v;
                /* 2 flops for scalar LJ+Coulomb force */
        }
#ifdef HALF_LJ
        else
        {
            fscalV.v = fcoulV.v;
        }
#endif
#else
        fscalV.v = frLJV.v*rinvsqV.v;
#endif
        fxV.v = fscalV.v*dxV.v;
        fyV.v = fscalV.v*dyV.v;
        fzV.v = fscalV.v*dzV.v;
        
        //DEVICE_CODE_FENCE();
        if(write_ci) {
            real fx = (fxV.p[0] + fxV.p[1] + fxV.p[2], fxV.p[3]);
            real fy = (fxV.p[0] + fxV.p[1] + fxV.p[2], fxV.p[3]);
            real fz = (fxV.p[0] + fxV.p[1] + fxV.p[2], fxV.p[3]);
            fi[i*FI_STRIDE+XX] += fx;
            fi[i*FI_STRIDE+YY] += fy;
            fi[i*FI_STRIDE+ZZ] += fz;
        }

        if(write_cj) {
            ldm_f[aj[0]*F_STRIDE+XX-start_f]  -= fxV.p[0];
            ldm_f[aj[0]*F_STRIDE+YY-start_f]  -= fyV.p[0];
            ldm_f[aj[0]*F_STRIDE+ZZ-start_f]  -= fzV.p[0];
            ldm_f[aj[1]*F_STRIDE+XX-start_f]  -= fxV.p[1];
            ldm_f[aj[1]*F_STRIDE+YY-start_f]  -= fyV.p[1];
            ldm_f[aj[1]*F_STRIDE+ZZ-start_f]  -= fzV.p[1];
            ldm_f[aj[2]*F_STRIDE+XX-start_f]  -= fxV.p[2];
            ldm_f[aj[2]*F_STRIDE+YY-start_f]  -= fyV.p[2];
            ldm_f[aj[2]*F_STRIDE+ZZ-start_f]  -= fzV.p[2];
            ldm_f[aj[3]*F_STRIDE+XX-start_f]  -= fxV.p[3];
            ldm_f[aj[3]*F_STRIDE+YY-start_f]  -= fyV.p[3];
            ldm_f[aj[3]*F_STRIDE+ZZ-start_f]  -= fzV.p[3];
        }
#undef interactV
#else // SIMD_INNER =======================================
        for (j = 0; j < UNROLLJ; j++)
        {
            int  aj;
            real dx, dy, dz;
            real rsq, rinv;
            real rinvsq, rinvsix;
            real c6, c12;
            real FrLJ6 = 0, FrLJ12 = 0, frLJ = 0, VLJ = 0;

#ifdef CALC_COULOMB
            real qq;
            real fcoul;
#ifdef CALC_COUL_TAB
            real rs, frac;
            int  ri;
            real fexcl;
#endif
#ifdef CALC_ENERGIES
            real vcoul;
#endif
#endif
            real fscal;
            real fx, fy, fz;

            /* A multiply mask used to zero an interaction
             * when either the distance cutoff is exceeded, or
             * (if appropriate) the i and j indices are
             * unsuitable for this kind of inner loop. */
            real skipmask;

#ifdef CHECK_EXCLS
            /* A multiply mask used to zero an interaction
             * when that interaction should be excluded
             * (e.g. because of bonding). */
            int interact;

            //interact = ((l_cj[cjind].excl>>(i*UNROLLI + j)) & 1);
            interact = ((Ccj_p->excl>>(i*UNROLLI + j)) & 1);
#ifndef EXCL_FORCES
            skipmask = interact;
#else   
            skipmask = !(cj == ci_sh && j <= i);
#endif  // EXCL_FORCES
#else
#define interact 1.0
            skipmask = 1.0;
#endif  // CHECK_EXCLS
            ////DEVICE_CODE_FENCE();
#ifdef DEBUG_SDLB
            TLOG("kaCHI 7.1.\n");
#endif
            aj = cj*UNROLLJ + j;

            // dx  = xi[i*XI_STRIDE+XX] - x[aj*X_STRIDE+XX];
            // dy  = xi[i*XI_STRIDE+YY] - x[aj*X_STRIDE+YY];
            // dz  = xi[i*XI_STRIDE+ZZ] - x[aj*X_STRIDE+ZZ];
            dx  = xi[i*XI_STRIDE+XX] - Cxj_p[j*X_STRIDE+XX];
            dy  = xi[i*XI_STRIDE+YY] - Cxj_p[j*X_STRIDE+YY];
            dz  = xi[i*XI_STRIDE+ZZ] - Cxj_p[j*X_STRIDE+ZZ];
#ifdef DEBUG_CACHE
            if(x[aj*X_STRIDE+XX] != Cxj_p[j*X_STRIDE+XX])
            {
                TLOG("KAAAA! Cache ERR: cj =%d\n", cj);
            }
            if(x[aj*X_STRIDE+YY] != Cxj_p[j*X_STRIDE+YY])
            {
                TLOG("KAAAA! Cache ERR: cj =%d\n", cj);
            }
            if(x[aj*X_STRIDE+ZZ] != Cxj_p[j*X_STRIDE+ZZ])
            {
                TLOG("KAAAA! Cache ERR: cj =%d\n", cj);
            }
#endif

            rsq = dx*dx + dy*dy + dz*dz;

#ifdef DEBUG_FPEX
            TLOG("xiX =%f, xiY =%f, xiZ =%f\n", xi[i*XI_STRIDE+XX], xi[i*XI_STRIDE+YY], xi[i*XI_STRIDE+ZZ]);
            TLOG("x_X =%f, x_Y =%f, x_Z =%f\n", x[aj*X_STRIDE+XX], x[aj*X_STRIDE+YY], x[aj*X_STRIDE+ZZ]);
            TLOG("dx =%f, dy =%f, dz =%f\n", dx, dy, dz);
            TLOG("rsq =%f, rcut2 =%f\n", rsq, rcut2);
#endif

            /* Prepare to enforce the cut-off. */
            skipmask = (rsq >= rcut2) ? 0 : skipmask;
            /* 9 flops for r^2 + cut-off check */
            ////DEVICE_CODE_FENCE();
#ifdef DEBUG_SDLB
            TLOG("kaCHI 7.2.\n");
            //wait_host(device_core_id);
#endif
#ifdef CHECK_EXCLS
            /* Excluded atoms are allowed to be on top of each other.
             * To avoid overflow of rinv, rinvsq and rinvsix
             * we add a small number to rsq for excluded pairs only.
             */
            rsq += (1 - interact)*NBNXN_AVOID_SING_R2_INC;
#endif

            rinv = gmx_invsqrt(rsq);
            /* 5 flops for invsqrt */

            /* Partially enforce the cut-off (and perhaps
             * exclusions) to avoid possible overflow of
             * rinvsix when computing LJ, and/or overflowing
             * the Coulomb table during lookup. */
            rinv = rinv * skipmask;

            rinvsq  = rinv*rinv;

#ifdef HALF_LJ
            if (i < UNROLLI/2)
#endif
            {
                //TODO: ldm load: nbfp
                //c6      = nbfp[type_i_off+type[aj]*2  ];
                //c12     = nbfp[type_i_off+type[aj]*2+1];
                /* SAMPLE 1:  NTYPE = 39, 39*2 is the MIN fetch size */
                /* SAMPLE 2:  NTYPE = 12, 12*2 is the MIN fetch size */
                /* OUR CACHE SIZE =  39*2*12*sizeof(int) = 3744 Byte */
                c6      = nbfp[type_i_off+Ctj_p[j]*2  ];
                c12     = nbfp[type_i_off+Ctj_p[j]*2+1];

#if defined LJ_CUT
                rinvsix = interact*rinvsq*rinvsq*rinvsq;
                FrLJ6   = c6*rinvsix;
                FrLJ12  = c12*rinvsix*rinvsix;
                frLJ    = FrLJ12 - FrLJ6;
                /* 7 flops for r^-2 + LJ force */
#if defined CALC_ENERGIES
                VLJ     = (FrLJ12 + c12*cpot)/12 -
                    (FrLJ6 + c6*cpot)/6;
                /* 7 flops for LJ energy */
#endif
#endif
                ////DEVICE_CODE_FENCE();
#ifdef DEBUG_SDLB
                TLOG("kaCHI 7.2.1.\n");
                //wait_host(device_core_id);
#endif

#if defined CALC_ENERGIES
                /* Masking should be done after force switching,
                 * but before potential switching.
                 */
                /* Need to zero the interaction if there should be exclusion. */
                VLJ     = VLJ * interact;
#endif

#if defined CALC_ENERGIES
                /* Need to zero the interaction if r >= rcut */
                VLJ     = VLJ * skipmask;
                /* 1 more flop for LJ energy */
#endif

                ////DEVICE_CODE_FENCE();
#ifdef CALC_ENERGIES
                Vvdw_ci += VLJ;
                /* 1 flop for LJ energy addition */
#endif
            }

#ifdef CALC_COULOMB
            /* Enforce the cut-off and perhaps exclusions. In
             * those cases, rinv is zero because of skipmask,
             * but fcoul and vcoul will later be non-zero (in
             * both RF and table cases) because of the
             * contributions that do not depend on rinv. These
             * contributions cannot be allowed to accumulate
             * to the force and potential, and the easiest way
             * to do this is to zero the charges in
             * advance. */
            //qq = skipmask * qi[i] * q[aj];
            qq = skipmask * qi[i] * Cqj_p[j];

#ifdef CALC_COUL_RF
            fcoul  = qq*(interact*rinv*rinvsq - k_rf2);
            /* 4 flops for RF force */
#ifdef CALC_ENERGIES
            vcoul  = qq*(interact*rinv + k_rf*rsq - c_rf);
            /* 4 flops for RF energy */
#endif // CALC_ENERGIES
#endif // CALC_COUL_RF

#ifdef CALC_COUL_TAB
            ////DEVICE_CODE_FENCE();
#ifdef DEBUG_FPEX
            TLOG("rsq =%f, tabq_scale =%f\n", rsq, ic.tabq_scale);
#endif //DEBUG_FPEX
            rs     = rsq*rinv*ic.tabq_scale;
            ri     = (int)rs;
            frac   = rs - ri;
#ifndef GMX_DOUBLE
            /* fexcl = F_i + frac * (F_(i+1)-F_i) */
            // ri =(int)(dx^2+dy^2+dz^2) / sqrt(dx^2+dy^2+dz^2) * ic.tabq_scale
            // (int)(dx^2+dy^2+dz^2) / sqrt(dx^2+dy^2+dz^2) = [0, +inf]??
            // SAMPLE1 tabq_scale=1071
            // SAMPLE2 WILL NOT INTO CALC_COUL_TAB
            // RI =[0,tabq_scale]
            // RI*4 =[0,tabq_scale*4]  =  [0, 4300]
            // 4500 * sizeof(real) = 17200 Byte
            // SAMPLE sizeof(tab_coul_FDV0) = 17168 Byte
            // like a shit! cache miss is very often!!!!
            //
            // Lets's LOAD all of it to LDM???? 17KB is acceptable
            // TLOG("RI =%d\t tabq_scale =%f\n", ri, ic.tabq_scale);
            fexcl  = tab_coul_FDV0[ri*4] + frac*tab_coul_FDV0[ri*4+1];
            //fexcl  = (1 - frac)*tab_coul_FDV0[ri*2] + frac*tab_coul_FDV0[(ri+1)*2];
#else
            /* fexcl = (1-frac) * F_i + frac * F_(i+1) */
            fexcl  = (1 - frac)*tab_coul_F[ri] + frac*tab_coul_F[ri+1];
#endif // GMX_DOUBLE
            fcoul  = interact*rinvsq - fexcl;
            /* 7 flops for float 1/r-table force */
#ifdef CALC_ENERGIES
            ////DEVICE_CODE_FENCE();
#ifndef GMX_DOUBLE
            //TODO: ldm load: tab_coul_FDV0, tab_coul_V, tab_coul_F
#ifdef DEBUG_FPEX
            TLOG("qq =%f, rinv =%f, interact =%f, sh_ewald =%f, halfsp =%f, frac =%f, fexcl =%f\n", qq, Vc_ci, interact, ic.sh_ewald, halfsp, frac, fexcl);
#endif // DEBUG_FPEX
            vcoul  = qq*(interact*(rinv - ic.sh_ewald)
                         -(tab_coul_FDV0[ri*4+2]
                           -halfsp*frac*(tab_coul_FDV0[ri*4] + fexcl)));
            // vcoul  = qq*(interact*(rinv - ic.sh_ewald)
            //              -(tab_coul_FDV0[ri*2+1]
            //                -halfsp*frac*(tab_coul_FDV0[ri*2] + fexcl)));
            /* 7 flops for float 1/r-table energy (8 with excls) */
#else
            vcoul  = qq*(interact*(rinv - ic.sh_ewald)
                         -(tab_coul_V[ri]
                           -halfsp*frac*(tab_coul_F[ri] + fexcl)));
#endif // GMX_DOUBLE
#endif // CALC_ENERGIES
            fcoul *= qq*rinv;
#endif // CALC_COUL_TAB
            ////DEVICE_CODE_FENCE();
#ifdef CALC_ENERGIES
#ifdef DEBUG_FPEX
            TLOG("Vc_ci =%f, vcoul =%f\n", Vc_ci, vcoul);
#endif // DEBUG_FPEX
            Vc_ci += vcoul;
            /* 1 flop for Coulomb energy addition */
#endif // CALC_ENERGIES
#endif // CALC_COULOMB

#ifdef DEBUG_SDLB
            TLOG("kaCHI 7.2.2.\n");
            //wait_host(device_core_id);
#endif

#ifdef CALC_COULOMB
#ifdef HALF_LJ
            if (i < UNROLLI/2)
#endif
            {
                fscal = frLJ*rinvsq + fcoul;
                /* 2 flops for scalar LJ+Coulomb force */
            }
#ifdef HALF_LJ
            else
            {
                fscal = fcoul;
            }
#endif
#else
            fscal = frLJ*rinvsq;
#endif
            fx = fscal*dx;
            fy = fscal*dy;
            fz = fscal*dz;

            /* Increment i-atom force */
            ////DEVICE_CODE_FENCE();
#ifdef SW_NEW_ALG
            // if(write_ci) {
#endif
                fi[i*FI_STRIDE+XX] += fx;
                fi[i*FI_STRIDE+YY] += fy;
                fi[i*FI_STRIDE+ZZ] += fz;
#ifdef SW_NEW_ALG
            // }
#endif
            /* Decrement j-atom force */
            ////DEVICE_CODE_FENCE();
#ifdef DEBUG_SDLB
            TLOG("kaCHI 7.3.\n");
            //wait_host(device_core_id);
#endif
#ifdef SW_NEW_ALG
            if(write_cj) {
#endif
                //TODO: REDUCE SUM
                ldm_f[aj*F_STRIDE+XX-start_f]  -= fx;
                ldm_f[aj*F_STRIDE+YY-start_f]  -= fy;
                ldm_f[aj*F_STRIDE+ZZ-start_f]  -= fz;
                /* 9 flops for force addition */
#ifdef SW_NEW_ALG
            }
#endif
#ifdef DEBUG_SDLB
            TLOG("kaCHI 7.4.\n");
            if(j == 3)
            {
                //wait_host(device_core_id);
            }
#endif
        }
#undef interact
#endif // SIMD_INNER
    }
    
    }
}

#undef EXCL_FORCES
