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

// #if defined CHECK_EXCLS && (defined CALC_COULOMB || defined LJ_EWALD)
// #define EXCL_FORCES
// #endif
if (macro_has(para_LJ_EWALD)) {
    #if defined CHECK_EXCLS
        macro_para |= 1 << para_EXCL_FORCES;
    #endif
}
else {
    #if defined CHECK_EXCLS && defined CALC_COULOMB
        macro_para |= 1 << para_EXCL_FORCES;
    #endif
}

{
    int cj;
// #ifdef ENERGY_GROUPS
    int egp_cj;
// #endif
    int i;

    cj = l_cj[cjind].cj;


    if (macro_has(para_ENERGY_GROUPS))
        egp_cj = nbat->energrp[cj];

    for (i = 0; i < UNROLLI; i++)
    {
        int ai;
        int type_i_off;
        int j;

        ai = ci*UNROLLI + i;

        type_i_off = type[ai]*ntype2;

        for (j = 0; j < UNROLLJ; j++)
        {
            int  aj;
            real dx, dy, dz;
            real rsq, rinv;
            real rinvsq, rinvsix;
            real c6, c12;
            real FrLJ6 = 0, FrLJ12 = 0, frLJ = 0, VLJ = 0;
            // #if defined LJ_FORCE_SWITCH || defined LJ_POT_SWITCH
            real r, rsw;
            // #endif

            #ifdef CALC_COULOMB
            real qq;
            real fcoul;
            // #ifdef CALC_COUL_TAB
            real rs, frac;
            int  ri;
            real fexcl;
            // #endif
            // #ifdef CALC_ENERGIES
            real vcoul;
            // #endif
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

            interact = ((l_cj[cjind].excl>>(i*UNROLLI + j)) & 1);
            // #ifndef EXCL_FORCES
            if (!macro_has(para_EXCL_FORCES))
                skipmask = interact;
            // #else
            else
                skipmask = !(cj == ci_sh && j <= i);
            // #endif
#else
#define interact 1.0
            skipmask = 1.0;
#endif

            aj = cj*UNROLLJ + j;

            dx  = xi[i*XI_STRIDE+XX] - x[aj*X_STRIDE+XX];
            dy  = xi[i*XI_STRIDE+YY] - x[aj*X_STRIDE+YY];
            dz  = xi[i*XI_STRIDE+ZZ] - x[aj*X_STRIDE+ZZ];

            rsq = dx*dx + dy*dy + dz*dz;

            /* Prepare to enforce the cut-off. */
            skipmask = (rsq >= rcut2) ? 0 : skipmask;
            /* 9 flops for r^2 + cut-off check */

#ifdef CHECK_EXCLS
            /* Excluded atoms are allowed to be on top of each other.
             * To avoid overflow of rinv, rinvsq and rinvsix
             * we add a small number to rsq for excluded pairs only.
             */
            rsq += (1 - interact)*NBNXN_AVOID_SING_R2_INC;
#endif

#ifdef COUNT_PAIRS
            npair++;
#endif

            rinv = (1.0/sqrt(rsq));//rsqrt(rsq);
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
                c6      = nbfp[type_i_off+type[aj]*2  ];
                c12     = nbfp[type_i_off+type[aj]*2+1];

                int mask_3 = macro_has(para_LJ_CUT) || macro_has(para_LJ_FORCE_SWITCH) || macro_has(para_LJ_POT_SWITCH);
                rinvsix = interact*rinvsq*rinvsq*rinvsq*mask_3 + rinvsix*!mask_3;
                FrLJ6   = c6*rinvsix*mask_3 + FrLJ6*!mask_3;
                FrLJ12  = c12*rinvsix*rinvsix*mask_3 + FrLJ12*!mask_3;
                frLJ    = (FrLJ12 - FrLJ6)*mask_3 + frLJ*!mask_3;
                /* 7 flops for r^-2 + LJ force */
                int mask_3_1 = mask_3&(macro_has(para_CALC_ENERGIES) || macro_has(para_LJ_POT_SWITCH))
                VLJ     = ((FrLJ12 + c12*ic->repulsion_shift.cpot)/12 -
                    (FrLJ6 + c6*ic->dispersion_shift.cpot)/6)*mask_3_1 + VLJ*!mask_3_1;
                /* 7 flops for LJ energy */

                /* Force or potential switching from ic->rvdw_switch */
                int mask_4 = macro_has(para_LJ_FORCE_SWITCH) || macro_has(para_LJ_POT_SWITCH;
                r       = rsq*rinv*mask_4 + r*!mask_4;
                rsw     = (r - ic->rvdw_switch)*mask_4 + rsw*!mask_4;
                rsw     = (rsw >= 0.0 ? rsw : 0.0)*mask_4 + rsw*!mask_4;

                int mask_5 = macro_has(para_LJ_FORCE_SWITCH);
                int mask_5_1 = macro_has(para_CALC_ENERGIES);
                frLJ   +=
                    (-c6*(ic->dispersion_shift.c2 + ic->dispersion_shift.c3*rsw)*rsw*rsw*r
                    + c12*(ic->repulsion_shift.c2 + ic->repulsion_shift.c3*rsw)*rsw*rsw*r)*mask_5;
                VLJ    +=
                    (-c6*(-ic->dispersion_shift.c2/3 - ic->dispersion_shift.c3/4*rsw)*rsw*rsw*rsw
                    + c12*(-ic->repulsion_shift.c2/3 - ic->repulsion_shift.c3/4*rsw)*rsw*rsw*rsw)*mask_5_1;

                int mask_6 = macro_has(para_CALC_ENERGIES) || macro_has(para_LJ_POT_SWITCH);
                /* Masking should be done after force switching,
                    * but before potential switching.
                    */
                /* Need to zero the interaction if there should be exclusion. */
                VLJ     = VLJ * interact*mask_6 + VLJ*!mask_6;

                int mask_7 = macro_has(para_LJ_POT_SWITCH);
                real sw, dsw;

                sw    = 1.0 + (swV3 + (swV4+ swV5*rsw)*rsw)*rsw*rsw*rsw;
                dsw   = (swF2 + (swF3 + swF4*rsw)*rsw)*rsw*rsw;

                frLJ  = (frLJ*sw - r*VLJ*dsw)*mask_7 + frLJ*!mask_7;
                VLJ  = VLJ*sw*mask_7 + VLJ*!mask_7;

                // #ifdef LJ_EWALD
                if (macro_has(para_LJ_EWALD))
                {
                    real c6grid, rinvsix_nm, cr2, expmcr2, poly, sh_mask;

                    int mask_8 = macro_has(para_LJ_EWALD_COMB_GEOM);
                    c6grid       = ljc[type[ai]*2]*ljc[type[aj]*2]*mask_8;
                    real sigma, sigma2, epsilon;

                    /* These sigma and epsilon are scaled to give 6*C6 */
                    sigma   = ljc[type[ai]*2] + ljc[type[aj]*2];
                    epsilon = ljc[type[ai]*2+1]*ljc[type[aj]*2+1];

                    sigma2  = sigma*sigma;
                    c6grid  = epsilon*sigma2*sigma2*sigma2*!mask_8 + c6grid;
                    //printf("No LJ Ewald combination rule defined");

                    #ifdef CHECK_EXCLS
                        /* Recalculate rinvsix without exclusion mask */
                        rinvsix_nm   = rinvsq*rinvsq*rinvsq;
                    #else
                        rinvsix_nm   = rinvsix;
                    #endif
                    cr2          = lje_coeff2*rsq;
                    #ifdef GMX_DOUBLE
                        expmcr2      = exp(-cr2);
                    #else
                        expmcr2      = expf(-cr2);
                    #endif
                    poly         = 1 + cr2 + 0.5*cr2*cr2;

                    /* Subtract the grid force from the total LJ force */
                    frLJ        += c6grid*(rinvsix_nm - expmcr2*(rinvsix_nm*poly + lje_coeff6_6));

                    int mask_9 = macro_has(para_CALC_ENERGIES);
                    /* Shift should only be applied to real LJ pairs */
                    sh_mask      = lje_vc*interact;

                    VLJ         += c6grid/6*(rinvsix_nm*(1 - expmcr2*poly) + sh_mask)*mask_9;

                }
                /* LJ_EWALD */

                int mask_10 = macro_has(para_VDW_CUTOFF_CHECK);
                int mask_10_1 = macro_has(para_CALC_ENERGIES);
                real skipmask_rvdw;

                skipmask_rvdw = (rsq < rvdw2);
                frLJ         = skipmask_rvdw*frLJ*mask_10 + frLJ*!mask_10;
                VLJ    = VLJ*skipmask_rvdw*mask_10*mask_10_1 + VLJ*!(mask_10*mask_10_1);
                VLJ   = VLJ * skipmask*(!mask_10*mask_10_1) + VLJ *!(!mask_10*mask_10_1)
                
                int mask1 = macro_has(para_CALC_ENERGIES);
                int mask2 = macro_has(para_ENERGY_GROUPS);
                Vvdw[egp_sh_i[i]+((egp_cj>>(nbat->neg_2log*j)) & egp_mask)] += VLJ*mask1*mask2;
                Vvdw_ci += VLJ*mask1*!mask2;
                /* 1 flop for LJ energy addition */
                }
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
            qq = skipmask * qi[i] * q[aj];

            int para_CALC_COUL_RFmask = macro_has(para_CALC_COUL_RF);
            int para_CALC_ENERGIESmask = macro_has(para_CALC_ENERGIES);
            fcoul  = qq*(interact*rinv*rinvsq - k_rf2)*para_CALC_COUL_RFmask + fcoul*para_CALC_COUL_RFmask+ fcoul*!para_CALC_COUL_RFmask;
            /* 4 flops for RF force */
            vcoul  = qq*(interact*rinv + k_rf*rsq - c_rf)*para_CALC_COUL_RFmask*!para_CALC_ENERGIESmask + vcoul*!para_CALC_COUL_RFmask*para_CALC_ENERGIESmask;
            /* 4 flops for RF energy */

            if (macro_has(para_CALC_COUL_TAB)) {
                rs     = rsq*rinv*ic->tabq_scale;
                ri     = (int)rs;
                frac   = rs - ri;
                #ifndef GMX_DOUBLE
                    /* fexcl = F_i + frac * (F_(i+1)-F_i) */
                    fexcl  = tab_coul_FDV0[ri*4] + frac*tab_coul_FDV0[ri*4+1];
                #else
                    /* fexcl = (1-frac) * F_i + frac * F_(i+1) */
                    fexcl  = (1 - frac)*tab_coul_F[ri] + frac*tab_coul_F[ri+1];
                #endif
                fcoul  = interact*rinvsq - fexcl;
                /* 7 flops for float 1/r-table force */

                int maskpara_CALC_ENERGIES = macro_has(para_CALC_ENERGIES);
                #ifndef GMX_DOUBLE
                vcoul  = (qq*(interact*(rinv - ic->sh_ewald)
                                -(tab_coul_FDV0[ri*4+2]
                                -halfsp*frac*(tab_coul_FDV0[ri*4] + fexcl))))*maskpara_CALC_ENERGIES + vcoul*!maskpara_CALC_ENERGIES;
                /* 7 flops for float 1/r-table energy (8 with excls) */
                #else
                vcoul  = (qq*(interact*(rinv - ic->sh_ewald)
                                -(tab_coul_V[ri]
                                -halfsp*frac*(tab_coul_F[ri] + fexcl))))*maskpara_CALC_ENERGIES + vcoul*!maskpara_CALC_ENERGIES;;
                #endif

                fcoul *= qq*rinv;
            }

            int mask1 = macro_has(para_CALC_ENERGIES);
            int mask2 = macro_has(para_ENERGY_GROUPS);
            Vc[egp_sh_i[i]+((egp_cj>>(nbat->neg_2log*j)) & egp_mask)] += vcoul*mask1*mask2;
            Vc_ci += vcoul*mask1*!mask2;
                /* 1 flop for Coulomb energy addition */
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
            fi[i*FI_STRIDE+XX] += fx;
            fi[i*FI_STRIDE+YY] += fy;
            fi[i*FI_STRIDE+ZZ] += fz;
            /* Decrement j-atom force */
            f[aj*F_STRIDE+XX]  -= fx;
            f[aj*F_STRIDE+YY]  -= fy;
            f[aj*F_STRIDE+ZZ]  -= fz;
            /* 9 flops for force addition */
        }
    }
}

#undef interact
// #undef EXCL_FORCES
macro_para &= ~(1 << para_EXCL_FORCES);