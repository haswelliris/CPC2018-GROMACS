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
            int  aj[4];
            realv4 dx, dy, dz;
            realv4 rsq, rinv;
            realv4 rinvsq, rinvsix;
            realv4 c6, c12;
            realv4 FrLJ6, FrLJ12, frLJ, VLJ;
            FrLJ6.v = 0.0, FrLJ12.v = 0.0, frLJ.v = 0.0, VLJ.v = 0.0;

#ifdef CALC_COULOMB
            realv4 qq;
            realv4 fcoul;
#ifdef CALC_COUL_TAB
            realv4 rs, frac;
            int  ri[4];
            realv4 fexcl;
#endif
#ifdef CALC_ENERGIES
            realv4 vcoul;
#endif
#endif
            realv4 fscal;
            realv4 fx, fy, fz;

            /* A multiply mask used to zero an interaction
             * when either the distance cutoff is exceeded, or
             * (if appropriate) the i and j indices are
             * unsuitable for this kind of inner loop. */
            realv4 skipmask;

#ifdef CHECK_EXCLS
            /* A multiply mask used to zero an interaction
             * when that interaction should be excluded
             * (e.g. because of bonding). */
            realv4 interact;

            //interact = ((l_cj[cjind].excl>>(i*UNROLLI + j)) & 1);
            //interact = ((Ccj_p->excl>>(i*UNROLLI + j)) & 1);
            interact.p[0] = (real)((Ccj_p->excl>>(i*UNROLLI + 0)) & 1);
            interact.p[1] = (real)((Ccj_p->excl>>(i*UNROLLI + 1)) & 1);
            interact.p[2] = (real)((Ccj_p->excl>>(i*UNROLLI + 2)) & 1);
            interact.p[3] = (real)((Ccj_p->excl>>(i*UNROLLI + 3)) & 1);
            // WHAT!!????
            // real it0 = (real)((Ccj_p->excl>>(i*UNROLLI + 0)) & 1);
            // real it1 = (real)((Ccj_p->excl>>(i*UNROLLI + 1)) & 1);
            // real it2 = (real)((Ccj_p->excl>>(i*UNROLLI + 2)) & 1);
            // real it3 = (real)((Ccj_p->excl>>(i*UNROLLI + 3)) & 1);
            // interact.v = simd_set_realv4(it0, it1, it2, it3);
#ifndef EXCL_FORCES
            skipmask.v = interact.v;
#else   
            //skipmask = !(cj == ci_sh && j <= i);
            skipmask.p[0] = !(cj == ci_sh && 0 <= i);
            skipmask.p[1] = !(cj == ci_sh && 1 <= i);
            skipmask.p[2] = !(cj == ci_sh && 2 <= i);
            skipmask.p[3] = !(cj == ci_sh && 3 <= i);
            // FaQ！
            // skipmask.v = simd_set_realv4(
            //     ((real)!(cj == ci_sh && 0 <= i)),
            //     ((real)!(cj == ci_sh && 1 <= i)),
            //     ((real)!(cj == ci_sh && 2 <= i)),
            //     ((real)!(cj == ci_sh && 3 <= i))
            // );
            
#endif  // EXCL_FORCES
#else
#define interact one
            skipmask.v = 1.0;
#endif  // CHECK_EXCLS
            //DEVICE_CODE_FENCE();
#ifdef DEBUG_SDLB
            TLOG("kaCHI 7.1.\n");
#endif
            //DEVICE_CODE_FENCE();
            //aj = cj*UNROLLJ + j;
            aj[0] = cj*UNROLLJ + 0;
            aj[1] = cj*UNROLLJ + 1;
            aj[2] = cj*UNROLLJ + 2;
            aj[3] = cj*UNROLLJ + 3;

            realv4 xiX, xiY, xiZ;
            realv4 xjX, xjY, xjZ;
            xiX.v = xi[i*XI_STRIDE+XX];
            xiY.v = xi[i*XI_STRIDE+YY];
            xiZ.v = xi[i*XI_STRIDE+ZZ];
            xjX.p[0] = Cxj_p[0*X_STRIDE+XX];
            xjY.p[0] = Cxj_p[0*X_STRIDE+YY];
            xjZ.p[0] = Cxj_p[0*X_STRIDE+ZZ];
            xjX.p[1] = Cxj_p[1*X_STRIDE+XX];
            xjY.p[1] = Cxj_p[1*X_STRIDE+YY];
            xjZ.p[1] = Cxj_p[1*X_STRIDE+ZZ];
            xjX.p[2] = Cxj_p[2*X_STRIDE+XX];
            xjY.p[2] = Cxj_p[2*X_STRIDE+YY];
            xjZ.p[2] = Cxj_p[2*X_STRIDE+ZZ];
            xjX.p[3] = Cxj_p[3*X_STRIDE+XX];
            xjY.p[3] = Cxj_p[3*X_STRIDE+YY];
            xjZ.p[3] = Cxj_p[3*X_STRIDE+ZZ];
            // dx  = xi[i*XI_STRIDE+XX] - x[aj*X_STRIDE+XX];
            // dy  = xi[i*XI_STRIDE+YY] - x[aj*X_STRIDE+YY];
            // dz  = xi[i*XI_STRIDE+ZZ] - x[aj*X_STRIDE+ZZ];
            // dx  = xi[i*XI_STRIDE+XX] - Cxj_p[j*X_STRIDE+XX];
            // dy  = xi[i*XI_STRIDE+YY] - Cxj_p[j*X_STRIDE+YY];
            // dz  = xi[i*XI_STRIDE+ZZ] - Cxj_p[j*X_STRIDE+ZZ];

            // dx.p[0]  = xiX.p[0] - xjX.p[0];
            // dy.p[0]  = xiY.p[0] - xjY.p[0];
            // dz.p[0]  = xiZ.p[0] - xjZ.p[0];
            // dx.p[1]  = xiX.p[1] - xjX.p[1];
            // dy.p[1]  = xiY.p[1] - xjY.p[1];
            // dz.p[1]  = xiZ.p[1] - xjZ.p[1];
            // dx.p[2]  = xiX.p[2] - xjX.p[2];
            // dy.p[2]  = xiY.p[2] - xjY.p[2];
            // dz.p[2]  = xiZ.p[2] - xjZ.p[2];
            // dx.p[3]  = xiX.p[3] - xjX.p[3];
            // dy.p[3]  = xiY.p[3] - xjY.p[3];
            // dz.p[3]  = xiZ.p[3] - xjZ.p[3];
            dx.v  = xiX.v - xjX.v;
            dy.v  = xiY.v - xjY.v;
            dz.v  = xiZ.v - xjZ.v;
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

            //rsq = dx*dx + dy*dy + dz*dz;
            // rsq.p[0] = dx.p[0]*dx.p[0] + dy.p[0]*dy.p[0] + dz.p[0]*dz.p[0];
            // rsq.p[1] = dx.p[1]*dx.p[1] + dy.p[1]*dy.p[1] + dz.p[1]*dz.p[1];
            // rsq.p[2] = dx.p[2]*dx.p[2] + dy.p[2]*dy.p[2] + dz.p[2]*dz.p[2];
            // rsq.p[3] = dx.p[3]*dx.p[3] + dy.p[3]*dy.p[3] + dz.p[3]*dz.p[3];
            rsq.v = dx.v*dx.v + dy.v*dy.v + dz.v*dz.v;

#ifdef DEBUG_FPEX
            TLOG("xiX =%f, xiY =%f, xiZ =%f\n", xi[i*XI_STRIDE+XX], xi[i*XI_STRIDE+YY], xi[i*XI_STRIDE+ZZ]);
            TLOG("x_X =%f, x_Y =%f, x_Z =%f\n", x[aj*X_STRIDE+XX], x[aj*X_STRIDE+YY], x[aj*X_STRIDE+ZZ]);
            TLOG("dx =%f, dy =%f, dz =%f\n", dx, dy, dz);
            TLOG("rsq =%f, rcut2 =%f\n", rsq, rcut2);
#endif

            /* Prepare to enforce the cut-off. */
            //skipmask = (rsq >= rcut2) ? 0 : skipmask;
            /* 9 flops for r^2 + cut-off check */
            // skipmask.p[0] = (rsq.p[0] >= rcut2) ? 0.0 : skipmask.p[0];
            // skipmask.p[1] = (rsq.p[1] >= rcut2) ? 0.0 : skipmask.p[1];
            // skipmask.p[2] = (rsq.p[2] >= rcut2) ? 0.0 : skipmask.p[2];
            // skipmask.p[3] = (rsq.p[3] >= rcut2) ? 0.0 : skipmask.p[3];
            skipmask.p[0] = (rsq.p[0] - rcut2 >= 0) ? 0.0 : skipmask.p[0];
            skipmask.p[1] = (rsq.p[1] - rcut2 >= 0) ? 0.0 : skipmask.p[1];
            skipmask.p[2] = (rsq.p[2] - rcut2 >= 0) ? 0.0 : skipmask.p[2];
            skipmask.p[3] = (rsq.p[3] - rcut2 >= 0) ? 0.0 : skipmask.p[3];
            // skipmask.p[0] = (rsq.p[0] - rcut2 < 0) ? skipmask.p[0] :0.0;
            // skipmask.p[1] = (rsq.p[1] - rcut2 < 0) ? skipmask.p[1] :0.0;
            // skipmask.p[2] = (rsq.p[2] - rcut2 < 0) ? skipmask.p[2] :0.0;
            // skipmask.p[3] = (rsq.p[3] - rcut2 < 0) ? skipmask.p[3] :0.0;
            // SIMD_SELLGE NOT SUPPORTED?! KIDDING ME!
            // DEVICE_CODE_FENCE();
            // realv4 rdiff;
            // rdiff.v = rsq.v - rcut2;
            // skipmask.v = __builtin_sw64_sllt(rdiff.v, skipmask.v, zero.v);
            //DEVICE_CODE_FENCE();
#ifdef DEBUG_SDLB
            TLOG("kaCHI 7.2.\n");
            //wait_host(device_core_id);
#endif
#ifdef CHECK_EXCLS
            /* Excluded atoms are allowed to be on top of each other.
             * To avoid overflow of rinv, rinvsq and rinvsix
             * we add a small number to rsq for excluded pairs only.
             */
            //rsq += (1 - interact)*NBNXN_AVOID_SING_R2_INC;
            // rsq.p[0] += (1.0 - interact.p[0])*NBNXN_AVOID_SING_R2_INC;
            // rsq.p[1] += (1.0 - interact.p[1])*NBNXN_AVOID_SING_R2_INC;
            // rsq.p[2] += (1.0 - interact.p[2])*NBNXN_AVOID_SING_R2_INC;
            // rsq.p[3] += (1.0 - interact.p[3])*NBNXN_AVOID_SING_R2_INC;
            rsq.v += (1.0 - interact.v)*NBNXN_AVOID_SING_R2_INC;
#endif

            //rinv = gmx_invsqrt(rsq);
            // rinv.p[0] = gmx_invsqrt(rsq.p[0]);
            // rinv.p[1] = gmx_invsqrt(rsq.p[1]);
            // rinv.p[2] = gmx_invsqrt(rsq.p[2]);
            // rinv.p[3] = gmx_invsqrt(rsq.p[3]);
            rinv.v = 1.0/simd_vsqrt(rsq.v);
            /* 5 flops for invsqrt */

            /* Partially enforce the cut-off (and perhaps
             * exclusions) to avoid possible overflow of
             * rinvsix when computing LJ, and/or overflowing
             * the Coulomb table during lookup. */
            //rinv = rinv * skipmask;
            // rinv.p[0] = rinv.p[0] * skipmask.p[0];
            // rinv.p[1] = rinv.p[1] * skipmask.p[1];
            // rinv.p[2] = rinv.p[2] * skipmask.p[2];
            // rinv.p[3] = rinv.p[3] * skipmask.p[3];
            // FIXME !!!!!!!!
            rinv.v  = rinv.v*skipmask.v;

            //rinvsq  = rinv*rinv;
            // rinvsq.p[0]  = rinv.p[0]*rinv.p[0];
            // rinvsq.p[1]  = rinv.p[1]*rinv.p[1];
            // rinvsq.p[2]  = rinv.p[2]*rinv.p[2];
            // rinvsq.p[3]  = rinv.p[3]*rinv.p[3];
            // FIXME !!!!!!!!
            rinvsq.v  = rinv.v*rinv.v;

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
                //c6      = nbfp[type_i_off+Ctj_p[j]*2  ];
                //c12     = nbfp[type_i_off+Ctj_p[j]*2+1];
                c6.p [0]      = nbfp[type_i_off+Ctj_p[0]*2  ];
                c12.p[0]      = nbfp[type_i_off+Ctj_p[0]*2+1];
                c6.p [1]      = nbfp[type_i_off+Ctj_p[1]*2  ];
                c12.p[1]      = nbfp[type_i_off+Ctj_p[1]*2+1];
                c6.p [2]      = nbfp[type_i_off+Ctj_p[2]*2  ];
                c12.p[2]      = nbfp[type_i_off+Ctj_p[2]*2+1];
                c6.p [3]      = nbfp[type_i_off+Ctj_p[3]*2  ];
                c12.p[3]      = nbfp[type_i_off+Ctj_p[3]*2+1];
                DEVICE_CODE_FENCE();
#if defined LJ_CUT
                //rinvsix = interact*rinvsq*rinvsq*rinvsq;
                // rinvsix.p[0] = interact.p[0]*rinvsq.p[0]*rinvsq.p[0]*rinvsq.p[0];
                // rinvsix.p[1] = interact.p[1]*rinvsq.p[1]*rinvsq.p[1]*rinvsq.p[1];
                // rinvsix.p[2] = interact.p[2]*rinvsq.p[2]*rinvsq.p[2]*rinvsq.p[2];
                // rinvsix.p[3] = interact.p[3]*rinvsq.p[3]*rinvsq.p[3]*rinvsq.p[3];
                // FIXME !!!!!!!!
                rinvsix.v = interact.v*rinvsq.v*rinvsq.v*rinvsq.v;
                //FrLJ6   = c6*rinvsix;
                // FrLJ6.p[0]   = c6.p[0]*rinvsix.p[0];
                // FrLJ6.p[1]   = c6.p[1]*rinvsix.p[1];
                // FrLJ6.p[2]   = c6.p[2]*rinvsix.p[2];
                // FrLJ6.p[3]   = c6.p[3]*rinvsix.p[3];
                // FIXME !!!!!!!!
                FrLJ6.v   = c6.v*rinvsix.v;
                // FrLJ12  = c12*rinvsix*rinvsix;
                // FrLJ12.p[0]  = c12.p[0]*rinvsix.p[0]*rinvsix.p[0];
                // FrLJ12.p[1]  = c12.p[1]*rinvsix.p[1]*rinvsix.p[1];
                // FrLJ12.p[2]  = c12.p[2]*rinvsix.p[2]*rinvsix.p[2];
                // FrLJ12.p[3]  = c12.p[3]*rinvsix.p[3]*rinvsix.p[3];
                FrLJ12.v  = c12.v*rinvsix.v*rinvsix.v;
                //frLJ    = FrLJ12 - FrLJ6;
                // frLJ.p[0]    = FrLJ12.p[0] - FrLJ6.p[0];
                // frLJ.p[1]    = FrLJ12.p[1] - FrLJ6.p[1];
                // frLJ.p[2]    = FrLJ12.p[2] - FrLJ6.p[2];
                // frLJ.p[3]    = FrLJ12.p[3] - FrLJ6.p[3];
                // FIXME !!!!!!!!
                frLJ.v    = FrLJ12.v - FrLJ6.v;
                /* 7 flops for r^-2 + LJ force */
#if defined CALC_ENERGIES
                //VLJ     = (FrLJ12 + c12*cpot)/12 - (FrLJ6 + c6*cpot)/6;
                // VLJ.p[0]     = (FrLJ12.p[0] + c12.p[0]*cpot)/12 - (FrLJ6.p[0] + c6.p[0]*cpot)/6;
                // VLJ.p[1]     = (FrLJ12.p[1] + c12.p[1]*cpot)/12 - (FrLJ6.p[1] + c6.p[1]*cpot)/6;
                // VLJ.p[2]     = (FrLJ12.p[2] + c12.p[2]*cpot)/12 - (FrLJ6.p[2] + c6.p[2]*cpot)/6;
                // VLJ.p[3]     = (FrLJ12.p[3] + c12.p[3]*cpot)/12 - (FrLJ6.p[3] + c6.p[3]*cpot)/6;
                VLJ.v     = (c12.v*cpot.v+FrLJ12.v)/12.0 - (c6.v*cpot.v+FrLJ6.v)/6.0;
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
                //VLJ     = VLJ * interact;
                // VLJ.p[0]     = VLJ.p[0] * interact.p[0];
                // VLJ.p[1]     = VLJ.p[1] * interact.p[1];
                // VLJ.p[2]     = VLJ.p[2] * interact.p[2];
                // VLJ.p[3]     = VLJ.p[3] * interact.p[3];
                VLJ.v     = VLJ.v * interact.v;
#endif

#if defined CALC_ENERGIES
                /* Need to zero the interaction if r >= rcut */
                //VLJ     = VLJ * skipmask;
                // VLJ.p[0]     = VLJ.p[0] * skipmask.p[0];
                // VLJ.p[1]     = VLJ.p[1] * skipmask.p[1];
                // VLJ.p[2]     = VLJ.p[2] * skipmask.p[2];
                // VLJ.p[3]     = VLJ.p[3] * skipmask.p[3];
                VLJ.v     = VLJ.v * skipmask.v;
                /* 1 more flop for LJ energy */
#endif

                ////DEVICE_CODE_FENCE();
#ifdef CALC_ENERGIES
                //Vvdw_ci += VLJ;
                Vvdw_ci += VLJ.p[0]+VLJ.p[1]+VLJ.p[2]+VLJ.p[3];
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
            //qq = skipmask * qi[i] * Cqj_p[j];
            // qq.p[0] = skipmask.p[0] * qi[i] * Cqj_p[0];
            // qq.p[1] = skipmask.p[1] * qi[i] * Cqj_p[1];
            // qq.p[2] = skipmask.p[2] * qi[i] * Cqj_p[2];
            // qq.p[3] = skipmask.p[3] * qi[i] * Cqj_p[3];
            realv4 qiV, qjV;
            qiV.v = qi[i];
            // qjV.p[0] = Cqj_p[0];
            // qjV.p[1] = Cqj_p[1];
            // qjV.p[2] = Cqj_p[2];
            // qjV.p[3] = Cqj_p[3];
            simd_load(qjV.v, &Cqj_p[0]);
            qq.v = skipmask.v * qiV.v * qjV.v;

#ifdef CALC_COUL_RF
            //fcoul  = qq*(interact*rinv*rinvsq - k_rf2);
            // fcoul.p[0]  = qq.p[0]*(interact.p[0]*rinv.p[0]*rinvsq.p[0] - k_rf2);
            // fcoul.p[1]  = qq.p[1]*(interact.p[1]*rinv.p[1]*rinvsq.p[1] - k_rf2);
            // fcoul.p[2]  = qq.p[2]*(interact.p[2]*rinv.p[2]*rinvsq.p[2] - k_rf2);
            // fcoul.p[3]  = qq.p[3]*(interact.p[3]*rinv.p[3]*rinvsq.p[3] - k_rf2);
            fcoul.v  = qq.v*(interact.v*rinv.v*rinvsq.v - k_rf2);
            /* 4 flops for RF force */
#ifdef CALC_ENERGIES
            //vcoul  = qq*(interact*rinv + k_rf*rsq - c_rf);
            // vcoul.p[0]  = qq.p[0]*(interact.p[0]*rinv.p[0] + k_rf*rsq.p[0] - c_rf);
            // vcoul.p[1]  = qq.p[1]*(interact.p[1]*rinv.p[1] + k_rf*rsq.p[1] - c_rf);
            // vcoul.p[2]  = qq.p[2]*(interact.p[2]*rinv.p[2] + k_rf*rsq.p[2] - c_rf);
            // vcoul.p[3]  = qq.p[3]*(interact.p[3]*rinv.p[3] + k_rf*rsq.p[3] - c_rf);
            vcoul.v  = qq.v*(interact.v*rinv.v + k_rf*rsq.v - c_rf);
            /* 4 flops for RF energy */
#endif // CALC_ENERGIES
#endif // CALC_COUL_RF

#ifdef CALC_COUL_TAB
            ////DEVICE_CODE_FENCE();
#ifdef DEBUG_FPEX
            TLOG("rsq =%f, tabq_scale =%f\n", rsq, ic.tabq_scale);
#endif //DEBUG_FPEX
            //rs     = rsq*rinv*ic.tabq_scale;
            // rs.p[0]     = rsq.p[0]*rinv.p[0]*ic.tabq_scale;
            // rs.p[1]     = rsq.p[1]*rinv.p[1]*ic.tabq_scale;
            // rs.p[2]     = rsq.p[2]*rinv.p[2]*ic.tabq_scale;
            // rs.p[3]     = rsq.p[3]*rinv.p[3]*ic.tabq_scale;
            rs.v     = rsq.v*rinv.v*ic.tabq_scale;
            //DEVICE_CODE_FENCE();
            //ri     = (int)rs;
            ri[0]     = (int)rs.p[0];
            ri[1]     = (int)rs.p[1];
            ri[2]     = (int)rs.p[2];
            ri[3]     = (int)rs.p[3];
            //frac   = rs - ri;
            frac.p[0]   = rs.p[0] - (real)ri[0];
            frac.p[1]   = rs.p[1] - (real)ri[1];
            frac.p[2]   = rs.p[2] - (real)ri[2];
            frac.p[3]   = rs.p[3] - (real)ri[3];
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
            //fexcl  = tab_coul_FDV0[ri*4] + frac*tab_coul_FDV0[ri*4+1];
            realv4 F, D, V;
            F.p[0] = tab_coul_FDV0[ri[0]*4];
            D.p[0] = tab_coul_FDV0[ri[0]*4+1];
#ifdef CALC_ENERGIES
            V.p[0] = tab_coul_FDV0[ri[0]*4+2];
#endif
            F.p[1] = tab_coul_FDV0[ri[1]*4];
            D.p[1] = tab_coul_FDV0[ri[1]*4+1];
#ifdef CALC_ENERGIES
            V.p[1] = tab_coul_FDV0[ri[1]*4+2];
#endif
            F.p[2] = tab_coul_FDV0[ri[2]*4];
            D.p[2] = tab_coul_FDV0[ri[2]*4+1];
#ifdef CALC_ENERGIES
            V.p[2] = tab_coul_FDV0[ri[2]*4+2];
#endif
            F.p[3] = tab_coul_FDV0[ri[3]*4];
            D.p[3] = tab_coul_FDV0[ri[3]*4+1];
#ifdef CALC_ENERGIES
            V.p[3] = tab_coul_FDV0[ri[3]*4+2];
#endif
            //DEVICE_CODE_FENCE();
            // fexcl.p[0]  = F.p[0] + frac.p[0]*D.p[0];
            // fexcl.p[1]  = F.p[1] + frac.p[1]*D.p[1];
            // fexcl.p[2]  = F.p[2] + frac.p[2]*D.p[2];
            // fexcl.p[3]  = F.p[3] + frac.p[3]*D.p[3];
            fexcl.v  = frac.v*D.v+F.v;
            //fexcl  = (1 - frac)*tab_coul_FDV0[ri*2] + frac*tab_coul_FDV0[(ri+1)*2];
#else
            /* fexcl = (1-frac) * F_i + frac * F_(i+1) */
            //fexcl  = (1 - frac)*tab_coul_F[ri] + frac*tab_coul_F[ri+1];
            fexcl.p[0]  = (1.0 - frac.p[0])*tab_coul_F[ri[0]] + frac.p[0]*tab_coul_F[ri[0]+1];
            fexcl.p[1]  = (1.0 - frac.p[1])*tab_coul_F[ri[1]] + frac.p[1]*tab_coul_F[ri[1]+1];
            fexcl.p[2]  = (1.0 - frac.p[2])*tab_coul_F[ri[2]] + frac.p[2]*tab_coul_F[ri[2]+1];
            fexcl.p[3]  = (1.0 - frac.p[3])*tab_coul_F[ri[3]] + frac.p[3]*tab_coul_F[ri[3]+1];
#endif // GMX_DOUBLE
            //fcoul  = interact*rinvsq - fexcl;
            // fcoul.p[0]  = interact.p[0]*rinvsq.p[0] - fexcl.p[0];
            // fcoul.p[1]  = interact.p[1]*rinvsq.p[1] - fexcl.p[1];
            // fcoul.p[2]  = interact.p[2]*rinvsq.p[2] - fexcl.p[2];
            // fcoul.p[3]  = interact.p[3]*rinvsq.p[3] - fexcl.p[3];
            fcoul.v  = interact.v*rinvsq.v-fexcl.v;
            /* 7 flops for float 1/r-table force */
#ifdef CALC_ENERGIES
            ////DEVICE_CODE_FENCE();
#ifndef GMX_DOUBLE
            //TODO: ldm load: tab_coul_FDV0, tab_coul_V, tab_coul_F
#ifdef DEBUG_FPEX
            TLOG("qq =%f, rinv =%f, interact =%f, sh_ewald =%f, halfsp =%f, frac =%f, fexcl =%f\n", qq, Vc_ci, interact, ic.sh_ewald, halfsp, frac, fexcl);
#endif // DEBUG_FPEX
            // vcoul  = qq*(interact*(rinv - ic.sh_ewald)
            //              -(tab_coul_FDV0[ri*4+2]
            //                -halfsp*frac*(tab_coul_FDV0[ri*4] + fexcl)));


            // vcoul.p[0]  = qq.p[0]*(interact.p[0]*(rinv.p[0] - ic.sh_ewald)
            //              -(V.p[0]
            //                -halfsp*frac.p[0]*(F.p[0] + fexcl.p[0])));
            // vcoul.p[1]  = qq.p[1]*(interact.p[1]*(rinv.p[1] - ic.sh_ewald)
            //              -(V.p[1]
            //                -halfsp*frac.p[1]*(F.p[1] + fexcl.p[1])));
            // vcoul.p[2]  = qq.p[2]*(interact.p[2]*(rinv.p[2] - ic.sh_ewald)
            //              -(V.p[2]
            //                -halfsp*frac.p[2]*(F.p[2] + fexcl.p[2])));
            // vcoul.p[3]  = qq.p[3]*(interact.p[3]*(rinv.p[3] - ic.sh_ewald)
            //              -(V.p[3]
            //                -halfsp*frac.p[3]*(F.p[3] + fexcl.p[3])));
            // Unaligned EX
            vcoul.v  = qq.v*(interact.v*(rinv.v - ic.sh_ewald)
                         -(V.v
                           -halfsp*frac.v*(F.v + fexcl.v)));

                           
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
            //fcoul *= qq*rinv;
            // fcoul.p[0] *= qq.p[0]*rinv.p[0];
            // fcoul.p[1] *= qq.p[1]*rinv.p[1];
            // fcoul.p[2] *= qq.p[2]*rinv.p[2];
            // fcoul.p[3] *= qq.p[3]*rinv.p[3];
            fcoul.v *= qq.v*rinv.v;
#endif // CALC_COUL_TAB
            ////DEVICE_CODE_FENCE();
#ifdef CALC_ENERGIES
#ifdef DEBUG_FPEX
            TLOG("Vc_ci =%f, vcoul =%f\n", Vc_ci, vcoul);
#endif // DEBUG_FPEX
            //Vc_ci += vcoul;
            Vc_ci += vcoul.p[0] + vcoul.p[1] + vcoul.p[2] + vcoul.p[3];
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
                //fscal = frLJ*rinvsq + fcoul;
                // fscal.p[0] = frLJ.p[0]*rinvsq.p[0] + fcoul.p[0];
                // fscal.p[1] = frLJ.p[1]*rinvsq.p[1] + fcoul.p[1];
                // fscal.p[2] = frLJ.p[2]*rinvsq.p[2] + fcoul.p[2];
                // fscal.p[3] = frLJ.p[3]*rinvsq.p[3] + fcoul.p[3];
                fscal.v = frLJ.v*rinvsq.v+fcoul.v;
                /* 2 flops for scalar LJ+Coulomb force */
            }
#ifdef HALF_LJ
            else
            {
                //fscal = fcoul;
                // fscal.p[0] = fcoul.p[0];
                // fscal.p[1] = fcoul.p[1];
                // fscal.p[2] = fcoul.p[2];
                // fscal.p[3] = fcoul.p[3];
                fscal.v = fcoul.v;
            }
#endif
#else
            //fscal = frLJ*rinvsq;
            // fscal.p[0] = frLJ.p[0]*rinvsq.p[0];
            // fscal.p[1] = frLJ.p[1]*rinvsq.p[1];
            // fscal.p[2] = frLJ.p[2]*rinvsq.p[2];
            // fscal.p[3] = frLJ.p[3]*rinvsq.p[3];
            fscal.v = frLJ.v*rinvsq.v;
#endif
            // fx = fscal*dx;
            // fy = fscal*dy;
            // fz = fscal*dz;
            // fx.p[0] = fscal.p[0]*dx.p[0];
            // fy.p[0] = fscal.p[0]*dy.p[0];
            // fz.p[0] = fscal.p[0]*dz.p[0];
            // fx.p[1] = fscal.p[1]*dx.p[1];
            // fy.p[1] = fscal.p[1]*dy.p[1];
            // fz.p[1] = fscal.p[1]*dz.p[1];
            // fx.p[2] = fscal.p[2]*dx.p[2];
            // fy.p[2] = fscal.p[2]*dy.p[2];
            // fz.p[2] = fscal.p[2]*dz.p[2];
            // fx.p[3] = fscal.p[3]*dx.p[3];
            // fy.p[3] = fscal.p[3]*dy.p[3];
            // fz.p[3] = fscal.p[3]*dz.p[3];
            fx.v = fscal.v*dx.v;
            fy.v = fscal.v*dy.v;
            fz.v = fscal.v*dz.v;

            /* Increment i-atom force */
            ////DEVICE_CODE_FENCE();
#ifdef SW_NEW_ALG
            if(write_ci) {
#endif
                // fi[i*FI_STRIDE+XX] += fx;
                // fi[i*FI_STRIDE+YY] += fy;
                // fi[i*FI_STRIDE+ZZ] += fz;
                fi[i*FI_STRIDE+XX] += fx.p[0] + fx.p[1] + fx.p[2] + fx.p[3];
                fi[i*FI_STRIDE+YY] += fy.p[0] + fy.p[1] + fy.p[2] + fy.p[3];
                fi[i*FI_STRIDE+ZZ] += fz.p[0] + fz.p[1] + fz.p[2] + fz.p[3];
#ifdef SW_NEW_ALG
            }
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
                ldm_f[aj[0]*F_STRIDE+XX-start_f]  -= fx.p[0];
                ldm_f[aj[0]*F_STRIDE+YY-start_f]  -= fy.p[0];
                ldm_f[aj[0]*F_STRIDE+ZZ-start_f]  -= fz.p[0];

                ldm_f[aj[1]*F_STRIDE+XX-start_f]  -= fx.p[1];
                ldm_f[aj[1]*F_STRIDE+YY-start_f]  -= fy.p[1];
                ldm_f[aj[1]*F_STRIDE+ZZ-start_f]  -= fz.p[1];

                ldm_f[aj[2]*F_STRIDE+XX-start_f]  -= fx.p[2];
                ldm_f[aj[2]*F_STRIDE+YY-start_f]  -= fy.p[2];
                ldm_f[aj[2]*F_STRIDE+ZZ-start_f]  -= fz.p[2];

                ldm_f[aj[3]*F_STRIDE+XX-start_f]  -= fx.p[3];
                ldm_f[aj[3]*F_STRIDE+YY-start_f]  -= fy.p[3];
                ldm_f[aj[3]*F_STRIDE+ZZ-start_f]  -= fz.p[3];
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
#undef interact
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
