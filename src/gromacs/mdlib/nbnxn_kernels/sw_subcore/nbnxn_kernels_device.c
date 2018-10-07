#include "sw/SwDevice.h"
#include "nbnxn_kernels_device.h"

#define M2_DIV(X,M2_LOG2) ((X)>>(M2_LOG2))
#define M2_MOD(X,MSK) ((X)&(MSK))

#define XC_GP   (12)
#ifndef GMX_DOUBLE
#define XC_SZ   (384) //(32*12)
#define XC_MSK  (31)
#define XC_LOG2 (5)
#else
#define XC_SZ   (192) //(16*12)
#define XC_MSK  (15)
#define XC_LOG2 (4)
#endif

// ===== CACHE =====
// -----   x   -----

typedef struct {
    real C[XC_SZ];
    int hd;
} _x_cache_t;

__thread_local real      *Hxi;      
__thread_local _x_cache_t Cxi; // 1536 B
__thread_local int        Sxi; // sizeof(x)/12
__thread_local real      *Hxj; 
__thread_local _x_cache_t Cxj; // 1536 B
__thread_local int        Sxj; // sizeof(x)/12

static inline real *xi_C(unsigned int gp_i)
{
    int Coffset = M2_MOD(gp_i,XC_MSK);
    int Chead   = M2_DIV(gp_i,XC_LOG2);
    int fetch_sz;
    if(Chead == Cxi.hd)
    {
        return &Cxi.C[Coffset*XC_GP];
    }
    DEVICE_CODE_FENCE();
    Cxi.hd = Chead;
#ifndef DEEP_DARK_FANTASY
    if((gp_i+XC_MSK+1) <= Sxi)
    {
        fetch_sz = XC_SZ;
    }
    else
    {
        fetch_sz = (Sxi-gp_i)*XC_GP;
    }
    sync_get(&Cxi.C[0], Hxi+XC_SZ*Chead, fetch_sz*sizeof(real));
#else
    sync_get(&Cxi.C[0], Hxi+XC_SZ*Chead, XC_SZ*sizeof(real));
#endif
    DEVICE_CODE_FENCE();
    return &Cxi.C[Coffset*XC_GP];
}

static inline real *xj_C(unsigned int gp_i)
{
    int Coffset = M2_MOD(gp_i,XC_MSK);
    int Chead   = M2_DIV(gp_i,XC_LOG2);
    int fetch_sz;
    if(Chead == Cxj.hd)
    {
        return &Cxj.C[Coffset*XC_GP];
    }
    DEVICE_CODE_FENCE();
    Cxj.hd = Chead;
#ifndef DEEP_DARK_FANTASY
    else if((gp_i+XC_MSK+1) <= Sxj)
    {
        fetch_sz = XC_SZ;
    }
    else
    {
        fetch_sz = (Sxj-gp_i)*XC_GP;   
    }
    sync_get(&Cxj.C[0], Hxj+XC_SZ*Chead, fetch_sz*sizeof(real));
#else
    sync_get(&Cxj.C[0], Hxj+XC_SZ*Chead, XC_SZ*sizeof(real));
#endif
    DEVICE_CODE_FENCE();
    return &Cxj.C[Coffset*XC_GP];
}

#define QC_GP   (4)
#ifndef GMX_DOUBLE
#define QC_SZ   (256) //(64*4)
#define QC_MSK  (63)
#define QC_LOG2 (6)
#else
#define QC_SZ   (128) //(32*4)
#define QC_MSK  (31)
#define QC_LOG2 (5)
#endif

// ===== CACHE =====
// -----   q   -----

typedef struct {
    real C[QC_SZ];
    int hd;
} _q_cache_t;

__thread_local real       *Hqi; 
__thread_local _q_cache_t  Cqi; // 1024 B
__thread_local int         Sqi;
__thread_local real       *Hqj; 
__thread_local _q_cache_t  Cqj; // 1024 B
__thread_local int         Sqj;

static inline real *qi_C(unsigned int gp_i)
{
    int Coffset = M2_MOD(gp_i,QC_MSK);
    int Chead   = M2_DIV(gp_i,QC_LOG2);
    int fetch_sz;
    if(Chead == Cqi.hd)
    {
        return &Cqi.C[Coffset*QC_GP];
    }
    DEVICE_CODE_FENCE();
    Cqi.hd = Chead;
#ifndef DEEP_DARK_FANTASY
    if((gp_i+QC_MSK+1) <= Sqi)
    {
        fetch_sz = QC_SZ;
    }
    else
    {
        fetch_sz = (Sqi-gp_i)*QC_GP;
    }
    sync_get(&Cqi.C[0], Hqi+QC_SZ*Chead, fetch_sz*sizeof(real));
#else
    sync_get(&Cqi.C[0], Hqi+QC_SZ*Chead, QC_SZ*sizeof(real));
#endif
    DEVICE_CODE_FENCE();
    return &Cqi.C[Coffset*QC_GP];
}

static inline real *qj_C(unsigned int gp_i)
{
    int Coffset = M2_MOD(gp_i,QC_MSK);
    int Chead   = M2_DIV(gp_i,QC_LOG2);
    int fetch_sz;
    if(Chead == Cqj.hd)
    {
        return &Cqj.C[Coffset*QC_GP];
    }
    DEVICE_CODE_FENCE();
    Cqj.hd = Chead;
#ifndef DEEP_DARK_FANTASY
    if((gp_i+QC_MSK+1) <= Sqj)
    {
        fetch_sz = QC_SZ;
    }
    else
    {
        fetch_sz = (Sqj-gp_i)*QC_GP;
    }
    sync_get(&Cqj.C[0], Hqj+QC_SZ*Chead, fetch_sz*sizeof(real));
#else
    sync_get(&Cqj.C[0], Hqj+QC_SZ*Chead, QC_SZ*sizeof(real));
#endif
    DEVICE_CODE_FENCE();
    return &Cqj.C[Coffset*QC_GP];
}

#define TC_GP   (4)
#define TC_SZ   (256) //(64*4)
#define TC_MSK  (63)
#define TC_LOG2 (6)

// ===== CACHE =====
// -----   t   -----

typedef struct {
    int C[TC_SZ];
    int hd;
} _t_cache_t;

__thread_local _t_cache_t  Cti; // 1024 B
__thread_local int         Sti; 
__thread_local _t_cache_t  Ctj; // 1024 B
__thread_local int         Stj;

static inline real *ti_C(int idx)
{

}

static inline real *tj_C(int idx)
{

}

#define CI_GP   (1)
#define CI_SZ   (64)
#define CI_MSK  (63)
#define CI_LOG2 (6)
typedef struct {
    nbnxn_ci_t C[CI_SZ];
    int hd;
} _ci_cache_t;

// ===== CACHE =====
// -----   ci  -----

__thread_local _ci_cache_t Cci; // 1024 B
__thread_local int         Sci;
static inline _ci_cache_t *ci_C(int idx)
{

}

#define CJ_GP   (1)
#define CJ_SZ   (128)
#define CJ_MSK  (127)
#define CJ_LOG2 (7)
typedef struct {
    nbnxn_cj_t C[CJ_SZ];
    int hd;
} _cj_cache_t;

// ===== CACHE =====
// -----   cj  -----

__thread_local int           Scj; // 1024 B
__thread_local _cj_cache_t   Ccj;
static inline _cj_cache_t *cj_C(int idx)
{

}

#define NF_GP   (2)
#ifndef GMX_DOUBLE
#define NF_SZ   (256) //(128*2)
#define NF_MSK  (127)
#define NF_LOG2 (7)
#else
#define NF_C_SZ (128) //(64*2)
#define NF_MSK  (63)
#define NF_LOG2 (6)
#endif

// ===== CACHE =====
// -----  nbfp -----

typedef struct {
    real C[NF_SZ];
    int hd;
} _nf_cache_t;

__thread_local _nf_cache_t Cnf; // 1024 B
__thread_local int         Snf;
static inline real *nf_C(int idx)
{

}

static inline void clear_C()
{
   Cxi.hd = -1;
   Cxj.hd = -1;

   Cqi.hd = -1;
   Cqj.hd = -1;
   Cti.hd = -1;
   Ctj.hd = -1;

   Cci.hd = -1;
   Ccj.hd = -1;

   Cnf.hd = -1;
}

typedef void (*p_nbk_func_noener)();

typedef void (*p_nbk_func_ener)();

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
} func_para_t;

__thread_local func_para_t device_func_para;

/* Analytical reaction-field kernels */
//++++++++++++++++++++++++++++++++++++
#define CALC_COUL_RF
    //===================================
    #define LJ_CUT
    #include "nbnxn_kernel_device_includes.h"
    #undef LJ_CUT

    #ifndef SW_TEST_FUNC /* in SwConfig.h */
        #define LJ_FORCE_SWITCH
        #include "nbnxn_kernel_device_includes.h"
        #undef LJ_FORCE_SWITCH

        #define LJ_POT_SWITCH
        #include "nbnxn_kernel_device_includes.h"
        #undef LJ_POT_SWITCH

        #define LJ_EWALD
            #define LJ_CUT
                #define LJ_EWALD_COMB_GEOM
                #include "nbnxn_kernel_device_includes.h"
                #undef LJ_EWALD_COMB_GEOM

                #define LJ_EWALD_COMB_LB
                #include "nbnxn_kernel_device_includes.h"
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
    #include "nbnxn_kernel_device_includes.h"
    #undef LJ_CUT

    #ifndef SW_TEST_FUNC /* in SwConfig.h */
        #define LJ_FORCE_SWITCH
        #include "nbnxn_kernel_device_includes.h"
        #undef LJ_FORCE_SWITCH
    
        #define LJ_POT_SWITCH
        #include "nbnxn_kernel_device_includes.h"
        #undef LJ_POT_SWITCH

        #define LJ_EWALD
            #define LJ_CUT
                #define LJ_EWALD_COMB_GEOM
                #include "nbnxn_kernel_device_includes.h"
                #undef LJ_EWALD_COMB_GEOM

                #define LJ_EWALD_COMB_LB
                #include "nbnxn_kernel_device_includes.h"
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
        #include "nbnxn_kernel_device_includes.h"
        #undef LJ_CUT

        #define LJ_FORCE_SWITCH
        #include "nbnxn_kernel_device_includes.h"
        #undef LJ_FORCE_SWITCH

        #define LJ_POT_SWITCH
        #include "nbnxn_kernel_device_includes.h"
        #undef LJ_POT_SWITCH

        #define LJ_EWALD
            #define LJ_CUT
                #define LJ_EWALD_COMB_GEOM
                #include "nbnxn_kernel_device_includes.h"
                #undef LJ_EWALD_COMB_GEOM

                #define LJ_EWALD_COMB_LB
                #include "nbnxn_kernel_device_includes.h"
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

__thread_local p_nbk_func_noener p_nbk_c_noener_device[coultNR][vdwtNR] =
{
#ifdef SW_TEST_FUNC /* in SwConfig.h */
    { nbnxn_kernel_ElecRF_VdwLJ_F_device,           NULL, NULL, NULL, NULL           },
    { nbnxn_kernel_ElecQSTab_VdwLJ_F_device,        NULL, NULL, NULL, NULL           },
    { NULL, NULL, NULL, NULL, NULL  }
#else
    { nbnxn_kernel_ElecRF_VdwLJ_F_device,           nbnxn_kernel_ElecRF_VdwLJFsw_F_device,           nbnxn_kernel_ElecRF_VdwLJPsw_F_device,           nbnxn_kernel_ElecRF_VdwLJEwCombGeom_F_device,           nbnxn_kernel_ElecRF_VdwLJEwCombLB_F_device           },
    { nbnxn_kernel_ElecQSTab_VdwLJ_F_device,        nbnxn_kernel_ElecQSTab_VdwLJFsw_F_device,        nbnxn_kernel_ElecQSTab_VdwLJPsw_F_device,        nbnxn_kernel_ElecQSTab_VdwLJEwCombGeom_F_device,        nbnxn_kernel_ElecQSTab_VdwLJEwCombLB_F_device        },
    { nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_F_device, nbnxn_kernel_ElecQSTabTwinCut_VdwLJFsw_F_device, nbnxn_kernel_ElecQSTabTwinCut_VdwLJPsw_F_device, nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombGeom_F_device, nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombLB_F_device }
#endif
};

__thread_local p_nbk_func_ener p_nbk_c_ener_device[coultNR][vdwtNR] =
{
#ifdef SW_TEST_FUNC /* in SwConfig.h */
    { nbnxn_kernel_ElecRF_VdwLJ_VF_device,           NULL, NULL, NULL, NULL           },
    { nbnxn_kernel_ElecQSTab_VdwLJ_VF_device,        NULL, NULL, NULL, NULL           },
    { NULL, NULL, NULL, NULL, NULL  }
#else
    { nbnxn_kernel_ElecRF_VdwLJ_VF_device,           nbnxn_kernel_ElecRF_VdwLJFsw_VF_device,           nbnxn_kernel_ElecRF_VdwLJPsw_VF_device,           nbnxn_kernel_ElecRF_VdwLJEwCombGeom_VF_device,           nbnxn_kernel_ElecRF_VdwLJEwCombLB_VF_device            },
    { nbnxn_kernel_ElecQSTab_VdwLJ_VF_device,        nbnxn_kernel_ElecQSTab_VdwLJFsw_VF_device,        nbnxn_kernel_ElecQSTab_VdwLJPsw_VF_device,        nbnxn_kernel_ElecQSTab_VdwLJEwCombGeom_VF_device,        nbnxn_kernel_ElecQSTab_VdwLJEwCombLB_VF_device         },
    { nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VF_device, nbnxn_kernel_ElecQSTabTwinCut_VdwLJFsw_VF_device, nbnxn_kernel_ElecQSTabTwinCut_VdwLJPsw_VF_device, nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombGeom_VF_device, nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombLB_VF_device  }
#endif
};
#ifdef SW_ENERGRP /* in SwConfig */

__thread_local p_nbk_func_ener p_nbk_c_energrp_device[coultNR][vdwtNR] =
{
#ifdef SW_TEST_FUNC /* in SwConfig.h */
    { nbnxn_kernel_ElecRF_VdwLJ_VgrpF_device,           NULL, NULL, NULL, NULL           },
    { nbnxn_kernel_ElecQSTab_VdwLJ_VgrpF_device,        NULL, NULL, NULL, NULL           },
    { NULL, NULL, NULL, NULL, NULL  }
#else
    { nbnxn_kernel_ElecRF_VdwLJ_VgrpF_device,           nbnxn_kernel_ElecRF_VdwLJFsw_VgrpF_device,           nbnxn_kernel_ElecRF_VdwLJPsw_VgrpF_device,           nbnxn_kernel_ElecRF_VdwLJEwCombGeom_VgrpF_device,           nbnxn_kernel_ElecRF_VdwLJEwCombLB_VgrpF_device           },
    { nbnxn_kernel_ElecQSTab_VdwLJ_VgrpF_device,        nbnxn_kernel_ElecQSTab_VdwLJFsw_VgrpF_device,        nbnxn_kernel_ElecQSTab_VdwLJPsw_VgrpF_device,        nbnxn_kernel_ElecQSTab_VdwLJEwCombGeom_VgrpF_device,        nbnxn_kernel_ElecQSTab_VdwLJEwCombLB_VgrpF_device        },
    { nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VgrpF_device, nbnxn_kernel_ElecQSTabTwinCut_VdwLJFsw_VgrpF_device, nbnxn_kernel_ElecQSTabTwinCut_VdwLJPsw_VgrpF_device, nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombGeom_VgrpF_device, nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombLB_VgrpF_device }
#endif
};

#endif


void kaCHI_func(int n)
{
    TLOG("kaCHI Func called %d.\n", n);
}


void device_run()
{
#ifdef DEBUG_SDLB
    TLOG("kaCHI 0.\n");
#endif
    int func_type = device_in_param[FUNC_TYPE];
    int func_i    = device_in_param[FUNC_I];
    int func_j    = device_in_param[FUNC_J];
    device_func_para     = *((func_para_t*)device_in_param[FUNC_PARAM_PTR]);
#ifdef SW_DEVICE_LOG /* in SwConfig */
    if(((device_notice_counter + device_param.host_rank - 1) % 64) == 0)
        OLOG("FuncType =%d, I =%d, J =%d\n", func_type, func_i, func_j);
#endif
    if(func_type > 1 || func_i > 1 || func_j > 0)
    {
        OLOG("UNKNOWN FUNC: FuncType =%d, I =%d, J =%d\n", func_type, func_i, func_j);
        return;
    }
#define SW_PRINT_PARASIZE
#undef SW_PRINT_PARASIZE
#ifdef SW_PRINT_PARASIZE
    {
        int natoms = device_func_para.nbat->natoms;
        int fstride = device_func_para.nbat->fstride;
        int sizeof_f = natoms*fstride;
        int sizeof_fshift = SHIFTS*DIM;
        if(((device_notice_counter + device_param.host_rank - 1) % 64) == 0)
            OLOG("natoms =%d, fstride =%d, sizeof_f =%d, sizeof_fshift =%d\n", natoms, fstride, sizeof_f, sizeof_fshift);
    }
#endif
#ifdef DEBUG_SDLB
    TLOG("kaCHI 0.1.\n");
    kaCHI_func(0);
    TLOG("kaCHI sizeof(nbnxn_pairlist_t) =%d\n", sizeof(nbnxn_pairlist_t));
    //wait_host(device_core_id);
#endif
    if(func_type == FUNC_NO_ENER)
    {
        p_nbk_c_noener_device[func_i][func_j]();
    }
    else if(func_type == FUNC_ENER)
    {
        p_nbk_c_ener_device[func_i][func_j]();
    }
#ifdef DEBUG_SDLB
    kaCHI_func(1);
    TLOG("kaCHI Final.\n");
#endif
}
