#ifndef __SW_CONFIG_H__
#define __SW_CONFIG_H__

#define PARAM_SIZE                          32
    #define PARAM_NOTICE                    0
    // 通知flag
    #define PARAM_MPI_RANK                  1
    #define PARAM_DEVICE_ACTION             2
        #define DEVICE_ACTION_EXIT                  999
        #define DEVICE_ACTION_RUN                   666
    // ---- CUSTOM PARAM
    #define WORKLOADPARA                    3
    #define OTHER_PTR_HEAD                  4
    #define SMALL_BLOCK_N                   5
    #define SMALL_BLOCK_X                   6
    #define SMALL_BLOCK_Y                   7
    #define LARGE_BLOCK_N                   8
    #define LARGE_BLOCK_X                   9
    #define LARGE_BLOCK_Y                   10
    // ----

#define CORE_SYNC_64 0x0000FFFF
#define CORE_SYNC_8  0x0000000F

#define DEVICE_SAFE_PAD 48

struct InitParam
{
    int host_rank;
    long *host_to_device, *device_to_host;

    // ---- CUSTOM VAR
    int x_sec, y_sec, Z;
    int *walls, *flags;
    float nu, omega, CSmago;
    // ----
};

struct WorkLoadPara {
    int                         macro_para;
    nbnxn_pairlist_t            *nbl;
    nbnxn_atomdata_t            *nbat;
    interaction_const_t         *ic;
    rvec                        *shift_vec;
    real                        *f;
    real                        *fshift;
    real                        *Vvdw;
    real                        *Vc;
};

// ---- CUSTOM CONFIG

#define BLOCK_SIZE(ID, ID_SZ, N) (((N) / (ID_SZ)) + (((N) % (ID_SZ) > (ID)) ? 1 : 0))
// 将N个任务分给ID_SZ个从核，根据从核ID返回每个从核的任务数量
#define BLOCK_HEAD(ID, ID_SZ, N) (((N) / (ID_SZ)) * (ID) + (((N) % (ID_SZ) > (ID)) ? (ID) : (N) % (ID_SZ)))
// 将N个任务分给ID_SZ个从核，根据从核ID返回第一个任务的编号

// ----

#endif
