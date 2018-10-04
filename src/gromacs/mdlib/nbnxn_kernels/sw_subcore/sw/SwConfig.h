#ifndef __SW_CONFIG_H__
#define __SW_CONFIG_H__

#define PARAM_SIZE                          32
    #define PARAM_NOTICE                    0
    #define PARAM_MPI_RANK                  1
    #define PARAM_DEVICE_ACTION             2
        #define DEVICE_ACTION_EXIT                  999
        #define DEVICE_ACTION_RUN                   666
    // ---- CUSTOM PARAM
    #define FUNC_TYPE                       3
        #define FUNC_NO_ENER                        0
        #define FUNC_ENER                           1
        #define FUNC_ENERGRP                        2
    #define FUNC_I                          4
    #define FUNC_J                          5
    #define FUNC_PARAM_PTR                  6
    // ----

#define CORE_SYNC_64 0x0000FFFF
#define CORE_SYNC_8  0x0000000F

#define DEVICE_SAFE_PAD 48

struct InitParam
{
    int host_rank;
    long *host_to_device, *device_to_host;

    // ---- CUSTOM VAR
    // ----
};

// ---- CUSTOM CONFIG
#define BLOCK_SIZE(ID, ID_SZ, N) (((N) / (ID_SZ)) + (((N) % (ID_SZ) > (ID)) ? 1 : 0))
#define BLOCK_HEAD(ID, ID_SZ, N) (((N) / (ID_SZ)) * (ID) + (((N) % (ID_SZ) > (ID)) ? (ID) : (N) % (ID_SZ)))

#define SW_HOST_LOG
#undef SW_HOST_LOG

#define SW_DEVICE_LOG
#undef SW_DEVICE_LOG

#define SW_ENERGRP
#undef SW_ENERGRP

#define SW_TEST_FUNC
// ----

#endif
