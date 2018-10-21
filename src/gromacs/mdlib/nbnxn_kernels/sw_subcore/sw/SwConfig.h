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

#define DEVICE_SAFE_PAD 64

#define DEVICE_CODE_FENCE() {asm volatile ("nop":::"memory");}

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

#define SW_NEW_ALG

#define DENUG_F
#undef DENUG_F

#define DEBUG_SDLB
#undef DEBUG_SDLB

#define DEBUG_FPEX
#undef DEBUG_FPEX

#define DEBUG_CACHE
#undef DEBUG_CACHE

//---------------------DEEP_DARK_FANTASY----------------------------
//open this to avoid memory boundary check, maybe cause MEMORY ACCESS ERROR
#define DEEP_DARK_FANTASY
//---------------------DEEP_DARK_FANTASY----------------------------

#define COUNT_CACHE_HITS
#undef COUNT_CACHE_HITS

#define SW_NOCACLU
#undef SW_NOCACLU

#define HAHAHAHAHAHAHAHA
#undef HAHAHAHAHAHAHAHA

#define DEBUG_LB
#undef DEBUG_LB

#define SIMD_INNER

#define DEBUG_SIMD
#undef DEBUG_SIMD

#define FORCE_OVERLAP

#define MAX_F_LDM_SIZE 360 // div 12
// ----

#endif
