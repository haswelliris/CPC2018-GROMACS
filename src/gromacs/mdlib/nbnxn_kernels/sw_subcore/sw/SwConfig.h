#ifndef __SW_CONFIG_H__
#define __SW_CONFIG_H__

#define PARAM_SIZE                          32
    #define PARAM_NOTICE                    0
    #define PARAM_MPI_RANK                  1
    #define PARAM_DEVICE_ACTION             2
        #define DEVICE_ACTION_EXIT                  999
        #define DEVICE_ACTION_RUN                   666
    // ---- CUSTOM PARAM
    #define CURRENT_PTR_HEAD                3
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

// ---- CUSTOM CONFIG
#define SMALL_BLOCK_Z 50
#define LARGE_BLOCK_Z 20

#define BLOCK_SIZE(ID, ID_SZ, N) (((N) / (ID_SZ)) + (((N) % (ID_SZ) > (ID)) ? 1 : 0))
#define BLOCK_HEAD(ID, ID_SZ, N) (((N) / (ID_SZ)) * (ID) + (((N) % (ID_SZ) > (ID)) ? (ID) : (N) % (ID_SZ)))

#define NODES_IDX(I1, I2, I3, I4, I5) (((((I1) * (x_sec + 2) + (I2)) * (y_sec + 2) + (I3)) * Z + (I4)) * 19 + (I5))
#define FLAGS_IDX(I2, I3, I4) (((I2) * (y_sec + 2) + (I3)) * Z + (I4))
#define WALLS_IDX(I2, I3, I4, I5) ((((I2) * y_sec + (I3)) * Z + (I4)) * 19 + (I5))

#define D_OTHER_IDX(I1, I2) (((I1) + 1) * 19 + (I2))

#define D_NODES_IDX(I2, I3, I4, I5) ((((I2) * (device_param.y_sec + 2) + (I3)) * device_param.Z + (I4)) * 19 + (I5))
#define D_FLAGS_IDX(I2, I3, I4) (((I2) * (device_param.y_sec + 2) + (I3)) * device_param.Z + (I4))
#define D_WALLS_IDX(I2, I3, I4, I5) ((((I2) * device_param.y_sec + (I3)) * device_param.Z + (I4)) * 19 + (I5))
// ----

#endif
