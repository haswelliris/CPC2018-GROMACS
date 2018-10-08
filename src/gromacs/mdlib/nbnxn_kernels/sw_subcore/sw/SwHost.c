#include "SwHost.h"

#include <stdlib.h>
#include "mpi.h"
#include "athread.h"

volatile long host_out_param[PARAM_SIZE];
volatile long host_in_param[PARAM_SIZE];
volatile long host_notice_counter;
struct InitParam host_param;
int load_balance_step = 1;


extern void SLAVE_FUN(device_main)();

void notice_device()
{
    asm volatile ("nop":::"memory");
    host_out_param[PARAM_NOTICE] = host_out_param[PARAM_NOTICE] + 1;
}

void wait_device()
{
    while(host_in_param[PARAM_NOTICE] < host_notice_counter);
    host_notice_counter++;
}

void init_device()
{
    int i;
    for(i = 0; i < PARAM_SIZE; i++)
    {
        host_out_param[i] = 0;
        host_in_param[i] = 0;
    }
    host_notice_counter = 1;

    MPI_Comm_rank(MPI_COMM_WORLD, &host_param.host_rank);
    host_param.host_to_device = (long*)&host_out_param[0];
    host_param.device_to_host = (long*)&host_in_param[0];

    host_out_param[PARAM_MPI_RANK] = host_param.host_rank;

    athread_init();
    athread_spawn(device_main, (void*)&host_param);

    char *load_balance_step_str = getenv("LOAD_BALANCE_STEP");
    if (load_balance_step_str != NULL)
        load_balance_step = atoi(load_balance_step_str);
    if (load_balance_step == 0)
        load_balance_step = 1;
    TLOG("LOAD BALANCE STEP =%d\n", load_balance_step);
}

void release_device()
{
    
    host_out_param[PARAM_DEVICE_ACTION] = DEVICE_ACTION_EXIT;
    notice_device();
    wait_device();

    athread_join();
}
