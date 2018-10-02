#include "SwDevice.h"

#include <memory.h>
#include <math.h>

__thread_local struct InitParam device_param;

__thread_local int device_core_id, device_core_x, device_core_y;
__thread_local long device_notice_counter;
__thread_local volatile long device_in_param[PARAM_SIZE]; 
__thread_local volatile long device_out_param[PARAM_SIZE]; 

__thread_local volatile unsigned long sync_get_reply, sync_put_reply; 
__thread_local volatile unsigned long async_get_reply, async_put_reply; 
__thread_local volatile unsigned long async_get_reply_counter, async_put_reply_counter; 

void* device_malloc(int sz)
// 申请内存
{
    long ret = ldm_malloc(sz);
    return (void*)ret;
}

void* device_align_malloc(int sz, int alignment, int right_shft)
// 从核内存申请+对齐
{
    void* ptr = ldm_malloc(sz);
    int rem = ((long)ptr) % alignment;
    if(rem > 0)
    {
        return ptr + (alignment - rem + right_shft);
    }
    else
    {
        return ptr + right_shft;
    }
}

__thread_local volatile intv8 bcast_buff; 
void device_bcast_32Byte(int device_core_id, volatile intv8* bcast_data)
// 全核广播
{
    bcast_buff = *bcast_data;

    if (device_core_id / 8 == 0)
    {
        __builtin_sw64_putC(bcast_buff, CORE_SYNC_8);
    } else {
        bcast_buff = __builtin_sw64_getxC(bcast_buff);
    }

    if (device_core_id % 8 == 0)
    {
        asm volatile ("nop":::"memory");
        __builtin_sw64_putR(bcast_buff, CORE_SYNC_8);
        *bcast_data = bcast_buff;
    } else {
        asm volatile ("nop":::"memory");
        *bcast_data = __builtin_sw64_getxR(bcast_buff);
    }
}

void wait_host(int device_core_id)
// 全核等待主核信号
{
    if(device_core_id == 0)
    {
        while(1)
        {
            sync_get_reply = 0;
            athread_get(
                PE_MODE, device_param.host_to_device, (void*)&device_in_param[0],
                sizeof(long) * PARAM_SIZE, (void*)(&sync_get_reply), 0, 0, 0
            );
            while(sync_get_reply != 1);
            asm volatile ("#nop":::"memory");
            if(device_in_param[PARAM_NOTICE] >= device_notice_counter)
                break;
        }

        device_notice_counter++;
    }

    int i;
    for(i = 0; i < PARAM_SIZE / 4; i++)
    {
        device_bcast_32Byte(device_core_id, (intv8*)&device_in_param[i * 4]);
    }
}

void notice_host(int device_core_id)
// 通知主核
{
    device_out_param[PARAM_NOTICE] = device_out_param[PARAM_NOTICE] + 1;

    athread_syn(ARRAY_SCOPE, CORE_SYNC_64);
    if(device_core_id == 0)
    {
        sync_put_reply = 0;
        athread_put(PE_MODE, (void*)&device_out_param[0], device_param.device_to_host,sizeof(long) * PARAM_SIZE, (void*)(&sync_put_reply), 0, 0);
        while(sync_put_reply != 1);

    }
    athread_syn(ARRAY_SCOPE, CORE_SYNC_64);
}

void sync_put(void* h_ptr, void* d_ptr, int size)
// 同步写
{
    sync_put_reply = 0;
    athread_put(PE_MODE, d_ptr, h_ptr, size, (void*)(&sync_put_reply), 0, 0);
    while (sync_put_reply != 1);
}

void sync_get(void* d_ptr, void* h_ptr, int size)
// 同步读
{
    sync_get_reply = 0;
    athread_get(PE_MODE, h_ptr, d_ptr, size, (void*)(&sync_get_reply), 0, 0, 0);
    while (sync_get_reply != 1);
}

void async_put(void* h_ptr, void* d_ptr, int size)
// 异步写
{
    athread_put(PE_MODE, d_ptr, h_ptr, size, (void*)(&async_put_reply), 0, 0);
    async_put_reply_counter++;
}

void wait_all_async_put()
// 等待异步写结束
{
    while(async_put_reply != async_put_reply_counter);
}

void async_get(void* d_ptr, void* h_ptr, int size)
// 异步读
{
    athread_get(PE_MODE, h_ptr, d_ptr, size, (void*)(&async_get_reply), 0, 0, 0);
    async_get_reply_counter++;
}

void wait_all_async_get()
// 等待异步读结束
{
    while(async_get_reply != async_get_reply_counter);
}

// ----------  CUSTOM FUNCTION ----------
extern void device_run();
// --------------------------------------

void device_main(void *_param)
{
    int i;

    device_core_id = athread_get_id(-1);
    device_core_x = device_core_id % 8;
    device_core_y = device_core_id / 8;

    device_param = *((struct InitParam*) _param);

    for (i = 0; i < PARAM_SIZE; i++)
    {
        device_in_param[i] = 0; device_out_param[i] = 0;
    }

    device_notice_counter = 1;
    async_get_reply = 0, async_put_reply = 0;
    async_get_reply_counter = 0, async_put_reply_counter = 0;


    int myrank = device_param.host_rank;
    TLOG("Device Ready\n");
    while(1)
    {
        wait_host(device_core_id);

        if (device_in_param[PARAM_DEVICE_ACTION] == DEVICE_ACTION_EXIT)
        {
            TLOG("Device Release\n");
            notice_host(device_core_id);
            break;
        }
        if (device_in_param[PARAM_DEVICE_ACTION] == DEVICE_ACTION_RUN)
        {
            device_run();
            notice_host(device_core_id);
        }
    }
}
