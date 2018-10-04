#ifndef __SW_DEVICE_H__
#define __SW_DEVICE_H__

#include "SwConfig.h"

#include "slave.h"
#include "dma.h"
#include "simd.h"
#include <stdio.h>
#include <stdlib.h>


#define TLOG(M, ...)  {if(device_rank == 0 && device_core_id == 0) { \
	                  fprintf(stderr, "[%s  %s] [DEV TLOG ] (%d): " M "", __DATE__, __TIME__, __LINE__, ##__VA_ARGS__); }}

#define OLOG(M, ...)  {if(device_core_id == 0) { \
	                  fprintf(stderr, "[%s  %s] [DEV OLOG ] (RANK%d  %d): " M "", __DATE__, __TIME__, device_rank, __LINE__, ##__VA_ARGS__); }}

#ifndef __thread_local
#define __thread_local const
#endif

extern __thread_local struct InitParam device_param;
extern __thread_local int device_core_id, device_core_x, device_core_y, device_rank;
extern __thread_local long device_notice_counter;     
extern __thread_local volatile long device_in_param[PARAM_SIZE]; 
extern __thread_local volatile long device_out_param[PARAM_SIZE];   

extern __thread_local volatile unsigned long sync_get_reply, sync_put_reply; 
extern __thread_local volatile unsigned long async_get_reply, async_put_reply; 
extern __thread_local volatile unsigned long async_get_reply_counter, async_put_reply_counter; 

void* device_malloc(int sz);
void* device_align_malloc(int sz, int alignment, int right_shft);

void device_main(void *_param);
void device_bcast_32Byte(int device_core_id, volatile intv8* bcast_data);
void wait_host(int device_core_id);
void notice_host(int device_core_id);

void sync_put(void* h_ptr, void* d_ptr, int size);
void sync_get(void* d_ptr, void* h_ptr, int size);

void async_put(void* h_ptr, void* d_ptr, int size);
void wait_all_async_put();

void async_get(void* d_ptr, void* h_ptr, int size);
void wait_all_async_get();

#endif
