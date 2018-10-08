#ifndef __SW_HOST_H__
#define __SW_HOST_H__
#include "SwConfig.h"

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "athread.h"

extern volatile long host_out_param[PARAM_SIZE];
extern volatile long host_in_param[PARAM_SIZE];
extern volatile long host_notice_counter;
extern struct InitParam host_param;

extern int load_balance_step;

void init_device();

void notice_device();

void wait_device();

void release_device();

#define OLOG(M, ...)  { \
	                  fprintf(stderr, "[HOST OLOG ] (RANK%d  %d): " M "", host_param.host_rank, __LINE__, ##__VA_ARGS__); }

#define TLOG(M, ...)  {if(host_param.host_rank == 0) { \
	                  fprintf(stderr, "[HOST TLOG ] (%d): " M "", __LINE__, ##__VA_ARGS__); }}

#endif
