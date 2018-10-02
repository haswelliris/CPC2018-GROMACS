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

void init_device();

void notice_device();

void wait_device();

void release_device();

#endif
