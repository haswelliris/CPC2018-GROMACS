#include "SwDevice.h"

extern void subcore_func(struct WorkLoadPara *);

void device_run() {
	// if (device_core_id == device_param.host_rank)
	// 	printf("hello from %d\n", device_core_id);
	// printf("%ld\n", device_in_param[WORKLOADPARA]);
	#ifndef HOST_RUN
		subcore_func((struct WorkLoadPara*)device_in_param[WORKLOADPARA]);
	#endif
}