#include "SwDevice.h"

extern __thread_local struct InitParam device_param;
extern __thread_local int device_core_id;

void device_run() {
	if (device_core_id == device_param.host_rank)
		printf("hello from %d\n", device_core_id);
}