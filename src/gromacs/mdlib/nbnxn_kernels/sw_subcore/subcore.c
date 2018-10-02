#include <stdio.h>
#include "slave.h"
#include "subcore.h"
#include <math.h>

__thread_local int subcore_id;

void subcoreController(void* paras) {
	subcore_id = athread_get_id(-1);
	if (subcore_id == 0)
		printf("hello world for %d!\n", subcore_id);
}
