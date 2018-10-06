#include "SwDevice.h"

// extern subcore_func();

void device_run() {
	// if (device_core_id == device_param.host_rank)
	// 	printf("hello from %d\n", device_core_id);
	// workLoadPara = *((struct WorkLoadPara*)device_in_param[WORKLOADPARA]);
	if (device_core_id == 0) {
		// printf("device_run=%ld\n", device_in_param[WORKLOADPARA]);
		// printf("%d %d\n",sizeof(workLoadPara),workLoadPara.macro_para);
		
	printf("----------2---------------\n");
	// printf("workLoadPara: %ld %d\n",(long)&workLoadPara,workLoadPara.macro_para);
	// printf("size %d\n", sizeof(struct WorkLoadPara));
	// printf("address work %ld  mod32 %ld\n",(long)workLoadPara,((long)workLoadPara)%32);
	// printf("in_param     %ld  mod16 %ld\n",(long)&device_in_param,((long)&device_in_param)%16);
	}
	        asm volatile ("nop":::"memory");
	sync_get(&workLoadPara, (struct WorkLoadPara*)device_in_param[WORKLOADPARA], sizeof(struct WorkLoadPara));
	// 暴力广播workloadpara的数据
	// int i;
    // for(i = 0; i < sizeof(struct WorkLoadPara)/8 / 4; i++)
    // {
    //     device_bcast_32Byte(device_core_id, (intv8*)(((long*)(p_workLoadPara))[i*4]));
    // }
	// device_bcast_32Byte(device_core_id, workLoadPara);
	        asm volatile ("nop":::"memory");
	if (device_core_id == 0) {
		printf("----------0-----------\n");
		printf("address work %ld  mod32 %ld\n",(long)&workLoadPara,((long)&workLoadPara)%32);
	}
	        asm volatile ("nop":::"memory");
	// asm volatile("halt");
	// printf("workLoadPara=%ld\n", workLoadPara.macro_para);
	subcore_func();
}
