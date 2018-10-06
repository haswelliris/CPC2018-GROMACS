CC=sw5cc

host :
	make clean
	make host_libsubcore.a
	make host_libfake.a

device :
	make clean
	make device_libsubcore.a
	make device_libfake.a


host_libsubcore.a : SwDevice.o nbnxn_kernel_ref_outer_subcore.o SwHost.o
	ar -crv libsubcore.a $^

device_libsubcore.a : SwDevice.o nbnxn_kernel_ref_outer_subcore.o SwHost.o nbnxn_kernel_ref_outer_fun.o
	ar -crv libsubcore.a $^

host_libfake.a : nbnxn_kernel_ref_outer_fun.o
	ar -crv libfake.a $^

device_libfake.a : nbnxn_kernel_ref_outer_fun.o
	touch libfake.a

SwDevice.o : SwDevice.c
	${CC} -slave -O3 -msimd -c $^

nbnxn_kernel_ref_outer_subcore.o : nbnxn_kernel_ref_outer_subcore.c
	${CC} -slave -O3 -msimd -c $^

SwHost.o : SwHost.c
	${CC} -host -O3 -I/usr/sw-mpp/mpi2/include/ -lm -msimd -c $^
	
nbnxn_kernel_ref_outer_fun.o : nbnxn_kernel_ref_outer_fun.c nbnxn_kernel_ref_outer_fun.h nbnxn_kernel_ref_inner.h
	${CC} -slave -O3 -I/usr/sw-mpp/mpi2/include/ -lm -msimd -c $<

clean :
	rm -rf *.o *.a