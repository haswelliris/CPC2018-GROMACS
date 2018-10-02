CC=sw5cc

all :
	make libsubcore.a
	make libswhost.a

libsubcore.a : SwDevice.o nbnxn_kernel_ref_outer_subcore.o SwHost.o
	ar -crv $@ $^

SwDevice.o : SwDevice.c
	${CC} -slave -O3 -msimd -c $^

nbnxn_kernel_ref_outer_subcore.o : nbnxn_kernel_ref_outer_subcore.c
	${CC} -slave -O3 -msimd -c $^

#libswhost.a : SwHost.o
#	ar -crv $@ $^

SwHost.o : SwHost.c
	${CC} -host -O3 -I/usr/sw-mpp/mpi2/include/ -lm -msimd -c $^
