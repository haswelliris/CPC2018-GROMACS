TARGET = libsubcore.a
CC = sw5cc

CFLAGS =  -O3 -host -I/usr/sw-mpp/mpi2/include/ -lm -msimd -OPT:ieee_arith=3
SCFLAGS = -O3 -slave -msimd -lmldm_1KB_slave -OPT:ieee_arith=3

COBJ = SwHost.o
SOBJ = SwDevice.o nbnxn_kernels_device.o

VPATH = ./sw/

$(TARGET): $(COBJ) $(SOBJ)
	ar -crv $@ $^
	
$(COBJ): %.o: %.c
	$(CC) $(CFLAGS) -c $<

$(SOBJ): %.o: %.c
	$(CC) $(SCFLAGS) -c $<

#-------------------------------------*
.PHONY : clean clear
clean:
	-rm -rf $(COBJ) $(SOBJ)
	
