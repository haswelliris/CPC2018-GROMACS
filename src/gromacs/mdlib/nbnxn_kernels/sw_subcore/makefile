CC=sw5cc

libsubcore.a : subcore.o
	ar -crv $@ $^

subcore.o : subcore.c
	${CC} -slave -c $^
