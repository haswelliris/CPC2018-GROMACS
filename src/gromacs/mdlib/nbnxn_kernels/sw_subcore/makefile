CC=sw5cc

all :
	make libsubcore.a

libsubcore.a : subcore.o
	ar -crv $@ $^

subcore.o : subcore.c
	${CC} -slave -c $^

