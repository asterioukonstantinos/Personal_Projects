#EXECS=test
MPICC?=/cm/shared/apps/openmpi/gcc/64/1.10.3/bin/mpicc

#all: ${EXECS}

par:nbody-par.c
	${MPICC} -Wall -O3 -o nbody-par nbody-par.c -lm -maxv2 -march=native

seq: nbody-seq.c
	${MPICC} -Wall -O3 -o nbody-seq nbody-seq.c -lm -mavx2 -march=native

#clean:
#	rm -f *.o ${EXECS} *~ *core
