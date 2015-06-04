.PHONY=clean

CFLAGS=-g -Wall -lgsl -lm -lgslcblas -fopenmp -std=c99 -Iinclude
CC=gcc

LIBDIR=lib
INCLUDEDIR=include

# enlaza proyecto
tight-binding: ${LIBDIR}/slater-koster.o ${LIBDIR}/hamiltonian.o 
	${CC} ${CFLAGS} src/tight-binding.c ${LIBDIR}/slater-koster.o ${LIBDIR}/hamiltonian.o -o tight-binding

# compila archivos lib
${LIBDIR}/slater-koster.o: ${LIBDIR}/slater-koster.c ${INCLUDEDIR}/slater-koster.h
	${CC} ${CFLAGS} ${LIBDIR}/slater-koster.c -c -o ${LIBDIR}/slater-koster.o
	

${LIBDIR}/hamiltonian.o: ${LIBDIR}/hamiltonian.c ${INCLUDEDIR}/hamiltonian.h
	${CC} ${CFLAGS} ${LIBDIR}/hamiltonian.c -c -o ${LIBDIR}/hamiltonian.o


clean:
	rm -f tight-binding ${LIBDIR}/*.o
