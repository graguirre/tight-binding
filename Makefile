.PHONY=clean

CFLAGS=-g -Wall -lgsl -lm -lgslcblas -fopenmp -std=c99 -Iinclude
CC=gcc

LIBDIR=lib
INCLUDEDIR=include

# enlaza proyecto
tight-binding: ${LIBDIR}/slater-koster.o 
	${CC} ${CFLAGS} src/tight-binding.c ${LIBDIR}/slater-koster.o -o tight-binding

# compila archivos lib
${LIBDIR}/%.o: ${LIBDIR}/%.c ${INCLUDEDIR}/%.h

clean:
	rm -f tight-binding ${LIBDIR}/*.o
