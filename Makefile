CFLAGS=-g -Wall -lgsl -lm -lgslcblas -fopenmp -std=c99 -DDEBUG
CC=gcc

tight-binding: tight-binding.c


clean:
	rm -f tight-binding
