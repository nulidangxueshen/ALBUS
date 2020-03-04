all:spmv_albus

.PHONY : spmv_albus
#compilers(g++ or icc)

CC=icpc
#CC = g++

#g++
#spmv_albus:
#	$(CC) -march=native -fopenmp -O3 -o spmv_albus main.cpp

#icc
spmv_albus:
	$(CC) -xCORE-AVX2 -qopt-prefetch=3 -fopenmp -O3 -o spmv_albus main.cpp

.PHONY:clean
	rm spmv_albus
