all:spmv_albus

.PHONY : spmv_albus
#compilers(g++ or icc)

CC=icpc
#CC = g++

#g++
#spmv_albus:
#	$(CC) -march=native -mavx -mavx2 -fopenmp -O3 -o spmv_albus main.cpp

#icc
spmv_albus:
	$(CC) -g -xCORE-AVX2 -fopenmp -O3 -o spmv_albus main.cpp

.PHONY:clean
	rm spmv_albus
