# ALBUS
A Method for efficiently processing SpMV using SIMD and load balancing  

## Manual:

(1) Enter the `make` command to compile the Makefile.    
`$ ls`  
`ALBUS(OpenMP)  ALBUS(OpenMP+AVX2)  LICENSE  README.md`  
`$ cd ALBUS'('OpenMP')'`  
`$ ls`  
`Albus_spmv.h  main.cpp  Makefile Storage_format.h`  
`$ make`  
`icpc -xCORE-AVX2 -qopt-prefetch=3 -fopenmp -O3 -o spmv_albus main.cpp`  
`$ ls`   
`Albus_spmv.h  main.cpp  Makefile  spmv_albus  Storage_format.h`  
  
(2) Run "spmv_albus" , Input:  
**`./spmv_albus <matrix name> <iterations times>`**  
`$ ls`  
`Albus_spmv.h  main.cpp  Makefile  spmv_albus  Storage_format.h dense.mtx`  
`$ ./spmv_albus dense.mtx 4096`  
**`------------------------------------------------------------------`  
`FileName   : dense.mtx`  
`Iterations : 4096`  
`--------------Matrix Information and Performance Data-------------`  
`(row_num : 2048) , (col_num : 2048) , (nzz_num : 4194304)`  
`( Matrix Kind : real ) , ( Data type : general )`  
`NZZ_NUM : 4194304      (Note: If the matrix is a symmetric, NZZ_NUM <= 2 * nzz_num)`  
`ALBUS parallel time        >>>>>> 0.110031 ms`  
`ALBUS parallel SPMV Gflops >>>>>> 76.238752 GFlops`  
`-------------------------ALBUS TEST ANSWER------------------------`  
`i = 0 , mtx_ans = 125136.773347`  
`i = 1 , mtx_ans = 157020.769696`  
`i = 2 , mtx_ans = 190806.910217`  
`i = 3 , mtx_ans = 190497.024711`  
`i = 4 , mtx_ans = 124703.888516`  
`i = 5 , mtx_ans = 659168.485721`  
`i = 6 , mtx_ans = 163069.766193`  
`i = 7 , mtx_ans = 332258.152177`  
`i = 8 , mtx_ans = 161110.817741`  
`i = 9 , mtx_ans = 187378.056863`  
` - - -`  
` - - -`  
`i = 2047 , mtx_ans = 125099.232999`  
`-----------------------Serial CSR TEST ANSWER---------------------`  
`i = 0 , mtx_ans = 125136.773347`  
`i = 1 , mtx_ans = 157020.769696`  
`i = 2 , mtx_ans = 190806.910217`  
`i = 3 , mtx_ans = 190497.024711`  
`i = 4 , mtx_ans = 124703.888516`  
`i = 5 , mtx_ans = 659168.485721`  
`i = 6 , mtx_ans = 163069.766193`  
`i = 7 , mtx_ans = 332258.152177`  
`i = 8 , mtx_ans = 161110.817741`  
`i = 9 , mtx_ans = 187378.056863`  
` - - -`  
` - - -`  
`i = 2047 , mtx_ans = 125099.232999`  
`--------------------------Check Answer----------------------------`  
`right.....PASS!`  
`------------------------------------------------------------------`**  
## Note:  
(1) The experimental data comes from formerly the University of Florida Sparse Matrix Collection.  

(2) This program supports `g++` and `icc` compilation. The Makefile file in this program is initially set to be compiled by the `icc` compiler. You can also comment the `icc` compiled code and use the `g++` compiler in the file to compile. For comparison experiments, please choose The same compilation environment.  

(3) This experiment was tested on two CPU platforms: `Dual Socket Intel(R) Xeon(R) CPU E5-2670 v3` and `Dual Socket Intel(R) Xeon(R) Silver 4110 CPU `.  
