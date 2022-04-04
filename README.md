# ALBUS
[A Method for efficiently processing SpMV using SIMD and load balancing](https://doi.org/10.1016/j.future.2020.10.036)  

our paper:  

H. Bian, J. Huang, L. Liu, D. Huang, X. Wang, ALBUS: A method for efficiently processing SpMV using SIMD and Load balancing, Future Generation Computer Systems, Volume 116, 2021, pp. 371-392.


## Performance
**Unit**：GFlops

**Equipment**：Intel(R) Xeon(R) Gold 6242 CPU (double socket)(16cores * 2)

**NOTE:**

CSR5(ICS'15)

Merge(SC'16)

CVR(CGO'18)

MKL(Intel oneAPI 2021)

|  Sparse Matrix  | ALBUS  | CSR5   | CVR  |  Merge   |   MKL  |
|      :----:     | :----: | :----: | :----: | :----: | :----: |
|       rma10     | **94.52** | 42.05 | 65.84 |  | |
| cant            | 65.38 | 52.66 | **68.66** |   |  |
| mac_econ_fwd500 | **68.56** | 25.33 | 40.79 |   |  |
| mc2depi         | **75.26** | 45.10 | 31.78 |   |  |
| shipsec1        | **33.03** | 31.14 | 31.85 |   |  |
| scircuit        | **43.32** | 20.25 | 34.35 |   |  |
| dc2             | **65.32** | 24.59 | 39.28 |   |  |
| amazon-2008     | **17.82** | 15.55 | 13.17 |   |  |
| consph          | **53.14** | 43.74 | 51.07 |   |  |
| cop20k_A        | **41.26** | 36.06 | 32.16 |   |  |
| crankseg_2      | **24.93** | 23.22 | 24.07 |   |  |
| ins2            | 23.06 | 18.92 | **34.16** |   |  |
| circuit5M       | **13.02** | 12.61 | 11.10 |   |  |
| FullChip        | **13.32** | 12.63 | 11.18 |   |  |
| mip1            | **26.21** | 25.80 | 26.10 |   |  |
| pdb1HYS         | 65.44 | 57.03 | **68.23** |   |  |
| pwtk            | **25.53** | 24.54 | 24.18 |   |  |
| LP              | **19.89** | 19.23 | 19.54 |   |  |
| in-2004         | **18.36** | 17.55 | 15.50 |   |  |
| eu-2005         | **19.98** | 18.76 | 16.37 |   |  |
| road_usa        |  **6.71** | 6.56 | 6.27 |   |  |
| roadNet-CA      | **14.60** | 12.23 | 10.86 |   |  |
| soc-LiveJournal1| **9.31**  | 7.95 | 8.39 |   |  |
| webbase-1M      | **34.14** | 21.92 | 22.13 |   |  |
| web-Google      | **11.23** | 10.40 | 8.09 |   |  |
| web-Stanford    | **26.62** | 18.28 | 13.42 |   |  |
| Stanford        | **19.65** | 15.93 | 11.24 |   |  |
| wiki-topcats    | **15.24** | 13.62 | 12.56 |   |  |


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

## Cite This:  
```
@article{BIAN2021371,
  title = {ALBUS: A method for efficiently processing SpMV using SIMD and Load balancing},
  journal = {Future Generation Computer Systems},
  volume = {116},
  pages = {371-392},
  year = {2021},
  issn = {0167-739X},
  doi = {https://doi.org/10.1016/j.future.2020.10.036},
  url = {https://www.sciencedirect.com/science/article/pii/S0167739X2033020X},
  author = {Haodong Bian and Jianqiang Huang and Lingbin Liu and Dongqiang Huang and Xiaoying Wang},
  keywords = {SpMV, ALBUS, CSR5, MKL, SIMD, Load balancing},
  abstract = {SpMV (Sparse matrix–vector multiplication) is widely used in many fields. Improving the performance of SpMV has been the pursuit of many researchers. Parallel SpMV using multi-core processors has been a standard parallel method used by researchers. In reality, the number of non-zero elements in many sparse matrices is not evenly distributed, so parallelism without preprocessing will cause a large amount of performance loss due to uneven load. In this paper, we propose ALBUS (Absolute Load Balancing Using SIMD (Single Instruction Multiple Data)), a method for efficiently processing SpMV using load balancing and SIMD vectorization. On the one hand, ALBUS can achieve multi-core balanced load processing; on the other hand, it gives full play to the ability of SIMD vectorization parallelism under the CPU. We selected 20 sets of regular matrices and 20 sets of irregular matrices to form the Benchmark suite. We performed SpMV performance comparison tests on ALBUS, CSR5 (Compressed Sparse Row 5), Merge(Merge-based SpMV), and MKL (Math Kernel Library) under the same conditions. On the E5-2670 v3 CPU platform, For 20 sets of regular matrices, ALBUS can achieve an average speedup of 1.59x, 1.32x, 1.48x (up to 2.53x, 2.22x, 2.31x) compared to CSR5, Merge, MKL, respectively. For 20 sets of irregular matrices, ALBUS can achieve an average speedup of 1.38x, 1.42x, 2.44x (up to 2.33x, 2.24x, 5.37x) compared to CSR5, Merge, MKL, respectively.}
}
```
