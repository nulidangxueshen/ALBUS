#include<stdio.h>
#include<math.h>
#include<time.h>
#include"Storage_format.h"
#include"Albus_spmv.h"
#include<immintrin.h>
#include<omp.h>
#include<sys/time.h>
#include<stdlib.h>
#define INT int
#define DOU double

int main(int argc , char ** argv)
{
        //------------------calculation iterations------------------//
        char * filename;
        char * iterations;
        if(argc>1)
        {
                filename   = argv[1];
                iterations = argv[2];
        }
        else
        {
                printf("Error!\n");
        }
        printf("------------------------------------------------------------------\n");
        printf("FileName   : %s\n",filename);
        //----------------------------------------------------------//
        //------------------calculation iterations------------------//
        INT ite=0;
        for(INT i=0;i<strlen(iterations);i++)
        {
                ite = ite * 10 + iterations[i]-'0';
        }
        printf("Iterations : %d\n",ite);
        //----------------------------------------------------------//
        FILE * fp_mtx;
        INT  * row_ptr;
        INT  * col_idx;
        DOU  * mtx_val;
        DOU  * mtx_ans;
        DOU  * mtx_ans_serial;
        DOU  * vec_val;
        INT  * par_set;
        INT  * start;
        INT  * end;
        DOU  * mid_ans;
        INT  * block_size;
        DOU    start_time,end_time;
        INT    thread_nums = omp_get_max_threads();
        fp_mtx = fopen(filename,"rb+");
        printf("--------------Matrix Information and Performance Data-------------\n");
        ReadFile(fp_mtx,row_ptr,col_idx,mtx_val,vec_val,par_set);
        //---------------------------------------------------------------//
        mtx_ans    = (DOU *)aligned_alloc(64,sizeof(DOU)*par_set[0]);
        //---------------------------------------------------------------//
        start      = (INT *)malloc(sizeof(INT)*thread_nums);
        end        = (INT *)malloc(sizeof(INT)*thread_nums);
        mid_ans    = (DOU *)malloc(sizeof(DOU)*thread_nums*2);
        block_size = (INT *)malloc(sizeof(INT)*(thread_nums+1));
	ite = par_set[8];
        //---------------------------------------------------------------//
        albus_balance(row_ptr,par_set,start,end,mid_ans,block_size,thread_nums);
        struct timeval startTime,endTime;
        //-----------------------------Warm Up---------------------------//
        for(INT i=0;i<ite;i++)
        {
                SPMV_DOU(row_ptr,col_idx,mtx_val,par_set,mtx_ans,vec_val,start,end,mid_ans,block_size,thread_nums);
        }
        //---------------------------------------------------------------//
        DOU Timeuse=0;
        gettimeofday(&startTime,NULL);
        for(INT i=0;i<ite;i++)
        {
                SPMV_DOU(row_ptr,col_idx,mtx_val,par_set,mtx_ans,vec_val,start,end,mid_ans,block_size,thread_nums);
        }
        gettimeofday(&endTime,NULL);
        Timeuse = (endTime.tv_sec - startTime.tv_sec) + (endTime.tv_usec - startTime.tv_usec)/1000000.0;
        //-----------------------Performance Data------------------------//
        printf("ALBUS parallel time        >>>>>> %lf ms \n",Timeuse*1000/ite*1.0);
        printf("ALBUS parallel SPMV Gflops >>>>>> %lf GFlops \n",(par_set[2]/((Timeuse)/ite*1.0)/1000000000)*2.0);
        //----------------------------check------------------------------//
        mtx_ans_serial = (DOU *)malloc(sizeof(DOU)*par_set[0]);
        INT error_cnt=0;
        printf("-------------------------ALBUS TEST ANSWER------------------------\n");
        for(INT i=0;i<10;i++)
        {
                printf("i = %d , mtx_ans = %lf\n",i,mtx_ans[i]);
        }
        printf(" - - - \n");
        printf(" - - - \n");
        printf("i = %d , mtx_ans = %lf\n",par_set[0]-1,mtx_ans[par_set[0]-1]);
        printf("-----------------------Serial CSR TEST ANSWER---------------------\n");
        for(INT i=0;i<par_set[0];i++)
        {
                mtx_ans_serial[i]=0.0;
                for(INT j=row_ptr[i];j<row_ptr[i+1];j++)
                {
                        mtx_ans_serial[i]+=mtx_val[j]*vec_val[col_idx[j]];
                }
                DOU res  = mtx_ans_serial[i]-mtx_ans[i]>=0.0?mtx_ans_serial[i]-mtx_ans[i]:-mtx_ans_serial[i]+mtx_ans[i];
                DOU res1 = mtx_ans_serial[i]>=0.0?mtx_ans_serial[i]:-mtx_ans_serial[i];
                if((res-0.01*res1)>0.0000001)
                {
                        printf("id = %d %.16lf %.16lf %.16lf\n",i,mtx_ans[i],mtx_ans_serial[i],res-0.01*res1);
                        error_cnt++;
                }
        }
        for(INT i=0;i<10;i++)
        {
                printf("i = %d , mtx_ans = %lf\n",i,mtx_ans_serial[i]);
        }
        printf(" - - - \n");
        printf(" - - - \n");
        printf("i = %d , mtx_ans = %lf\n",par_set[0]-1,mtx_ans_serial[par_set[0]-1]);
        printf("--------------------------Check Answer----------------------------\n");
        if(error_cnt==0)
        {
                printf("right.....PASS!\n");
        }
        else
        {
                printf(" %d error.......please you check your program\n",error_cnt);
        }
        printf("------------------------------------------------------------------\n");
        //---------------------------------------------------------------------//
        free(row_ptr);
        free(col_idx);
        free(mtx_val);
        free(mtx_ans);
        free(mtx_ans_serial);
        free(vec_val);
        free(par_set);
        free(start);
        free(end);
        free(mid_ans);
        free(block_size);
        return 0;
}
