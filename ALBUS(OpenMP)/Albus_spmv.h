#include<iostream>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<omp.h>
#include<cstring>
#include<sys/time.h>
#include<stdlib.h>
using namespace std;

#define INT int
#define DOU double
#define AVX_DOU __m256d

void thread_block(INT thread_id,INT start,INT end,INT num2,INT num1,INT *&row_ptr,INT *&col_idx,DOU *&mtx_val,DOU *&mtx_ans,DOU *&mid_ans,DOU *&vec_val)
{
        INT start1,end1,num,Thread,i,j;
        DOU sum=0.0;
        Thread = thread_id<<1;
        end1   = row_ptr[start];
        sum    = 0.0;
        for(j=num2;j<end1;j++)
        {
                sum+=mtx_val[j]*vec_val[col_idx[j]];
        }
        mid_ans[Thread] = sum;
        #pragma prefetch
        for(i=start;i<end;++i)
        {
                j    = row_ptr[i];
                end1 = row_ptr[i+1];
                sum  = 0.0;
                for(;j<end1;j++)
                {
                        sum+=mtx_val[j]*vec_val[col_idx[j]];
                }
                mtx_ans[i] = sum;
        }
        start1 = row_ptr[end];
        sum    = 0.0;
        for(j=start1;j<num1;j++)
        {
                sum+=mtx_val[j]*vec_val[col_idx[j]];
        }
        mid_ans[Thread+1] = sum;
}

INT binary_search(INT *&row_ptr,INT num,INT end)
{
        INT l,r,h,t=0;
        l=0,r=end;
        while(l<=r)
        {
                h = (l+r)>>1;
                if(row_ptr[h]>=num)
                {
                        r=h-1;
                }
                else
                {
                        l=h+1;
                        t=h;
                }
        }
        return t;
}

void albus_balance(INT *&row_ptr,INT *&par_set,INT *&start,INT *&end,DOU *&mid_ans,INT *&block_size,INT thread_nums)
{
        block_size[thread_nums] = par_set[2];
        start[0] = 0;
        end[thread_nums-1] = par_set[0];
        INT tt=par_set[2]/thread_nums;
        for(INT i=0;i<thread_nums;i++)
        {
                block_size[i]=tt*i;
                start[i] = binary_search(row_ptr,block_size[i],par_set[0]);
                end[i-1] = start[i];
        }
}

void SPMV_DOU(INT *&row_ptr,INT *&col_idx,DOU *&mtx_val,INT *&par_set,DOU *&mtx_ans,DOU *& vec_val,INT * start,INT * end,DOU * mid_ans,INT * block_size,INT thread_nums)
{
        INT i;
        #pragma omp parallel private(i)
        {
                #pragma omp for schedule(static) nowait
                for(i=0;i<thread_nums;++i)
                {
                        thread_block(i,start[i]+1,end[i],block_size[i],block_size[i+1],row_ptr,col_idx,mtx_val,mtx_ans,mid_ans,vec_val);
                }
        }
        mtx_ans[0] = mid_ans[0];
        INT sub;
        for(i=1;i<thread_nums;++i)
        {
                sub = i<<1;
                if(end[i-1] == start[i])
                {
                        mtx_ans[start[i]] = mid_ans[sub-1] + mid_ans[sub];
                }
                else
                {
                        mtx_ans[start[i]] = mid_ans[sub];
                        mtx_ans[end[i-1]] = mid_ans[sub-1];
                }
        }
}
