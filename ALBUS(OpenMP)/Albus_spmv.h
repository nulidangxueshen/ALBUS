#include<iostream>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<omp.h>
#include<immintrin.h>
#include<cstring>
#include<sys/time.h>
#include<stdlib.h>
using namespace std;

#define INT int
#define DOU double
#define AVX_DOU __m256d
#define SSE_DOU __m128d


inline void thread_block(INT thread_id,INT start,INT end,INT start2,INT end2,INT * __restrict row_ptr,INT * __restrict col_idx,DOU * __restrict mtx_val,DOU * __restrict mtx_ans,DOU * __restrict mid_ans,DOU * __restrict vec_val)
{

        register INT start1,end1,num,Thread,i,j;
	register DOU sum;
	switch(start < end)
	{
		case true: {	
			mtx_ans[start]  = 0.0;
                        mtx_ans[end]    = 0.0;
                        start1          = row_ptr[start] + start2;
                        start++;
                        end1            = row_ptr[start];
                        Thread          = thread_id<<1;
                        sum             = 0.0;
                        #pragma unroll(8)
                        for(j=start1;j<end1;j++)
                        {
                                sum+=mtx_val[j]*vec_val[col_idx[j]];
                        }

                        mid_ans[Thread] = sum;
                        start1 = end1;

                        for(i=start;i<end;++i)
                        {
                                end1       = row_ptr[i+1];
                                sum        = 0.0;
				#pragma simd
                                for(j=start1;j<end1;j++)
                                {
                                        sum+=mtx_val[j]*vec_val[col_idx[j]];
                                }
                                mtx_ans[i] = sum;
                                start1 = end1;
                        }
                        start1 = row_ptr[end];
                        end1   = start1 + end2;
                        sum    = 0.0;
                        #pragma unroll(8)
                        for(j=start1;j<end1;j++)
                        {
                                sum+=mtx_val[j]*vec_val[col_idx[j]];
                        }
                        mid_ans[Thread | 1] = sum;
			return ;
		}
		default : {
			mtx_ans[start] = 0.0;
                        sum            = 0.0;
                        Thread         = thread_id<<1;
                        start1         = row_ptr[start] + start2;
                        end1           = row_ptr[end] + end2;
                        #pragma unroll(8)
                        for(j=start1;j<end1;j++)
                        {
                                sum+=mtx_val[j]*vec_val[col_idx[j]];
                        }
                        mid_ans[Thread]     = sum;
                        mid_ans[Thread | 1] = 0.0;
			return ;
		}
	}
}

inline INT binary_search(INT *&row_ptr,INT num,INT end)
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

inline void albus_balance(INT *&row_ptr,INT *&par_set,INT *&start,INT *&end,INT *&start1,INT *&end1,DOU *&mid_ans,INT thread_nums)
{
        register int tmp;
	start[0]            = 0;
	start1[0]           = 0;
        end[thread_nums-1]  = par_set[0];
	end1[thread_nums-1] = 0;
        INT tt=par_set[2]/thread_nums;
        for(INT i=1;i<thread_nums;i++)
        {
                tmp=tt*i;
                start[i]  = binary_search(row_ptr,tmp,par_set[0]);
		start1[i] = tmp - row_ptr[start[i]];
                end[i-1]  = start[i];
		end1[i-1] = start1[i];
        }
}

inline void SPMV_DOU(INT * __restrict row_ptr,INT * __restrict col_idx,DOU * __restrict mtx_val,INT * __restrict par_set,DOU * __restrict mtx_ans,DOU * __restrict vec_val,INT * __restrict start,INT * __restrict end,INT * __restrict start1,INT * __restrict end1,DOU * __restrict mid_ans, INT thread_nums)
{
        register INT i;
        #pragma omp parallel private(i)
        {
                #pragma omp for schedule(static) nowait
		for(i=0;i<thread_nums;++i)
                {
                        thread_block(i,start[i],end[i],start1[i],end1[i],row_ptr,col_idx,mtx_val,mtx_ans,mid_ans,vec_val);
                }
        }
        mtx_ans[0] = mid_ans[0];
        INT sub;
	#pragma unroll(32)
        for(i=1;i<thread_nums;++i)
        {
		sub = i<<1;
		register INT tmp1 = start[i];
		register INT tmp2 = end[i-1];
                if(tmp1 == tmp2)
                {
                        mtx_ans[tmp1] += (mid_ans[sub-1] + mid_ans[sub]);
                }
                else
                {
                        mtx_ans[tmp1] += mid_ans[sub];
                        mtx_ans[tmp2] += mid_ans[sub-1];
                }
        }
}
