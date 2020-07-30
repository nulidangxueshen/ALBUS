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

inline DOU SIMD_fast1(INT start1,INT end1,INT num,INT * __restrict row_ptr,INT * __restrict col_idx,DOU * __restrict mtx_val,DOU * __restrict vec_val)
{
	DOU answer;
	switch(num)
                {
                        case 4 : {
                                SSE_DOU mtx_3 , vec_3 , ans_3 , mtx_3_1 , vec_3_1;
                                INT s1,s2,s3;
                                s1=start1+1;
                                s2=start1+2;
                                s3=start1+3;
                                mtx_3   = _mm_load_pd(mtx_val+start1);
                                mtx_3_1 = _mm_load_pd(mtx_val+s2);
                                vec_3   = _mm_set_pd(vec_val[col_idx[s1]],vec_val[col_idx[start1]]);
                                vec_3_1 = _mm_set_pd(vec_val[col_idx[s3]],vec_val[col_idx[s2]]);
                                ans_3   = _mm_fmadd_pd(mtx_3_1,vec_3_1,_mm_mul_pd(mtx_3,vec_3));
                                answer  = ans_3[0]+ans_3[1];
                                return answer;
                        }
                        default : {
                                AVX_DOU mtx_ans_1,mtx_3,vec_3;
                                register INT s1,s2,s3;
                                register INT t      = (num>>2)<<2;
                                register INT j_end  = start1 + t;
                                register INT num_1  = num&3;
                                register INT start2 = j_end;
                                s1=start1+1;
                                s2=start1+2;
                                s3=start1+3;
                                _mm_prefetch((DOU *)&mtx_val[start1+16],_MM_HINT_T0);
                                _mm_prefetch((DOU *)&col_idx[start1+16],_MM_HINT_T0);
                                mtx_3     = _mm256_load_pd(mtx_val+start1);
                                vec_3     = _mm256_set_pd(vec_val[col_idx[s3]],vec_val[col_idx[s2]],vec_val[col_idx[s1]],vec_val[col_idx[start1]]);
                                mtx_ans_1 = _mm256_mul_pd(mtx_3,vec_3);
                                start1+=4;
                                #pragma unroll(8)
                                for(;start1<j_end;start1+=4)
                                {
                                        s1=start1+1;
                                        s2=start1+2;
                                        s3=start1+3;
                                        mtx_3     = _mm256_load_pd(mtx_val+start1);
                                        vec_3     = _mm256_set_pd(vec_val[col_idx[s3]],vec_val[col_idx[s2]],vec_val[col_idx[s1]],vec_val[col_idx[start1]]);
                                        mtx_ans_1 = _mm256_fmadd_pd(mtx_3,vec_3,mtx_ans_1);
                                }
				switch (num_1)
                                {
                                        case 0 : {
                                                mtx_ans_1 = _mm256_hadd_pd(mtx_ans_1,mtx_ans_1);
                                                answer    = mtx_ans_1[0]+mtx_ans_1[2];
                                                break;
                                        }
                                        case 1 : {
                                                mtx_ans_1 = _mm256_hadd_pd(mtx_ans_1,mtx_ans_1);
                                                answer    = mtx_ans_1[0]+mtx_ans_1[2];
                                                answer    = answer + (mtx_val[start2]*vec_val[col_idx[start2]]);
                                                break;
                                        }
                                        case 2 : {
                                                mtx_ans_1 = _mm256_hadd_pd(mtx_ans_1,mtx_ans_1);
                                                answer    = mtx_ans_1[0]+mtx_ans_1[2];
                                                s1        = start2+1;
                                                answer    = answer+(mtx_val[start2]*vec_val[col_idx[start2]]+mtx_val[s1]*vec_val[col_idx[s1]]);
                                                break;
                                        }
                                        case 3 : {
                                                s1  = start2+1;
                                                s2  = start2+2;
                                                mtx_3     = _mm256_load_pd(mtx_val+start2);
                                                vec_3     = _mm256_set_pd(0,vec_val[col_idx[s2]],vec_val[col_idx[s1]],vec_val[col_idx[start2]]);
                                                mtx_ans_1 = _mm256_fmadd_pd(mtx_3,vec_3,mtx_ans_1);
                                                mtx_ans_1 = _mm256_hadd_pd(mtx_ans_1,mtx_ans_1);
                                                answer = mtx_ans_1[0]+mtx_ans_1[2];
                                        }
                                }
                        }
                }
	return answer;
}

inline DOU SIMD_fast2(INT start1,INT end1,INT num,INT * __restrict__ row_ptr,INT * __restrict__ col_idx,DOU * __restrict__ mtx_val,DOU * __restrict__ vec_val)
{
	DOU answer;
	switch(num)
                {
                        case 0 : break;
                        case 1 : {
                                answer=mtx_val[start1]*vec_val[col_idx[start1]];
                                break;
                        }
                        case 2 : {
                                INT s1;
                                s1 = start1+1;
                                answer=mtx_val[start1]*vec_val[col_idx[start1]]+mtx_val[s1]*vec_val[col_idx[s1]];
                                break;
                        }
                        case 3 : {
                                SSE_DOU mtx_3 , vec_3 , ans_3 , mtx_3_1 , vec_3_1;
                                INT s1,s2;
                                s1=start1+1;
                                s2=start1+2;
                                mtx_3   = _mm_load_pd(mtx_val+start1);
                                mtx_3_1 = _mm_load_pd(mtx_val+s2);
                                vec_3   = _mm_set_pd(vec_val[col_idx[s1]],vec_val[col_idx[start1]]);
                                vec_3_1 = _mm_set_pd(0,vec_val[col_idx[s2]]);
                                ans_3   = _mm_fmadd_pd(mtx_3_1,vec_3_1,_mm_mul_pd(mtx_3,vec_3));
                                answer  = ans_3[0]+ans_3[1];
                                return answer;
                        }
                }
	return answer;
}

inline DOU calculation(INT start1,INT end1,INT num,INT * __restrict row_ptr,INT * __restrict col_idx,DOU * __restrict mtx_val,DOU * __restrict vec_val)
{
	if(num>=4)
	{
		return SIMD_fast1(start1,end1,num,row_ptr,col_idx,mtx_val,vec_val);
	}
	else
	{
		return SIMD_fast2(start1,end1,num,row_ptr,col_idx,mtx_val,vec_val);
	}
}

inline void thread_block(INT thread_id,INT start,INT end,INT num2,INT num1,INT * __restrict row_ptr,INT * __restrict col_idx,DOU * __restrict mtx_val,DOU * __restrict mtx_ans,DOU * __restrict mid_ans,DOU * __restrict vec_val)
{
        register INT start1,end1,num,Thread,i;
	DOU sum;
        start1 = num2;
        end1   = row_ptr[start];
        num    = end1 - start1;
        Thread = thread_id<<1;
        mid_ans[Thread] = calculation(start1,end1,num,row_ptr,col_idx,mtx_val,vec_val);
	start1 = end1;
	//#pragma unroll(8)
        for(i=start;i<end;++i)
        {
		//_mm_prefetch((char *)&mtx_val[i+32],_MM_HINT_T2);
		//_mm_prefetch((char *)&vec_val[i+64],_MM_HINT_T0);
		//_mm_prefetch((char *)&col_idx[i+32],_MM_HINT_T2);
                end1   = row_ptr[i+1];
                num    = end1 - start1;
                sum        = calculation(start1,end1,num,row_ptr,col_idx,mtx_val,vec_val);
		mtx_ans[i] = sum;
		start1 = end1;
        }
        start1 = row_ptr[end];
        end1   = num1;
        num    = end1 - start1;
        mid_ans[Thread + 1] = calculation(start1,end1,num,row_ptr,col_idx,mtx_val,vec_val);
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

inline void albus_balance(INT *&row_ptr,INT *&par_set,INT *&start,INT *&end,DOU *&mid_ans,INT *&block_size,INT thread_nums)
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

inline void SPMV_DOU(INT * __restrict row_ptr,INT * __restrict col_idx,DOU * __restrict mtx_val,INT * __restrict par_set,DOU * __restrict mtx_ans,DOU * __restrict vec_val,INT * __restrict start,INT * __restrict end,DOU * __restrict mid_ans,INT * __restrict block_size,INT thread_nums)
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
