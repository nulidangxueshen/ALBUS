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

inline DOU SIMD_fast1(INT start1,INT num,INT * __restrict row_ptr,INT * __restrict col_idx,DOU * __restrict mtx_val,DOU * __restrict vec_val)
{
	DOU answer;
	switch(num)
	{
		case 4 : {
			register SSE_DOU mtx_3 , vec_3 , ans_3 , mtx_3_1 , vec_3_1;
			register INT s1,s2,s3;
			s1      = start1 + 1;
			s2      = start1 + 2;
			s3      = start1 + 3;
			mtx_3   = _mm_load_pd(mtx_val+start1);
			mtx_3_1 = _mm_load_pd(mtx_val+s2);
			vec_3   = _mm_set_pd(vec_val[col_idx[s1]],vec_val[col_idx[start1]]);
			vec_3_1 = _mm_set_pd(vec_val[col_idx[s3]],vec_val[col_idx[s2]]);
			ans_3   = _mm_fmadd_pd(mtx_3_1,vec_3_1,_mm_mul_pd(mtx_3,vec_3));
			answer  = ans_3[0]+ans_3[1];
			return answer;
		}
		default : {
			register AVX_DOU mtx_ans_1,mtx_3,vec_3;
			register INT s1,s2,s3;
			register INT t      = num & (~3);
			register INT start2 = start1 + t;
			register INT num_1  = num & 3;
			s1 = start1 + 1;
			s2 = start1 + 2;
			s3 = start1 + 3;
			_mm_prefetch((DOU *)&mtx_val[start1+16],_MM_HINT_T0);
			_mm_prefetch((DOU *)&col_idx[start1+16],_MM_HINT_T0);
			mtx_3     = _mm256_load_pd(mtx_val+start1);
			vec_3     = _mm256_set_pd(vec_val[col_idx[s3]],vec_val[col_idx[s2]],vec_val[col_idx[s1]],vec_val[col_idx[start1]]);
			mtx_ans_1 = _mm256_mul_pd(mtx_3,vec_3);
			start1 += 4;
			#pragma unroll(32)
			for(;start1<start2;start1+=4)
			{
				s1        = start1 + 1;
				s2        = start1 + 2;
				s3        = start1 + 3;
				mtx_3     = _mm256_load_pd(mtx_val+start1);
				vec_3     = _mm256_setr_pd(vec_val[col_idx[start1]],vec_val[col_idx[s1]],vec_val[col_idx[s2]],vec_val[col_idx[s3]]);
				mtx_ans_1 = _mm256_fmadd_pd(mtx_3,vec_3,mtx_ans_1);
			}
			switch (num_1)
			{
				case 0 : {
					mtx_ans_1 = _mm256_hadd_pd(mtx_ans_1,mtx_ans_1);
					answer    = mtx_ans_1[0] + mtx_ans_1[2];
					return answer;
				}
				case 1 : {
					mtx_ans_1 = _mm256_hadd_pd(mtx_ans_1,mtx_ans_1);
					answer    = mtx_ans_1[0] + mtx_ans_1[2];
					answer    = answer + (mtx_val[start2]*vec_val[col_idx[start2]]);
					return answer;
				}
				case 2 : {
					mtx_ans_1 = _mm256_hadd_pd(mtx_ans_1,mtx_ans_1);
					answer    = mtx_ans_1[0] + mtx_ans_1[2];
					s1        = start2 + 1;
					answer    = answer+(mtx_val[start2]*vec_val[col_idx[start2]]+mtx_val[s1]*vec_val[col_idx[s1]]);
					return answer;
				}
				default : {
					s1        = start2 + 1;
					s2        = start2 + 2;
					mtx_3     = _mm256_load_pd(mtx_val+start2);
					vec_3     = _mm256_set_pd(0,vec_val[col_idx[s2]],vec_val[col_idx[s1]],vec_val[col_idx[start2]]);
					mtx_ans_1 = _mm256_fmadd_pd(mtx_3,vec_3,mtx_ans_1);
					mtx_ans_1 = _mm256_hadd_pd(mtx_ans_1,mtx_ans_1);
					answer    = mtx_ans_1[0] + mtx_ans_1[2];
					return answer;
				}
			}
		}
	}
}

inline DOU SIMD_fast2(INT start1,INT num,INT * __restrict__ row_ptr,INT * __restrict__ col_idx,DOU * __restrict__ mtx_val,DOU * __restrict__ vec_val)
{
	register DOU answer;
	switch(num)
        {
		case 0 : return 0;
		case 1 : {
			answer = mtx_val[start1] * vec_val[col_idx[start1]];
			return answer;
		}
		case 2 : {
			register INT s1;
			s1     = start1 + 1;
			answer = mtx_val[start1] * vec_val[col_idx[start1]] + mtx_val[s1] * vec_val[col_idx[s1]];
			return answer;
		}
                case 3 : {
			register SSE_DOU mtx_3 , vec_3 , ans_3 , mtx_3_1 , vec_3_1;
                        register INT s1,s2;
                        s1      = start1 + 1;
                        s2      = start1 + 2;
                        mtx_3   = _mm_load_pd(mtx_val+start1);
                        mtx_3_1 = _mm_load_pd(mtx_val+s2);
                        vec_3   = _mm_set_pd(vec_val[col_idx[s1]],vec_val[col_idx[start1]]);
                        vec_3_1 = _mm_set_pd(0,vec_val[col_idx[s2]]);
                        ans_3   = _mm_fmadd_pd(mtx_3_1,vec_3_1,_mm_mul_pd(mtx_3,vec_3));
                        answer  = ans_3[0] + ans_3[1];
                        return answer;
		}
	}
}

inline DOU calculation(INT start1,INT num,INT * __restrict row_ptr,INT * __restrict col_idx,DOU * __restrict mtx_val,DOU * __restrict vec_val)
{
	if(num>=4)
	{
		return SIMD_fast1(start1,num,row_ptr,col_idx,mtx_val,vec_val);
	}
	else
	{
		return SIMD_fast2(start1,num,row_ptr,col_idx,mtx_val,vec_val);
	}
}

inline void thread_block(INT thread_id,INT start,INT end,INT start2,INT end2,INT * __restrict row_ptr,INT * __restrict col_idx,DOU * __restrict mtx_val,DOU * __restrict mtx_ans,DOU * __restrict mid_ans,DOU * __restrict vec_val)
{

        register INT start1,end1,num,Thread,i;
	register DOU sum;
	switch(start < end)
	{
		case true: {
			mtx_ans[start]  = 0.0;
			mtx_ans[end]    = 0.0;
        		start1          = row_ptr[start] + start2;
			start++;
        		end1            = row_ptr[start];
        		num             = end1 - start1;
        		Thread          = thread_id<<1;
        		mid_ans[Thread] = calculation(start1,num,row_ptr,col_idx,mtx_val,vec_val);
			start1 = end1;
	
			#pragma simd
        		for(i=start;i<end;++i)
        		{
               		 	end1       = row_ptr[i+1];
                		num        = end1 - start1;
                		sum        = calculation(start1,num,row_ptr,col_idx,mtx_val,vec_val);
				mtx_ans[i] = sum;
				start1 = end1;
        		}
			start1 = row_ptr[end];
        		end1   = start1 + end2;
        		mid_ans[Thread | 1] = calculation(start1,end2,row_ptr,col_idx,mtx_val,vec_val);
			return ;
		}
		default : {
			mtx_ans[start]      = 0.0;
			Thread              = thread_id<<1;
			start1              = row_ptr[start] + start2;
			num                 = end2 - start2;
			mid_ans[Thread]     = calculation(start1,num,row_ptr,col_idx,mtx_val,vec_val);
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
