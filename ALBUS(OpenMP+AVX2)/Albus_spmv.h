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

DOU calculation(INT start1,INT end1,INT num,INT *&row_ptr,INT *&col_idx,DOU *&mtx_val,DOU *&vec_val)
{
        DOU answer=0;
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
                        AVX_DOU mtx_3;
                        AVX_DOU vec_3;
                        AVX_DOU ans_3;
                        INT s1,s2;
                        s1=start1+1;
                        s2=start1+2;
                        mtx_3 = _mm256_load_pd(mtx_val+start1);
                        vec_3 = _mm256_set_pd(0,vec_val[col_idx[s2]],vec_val[col_idx[s1]],vec_val[col_idx[start1]]);
                        ans_3 = _mm256_mul_pd(mtx_3,vec_3);
                        ans_3 = _mm256_hadd_pd(ans_3,ans_3);
                        answer = ans_3[0]+ans_3[2];
                        break;
                }
                case 4 : {
                        AVX_DOU mtx_3 , vec_3 , ans_3;
                        INT s1,s2,s3;
                        s1=start1+1;
                        s2=start1+2;
                        s3=start1+3;
                        mtx_3 = _mm256_load_pd(mtx_val+start1);
                        vec_3 = _mm256_set_pd(vec_val[col_idx[s3]],vec_val[col_idx[s2]],vec_val[col_idx[s1]],vec_val[col_idx[start1]]);
                        ans_3 = _mm256_mul_pd(mtx_3,vec_3);
                        ans_3 = _mm256_hadd_pd(ans_3,ans_3);
                        answer = ans_3[0]+ans_3[2];
                        break;
                }
                default : {
                        INT s1,s2,s3,j;
                        INT j_end = start1 + num/4*4;
                        INT num_1 = num % 4;
                        INT start2 = j_end;
                        AVX_DOU mtx_ans_1,mtx_3,vec_3;
                        mtx_ans_1 = _mm256_setzero_pd();
                        j=start1;
                        for(;j<j_end;j+=4)
                        {
                                s1=j+1;
                                s2=j+2;
                                s3=j+3;
                                mtx_3 = _mm256_load_pd(mtx_val+j);
                        vec_3 = _mm256_set_pd(vec_val[col_idx[s3]],vec_val[col_idx[s2]],vec_val[col_idx[s1]],vec_val[col_idx[j]]);
                                mtx_ans_1 = _mm256_fmadd_pd(mtx_3,vec_3,mtx_ans_1);
                        }
                        switch (num_1)
                        {
                                case 0 : {
                                        mtx_ans_1 = _mm256_hadd_pd(mtx_ans_1,mtx_ans_1);
                                        answer = mtx_ans_1[0]+mtx_ans_1[2];
                                        break;
                                }
                                case 1 : {
                                        mtx_ans_1 = _mm256_hadd_pd(mtx_ans_1,mtx_ans_1);
                                        answer = mtx_ans_1[0]+mtx_ans_1[2];
                                        answer=answer+(mtx_val[start2]*vec_val[col_idx[start2]]);
                                        break;
                                }
                                case 2 : {
                                        mtx_ans_1 = _mm256_hadd_pd(mtx_ans_1,mtx_ans_1);
                                        answer = mtx_ans_1[0]+mtx_ans_1[2];
                                        s1 = start2+1;
                                        answer=answer+(mtx_val[start2]*vec_val[col_idx[start2]]+mtx_val[s1]*vec_val[col_idx[s1]]);
                                        break;
                                }
                                case 3 : {
                                        s1=start2+1;
                                        s2=start2+2;
                                        mtx_3 = _mm256_load_pd(mtx_val+start2);
                                        vec_3 = _mm256_set_pd(0,vec_val[col_idx[s2]],vec_val[col_idx[s1]],vec_val[col_idx[start2]]);
                                        mtx_ans_1 = _mm256_fmadd_pd(mtx_3,vec_3,mtx_ans_1);
                                        mtx_ans_1 = _mm256_hadd_pd(mtx_ans_1,mtx_ans_1);
                                        answer = mtx_ans_1[0]+mtx_ans_1[2];
                                }
                        }
                }
        }
        return answer;
}

void thread_block(INT thread_id,INT start,INT end,INT num2,INT num1,INT *&row_ptr,INT *&col_idx,DOU *&mtx_val,DOU *&mtx_ans,DOU *&mid_ans,DOU *&vec_val)
{
        INT start1,end1,num,Thread,i;
        start1 = num2;
        end1   = row_ptr[start];
        num    = end1 - start1;
        Thread = thread_id<<1;
        mid_ans[Thread] = calculation(start1,end1,num,row_ptr,col_idx,mtx_val,vec_val);
        #pragma prefetch
        for(i=start;i<end;++i)
        {
                start1 = row_ptr[i];
                end1   = row_ptr[i+1];
                num    = end1 - start1;
                mtx_ans[i] = calculation(start1,end1,num,row_ptr,col_idx,mtx_val,vec_val);
        }
        start1 = row_ptr[end];
        end1   = num1;
        num    = end1 - start1;
        mid_ans[Thread + 1] = calculation(start1,end1,num,row_ptr,col_idx,mtx_val,vec_val);
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
