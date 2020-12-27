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
#define FLO float

void ReadFile(FILE * fp_mtx,INT *&row_ptr,INT *&col_idx,DOU *&mtx_val,DOU *&vec_val,INT *&par_set)
{
        INT i,row_num=0,col_num=0,nzz_num=0,v=-1,cnt=0,row,col,cnt1=0,nzz_int;
        char s[1024];
        string str[6];
        while((fgets(s,1024,fp_mtx))!=NULL)
        {
                if(cnt1==0)
                {
                        cnt=0;
                        for(i=0;i<strlen(s)-1;i++)
                        {
                                if(s[i]==' ')
                                {
                                        cnt++;
                                        continue;
                                }
                                str[cnt].push_back(s[i]);
                        }
                        cnt1=1;
                }
                if(s[0]>='0'&&s[0]<='9')
                {
                        cnt=0;
                        for(i=0;i<strlen(s)-1;i++)
                        {
                                if(s[i]==' ')
                                {
                                        cnt++;
                                        continue;
                                }
                                if(cnt==0)
                                {
                                        row_num=row_num*10+s[i]-'0';
                                }
                                if(cnt==1)
                                {
                                        col_num=col_num*10+s[i]-'0';
                                }
                                if(cnt==2)
                                {
                                        nzz_num=nzz_num*10+s[i]-'0';
                                }
                        }
                        break;
                }
        }
        cnt=0;
        printf("(row_num : %d) , (col_num : %d) , (nzz_num : %d)\n",row_num,col_num,nzz_num);
        cout<<"( Matrix Kind : "<<str[3]<<" ) , ( Data type : "<<str[4]<<" )"<<endl;
        DOU nzz,nzz1;
        par_set = (INT *)malloc(sizeof(INT)*10);
        INT * row_ptr1 = (INT *)aligned_alloc(64,sizeof(INT)*(row_num+1));
        INT * col_idx1 = (INT *)aligned_alloc(64,sizeof(INT)*nzz_num);
        INT * row_idx1 = (INT *)aligned_alloc(64,sizeof(INT)*nzz_num);
        DOU * mtx_val1 = (DOU *)aligned_alloc(64,sizeof(DOU)*nzz_num);
        row_ptr = (INT *)aligned_alloc(64,sizeof(INT)*(row_num+1));
        memset(row_ptr1,0,sizeof(INT)*(row_num+1));
        memset(row_ptr ,0,sizeof(INT)*(row_num+1));
        memset(mtx_val1,0,sizeof(DOU)*nzz_num);
        memset(col_idx1,0,sizeof(INT)*nzz_num);
        for(i=0;i<nzz_num;i++)
        {
                if(str[3]=="real")
                {
                        fscanf(fp_mtx,"%d %d %lg",&row,&col,&nzz);
                }
                if(str[3]=="pattern")
                {
                        fscanf(fp_mtx,"%d %d",&row,&col);
                        nzz=1.0;
                }
                if(str[3]=="integer")
                {
                        fscanf(fp_mtx,"%d %d %d",&row,&col,&nzz_int);
                        nzz = nzz_int * 1.0;
                }
                row--;
                col--;
                row_ptr1[row]++;
                row_idx1[i] = row;
                col_idx1[i] = col;
                mtx_val1[i] = nzz;
        }
        if(str[4]=="symmetric")
        {
                for(i=0;i<nzz_num;i++)
                {
                        if(row_idx1[i]!=col_idx1[i])
                        {
                                row_ptr1[col_idx1[i]]++;
                        }
                }
        }
        INT NZZ_NUM,num1,num2;
        num1 = row_ptr1[0];
        row_ptr1[0] = 0;
        for(i=1;i<=row_num;i++)
        {
                num2 = num1;
                num1 += row_ptr1[i];
                row_ptr1[i] = num2;
        }
        NZZ_NUM = num1;
        memcpy(row_ptr,row_ptr1,(row_num+1)*sizeof(INT));
        memset(row_ptr1,0,sizeof(INT)*(row_num+1));
        col_idx = (INT *)aligned_alloc(64,sizeof(INT)*NZZ_NUM);
        mtx_val = (DOU *)aligned_alloc(64,sizeof(DOU)*NZZ_NUM);
	vec_val = (DOU *)aligned_alloc(64,sizeof(DOU)*col_num);
        memset(mtx_val,0,sizeof(DOU)*NZZ_NUM);
        memset(col_idx,0,sizeof(INT)*NZZ_NUM);
        if(str[4]=="symmetric")
        {
                for(i=0;i<nzz_num;i++)
                {
                        if(row_idx1[i]!=col_idx1[i])
                        {
                                INT offset      = row_ptr[row_idx1[i]] + row_ptr1[row_idx1[i]];
                                col_idx[offset] = col_idx1[i];
                        	mtx_val[offset] = mtx_val1[i];
                        	row_ptr1[row_idx1[i]]++;

                                offset          = row_ptr[col_idx1[i]] + row_ptr1[col_idx1[i]];
                                col_idx[offset] = row_idx1[i];
                                mtx_val[offset] = mtx_val1[i];
                                row_ptr1[col_idx1[i]]++;
                        }
                        else
                        {
                                INT offset      = row_ptr[row_idx1[i]] + row_ptr1[row_idx1[i]];
                                col_idx[offset] = col_idx1[i];
                        	mtx_val[offset] = mtx_val1[i];
                        	row_ptr1[row_idx1[i]]++;
                        }
                }
        }
        else
        {
                for(i=0;i<nzz_num;i++)
                {
                        INT offset      = row_ptr[row_idx1[i]] + row_ptr1[row_idx1[i]];
			col_idx[offset] = col_idx1[i];
			mtx_val[offset] = mtx_val1[i];
                        row_ptr1[row_idx1[i]]++;
                }
        }
        free(row_idx1);
        free(row_ptr1);
        free(mtx_val1);
        free(col_idx1);
        nzz_num = NZZ_NUM;
        srand(time(NULL));
        for(i=0;i<col_num;i++)
        {
                vec_val[i] = ((rand()%10)+12)*1.0;
        }
        //----------------------------------------------------------------//
        cout<<"NZZ_NUM : "<<nzz_num<<"      (Note: If the matrix is a symmetric, NZZ_NUM <= 2 * nzz_num)"<<endl;
        par_set[0] = row_num;
        par_set[1] = col_num;
        par_set[2] = nzz_num;
	par_set[8] = 4 * min(200000ull, max(100ull, ((16ull << 30) / nzz_num)));
	cout<<par_set[8]<<endl;
}
