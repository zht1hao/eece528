#include<stdio.h>
#include<math.h>
#include<time.h>
#include<immintrin.h>

//passed from compile time argument
//#define DATA_TYPE_SELECTED 1
//#define MATRIX_SIZE  4096 
//#define PRINT_OUTPUT
//#define DEBUG

#define ROW 	MATRIX_SIZE 
#define COL	MATRIX_SIZE

#define RANDOM 1
#define MAX_SIZE 4096
#define FROM_MASTER 1
#define FROM_WORKER 2 
#define MASTER 0

#if DATA_TYPE_SELECTED ==0
#define INTEGER int
#define DATA_TYPE INTEGER 
#define FORMAT_SPECIFIER "%d "
#elif DATA_TYPE_SELECTED==1
#define FLOAT float
#define DATA_TYPE FLOAT 
#define FORMAT_SPECIFIER "%f "
#elif DATA_TYPE_SELECTED==2
#define DOUBLE double
#define DATA_TYPE DOUBLE
#define FORMAT_SPECIFIER "%f "
#endif

DATA_TYPE *alloc_2d(int rows, int cols) {
   int i;
   DATA_TYPE *array= (DATA_TYPE *)_mm_malloc(cols*rows*sizeof(DATA_TYPE),256);
   return array;
}

//DATA_TYPE **array= (DATA_TYPE **)malloc(rows*sizeof(DATA_TYPE*));
//for (i=0; i<rows; i++)
//    array[i] = (DATA_TYPE *) _mm_malloc(sizeof(DATA_TYPE)*cols,256); 
//DATA_TYPE **alloc_2d(int rows, int cols) {
//    int i;
//    DATA_TYPE *data = (DATA_TYPE *)_mm_malloc(rows*cols*sizeof(DATA_TYPE),256);
//    DATA_TYPE **array= (DATA_TYPE **)malloc(rows*sizeof(DATA_TYPE*));
//    for (i=0; i<rows; i++)
//        array[i] = &(data[cols*i]);
//
//    return array;
//}

 int compare_float(double f1, double f2)
 {
  double precision = 0.00001;
  if (((f1 - precision) < f2) && 
      ((f1 + precision) > f2))
   {
    return 1;
   }
  else
   {
    return 0;
   }
 }


void copy_matrix(DATA_TYPE *A,DATA_TYPE *B,int row,int col)
{
	int i,j;
	for(i=0;i<row;i++)
	for(j=0;j<col;j++)
	B[i*col+j]=A[i*col+j];
}
inline void baseline_algorithm_without_vectorization( DATA_TYPE *A,int row,int col)
{
	int i,j,k;
	DATA_TYPE factor;
	for(i=0;i<row-1;i++)
	{
		for(j=i+1;j<row;j++)
		{
			factor=A[j*col+i]/A[i*col+i];
			
			for(k=0;k<col;k++)
				A[j*col+k]=A[j*col+k]-factor*A[i*col+k];
		}
	}
}
//void baseline_algorithm_with_vectorization(DATA_TYPE A[][MATRIX_SIZE+1],int row,int col)
inline void baseline_algorithm_with_vectorization(DATA_TYPE *A,int row,int col)
{
	int i,j,k;
	register DATA_TYPE factor;
	__m256 temp_load;
	__m256 temp_load2;
	__m256 temp_result;
	__m256 temp_result_inter;
	__m256 temp_factor;
 	const int adj=col%8;
	
	for(i=0;i<row-1;i++)
	{
		for(j=i+1;j<row;j++)
		{
			factor=A[j*col+i]/A[i*col+i];
			#ifdef FLOAT
			temp_factor=_mm256_set_ps(factor,factor,factor,factor,factor,factor,factor,factor);
			for(k=0;k<col;k=k+16)
			{
				temp_load= _mm256_load_ps(&A[j*col+k]);
				temp_load2=_mm256_load_ps(&A[i*col+k]);
				temp_result_inter=_mm256_mul_ps(temp_load2,temp_factor);	
				temp_result=_mm256_sub_ps(temp_load,temp_result_inter);
				_mm256_store_ps(&A[j*col+k],temp_result);

				temp_load= _mm256_load_ps(&A[j*col+k+8]);
				temp_load2=_mm256_load_ps(&A[i*col+k+8]);
				temp_result_inter=_mm256_mul_ps(temp_load2,temp_factor);	
				temp_result=_mm256_sub_ps(temp_load,temp_result_inter);
				_mm256_store_ps(&A[j*col+k+8],temp_result);
			}
			if(adj!=0)
			{
				for(k=col-adj;k<col;k++)
		      		{
		      		 A[j*col+k]=A[j*col+k]-factor*A[i*col+k];
		      		}
			}
			#endif
			#ifdef DOUBLE 
			for(k=0;k<col;k++)
			{
		      		A[j*col+k]=A[j*col+k]-factor*A[i*col+k];
			}
			#endif
		}
	}
	
}


//void initialize_matrix( DATA_TYPE A[][MATRIX_SIZE+1],int row,int col,int random)
void initialize_matrix( DATA_TYPE *A,int row,int col,int random)
{
	int i,j;
	#ifdef DEBUG
	printf("initializing matrix\n");
	#endif
	for(i=0;i<row;i++)
	for(j=0;j<col;j++)
	if(random)
		A[i*col+j]=(DATA_TYPE)(rand()%5000);//initializing within the random range of [0 5000]
	else 
		A[i*col+j]=0;
}

//void print_matrix(DATA_TYPE A[][MATRIX_SIZE+1],int row,int col,int id)
void print_matrix(DATA_TYPE *A,int row,int col)
{
	printf("\n");
	int i,j;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
			printf(FORMAT_SPECIFIER,A[i*col+j]);
		
		printf("\n");
	}	
}

int main(int argc)
{
	int i,j,k,index;
	int m1,n1;
	int averow,row_extra,rows,row_offset,tm_rows,tm_row_offset;
	int mtype,dest,source;
	DATA_TYPE factor_array[MATRIX_SIZE+1];
	int broadcasting_thread,index_consumed,row_consumed;
	int range;
	DATA_TYPE factor;
    	DATA_TYPE *A;
    	DATA_TYPE *B;
	clock_t start,end;

	
	m1=n1=MATRIX_SIZE;	
	A=alloc_2d(m1,n1);	
	B=alloc_2d(m1,n1);	
	initialize_matrix(A,m1,n1,RANDOM);
	#ifdef PRINT_OUTPUT
	print_matrix(A,m1,n1);
	#endif
	copy_matrix(A,B,m1,n1);
	if(MATRIX_SIZE%8==0)
	{
		start=clock();	
		baseline_algorithm_with_vectorization(A,m1,n1);
		end=clock();	
		printf("baseline algorithm with vectorization time =%ld milliseconds\n",(end-start)*1000/CLOCKS_PER_SEC);
	}
	#ifdef PRINT_OUTPUT
		print_matrix(A,m1,n1);
	#endif

	start=clock();	
	baseline_algorithm_without_vectorization(B,m1,n1);
	end=clock();	
	printf("baseline algorithm without vectorization time =%ld milliseconds\n",(end-start)*1000/CLOCKS_PER_SEC);
	#ifdef PRINT_OUTPUT
	print_matrix(B,m1,n1);
	#endif

}
