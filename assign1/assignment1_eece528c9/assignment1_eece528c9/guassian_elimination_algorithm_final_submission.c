#include <mpi.h>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<immintrin.h>

//passed from compile time argument
//DATA_TYPE_SELECTED=1 will select float matrix
//DATA_TYPE_SELECTED=2 will select doublematrix
//MATRIX_SIZE will select the matrix size
//PRINT_OUTPUT will print the output
//SANITY_CHECK will perform the sanity check

//passed from compile time argument
//#define DATA_TYPE_SELECTED 1
//#define MATRIX_SIZE 4096 
//#define DEBUG
//#define PRINT_OUTPUT
//#define SANITY_CHECK

#define RANDOM 1
#define MAX_SIZE 4096
#define FROM_MASTER 1
#define FROM_WORKER 2 
#define MASTER 0
#define ROW MATRIX_SIZE
#define COL MATRIX_SIZE+1

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
#define ALIGN



DATA_TYPE **alloc_2d(int rows, int cols) {
    int i;
    DATA_TYPE *data = (DATA_TYPE *)_mm_malloc(rows*cols*sizeof(DATA_TYPE),256);
    DATA_TYPE **array= (DATA_TYPE **)malloc(rows*sizeof(DATA_TYPE*));
    for (i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}

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


void copy_matrix(DATA_TYPE **A,DATA_TYPE **B,int row,int col)
{
	int i,j;
	for(i=0;i<row;i++)
	for(j=0;j<col;j++)
	B[i][j]=A[i][j];
}
void baseline_algorithm(DATA_TYPE **A,int row,int col)
{
	int i,j,k;
	DATA_TYPE factor;
	for(i=0;i<row-1;i++)
	{
		for(j=i+1;j<row;j++)
		{
			factor=A[j][i]/A[i][i];
			for(k=0;k<col;k++)
				A[j][k]=A[j][k]-factor*A[i][k];
		}
	}
	
}


//void initialize_matrix( DATA_TYPE A[][COL],int row,int col,int random)
void initialize_matrix( DATA_TYPE **A,int row,int col,int random)
{
	int i,j;
	#ifdef DEBUG
	printf("initializing matrix\n");
	#endif
	for(i=0;i<row;i++)
	for(j=0;j<col;j++)
	if(random)
			A[i][j]=(DATA_TYPE)(rand()%5000);//initializing within the random range of [0 5000]
	else 
		A[i][j]=0;
}

void print_matrix(DATA_TYPE A[][COL],int row,int col,int id)
//void print_matrix(DATA_TYPE **A,int row,int col,int id)
{
	printf("\n");
	int i,j;
	for(i=0;i<row;i++)
	{
		if(id==0)
			printf("MASTER:");
		else
			printf("SLAVE%d:",id);
		for(j=0;j<col;j++)
			printf(FORMAT_SPECIFIER,A[i][j]);
		
		printf("\n");
	}	
}
//void print_matrix(DATA_TYPE A[][COL],int row,int col,int id)
void print_matrix_2(DATA_TYPE **A,int row,int col,int id)
{
	printf("\n");
	int i,j;
	for(i=0;i<row;i++)
	{
		if(id==0)
			printf("MASTER:");
		else
			printf("SLAVE%d:",id);
		for(j=0;j<col;j++)
			printf(FORMAT_SPECIFIER,A[i][j]);
		
		printf("\n");
	}	
}
int main(int argc)
{
//	int source_completion_row;
	MPI_Status status;
	int i,j,k,index;
	int m1,n1;
	int averow,row_extra,rows,row_offset,tm_rows,tm_row_offset;
	int mtype,dest,source;
	DATA_TYPE factor_array[COL];
	int broadcasting_thread,index_consumed,row_consumed;
	int range;
	DATA_TYPE factor;
    	DATA_TYPE **MPI_OUTPUT;
    	DATA_TYPE **A;
	#ifdef SANITY_CHECK    	
	DATA_TYPE **baselineout;
	#endif
//	int barrier,BARRIER;
	clock_t start,end;

	//Initialize the MPI environment
    	MPI_Init(NULL, NULL);
	
    	// Get the number of processes
    	int num_tasks,numworkers;
    	MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);

    	// Get the rank of the process
    	int task_id;
    	MPI_Comm_rank(MPI_COMM_WORLD, &task_id);

    	// Get the name of the processor
    	char processor_name[MPI_MAX_PROCESSOR_NAME];
    	int name_len;
    	MPI_Get_processor_name(processor_name, &name_len);

	#ifdef SANITY_CHECK
	if(task_id==MASTER)
		printf("Time taken will be more since SANITY_CHECK define has been passed. To get the actual execution time run without passing sanity check\n");
	#endif
	numworkers=num_tasks;

	if(task_id==MASTER)
		start=clock();	
	
	m1=ROW;	
	//Note the last column will represent the right hand side of equation
	//Adding the last column in the matrix
	n1=COL;

	averow=m1/numworkers;
	row_extra=m1%numworkers;
	
	if(task_id<row_extra)
		rows=averow+1;
	else
		rows=averow;
	
	row_offset=task_id*averow;
	for(i=0;i<=task_id-1;i++)
	{
		if(i<row_extra)
		row_offset=row_offset+1;
	}
	
    	DATA_TYPE **B;
	B=alloc_2d(rows,COL);	
    	//DATA_TYPE B[rows][COL];

	if(task_id==MASTER)
	{
//		DATA_TYPE A[ROW][COL];
		#ifdef DEBUG
			printf("MASTER: m1=%d n1=%d \n ",m1,n1);
		#endif
		A=alloc_2d(m1,n1);	
		MPI_OUTPUT=alloc_2d(m1,n1);	
		initialize_matrix(A,m1,n1,RANDOM);//rand
		#ifdef PRINT_OUTPUT 
		print_matrix_2(A,m1,n1,task_id);
		#endif
		#ifdef SANITY_CHECK
		baselineout=alloc_2d(m1,n1);	
		copy_matrix(A,baselineout,m1,n1);
		//print_matrix_2(baselineout,m1,n1,task_id);
		#endif
		//Allocating work for master thread
		//Can be parallelize
		for(i=0;i<rows;i++)
		for(j=0;j<n1;j++)
		{
		   B[i][j]=A[i][j];
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(task_id==MASTER)
	{
		mtype=FROM_MASTER;

      		//Sending work to different thread
		for (dest=1; dest<numworkers; dest++)
      		{
      		   tm_rows = (dest+1 <= row_extra) ? averow+1 : averow;   	
		   tm_row_offset=dest*averow;
		   for(i=0;i<=dest-1;i++)
		   {
		   	if(i<row_extra)
		   	tm_row_offset=tm_row_offset+1;
		   }
		   #ifdef DOUBLE 
      		   	MPI_Send(&A[tm_row_offset][0],tm_rows*n1, MPI_DOUBLE, dest, mtype,  MPI_COMM_WORLD);
		   #endif
		   #ifdef  FLOAT
      		   	MPI_Send(&A[tm_row_offset][0],tm_rows*n1, MPI_FLOAT, dest, mtype,  MPI_COMM_WORLD);
		   #endif
		   #ifdef DEBUG
      		   	printf("MASTER: Sending %d row to task %d tm_row_offset=%d\n",tm_rows,dest,tm_row_offset);
		   #endif
		
		}
		free(A);

	}
	
	
	if(task_id>MASTER)
	{
     	 	mtype = FROM_MASTER;
		#ifdef DOUBLE
        	MPI_Recv(*B,rows*n1, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
		#endif
		#ifdef FLOAT 
        	MPI_Recv(*B,rows*n1, MPI_FLOAT, MASTER, mtype, MPI_COMM_WORLD, &status);
		#endif
	}



   	#ifdef DEBUG
		if(task_id==MASTER)
			printf("\nMASTER%d: Received %d row ,row_offset=%d",task_id,rows,row_offset);
		else
			printf("\nSLAVE%d: Received %d row ,row_offset=%d",task_id,rows,row_offset);
		print_matrix_2(B,rows,n1,task_id);
	#endif

	broadcasting_thread=0;
	index_consumed=0;
	if(broadcasting_thread<row_extra)
		row_consumed=averow+1;
	else 
		row_consumed=averow;
	
	for(index=0;index<m1-1;index++)
	{
	
		if(index>=index_consumed+row_consumed)
		{
				
			index_consumed+=row_consumed;
			broadcasting_thread++;
		
			if(broadcasting_thread<row_extra)
				row_consumed=averow+1;
			else 
				row_consumed=averow;
		}
	
		if(task_id==broadcasting_thread)
		{
			for(i=0;i<n1;i++)
			factor_array[i]=B[index-row_offset][i];
			
		#ifdef DEBUG
			if(task_id==MASTER)
				printf("MASTER=%d: row_offset=%d rows=%d index=%d broadcasting thread=%d\n",task_id,row_offset,rows,index,broadcasting_thread);	
			else			
				printf("SLAVE=%d: row_offset=%d rows=%d index=%d broadcasting thread=%d\n",task_id,row_offset,rows,index,broadcasting_thread);
		#endif	
		}
		#ifdef DOUBLE
			MPI_Bcast(&factor_array[0], n1, MPI_DOUBLE,broadcasting_thread,MPI_COMM_WORLD);
		#endif
		#ifdef FLOAT 
			MPI_Bcast(&factor_array[0], n1, MPI_FLOAT,broadcasting_thread,MPI_COMM_WORLD);
		#endif
		MPI_Barrier(MPI_COMM_WORLD);

		if(task_id==MASTER)
		{
			for(k=0;k<n1;k++)
				MPI_OUTPUT[index][k]=factor_array[k];
		}


		#ifdef DEBUG
			if(task_id==MASTER)
				printf("MASTER%d:INDEX=%dFACTOR_ARRAY ",task_id,index);
			else
				printf("SLAVE%d:INDEX=%dFACTOR_ARRAY ",task_id,index);
			for(i=0;i<n1;i++)
			printf(" %f",factor_array[i]); 	
			printf("\n");
		#endif
		



		if(index<row_offset+rows)
		{	
			if(row_offset<=index)
				range=index-row_offset+1;	
			else
				range=0;
			
			for(i=range;i<rows;i++)
			{
				factor=B[i][index]/factor_array[index];
		   	#ifdef DEBUG
			if(task_id==MASTER)
				printf("MASTER%d:INDEX=%d factor=%f\n",task_id,index,factor);
			else
				printf("SLAVE%d:INDEX=%d factor=%f\n",task_id,index,factor);
			#endif
					for(j=0;j<n1;j++)
					{
						B[i][j]=B[i][j]-factor*factor_array[j];	
					}
			}

		}
	
	  	#ifdef DEBUG
			print_matrix_2(B,rows,n1,task_id);
		#endif

			
	}
	mtype=FROM_WORKER;
	if(task_id==numworkers-1){
		#ifdef DOUBLE
			MPI_Send(&B[averow-1][0], n1, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
		#endif
		#ifdef FLOAT 
			MPI_Send(&B[averow-1][0], n1, MPI_FLOAT, MASTER, mtype, MPI_COMM_WORLD);
		#endif
	}
	if(task_id==MASTER)
	{
		#ifdef DOUBLE
        		MPI_Recv(&MPI_OUTPUT[m1-1][0],n1, MPI_DOUBLE, numworkers-1, mtype, MPI_COMM_WORLD, &status);
		#endif
		#ifdef FLOAT 
        		MPI_Recv(&MPI_OUTPUT[m1-1][0],n1, MPI_FLOAT, numworkers-1, mtype, MPI_COMM_WORLD, &status);
		#endif
	
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(task_id==MASTER)
	{
		#ifdef PRINT_OUTPUT 
		print_matrix_2(MPI_OUTPUT,m1,n1,task_id);
		#endif
  		#ifdef DEBUG
		for(i=0;i<m1;i++)
		{
			for(j=0;j<n1;j++)
			printf("%f ",MPI_OUTPUT[i][j]);
			printf("\n");
		}
		#endif
		end=clock();
		printf("MPI time with %d thread =%ld miliseconds\n",num_tasks,(end-start)*1000/CLOCKS_PER_SEC);
		
		#ifdef SANITY_CHECK	
		baseline_algorithm(baselineout,m1,n1);
		int error=0;
		//print_matrix_2(baselineout,m1,n1,0);
		for(i=0;i<m1;i++)
		for(j=0;j<n1;j++)
		if(!compare_float(baselineout[i][j],MPI_OUTPUT[i][j]))	
		{
			printf("ERROR baselineout[%d][%d]=%f MPI_OUTPUT[%d][%d]=%f\n",i,j,baselineout[i][j],i,j,MPI_OUTPUT[i][j]);
			error=1;
		}
		if(!error)
		printf("SANITY_CHECK PASSED\n");
		free(MPI_OUTPUT);
		#endif
	}
  	#ifdef DEBUG
		printf("TASK_ID%d Completed Successfully\n",task_id);
	#endif
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	

}
