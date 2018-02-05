/***************************************************************************************/
/*Assignment No: 03                                                                    */
/*Program      : Gaussian Elimination using CUDA                                       */
/*Author       : Krupa Harthi, EECE528C1                                               */
/***************************************************************************************/

/**********Include Libraries************************************************************/
#include <stdio.h>
#include <cuda.h>
#include <time.h>
#include <stdlib.h>


/**********Kernel Logic for CUDA GPU***************************************************/
__global__ void gaus_row(float **d_in,int j,int n, int thread_slab)
{
int i = threadIdx.x;
float C_Mul;
int k;

i = i + thread_slab;
	
	if(i>j)
	{
		C_Mul = d_in[i][j]/d_in[j][j];
		
		for (k=0;k<n+1;k++)
		{
			d_in[i][k]=d_in[i][k]-(C_Mul*d_in[j][k]);
		}
	}
}

/**********Kernel Logic for CUDA GPU***************************************************/
int main (int argc, char **argv)
{
int ROW_SIZE, COL_SIZE, i, j, k=0, n=0;
float **h_A=NULL;
float *h_x=NULL;
float **h_A1=NULL;
float sum;
srand(time(NULL));
clock_t t1,t2;
int nDevices, GridSize, Max_Block = 512, blocksize = 0;

/**This Sample Data can be uncommented for testing results for small matrix size*******/
//float h_input[20]= {2,1,-1,2,5,4,5,-3,6,9,-2,5,-2,6,4,4,11,-4,8,2};  
//expect result  x(3) = 3, x(2) = 1, x(1) = -1, x(0) = 1 */
//float h_input[16]= {2,1,-1,2,5,4,5,-3,6,9,-2,5,-2,6,4,4};



n = atoi(argv[1]);

/***Set up Matrix size based on the input value n***************************************/

ROW_SIZE = n * sizeof(float);
COL_SIZE = n * sizeof(float);

/***Grid Size Calculation***************************************************************/

GridSize = (n + Max_Block - 1) / Max_Block;
printf("grid size = %d\n", GridSize);

/***cuda verfication and resource details***********************************************/

if (n >= Max_Block)
	blocksize = Max_Block;
else
	blocksize = n;

cudaError_t err = cudaGetDeviceCount(&nDevices);
  if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));


h_A = (float**)malloc(sizeof(float*)*n);
h_x = (float*)malloc(COL_SIZE);
h_A1= (float**)malloc(sizeof(float*)*n);

/***Allocating the array****************************************************************/
for(i=0; i<n; i++)
{
	h_A[i] = (float*)malloc(ROW_SIZE);
    h_A1[i]= (float*)malloc(ROW_SIZE);
}

/***Assignment and prtint Matrix********************************************************/
printf("\n ----------------------------------------------");
for(i=0; i<n; i++)
{
	if (n<10)
	printf("\n");
    for(j=0; j<n; j++)
    {
		h_A[i][j] = ((float)rand()/(float)(RAND_MAX)) * 5.0;

		//Comment above line and Uncomment this for testing with small matrix
		//h_A[i][j] = h_input[k];
		
		k++;
		
		if(n<10)
		printf("%.1f\t",h_A[i][j]);
	}
}

k=0;
printf("\n ----------------------------------------------");

/***declare GPU memory pointers************************************************************/
float **d_array_in;
float *d_in[n];

cudaMalloc((void **)&d_array_in, (n * sizeof(float *)));

/***Allocate memory on CPU****************************************************************/
for (i=0;i<n;i++)
{
	cudaMalloc((void **)&d_in[i], ROW_SIZE);
}

/***Transfer the array to the GPU*********************************************************/
for (i=0;i<n;i++)
{
	cudaMemcpy(d_in[i],h_A[i], ROW_SIZE, cudaMemcpyHostToDevice);
}

cudaMemcpy(d_array_in, d_in, (n * sizeof(float *)),cudaMemcpyHostToDevice);

//Starting the clock
t1= clock();

/***GE loop for the generation of upper triangular matrix*********************************/
for(j=0; j<n; j++)         
{
    // launch the kernel
    gaus_row<<<GridSize, blocksize >>>(d_array_in,j,n,0);
    cudaDeviceSynchronize();
}

//Ending the clock
t2 = clock();

/* cuda errors */
cudaError_t errSync  = cudaGetLastError();
cudaError_t errAsync = cudaDeviceSynchronize();

if (errSync != cudaSuccess)
  printf("Sync kernel error: %s\n", cudaGetErrorString(errSync));
if (errAsync != cudaSuccess)
  printf("Async kernel error: %s\n", cudaGetErrorString(errAsync));

/* time taken for parallel process */
printf("\nThe time taken is %f seconds\n",(float)(t2-t1)/CLOCKS_PER_SEC);



/***copy back the result array to the CPU***********************************************/
for(i=0;i<n;i++)
{
	cudaMemcpy(h_A1[i],d_in[i], ROW_SIZE, cudaMemcpyDeviceToHost);
}

/***Print the output matrix*************************************************************/
printf("\n ----------------------------------------------");
if(n<10)
{
	for(i=0; i<n; i++) 
	{
		printf("\n");
   		for(j=0; j<n; j++)
   		{
   			printf("%.1f\t",h_A1[i][j]);
   		}
	}
}

printf("\n ----------------------------------------------");

/***free GPU memory allocation*********************************************************/
cudaFree(d_array_in);
cudaFree(d_in);

/***Back substitution******************************************************************/
/*
h_x[n-1]=h_A1[n-1][n]/h_A1[n-1][n-1];
if(n<10)
	printf("\n Value of X(%d) = %f", (n-1),h_x[(n-1)]);

for(i=n-2;i>=0;i--)
{
	sum = 0;

	for(j=i;j<n;j++)
	{
		sum +=h_A1[i][j]*h_x[j];
	}

	h_x[i] = (h_A1[i][n]-sum)/h_A1[i][i];

	if(n<10)
		printf("\n value of X(%d)=%f",i,h_x[i]);
}
*/

free(h_A);
free(h_A1);
free(h_x);

printf("\n ---------JOB OVER -------------------------");
return 0;
}
