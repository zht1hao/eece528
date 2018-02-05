/*
NORMAL RUN
Run : ./exe_name matrix-size

TEST for validation
Enable the validation witht Setting the test_flag to 1
run : ./exe-name
     without any parameters
*/
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "cuda.h"
#include <string.h>

#define BLOCK_SIZE 512
#define MIN_SIZE 64
#define MAX_SIZE 4096
#define MIN(x, y) ((x) < (y) ? x : y)
#define MAX(x, y) ((x) > (y) ? x : y)
#define DEBUG 0
#define RESULT_COLUMN 1
#define TEST_ELEMENT 3
#define TEST_DATA 0
#define TEST_SIZE 4


int g_iSizeParameter = -1;
float *g_fGuassionInput;
float *g_fGuassionResult;


unsigned int totalKernelTime = 0;

__global__ void performOuterLoop(float *pr_fCudaStr, float *pr_fGuassMatrix, int pr_iOffset, int pr_iSizeParameter)
{
	// getting the coffecient for the respective row block as per the row offset	
        if(threadIdx.x + blockIdx.x * blockDim.x >= pr_iSizeParameter-1-pr_iOffset) return;
        *(pr_fCudaStr+pr_iSizeParameter*(blockDim.x*blockIdx.x+threadIdx.x+pr_iOffset+1)+pr_iOffset) = 
		*(pr_fGuassMatrix+pr_iSizeParameter*(blockDim.x*blockIdx.x+threadIdx.x+pr_iOffset+1)+pr_iOffset) / 
		*(pr_fGuassMatrix+pr_iSizeParameter*pr_iOffset+pr_iOffset);
}


__global__ void performInnerLoop(float *pr_fCudaStr, float *pr_fGuassMatrix, int pr_iOffset, int pr_iSizeParameter)
{       
        if(threadIdx.x + blockIdx.x * blockDim.x >= pr_iSizeParameter-1-pr_iOffset) return;
        if(threadIdx.y + blockIdx.y * blockDim.y >= pr_iSizeParameter-pr_iOffset) return;
        
        int l_rowValue = blockIdx.x * blockDim.x + threadIdx.x;
        int l_colValue = blockIdx.y * blockDim.y + threadIdx.y;
        // applying the above coffecients through out the block as per the row offset 
        pr_fGuassMatrix[ pr_iSizeParameter*(l_rowValue+1+pr_iOffset)+(l_colValue+pr_iOffset)] -= 
		pr_fCudaStr[ pr_iSizeParameter*(l_rowValue+1+pr_iOffset)+pr_iOffset] * 
		pr_fGuassMatrix[ pr_iSizeParameter*pr_iOffset+(l_colValue+pr_iOffset)];
}


void defineFloatMultiDimensionalArrayFunc( float **pr_fMultiDimArray , int pr_iSizeParameter){
        // Defination of multi dimensional array 
        *pr_fMultiDimArray = (float *) calloc(  pr_iSizeParameter ,  sizeof(float));
        if(pr_fMultiDimArray == NULL){
                printf(" Error : Memory not allocated ");
                exit(0);
        }
}//defineFloatMultiDimensionalArray Ends

void checkAnswer( float  **prStrSharedData ){

        FILE *l_ptrFileOutput ;
        l_ptrFileOutput = fopen("result.txt","r+");
        if( l_ptrFileOutput == NULL ){
                printf("\n Error : In opening file ");
        }else{
                float* result_vector = (float*) malloc(sizeof(float)*TEST_ELEMENT);
                float row_sum;

                for (int j=0; j<TEST_ELEMENT; j++){
                        row_sum = 0;
                        for (int k=0; k<=TEST_ELEMENT; k++){
                                row_sum += (*prStrSharedData)[j*TEST_ELEMENT+k];
                        }
                        result_vector[j] = row_sum;
                        printf(" \n result_vector[ %f ]",result_vector[j]);
                }

                float sumOfSquares = 0;
                float entryOfResidual;
                float element = 0 ;
                for (int i=0; i<TEST_ELEMENT; i++){
                        int scn = fscanf(l_ptrFileOutput,"%f", &element) ;
                        entryOfResidual = result_vector[i] - element;
                        sumOfSquares += entryOfResidual*entryOfResidual;
                }
                sumOfSquares = sqrt(sumOfSquares);

                if( (int)sumOfSquares   )
			printf("\n FAIL :: The Bench Mark for L2-Norm  is 0 &  result is %.2f \n",  sumOfSquares);
                else
                        printf("\n SUCESSS :: The L2-Norm of the result vector from Ax-b is: %.2f\n", sumOfSquares);

                free(result_vector);
                fclose(l_ptrFileOutput);
        }
}


///////////////////////////////// DISPLAY FUNCTION ////////////////////////////////////
void displayMatrixFunc(int pr_iSizeParameter  , int **pr_fMultiDimArray){

        int l_iOuterBound , l_iInnerBound , l_iSize = sqrt(pr_iSizeParameter) ;
        printf("--------------------------- \n");
        for ( l_iOuterBound = 0 ; l_iOuterBound < l_iSize ; l_iOuterBound++ ){
                printf("\n");
                for ( l_iInnerBound = 0 ; l_iInnerBound < l_iSize ; l_iInnerBound++ ){
                                printf( "%d ", (int) (*pr_fMultiDimArray)[l_iOuterBound * l_iSize + l_iInnerBound]);
                }//innerBound Ends
        }//outerBound ends
        printf("\n");
}//displayMatrix ends

//////////////////////// INITALIZE THE REQUIRED MATRIX //////////////////////////////////////////
void initializeMatrixFunc(int pr_iSizeParameter , float *pr_fMultiDimArray ){
        int l_iOuterBound   ;
        for ( l_iOuterBound = 0 ; l_iOuterBound < (pr_iSizeParameter) ; l_iOuterBound++ ){
                                (pr_fMultiDimArray)[l_iOuterBound] = ( (int)rand()/(int) (RAND_MAX/10)  );
        }//outerBound ends
}//intializeMatrix ends


void displayMatrix(float **ary, int nrow, int ncol)
{
	int i, j;
	
	for (i=0; i<nrow; i++) {
		for (j=0; j<ncol; j++) {
			printf("%3.5f ", (*ary)[g_iSizeParameter*i+j]);
		}
		printf("\n");
	}
	printf("\n");
}

void checkError(const char *pr_cMssg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err) 
    {
        printf("Cuda error: %s: %s.\n", pr_cMssg,cudaGetErrorString( err) );
        exit(EXIT_FAILURE);
    }                         
}

 
void performUpperTraingleOPerations( int pr_iSizeParameter)
{
	float *l_fCoffCudaPtr,*l_fInputCudaPtr;
	
	// allocate memory on GPU
	cudaMalloc((void **) &l_fCoffCudaPtr, pr_iSizeParameter * pr_iSizeParameter * sizeof(float));
	cudaMalloc((void **) &l_fInputCudaPtr, pr_iSizeParameter * pr_iSizeParameter * sizeof(float));
	
	// copy memory to GPU
	cudaMemcpy(l_fInputCudaPtr, g_fGuassionInput, pr_iSizeParameter * pr_iSizeParameter * sizeof(float),cudaMemcpyHostToDevice );
	cudaMemcpy(l_fCoffCudaPtr, g_fGuassionResult, pr_iSizeParameter * pr_iSizeParameter * sizeof(float),cudaMemcpyHostToDevice );
#if DEBUG
	displayMatrix(&g_fGuassionInput,pr_iSizeParameter,pr_iSizeParameter);	
	displayMatrix(&g_fGuassionResult,pr_iSizeParameter,pr_iSizeParameter);	
#endif	
	//SETTING THE GRID SIZE OF THE OUTER & INNER LOOP THING
	int l_iGridSize = (pr_iSizeParameter/BLOCK_SIZE) + (!(pr_iSizeParameter%BLOCK_SIZE)? 0:1);
	dim3 dimBlock(BLOCK_SIZE);
	dim3 dimGrid(l_iGridSize);
	//dim3 dimGrid( (N/dimBlock.x) + (!(N%dimBlock.x)?0:1) );
	int m_i2DBlockSize, m_i2DGridSize;
	if(  TEST_DATA )
                m_i2DBlockSize = TEST_SIZE;
        else
                m_i2DBlockSize = MIN_SIZE ;
	//using the size parameter & block_size to wrap around the block size
	m_i2DGridSize = (pr_iSizeParameter/m_i2DBlockSize) + (!(pr_iSizeParameter%m_i2DBlockSize?0:1)); 
	
	dim3 dimBlockXY(m_i2DBlockSize,m_i2DBlockSize);
	dim3 dimGridXY(m_i2DGridSize,m_i2DGridSize);
	// begin timing kernels
	int l_iOffSet;
	 cudaEvent_t l_cStart, l_cStop;
        cudaEventCreate(&l_cStart);
        cudaEventCreate(&l_cStop);
        cudaEventRecord(l_cStart, 0);
	for (l_iOffSet=0; l_iOffSet<(pr_iSizeParameter-1); l_iOffSet++) {
		// doing the call to the outer loop for performing divide operarions i.e. getting coffecients 
		performOuterLoop<<<dimGrid,dimBlock>>>(l_fCoffCudaPtr,l_fInputCudaPtr,l_iOffSet,pr_iSizeParameter);
		cudaThreadSynchronize();
		
		// doing the call to the inner loop for performing divide operarions i.e applying coffecient through out the row
		performInnerLoop<<<dimGridXY,dimBlockXY>>>(l_fCoffCudaPtr,l_fInputCudaPtr,l_iOffSet,pr_iSizeParameter);
		cudaThreadSynchronize();
	}
	cudaEventRecord(l_cStop, 0);
	cudaEventSynchronize(l_cStop);
	float l_fTime = 0 ;
        cudaEventElapsedTime(&l_fTime, l_cStart, l_cStop);
        printf("Guassion GPU time:  %f ms \n", l_fTime/1000.0);
	// copy memory back to CPU
	cudaMemcpy(g_fGuassionInput, l_fInputCudaPtr, pr_iSizeParameter * pr_iSizeParameter * sizeof(float),cudaMemcpyDeviceToHost );
	cudaMemcpy(g_fGuassionResult, l_fCoffCudaPtr, pr_iSizeParameter * pr_iSizeParameter * sizeof(float),cudaMemcpyDeviceToHost );
#if DEBUG
	displayMatrix(&g_fGuassionInput,pr_iSizeParameter,pr_iSizeParameter); 
#endif
	if( TEST_DATA )
		checkAnswer(&g_fGuassionInput);
	
	cudaFree(l_fCoffCudaPtr);
	cudaFree(l_fInputCudaPtr);
}

void performGaussElimination( int pr_iSizeParameter)
{

	if(TEST_DATA)
		g_iSizeParameter = 3; 

	defineFloatMultiDimensionalArrayFunc( &g_fGuassionInput , g_iSizeParameter * g_iSizeParameter );
	defineFloatMultiDimensionalArrayFunc( &g_fGuassionResult , g_iSizeParameter * g_iSizeParameter );

	if(TEST_DATA){
		g_fGuassionInput[0] = 10;
		g_fGuassionInput[1] = -7;
		g_fGuassionInput[2] = 3;
		g_fGuassionInput[3] = -6;
		g_fGuassionInput[4] = 8;
		g_fGuassionInput[5] = 4;
		g_fGuassionInput[6] = 2;
		g_fGuassionInput[7] = 6;
		g_fGuassionInput[8] = 9;
	}else{
		initializeMatrixFunc(g_iSizeParameter * g_iSizeParameter, g_fGuassionInput);
	}	


	performUpperTraingleOPerations(pr_iSizeParameter);
	free(g_fGuassionInput);
	free(g_fGuassionResult);

}//performGaussElimination ends

//FUNCTIONALITY START
// main funtionality methods
int checkInputLimits(int pr_iSizeParameter)
{       
        return (MIN_SIZE <= pr_iSizeParameter && pr_iSizeParameter <= MAX_SIZE);
} //checkInputLimits ends

int main(int argc, char *argv[])
{

        if( argc == 2 ){
                g_iSizeParameter = atoi(argv[1]);

                if( checkInputLimits( g_iSizeParameter )){
                        // Defination of multi dimensional array
                        performGaussElimination( g_iSizeParameter );
                        printf("\n sucess ");
                }
                else{

                        printf("\n Matrix size range should be 64 - 4096 ");
                        exit(-1);
                }
        }//argc if ends
        else{
                if (TEST_DATA)
                        performGaussElimination( TEST_ELEMENT );
                else
                {
                        printf("\n Need Matrix size (Two argument expected i.e size & dataType)!");
                        exit(-1);
                }
        }
        return -1;

}
