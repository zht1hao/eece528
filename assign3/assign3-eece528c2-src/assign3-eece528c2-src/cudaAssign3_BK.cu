//Author: 	Navdeep 79776167
//About	:	Guassion Elimination
//Date	:	OCt 2017
/*Decription :
Input : 
Output:
RUn:		./object_file number_of_parameters 
Description:	Implement the CUDA using the 1D array 

 */


//Header files
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>

//global const
#define MIN_SIZE 64
#define MAX_SIZE 4096
#define DEBUG 0
#define TEST_DATA 1
#define TEST_ELEMENT 3
#define RESULT_COLUMN 0
#define MIN(x,y)((x)<(y)? x:y )
#define MAX(x,y)((x)>(y)? x:y )
#define BLOCK_SIZE 3


// Structure declarationi for sahring data among threads

struct STR_SHARED_DATA
{
	unsigned int th_uiColumns;
	unsigned int th_uiRows;
	unsigned int th_iSizeParameter;
	float *th_fGuassionElements;
};
typedef STR_SHARED_DATA StrSharedData;

//global vraibles
double g_dFirstTimeTaken = 0 , g_dSecondTimeTaken = 0;

// main funtionality methods
int checkInputLimits(int pr_iSizeParameter){
	return ( MIN_SIZE <= pr_iSizeParameter && pr_iSizeParameter <= MAX_SIZE) ;
}//checkInputLimits ends


void defineFloatMultiDimensionalArrayFunc( float **pr_fMultiDimArray , int pr_iSizeParameter){
	// Defination of multi dimensional array 
	*pr_fMultiDimArray = (float *) calloc( pr_iSizeParameter ,  sizeof(float));
	if(pr_fMultiDimArray == NULL){
		printf(" Error : Memory not allocated ");
		exit(0);
	}
}//defineFloatMultiDimensionalArray Ends

void testData(  float **pr_vMatrix ){
        (*pr_vMatrix)[0] = 10;
        (*pr_vMatrix)[1] = -7;
        (*pr_vMatrix)[2] = 3;
        (*pr_vMatrix)[3] = -6;
        (*pr_vMatrix)[4] = 8;
        (*pr_vMatrix)[5] = 4;
        (*pr_vMatrix)[6] = 2;
        (*pr_vMatrix)[7] = 6;
        (*pr_vMatrix)[8] = 9;
}
// checking the answwer of the Test data
void checkAnswer(const StrSharedData prStrSharedDataOnGPU ){

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
                                row_sum += prStrSharedDataOnGPU.th_fGuassionElements[j*TEST_ELEMENT+k];
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
                float squares = 0 ;
                int scn = fscanf(l_ptrFileOutput,"%f", &squares) ;
                sumOfSquares = sqrt(sumOfSquares);
                printf("squares === %.20f \n", squares);

                if( squares == sumOfSquares)
                 printf("\n SUCESSS :: The L2-Norm of the result vector from Ax-b is: %.20f\n", sumOfSquares);
                else
                 printf("\n FAIL :: The Bench Mark for L2-Norm  %.20f & result is %.20f \n", squares, sumOfSquares);

                free(result_vector);
                fclose(l_ptrFileOutput);
        }
}


/////////////////// CUDA FUNCTIONS ///////////////////////
/* DATA TO BE COPY BACK & FORTH FROM GPPU TO CPU
   th_uiColumns;
   th_uiRows;
   th_iSizeParameter;
   th_fGuassionElements
 */
StrSharedData defineDataOnGPU(const StrSharedData prStrSharedData){
	StrSharedData l_StrSharedData = prStrSharedData;
	int size = l_StrSharedData.th_uiRows * l_StrSharedData.th_uiColumns * sizeof(float);
	// Alloc space for device copies of l_StrSharedData.th_fGuassionElements array 
	cudaMalloc((void**)&l_StrSharedData.th_fGuassionElements, size);
	cudaMemset(l_StrSharedData.th_fGuassionElements, 0, size );
	return l_StrSharedData;
}// define data on GPU ends

//function to copy matrix to GPU
void sendMatrixToGPU( const StrSharedData prStrSharedDataOnCPU, StrSharedData prStrSharedDataOnGPU)
{
	int size = prStrSharedDataOnCPU.th_uiRows * prStrSharedDataOnCPU.th_uiColumns * sizeof(float);
	prStrSharedDataOnGPU.th_uiRows = prStrSharedDataOnCPU.th_uiRows;
	prStrSharedDataOnGPU.th_uiColumns = prStrSharedDataOnCPU.th_uiColumns;
	prStrSharedDataOnGPU.th_iSizeParameter = prStrSharedDataOnCPU.th_iSizeParameter;
	// Copy inputs to device
	cudaMemcpy(prStrSharedDataOnGPU.th_fGuassionElements, prStrSharedDataOnCPU.th_fGuassionElements, size, cudaMemcpyHostToDevice);
}//sendMatrixToGPU ends

// Copy a device matrix to a host matrix.
void recvMatrixFromGPU(const StrSharedData prStrSharedDataOnGPU, StrSharedData prStrSharedDataOnCPU){
	int size = prStrSharedDataOnGPU.th_uiRows * prStrSharedDataOnGPU.th_uiColumns * sizeof(float);
	// Copy input from the device
	// Copy result back to host
	cudaMemcpy(prStrSharedDataOnCPU.th_fGuassionElements, prStrSharedDataOnGPU.th_fGuassionElements, size, cudaMemcpyDeviceToHost);
	
}//recvMatrixFromGPU ends

//function to check the cuda errors
void checkError(const char *msg){
	cudaError_t err = cudaGetLastError();
	if( cudaSuccess != err){
		printf("CUDA ERROR: %s (%s).\n", msg, cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
}//checkError Ends

//function to work on GPU to caculate the intermerdiate data
/*
blockIdx.x, blockIdx.y, blockIdx.z are built-in variables that returns the block ID
in the x-axis, y-axis, and z-axis of the block that is executing the given block of code.
• threadIdx.x, threadIdx.y, threadIdx.z are built-in variables that return the
thread ID in the x-axis, y-axis, and z-axis of the thread that is being executed by this
stream processor in this particular block.
• blockDim.x, blockDim.y, blockDim.z are built-in variables that return the “block
dimension” (i.e., the number of threads in a block in the x-axis, y-axis, and z-axis).

*/
__global__ void performGuassEliminationOuterLoop(float *pr_fCudaElements, float *pr_fGuassionElements, int pr_iSizeParameter, int pr_iOffest)
{

	//Used blockIdx.x,y to access block index within grid
	//Used threadIdx.x,y to access thread index within block
	//Used blockDim.x,y gives the number of threads in a block
	//USed gridDim.x,y gives the number of blocks in a grid
        if(threadIdx.x + blockIdx.x * blockDim.x >= pr_iSizeParameter-1-pr_iOffest) return;
        *(pr_fCudaElements+pr_iSizeParameter*(blockDim.x*blockIdx.x+threadIdx.x+pr_iOffest+1)+pr_iOffest) = 
	*(pr_fGuassionElements+pr_iSizeParameter*(blockDim.x*blockIdx.x+threadIdx.x+pr_iOffest+1)+pr_iOffest) /
	 *(pr_fGuassionElements+pr_iSizeParameter*pr_iOffest+pr_iOffest);
}

__global__ void performGuassEliminationInnerLoop(float *pr_fCudaElements, float *pr_fGuassionElements,  int pr_iSizeParameter, int pr_iOffest)
{       
	//Used blockIdx.x,y to access block index within grid
	//Used threadIdx.x,y to access thread index within block
	//Used blockDim.x,y gives the number of threads in a block
	//USed gridDim.x,y gives the number of blocks in a grid
        if(threadIdx.x + blockIdx.x * blockDim.x >= pr_iSizeParameter-1-pr_iOffest) return;
        if(threadIdx.y + blockIdx.y * blockDim.y >= pr_iSizeParameter-pr_iOffest) return;
        
        int l_iRowID = blockIdx.x * blockDim.x + threadIdx.x;
        int l_iColID = blockIdx.y * blockDim.y + threadIdx.y;
        
	__syncthreads();
        pr_fGuassionElements[pr_iSizeParameter*(l_iRowID+1+pr_iOffest)+(l_iColID+pr_iOffest)] -= pr_fCudaElements[pr_iSizeParameter*(l_iRowID+1+pr_iOffest)+pr_iOffest] * 
											pr_fGuassionElements[pr_iSizeParameter*pr_iOffest+(l_iColID+pr_iOffest)];
}


///////////////////////////////// DISPLAY FUNCTION ////////////////////////////////////
void displayMatrixFunc(const StrSharedData pr_MatrixData){
	printf("\n -------- print final value -------- \n");
	unsigned int l_iOuterLoop, l_iInnerLoop;
        for(l_iOuterLoop = 0; l_iOuterLoop < pr_MatrixData.th_uiRows; l_iOuterLoop++){
                for(l_iInnerLoop = 0; l_iInnerLoop < pr_MatrixData.th_uiColumns; l_iInnerLoop++)
                        printf("%f ", pr_MatrixData.th_fGuassionElements[l_iOuterLoop*pr_MatrixData.th_uiRows + l_iInnerLoop]);
                printf("\n");
	}//outerBound ends
	printf("\n");
}//displayMatrix ends

//////////////////////// INITALIZE THE REQUIRED MATRIX //////////////////////////////////////////
void initializeMatrixFunc(int pr_iSizeParameter , float **pr_fMultiDimArray ){
	int l_iOuterBound   ;
	for ( l_iOuterBound = 0 ; l_iOuterBound < pr_iSizeParameter ; l_iOuterBound++ ){
		(*pr_fMultiDimArray)[l_iOuterBound] = ( (float)rand()/(float) (RAND_MAX/10)  );
	}//outerBound ends
}//intializeMatrix ends

///////////////////////////// GUASSIAN ELIMINATION /////////////////////////////////////////////

/*
   th_uiColumns;
   th_uiRows;
   th_iSizeParameter;
   th_fGuassionElements
 */
void performGaussElimination(int pr_iSizeParameter)
{

	// Matrices for the program
	StrSharedData  l_iGuassionMatrix; // The NxN input matrix
	// Allocate and initialize the matrices
	l_iGuassionMatrix.th_iSizeParameter =pr_iSizeParameter; 
	l_iGuassionMatrix.th_uiColumns = pr_iSizeParameter;
        l_iGuassionMatrix.th_uiRows = pr_iSizeParameter;   

	// Alloc space for host copies of a, b, c and setup input values
	defineFloatMultiDimensionalArrayFunc(&l_iGuassionMatrix.th_fGuassionElements, l_iGuassionMatrix.th_uiColumns*l_iGuassionMatrix.th_uiRows );
	if( TEST_DATA )
	{
		testData( &l_iGuassionMatrix.th_fGuassionElements);
		//displayMatrixFunc(l_iGuassionMatrix);

	}
	else
	{
		initializeMatrixFunc(pr_iSizeParameter*pr_iSizeParameter, &l_iGuassionMatrix.th_fGuassionElements);
	}

	/// PErform the pivoting operations ///
#if 0 
	int i,k,j,p;
	float app, *temp;
	temp = (float*)malloc(l_iGuassionMatrix.th_iSizeParameter*sizeof(float));;
 	
	for(i=0;i<(l_iGuassionMatrix.th_iSizeParameter);i++){
                app = l_iGuassionMatrix.th_fGuassionElements[i* l_iGuassionMatrix.th_uiColumns+i];
                //initialization of p
                p = i;
                printf("\n  l_iMax = %f , l_iFlag = %d \n ",app,p);
                //find largest no of the columns and row no. of largest no.
                for(k = i+1; k < l_iGuassionMatrix.th_iSizeParameter ; k++)
                        if(fabs(app) < fabs(l_iGuassionMatrix.th_fGuassionElements[k+ i *l_iGuassionMatrix.th_uiColumns ])){
                                app = l_iGuassionMatrix.th_fGuassionElements[k+ i *l_iGuassionMatrix.th_uiColumns] ;
                                p = k;
                                printf("\n ==>  l_iMax = %f , l_iFlag = %d \n ",app,p);
                        }
                //swaping the elements of diagonal row and row containing largest no
                for(j = 0; j < l_iGuassionMatrix.th_iSizeParameter; j++)
                {
                        temp[j] = l_iGuassionMatrix.th_fGuassionElements[p*l_iGuassionMatrix.th_uiColumns +j ];
                        l_iGuassionMatrix.th_fGuassionElements[j+l_iGuassionMatrix.th_uiColumns*p] = l_iGuassionMatrix.th_fGuassionElements[j*l_iGuassionMatrix.th_uiColumns+i];
                        l_iGuassionMatrix.th_fGuassionElements[j*l_iGuassionMatrix.th_uiColumns+i] = temp[j];
                }


        }
#endif
	displayMatrixFunc(l_iGuassionMatrix);
	//////////////////////////////////////
	
	//Define the matrix data & send it to the GPU  
	StrSharedData l_strGuassionMatrixGPU = defineDataOnGPU(l_iGuassionMatrix);
	sendMatrixToGPU( l_iGuassionMatrix , l_strGuassionMatrixGPU );

	dim3 gpuGridOuter(BLOCK_SIZE);
        dim3 gpuGridBlock(BLOCK_SIZE);


	// Setting the thread block to the BLOCK_SIZE
	dim3 gpuThreadBlock(BLOCK_SIZE, BLOCK_SIZE);
	//Setting the thread the GPU grid accordingly
	dim3 gpuGrid(1, ceil((float)l_iGuassionMatrix.th_uiRows / BLOCK_SIZE));

	printf("Performing gaussian elimination on the GPU\n");

	//Intra-block sync is implemented with __syncthreads()	
	int l_iLoop;
	float l_fTime; 
	cudaDeviceSynchronize();
	//Creating event to find the elapsed time in the GPU while performing the event
	cudaEvent_t l_cStart, l_cStop;
	cudaEventCreate(&l_cStart);
	cudaEventCreate(&l_cStop);
	cudaEventRecord(l_cStart, 0);

//	printf("\n %d", l_strGuassionMatrixGPU.th_fGuassionElements);
	for (l_iLoop = 0; l_iLoop < pr_iSizeParameter - 1; l_iLoop++ )
	{
		// reduce the row elemintaion by splitting the basic logic into two parts.
		// doing the outer loop of guassion to compute the constant the Get a 1 in the row & column. we need to multiply by a constant in the second call i.e inner elimination loop 
		performGuassEliminationOuterLoop <<<gpuGridOuter,gpuGridOuter>>> (l_strGuassionMatrixGPU.th_fGuassionElements,l_iGuassionMatrix.th_fGuassionElements,l_strGuassionMatrixGPU.th_uiColumns,l_iLoop );
		cudaDeviceSynchronize();
		
		performGuassEliminationInnerLoop <<<gpuGrid,gpuThreadBlock>>> (l_strGuassionMatrixGPU.th_fGuassionElements,l_iGuassionMatrix.th_fGuassionElements,l_strGuassionMatrixGPU.th_uiColumns-l_iLoop,l_iLoop );
		cudaDeviceSynchronize();
	}
	cudaEventRecord(l_cStop, 0);
	cudaEventSynchronize(l_cStop);
	cudaEventElapsedTime(&l_fTime, l_cStart, l_cStop);
	printf("Guassion GPU time:  %f ms \n", l_fTime/1000);

	//reset the guassion matrix 
#if DEBUG
	displayMatrixFunc(l_iGuassionMatrix);
#endif	
	memset(l_iGuassionMatrix.th_fGuassionElements, 0, (l_iGuassionMatrix.th_uiRows * l_iGuassionMatrix.th_uiColumns * sizeof(float)));
#if DEBUG
	displayMatrixFunc(l_iGuassionMatrix);
#endif	
	recvMatrixFromGPU( l_strGuassionMatrixGPU, l_iGuassionMatrix  );
#if DEBUG
#endif
	displayMatrixFunc(l_iGuassionMatrix);
	if( TEST_DATA )
	 checkAnswer(l_iGuassionMatrix);
		
	cudaFree(l_strGuassionMatrixGPU.th_fGuassionElements);
	free(l_iGuassionMatrix.th_fGuassionElements);
	l_iGuassionMatrix.th_fGuassionElements = NULL;
}//peformGuassion ends 

///////////////////////// CONTROLLER FUNCTIONALITY ////////////////////////////////////////////

int main(int argc, char *argv[]){

	int l_iSizeParameter = -1 ;
	if( argc == 2 ){
		l_iSizeParameter = atoi(argv[1]);

		if( checkInputLimits( l_iSizeParameter )){
			// Defination of multi dimensional array 
			performGaussElimination( l_iSizeParameter );
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
}//main ends

