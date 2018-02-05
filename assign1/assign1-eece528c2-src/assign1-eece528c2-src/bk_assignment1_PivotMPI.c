
/*
Arguments ==    1st argument is matrix size
 RUN :  mpirun -n number_of_process objectfile matrix_size
FOR checking answer = ENable the TEST_DATA flag 
                This will enable the check answer functionality with respect to the becnh mark result using pivot feeded from result.txt 
*/


#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "mpi.h"

#define MASTER 0               /*   first task */
#define FROM_MASTER 1          /* setting a message type */
#define FROM_WORKER 2          /* setting a message type */
#define MIN_SIZE 64
#define MAX_SIZE 4096
#define DEBUG 0
#define TEST_DATA 1
#define RESULT_COLUMN 1
#define MIN(x,y)((x)<(y)? x:y )
#define MAX(x,y)((x)>(y)? x:y )
#define TEST_ELEMENT 3

double g_dStartTime, g_dEndTime;


//function to free memory
int free2dfloat(float ***array) {
	free(&((*array)[0][0]));
	free(*array);
	return 0;
}//free2dfloat ends

//function to allocate the array 
int define2DFloatMatrix(float ***array, int pr_iRowParameter, int pr_iColParameter) {
	float *l_ptrArray = (float *)malloc(pr_iRowParameter*pr_iColParameter*sizeof(float));
	if (!l_ptrArray) return -1;

	(*array) = (float **)malloc((pr_iRowParameter)*sizeof(float*));
	if (!(*array)) {
		free(l_ptrArray);
		return -1;
	}
	int i ;
	for ( i=0; i<(pr_iRowParameter); i++) 
		(*array)[i] = &(l_ptrArray[i*pr_iColParameter]);

	int l_iOuterLoop,l_iInnerLoop;
	 for (l_iOuterLoop=0; l_iOuterLoop<pr_iRowParameter; l_iOuterLoop++)
                {       
                        for (l_iInnerLoop=0; l_iInnerLoop<pr_iColParameter; l_iInnerLoop++)
                        {               (*array)[l_iOuterLoop][l_iInnerLoop] = ( (float)rand()/(float) (RAND_MAX/10)  );
                        }
                }

	return 0;
}//define2DFloatMatrix ends

int defineTestData(float ***array, int pr_iRowParameter, int pr_iColParameter) {
        float *l_ptrArray = (float *)malloc(pr_iRowParameter*pr_iColParameter*sizeof(float));
        if (!l_ptrArray) return -1;

        (*array) = (float **)malloc((pr_iRowParameter)*sizeof(float*));
        if (!(*array)) {
                free(l_ptrArray);
                return -1;
        }
        int i ;
        for ( i=0; i<(pr_iRowParameter); i++)
                (*array)[i] = &(l_ptrArray[i*pr_iColParameter]);

       
	 (*array)[0][0] = 10;
        (*array)[0][1] = -7;
        (*array)[0][2] = 3;
        (*array)[0][3] = 5;
        (*array)[1][0] = -6;
        (*array)[1][1] = 8;
        (*array)[1][2] = 4;
        (*array)[1][3] = 7;
        (*array)[2][0] = 2;
        (*array)[2][1] = 6;
        (*array)[2][2] = 9;
        (*array)[2][3] = -1;

	 return 0;
	
		

}//define2DFloatMatrix ends


//performBackwardSubstitution for guassion elemination
void performBackwardSubstitution(int pr_iSizeParameter , float ***pr_vMatrix , float **pr_vResultMatrix ){

	int l_iOuterLoop, l_iInnerLoop ;
	float l_fCompute;
	// CODE TO CHECKED FOR THE RUNTIME MATRIX PART
	(*pr_vResultMatrix)[pr_iSizeParameter-1]= (*pr_vMatrix) [pr_iSizeParameter-1][pr_iSizeParameter]/(*pr_vMatrix)[pr_iSizeParameter-1][pr_iSizeParameter-1];

	for( l_iOuterLoop = pr_iSizeParameter-2 ; l_iOuterLoop >= 0 ; l_iOuterLoop-- ){
		l_fCompute = 0 ;
		for( l_iInnerLoop = l_iOuterLoop+1 ; l_iInnerLoop < pr_iSizeParameter ; l_iInnerLoop++){
			l_fCompute += (*pr_vMatrix)[l_iOuterLoop][l_iInnerLoop] * (*pr_vResultMatrix)[l_iInnerLoop];
		}//Inner loop ends
		(*pr_vResultMatrix)[l_iOuterLoop] = ((*pr_vMatrix)[l_iOuterLoop][pr_iSizeParameter]- l_fCompute)/ (*pr_vMatrix)[l_iOuterLoop][l_iOuterLoop];
	}//Outer loop ends 
}// performBackwardSubs ends


//Checking the inputs for inner & outer bound
int checkInputLimits(int pr_iSizeParameter){
	return ( MIN_SIZE <= pr_iSizeParameter && pr_iSizeParameter <= MAX_SIZE) ;
}//checkInputLimits ends

//function to display matrix
void display1DimMatrixFunc(int pr_iSizeParameter , float **pr_fMultiDimArray){
	printf("\n -------- print final value -------- \n");
	int l_iOuterBound  ;
	for ( l_iOuterBound = 0 ; l_iOuterBound < pr_iSizeParameter ; l_iOuterBound++ ){
		printf( "%.3f ", (float) (*pr_fMultiDimArray)[l_iOuterBound]);
	}//outerBound ends
}//displayMatrix ends

//function to lof the time process
void logTimePeocessData(int pr_iSizeParameter, int pr_iProcessID){

	FILE *l_ptrFileOutput ;
	l_ptrFileOutput = fopen("output.txt","a+");
	if( l_ptrFileOutput == NULL ){
		printf("\n Error : In opening file ");
	}else{
		fprintf(l_ptrFileOutput,"%d ", pr_iSizeParameter);
		fprintf(l_ptrFileOutput,"%d ", pr_iProcessID);
		fprintf(l_ptrFileOutput,"%f", (double)(g_dEndTime - g_dStartTime));
		fprintf(l_ptrFileOutput,"\n");
		printf("\n Time taken to process guassion on matrix: %f seconds\n", (g_dEndTime - g_dStartTime));
		fclose(l_ptrFileOutput);
	}
}//logTimePeocessData ends

//////////////////////// INITALIZE THE REQUIRED MATRIX //////////////////////////////////////////
void initializeMatrixFunc(int pr_iSizeParameter , float ***pr_fMultiDimArray ){
	int l_iOuterBound , l_iInnerBound ;
	for ( l_iOuterBound = 0 ; l_iOuterBound < pr_iSizeParameter ; l_iOuterBound++ ){
		for ( l_iInnerBound = 0 ; l_iInnerBound < pr_iSizeParameter +RESULT_COLUMN; l_iInnerBound++ ){
			(*pr_fMultiDimArray)[l_iOuterBound][l_iInnerBound] = ( (float)rand()/(float) (RAND_MAX/10)  );
		}//innerBound Ends
		//INITILZING THE 
	}//outerBound ends
}//intializeMatrix ends


//perform the upper traingular operations
void pperformUpperTraingular( int pr_iRowParameter , int pr_iColParameter  , float ***pr_vMatrix  ){
	int i , k , j ,  l_iFlag ;
	float l_iMax , l_fCompute;

	for( i = 0 ; i < pr_iRowParameter-1 ; i++){
		for( j = i+1  ; j < pr_iRowParameter ; j++){
			l_fCompute = (*pr_vMatrix)[j][i] / (*pr_vMatrix)[i][i];

			for( k = 0 ; k < pr_iColParameter ; k++){
				(*pr_vMatrix)[j][k] -=  l_fCompute * (*pr_vMatrix)[i][k];
#if DEBUG
				printf("\n[ %.3f ] %.3f ",l_fCompute, (*pr_vMatrix)[j][k]);
#endif
			}//Inner loop ends
		}//Intermediate Loop ends
	}//Outer loop ends 
}//pperformUpperTraingular ends

//perform the pivot transforamtions 
void pivotPerformUpperTraingular( int pr_iRowParameter , int pr_iColParameter  , float ***pr_vMatrix  ){
	int i , k , j ,  l_iFlag ;
	float l_iMax , l_fCompute;
	for( i = 0 ; i < pr_iRowParameter ; i++){
		l_iMax = (*pr_vMatrix)[i][i];
		l_iFlag = i;
#if DEBUG	
		printf("\n l_iMax = %f , l_iFlag = %d \n ", l_iMax,l_iFlag);
#endif	
		for( k = i+1  ; k < pr_iRowParameter ; k++){
			if( fabs(l_iMax) < (fabs((*pr_vMatrix)[k][i]))){
                                l_iMax = (*pr_vMatrix)[k][i];
                                l_iFlag = k;
#if DEBUG	
			printf("\n ==> l_iMax = %f , l_iFlag = %d \n ", l_iMax,l_iFlag);
#endif
                        }//if ends
                }//for loop ends

		for( j = 0 ; j <= pr_iRowParameter ; j++){
			l_fCompute = (*pr_vMatrix)[l_iFlag][j];
			(*pr_vMatrix)[l_iFlag][j] = (*pr_vMatrix)[i][j];
			(*pr_vMatrix)[i][j] = l_fCompute ;
		}//for ends
        }//Outer loop ends

#if DEBUG	
        printf("\n pperform pivot"); 
#endif	
}//pivotPerformUpperTraingular ends



void checkAnswer(float ***pr_fMatrix ){
        
        FILE *l_ptrFileOutput ;
        l_ptrFileOutput = fopen("result.txt","r+");
        if( l_ptrFileOutput == NULL ){
                printf("\n Error : In opening file ");
        }else{  
                float* result_vector = (float*) malloc(sizeof(float)*TEST_ELEMENT);
                float row_sum;
               printf("\n"); 
                for (int j=0; j<TEST_ELEMENT; j++){
                        row_sum = 0;
                        for (int k=0; k<=TEST_ELEMENT; k++){
                                row_sum += (*pr_fMatrix)[j][k];
#if DEBUG
                	printf("\t%f",(*pr_fMatrix)[j][k]); 
#endif	
		        }
                        result_vector[j] = row_sum;
#if DEBUG
                        printf(" \n result_vector[ %f ] \n",result_vector[j]);
#endif	
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
#if DEBUG
                printf("squares === %.20f \n", squares);
#endif                
                if( squares == sumOfSquares)
                 printf("\n SUCESSS :: The L2-Norm of the result vector is: %.20f\n", sumOfSquares);
                else
                 printf("\n FAIL :: The Bench Mark for L2-Norm  %.20f & result is %.20f \n", squares, sumOfSquares);
                
                free(result_vector);
                fclose(l_ptrFileOutput);
        }
}


// MAIN functions 
int main (int argc, char *argv[]){

	int	l_iNoOfTask, l_iTaskId, l_iNoOfProcess, l_iSourceId, l_iAbortSattus, l_iDestId, l_iMessageSource;
	int	l_iNoOfRows, l_iRowOffset, l_iSizeParameter ; 
	int 	l_iOuterLoop, l_iInnerLoop ;
	float 	*l_fGuassionFinal = NULL; //SOLUTION
	float	**l_fMatrixPtr;           /* matrix A to be multiplied */

	l_iSizeParameter = -1;

	if(TEST_DATA){
		l_iSizeParameter = TEST_ELEMENT;
		printf("\n TESTE DATA ");
	}
	else
	{
		if( argc == 2 ){
			l_iSizeParameter = atoi(argv[1]);
			if( checkInputLimits( l_iSizeParameter )){
#if DEBUG
				printf("\n Continue the Process .... ");
#endif	
			}
			else{
				printf("\n Matrix size range should be 64 - 4096 ");
				exit(-1);
			}
		}//argc if ends
		else{
			printf("\n Need Matrix size (one argument expected i.e size )!");
			exit(-1);
		}
	}

	if( TEST_DATA)
	{
		printf("\n Initlize the Test data");
		l_iSizeParameter = TEST_ELEMENT;
		defineTestData(&l_fMatrixPtr,l_iSizeParameter,l_iSizeParameter+RESULT_COLUMN);
	}
	else
	{
		define2DFloatMatrix(&l_fMatrixPtr,l_iSizeParameter,l_iSizeParameter+RESULT_COLUMN);
	}


	//Define the SOLUTION SIngle dimensional matrix
	l_fGuassionFinal = (float *) calloc( l_iSizeParameter ,  sizeof(float));
	if(l_fGuassionFinal == NULL){
		printf(" Error : Memory not allocated ");
		exit(0);
	}

	MPI_Status status;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&l_iTaskId);
	MPI_Comm_size(MPI_COMM_WORLD,&l_iNoOfTask);

	if (l_iNoOfTask < 2 ) {
		printf("Need at least two MPI tasks. Quitting...\n");
		MPI_Abort(MPI_COMM_WORLD, l_iAbortSattus);
		exit(1);
	}
	l_iNoOfProcess = l_iNoOfTask-1;
//	printf("l_iNoOfProcess - %d",l_iNoOfProcess);

	/**************************** master task ************************************/
	if (l_iTaskId == MASTER)
	{
#if DEBUG
		printf("mpi started for %d tasks.\n",l_iNoOfTask);
		printf("Initializing arrays...\n");
#endif	
		//perform pivoting 
		pivotPerformUpperTraingular(l_iSizeParameter,l_iSizeParameter+RESULT_COLUMN , &l_fMatrixPtr  );
		/* Send matrix data to the worker tasks */
		int l_iTempRow, l_iTempCompute; 
		l_iTempRow = l_iSizeParameter/l_iNoOfProcess;
		l_iTempCompute = l_iSizeParameter%l_iNoOfProcess;
		l_iRowOffset = 0;
		l_iMessageSource = FROM_MASTER;
		g_dStartTime = MPI_Wtime();
		for (l_iDestId=1; l_iDestId<=l_iNoOfProcess; l_iDestId++)
		{
			l_iNoOfRows = (l_iDestId <= l_iTempCompute) ? l_iTempRow+1 : l_iTempRow;   	
#if DEBUG		
			printf("\n Sending %d l_iNoOfRows to task %d l_iRowOffset=%d\n",l_iNoOfRows,l_iDestId,l_iRowOffset);
#endif	
			MPI_Send(&l_iRowOffset, 1, MPI_INT, l_iDestId, l_iMessageSource, MPI_COMM_WORLD);
			MPI_Send(&l_iNoOfRows, 1, MPI_INT, l_iDestId, l_iMessageSource, MPI_COMM_WORLD);
			MPI_Send(&(l_fMatrixPtr[l_iRowOffset][0]), l_iNoOfRows*(l_iSizeParameter+RESULT_COLUMN), MPI_FLOAT, l_iDestId, l_iMessageSource,
					MPI_COMM_WORLD);
			l_iRowOffset = l_iRowOffset + l_iNoOfRows;
		}

		/* Receive results from worker tasks */
		l_iMessageSource = FROM_WORKER;
		for (l_iOuterLoop=1; l_iOuterLoop<=l_iNoOfProcess; l_iOuterLoop++)
		{
			l_iSourceId = l_iOuterLoop;
			MPI_Recv(&l_iRowOffset, 1, MPI_INT, l_iSourceId, l_iMessageSource, MPI_COMM_WORLD, &status);
			MPI_Recv(&l_iNoOfRows, 1, MPI_INT, l_iSourceId, l_iMessageSource, MPI_COMM_WORLD, &status);
			MPI_Recv(&(l_fMatrixPtr[l_iRowOffset][0]), l_iNoOfRows*(l_iSizeParameter+RESULT_COLUMN), MPI_FLOAT, l_iSourceId, l_iMessageSource, 
					MPI_COMM_WORLD, &status);
#if DEBUG			
			printf("Received results from task %d\n",l_iSourceId);
#endif	
		}
		g_dEndTime = MPI_Wtime();
		performBackwardSubstitution(l_iNoOfRows,&l_fMatrixPtr,&l_fGuassionFinal);
		logTimePeocessData( l_iSizeParameter,  l_iNoOfProcess) ;
		if( TEST_DATA)
		{
			checkAnswer(&l_fMatrixPtr);
		}
		//display1DimMatrixFunc(l_iNoOfRows, &l_fGuassionFinal);
	}//MASTER IF ends


	if (l_iTaskId > MASTER)
	{
		l_iMessageSource = FROM_MASTER;
		MPI_Recv(&l_iRowOffset, 1, MPI_INT, MASTER, l_iMessageSource, MPI_COMM_WORLD, &status);
		MPI_Recv(&l_iNoOfRows, 1, MPI_INT, MASTER, l_iMessageSource, MPI_COMM_WORLD, &status);
		MPI_Recv(&l_fMatrixPtr[0][0], l_iNoOfRows*(l_iSizeParameter+RESULT_COLUMN), MPI_FLOAT, MASTER, l_iMessageSource, MPI_COMM_WORLD, &status);
#if DEBUG
		printf("\n l_iNoOfRows -%d Col -%d ",l_iNoOfRows,l_iSizeParameter+RESULT_COLUMN);
#endif	
		pperformUpperTraingular(l_iNoOfRows,l_iSizeParameter+RESULT_COLUMN , &l_fMatrixPtr  );

	
		l_iMessageSource = FROM_WORKER;
		MPI_Send(&l_iRowOffset, 1, MPI_INT, MASTER, l_iMessageSource, MPI_COMM_WORLD);
		MPI_Send(&l_iNoOfRows, 1, MPI_INT, MASTER, l_iMessageSource, MPI_COMM_WORLD);
		MPI_Send(&l_fMatrixPtr[0][0], l_iNoOfRows*(l_iSizeParameter+RESULT_COLUMN), MPI_FLOAT, MASTER, l_iMessageSource, MPI_COMM_WORLD);
	}//SOURCE if ends

	MPI_Finalize();
	free2dfloat(&l_fMatrixPtr);
	return 0;	
}


