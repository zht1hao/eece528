//Author: 	Navdeep 79776167
//About	:	Vanilla matrix multiply
//Date	:	Sep 28 2017
/*Decription :
Input : 
Output: 
 */


//Header files
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


//global const
#define MIN_SIZE 64
#define MAX_SIZE 4096
#define DEBUG 0
#define TEST_DATA 0
#define RESULT_COLUMN 1
#define MIN(x,y)((x)<(y)? x:y )
#define MAX(x,y)((x)>(y)? x:y )
#define TEST_ELEMENT 3



//global vraibles
float *g_fGuassionFinal = NULL; //SOLUTION
double g_dFirstTimeTaken = 0 , g_dSecondTimeTaken = 0;

// main funtionality methods
int checkInputLimits(int pr_iSizeParameter){
	return ( MIN_SIZE <= pr_iSizeParameter && pr_iSizeParameter <= MAX_SIZE) ;
}//checkInputLimits ends

void defineFloatMultiDimensionalArrayFunc( float ***pr_fMultiDimArray , int pr_iSizeParameter){
	int l_iBound;
	// Defination of multi dimensional array 
	*pr_fMultiDimArray = (float **) calloc( pr_iSizeParameter ,  sizeof(float*));
	if(pr_fMultiDimArray == NULL){
		printf(" Error : Memory not allocated ");
		exit(0);	
	} 

	for ( l_iBound = 0 ; l_iBound < pr_iSizeParameter ; l_iBound++ )
	{	
		(*pr_fMultiDimArray)[l_iBound] = (float *) calloc( pr_iSizeParameter+RESULT_COLUMN ,  sizeof( float));
	}//for loop ends
}//defineFloatMultiDimensionalArray Ends

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

void setAnswer(float **A, float* b){

        FILE *l_ptrFileOutput ;
        l_ptrFileOutput = fopen("result.txt","w+");
        if( l_ptrFileOutput == NULL ){
                printf("\n Error : In opening file ");
        }else{
                float* result_vector = (float*) malloc(sizeof(float)*TEST_ELEMENT);
                float row_sum;

                for (int j=0; j<TEST_ELEMENT; j++){
                        row_sum = 0;
                        for (int k=0; k<=TEST_ELEMENT; k++){
                                row_sum += A[j][k];
                        }
                        result_vector[j] = row_sum;
                        printf(" \n result_vector[ %f ]",result_vector[j]);
                }

                float sumOfSquares = 0;
                float entryOfResidual;
                for (int i=0; i<TEST_ELEMENT; i++){
                        entryOfResidual = result_vector[i] - b[i];
                        sumOfSquares += entryOfResidual*entryOfResidual;
                        fprintf(l_ptrFileOutput,"%.20f\n", b[i]);
                }
                sumOfSquares = sqrt(sumOfSquares);
                printf("\nThe L2-Norm of the result vector from Ax-b is: %.20f\n", sumOfSquares);
                fprintf(l_ptrFileOutput,"%.20f\n",sumOfSquares);

                free(result_vector);
                fclose(l_ptrFileOutput);
        }
}

void checkAnswer(float** A ){

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
                                row_sum += A[j][k];
                        }
                        result_vector[j] = row_sum;
                        printf(" \n result_vector[ %f ]",result_vector[j]);
                }

                float sumOfSquares = 0;
                float entryOfResidual;
                float element = 0 ;
                for (int i=0; i<TEST_ELEMENT; i++){
                        fscanf(l_ptrFileOutput,"%f", &element) ;
                        entryOfResidual = result_vector[i] - element;
                        sumOfSquares += entryOfResidual*entryOfResidual;
                }
                float squares = 0 ;
                fscanf(l_ptrFileOutput,"%f", &squares) ;
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


///////////////////////////////// DISPLAY FUNCTION ////////////////////////////////////
void displayMultiDimMatrixFunc(int pr_iSizeParameter , float ***pr_fMultiDimArray){

	printf("\n -------- print matrix value -------- \n");
	int l_iOuterBound , l_iInnerBound ;
	for ( l_iOuterBound = 0 ; l_iOuterBound < pr_iSizeParameter ; l_iOuterBound++ ){
		printf("\n");
		for ( l_iInnerBound = 0 ; l_iInnerBound < pr_iSizeParameter + RESULT_COLUMN ; l_iInnerBound++ ){
			printf( "%.3f ", (float) (*pr_fMultiDimArray)[l_iOuterBound][l_iInnerBound]);
		}//innerBound Ends
	}//outerBound ends
}//displayMatrix ends

void display1DimMatrixFunc(int pr_iSizeParameter , float **pr_fMultiDimArray){
	printf("\n -------- print final value -------- \n");
	int l_iOuterBound  ;
	for ( l_iOuterBound = 0 ; l_iOuterBound < pr_iSizeParameter ; l_iOuterBound++ ){
		printf( "%.3f ", (float) (*pr_fMultiDimArray)[l_iOuterBound]);
	}//outerBound ends
}//displayMatrix ends

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

///////////////////////////// GUASSIAN ELIMINATION /////////////////////////////////////////////
void performUpperTraingular(int pr_iSizeParameter , float ***pr_vMatrix  ){
	int i , j , k ;
        float l_fCompute;

	clock_t l_cBegin = clock();
        for( i = 0 ; i < pr_iSizeParameter-1 ; i++){
                for( j = i+1  ; j < pr_iSizeParameter ; j++){
                        l_fCompute = (*pr_vMatrix)[j][i] / (*pr_vMatrix)[i][i];

                        for( k = 0 ; k < pr_iSizeParameter ; k++){
                                (*pr_vMatrix)[j][k] -=  l_fCompute * (*pr_vMatrix)[i][k];
#if DEBUG
                                printf("\n[ %.3f ] %.3f ",l_fCompute, (*pr_vMatrix)[j][k]);
#endif
                        }//Inner loop ends
                }//Intermediate Loop ends
        }//Outer loop ends 



	clock_t l_cEnd = clock();
	g_dFirstTimeTaken = ((double)(l_cEnd-l_cBegin)/ CLOCKS_PER_SEC);
#if DEBUG

	displayMultiDimMatrixFunc(pr_iSizeParameter,  pr_vMatrix);
#endif	
}

void performPivotUpperTraingular(int pr_iSizeParameter , float ***pr_vMatrix ){
	int i , k , j ,  l_iFlag ;
	float l_iMax , l_fCompute;
	for( i = 0 ; i < pr_iSizeParameter ; i++){
		l_iMax = (*pr_vMatrix)[i][i];
		l_iFlag = i;
#if DEBUG       
		printf("\n l_iMax = %f , l_iFlag = %d \n ", l_iMax,l_iFlag);
#endif
		for( k = i+1  ; k < pr_iSizeParameter ; k++){
			if( fabs(l_iMax) < (fabs((*pr_vMatrix)[k][i]))){
				l_iMax = (*pr_vMatrix)[k][i];
				l_iFlag = k;
#if DEBUG       
				printf("\n ==> l_iMax = %f , l_iFlag = %d \n ", l_iMax,l_iFlag);
#endif
			}//if ends
		}//for loop ends

		for( j = 0 ; j <= pr_iSizeParameter ; j++){
			l_fCompute = (*pr_vMatrix)[l_iFlag][j];
			(*pr_vMatrix)[l_iFlag][j] = (*pr_vMatrix)[i][j];
			(*pr_vMatrix)[i][j] = l_fCompute ;
		}//for ends
	}//Outer loop ends

#if DEBUG       
	printf("\n pperform pivot");
#endif

}

void performBackwardSubstitution(int pr_iSizeParameter , float ***pr_vMatrix , float **pr_vResultMatrix  ){

	int l_iOuterLoop, l_iInnerLoop ;
	float l_fCompute;
	// CODE TO CHECKED FOR THE RUNTIME MATRIX PART
	clock_t l_cBegin = clock();
	(*pr_vResultMatrix)[pr_iSizeParameter-1]= (*pr_vMatrix) [pr_iSizeParameter-1][pr_iSizeParameter]/(*pr_vMatrix)[pr_iSizeParameter-1][pr_iSizeParameter-1];

	for( l_iOuterLoop = pr_iSizeParameter-2 ; l_iOuterLoop >= 0 ; l_iOuterLoop-- ){
		l_fCompute = 0 ;
		for( l_iInnerLoop = l_iOuterLoop+1 ; l_iInnerLoop < pr_iSizeParameter ; l_iInnerLoop++){
			l_fCompute += (*pr_vMatrix)[l_iOuterLoop][l_iInnerLoop] * (*pr_vResultMatrix)[l_iInnerLoop];
		}//Inner loop ends
		(*pr_vResultMatrix)[l_iOuterLoop] = ((*pr_vMatrix)[l_iOuterLoop][pr_iSizeParameter]- l_fCompute)/ (*pr_vMatrix)[l_iOuterLoop][l_iOuterLoop];
	}//Outer loop ends 
	clock_t l_cEnd = clock();
	g_dSecondTimeTaken = ((double)(l_cEnd-l_cBegin)/ CLOCKS_PER_SEC);
#if DEBUG
	printf("\n Runtime For performBackwardSubstitution : %.3f sec \n ",g_dSecondTimeTaken);
#endif	

}// performBackwardSubs ends


void performGuassionElimination(int pr_iSizeParameter , float ***pr_vFirstMatrix ){

	performPivotUpperTraingular(pr_iSizeParameter , pr_vFirstMatrix );
	performUpperTraingular(pr_iSizeParameter , pr_vFirstMatrix  );
	performBackwardSubstitution(pr_iSizeParameter , pr_vFirstMatrix, &g_fGuassionFinal );
	printf("\n Runtime For Forward Guassion Elimination : %.3f sec \n ",(g_dFirstTimeTaken  ));
	printf("\n Runtime For Subs Guassion Elimination : %.3f sec \n ",(g_dSecondTimeTaken ));
}// performGuassionElimination ends




//////////////////////// CUSTOMIZE FREE FUNCTION ////////////////////////////////////////////////
void trashMemory( int pr_iSizeParameter, void ***pr_vMultiDimFreeMemory){
	int l_iBound ;	
	for ( l_iBound = 0 ; l_iBound < pr_iSizeParameter ; l_iBound++ ){
		free((*pr_vMultiDimFreeMemory)[l_iBound]);		
	}
	free((*pr_vMultiDimFreeMemory));		
}//trashMemory Ends


void testData(int pr_iSizeParameter, float ***pr_vMatrix ){

	(*pr_vMatrix)[0][0] = 10;
	(*pr_vMatrix)[0][1] = -7;
	(*pr_vMatrix)[0][2] = 3;
	(*pr_vMatrix)[0][3] = 5;
	(*pr_vMatrix)[1][0] = -6;
	(*pr_vMatrix)[1][1] = 8;
	(*pr_vMatrix)[1][2] = 4;
	(*pr_vMatrix)[1][3] = 7;
	(*pr_vMatrix)[2][0] = 2;
	(*pr_vMatrix)[2][1] = 6;
	(*pr_vMatrix)[2][2] = 9;
	(*pr_vMatrix)[2][3] = -1;
}

///////////////////////// CONTROLLER FUNCTIONALITY ////////////////////////////////////////////
void performDataTypeOperation(int pr_iSizeParameter ){

	//Define the SOLUTION SIngle dimensional matrix
	g_fGuassionFinal = (float *) calloc( pr_iSizeParameter ,  sizeof(float));
	if(g_fGuassionFinal == NULL){
		printf(" Error : Memory not allocated ");
		exit(0);	
	} 

#if TEST_DATA
	pr_iSizeParameter = 3;
#endif

#if DEBUG	
	display1DimMatrixFunc(pr_iSizeParameter, &g_fGuassionFinal);
#endif	
	//Initialize the other operands for multiplication	
	float **l_fFirstMultiDimArray = NULL;
	// Frist 2d Array
	defineFloatMultiDimensionalArrayFunc( &l_fFirstMultiDimArray ,pr_iSizeParameter );
	//intializing the array with random values
	initializeMatrixFunc(pr_iSizeParameter,  &l_fFirstMultiDimArray);
#if TEST_DATA
	testData(pr_iSizeParameter,  &l_fFirstMultiDimArray);
#endif

#if DEBUG
	displayMultiDimMatrixFunc(pr_iSizeParameter,  &l_fFirstMultiDimArray);
#endif

	//perform GUASSION ELIMINATIONS 
	performGuassionElimination( pr_iSizeParameter , &l_fFirstMultiDimArray );

#if DEBUG		
	display1DimMatrixFunc(pr_iSizeParameter, &g_fGuassionFinal);
#endif
	//free the used memory
	free(g_fGuassionFinal);
	trashMemory( pr_iSizeParameter ,  (void *) &l_fFirstMultiDimArray);


}//perform Functionality ends


int main(int argc, char *argv[]){

	int l_iSizeParameter = -1 ;
	if( argc == 2 ){
		l_iSizeParameter = atoi(argv[1]);

#if DEBUG
		performDataTypeOperation( l_iSizeParameter );
#endif
		if( checkInputLimits( l_iSizeParameter )){
			// Defination of multi dimensional array 
			performDataTypeOperation( l_iSizeParameter );
			printf("\n sucess ");
		}
		else{

			printf("\n Matrix size range should be 64 - 4096 ");
			exit(-1);	
		}
	}//argc if ends
	else{
		printf("\n Need Matrix size (Two argument expected i.e size & dataType)!");
		exit(-1);	
	}	
	return -1;
}//main ends

