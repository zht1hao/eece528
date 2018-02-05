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


//global const
#define MIN_SIZE 64
#define MAX_SIZE 4096
#define DEBUG 0
#define TEST_DATA 0
#define RESULT_COLUMN 1
#define MIN(x,y)((x)<(y)? x:y )
#define MAX(x,y)((x)>(y)? x:y )




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
void performUpperTraingular(int pr_iSizeParameter , float ***pr_vMatrix ){

	int l_iOuterLoop , l_iIntermediateLoop , l_iInnerLoop ;
	float l_fCompute;
	// CODE TO CHECKED FOR THE RUNTIME MATRIX PART
	clock_t l_cBegin = clock();
	for( l_iOuterLoop = 0 ; l_iOuterLoop < pr_iSizeParameter ; l_iOuterLoop++){
		for( l_iIntermediateLoop = 0 ; l_iIntermediateLoop < pr_iSizeParameter ; l_iIntermediateLoop++){
			if( l_iIntermediateLoop > l_iOuterLoop ){
				l_fCompute = (*pr_vMatrix)[l_iIntermediateLoop][l_iOuterLoop] / (*pr_vMatrix)[l_iOuterLoop][l_iOuterLoop];

				for( l_iInnerLoop = 0 ; l_iInnerLoop < pr_iSizeParameter+RESULT_COLUMN ; l_iInnerLoop++){
#if DEBUG
					printf(" %d %d %d ",l_iOuterLoop,l_iIntermediateLoop,l_iInnerLoop);	
#endif			
					(*pr_vMatrix)[l_iIntermediateLoop][l_iInnerLoop] -=  l_fCompute * (*pr_vMatrix)[l_iOuterLoop][l_iInnerLoop];
				}//Inner loop ends
#if DEBUG
				printf("\n");
#endif	
			}//if ends

		}//Intermediate Loop ends
#if DEBUG
		printf("-----\n");
#endif	
	}//Outer loop ends 
	clock_t l_cEnd = clock();
	g_dFirstTimeTaken = ((double)(l_cEnd-l_cBegin)/ CLOCKS_PER_SEC);
#if DEBUG
	printf("\n Runtime For performUpperTraingular : %.3f sec \n ",(g_dFirstTimeTaken);
#endif	
}

void performBackwardSubstitution(int pr_iSizeParameter , float ***pr_vMatrix , float **pr_vResultMatrix ){

	int l_iOuterLoop, l_iInnerLoop ;
	float l_fCompute;
	// CODE TO CHECKED FOR THE RUNTIME MATRIX PART
	 clock_t l_cBegin = clock();

	(*pr_vResultMatrix)[pr_iSizeParameter-1]= (*pr_vMatrix) [pr_iSizeParameter-1][pr_iSizeParameter]/(*pr_vMatrix)[pr_iSizeParameter-1][pr_iSizeParameter-1];

	for( l_iOuterLoop = pr_iSizeParameter-2 ; l_iOuterLoop >= 0 ; l_iOuterLoop-- ){
		l_fCompute = 0 ;
		for( l_iInnerLoop = l_iOuterLoop+1 ; l_iInnerLoop < pr_iSizeParameter ; l_iInnerLoop++){
#if DEBUG
			printf(" %d %d ",l_iOuterLoop,l_iInnerLoop);	
#endif			
			l_fCompute += (*pr_vMatrix)[l_iOuterLoop][l_iInnerLoop] * (*pr_vResultMatrix)[l_iInnerLoop];
		}//Inner loop ends
#if DEBUG
		printf("\n");
#endif	
		(*pr_vResultMatrix)[l_iOuterLoop] = ((*pr_vMatrix)[l_iOuterLoop][pr_iSizeParameter]- l_fCompute)/ (*pr_vMatrix)[l_iOuterLoop][l_iOuterLoop];
	}//Outer loop ends 
	clock_t l_cEnd = clock();
	g_dSecondTimeTaken = ((double)(l_cEnd-l_cBegin)/ CLOCKS_PER_SEC);
#if DEBUG
	printf("\n Runtime For performBackwardSubstitution : %.3f sec \n ",(g_dSecondTimeTaken);
#endif	

}// performGuassionElimination ends


void performGuassionElimination(int pr_iSizeParameter , float ***pr_vFirstMatrix ){
	performUpperTraingular(pr_iSizeParameter , pr_vFirstMatrix );
	performBackwardSubstitution(pr_iSizeParameter , pr_vFirstMatrix , &g_fGuassionFinal );
	//printf("\n Runtime For Guassion Elimination : %.3f sec \n ",(g_dFirstTimeTaken + g_dSecondTimeTaken ));
	printf("\n Runtime For Guassion Elimination : %.3f sec \n ",(g_dFirstTimeTaken));
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

