//Author: 	Navdeep 79776167
//About	:	Vanilla matrix multiply
//Date	:	Sep 20 2017
/*Decription :
Input : An m*n matrix A and an n*p matrix B
Output: The product of AB
 */



//Header files
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


//global const
#define MIN_SIZE 64
#define MAX_SIZE 4096
#define DEBUG 0
#define NO_OPTIMIZATION 0 
#define TILED_OPTIMIZATION 1 
#define MISC_OPTIMIZATION 0 
#define TILE_SIZE 32 
#define MIN(x,y)((x)<(y)? x:y )
#define MAX(x,y)((x)>(y)? x:y )




//global vraibles
float **g_fMultiDimArray = NULL;

// main funtionality methods
int checkInputLimits(int pr_iSizeParameter){
	return ( MIN_SIZE <= pr_iSizeParameter && pr_iSizeParameter <= MAX_SIZE) ;
}//checkInputLimits ends

void defineFloatMultiDimensionalArrayFunc( float ***pr_fMultiDimArray , int pr_iSizeParameter){
	int l_iBound;
	// Defination of multi dimensional array 
	*pr_fMultiDimArray = (float **) calloc( pr_iSizeParameter ,  sizeof(double*));
	if(pr_fMultiDimArray == NULL){
		printf(" Error : Memory not allocated ");
		exit(0);	
	} 

	for ( l_iBound = 0 ; l_iBound < pr_iSizeParameter ; l_iBound++ )
	{	
		(*pr_fMultiDimArray)[l_iBound] = (float *) calloc( pr_iSizeParameter ,  sizeof( double));
	}//for loop ends
}//defineFloatMultiDimensionalArray Ends

///////////////////////////////// DISPLAY FUNCTION ////////////////////////////////////
void displayMatrixFunc(int pr_iSizeParameter , char *pr_cDataType , float ***pr_fMultiDimArray, int ***pr_iMultiDimArray , double ***pr_dMultiDimArray){

	int l_iOuterBound , l_iInnerBound ;

	srand(time(NULL)); // intialize the random seed
	for ( l_iOuterBound = 0 ; l_iOuterBound < pr_iSizeParameter ; l_iOuterBound++ ){
		printf("\n");
		for ( l_iInnerBound = 0 ; l_iInnerBound < pr_iSizeParameter ; l_iInnerBound++ ){
				printf( "%.3f ", (float) (*pr_fMultiDimArray)[l_iOuterBound][l_iInnerBound]);
		}//innerBound Ends
	}//outerBound ends
}//displayMatrix ends

//////////////////////// INITALIZE THE REQUIRED MATRIX //////////////////////////////////////////
void initializeMatrixFunc(int pr_iSizeParameter, char *pr_cDataType , float ***pr_fMultiDimArray , int ***pr_iMultiDimArray , double ***pr_dMultiDimArray){
	int l_iOuterBound , l_iInnerBound ;
	for ( l_iOuterBound = 0 ; l_iOuterBound < pr_iSizeParameter ; l_iOuterBound++ ){
		for ( l_iInnerBound = 0 ; l_iInnerBound < pr_iSizeParameter ; l_iInnerBound++ ){
				(*pr_fMultiDimArray)[l_iOuterBound][l_iInnerBound] = ( (float)rand()/(float) (RAND_MAX/10)  );
		}//innerBound Ends
	}//outerBound ends
}//intializeMatrix ends
///////////////////////////// VANILLA MULTIPLICATION /////////////////////////////////////////////
void  performFloatVanillaMultiplication(int pr_iSizeParameter , float ***pr_vFirstMatrix , float ***pr_vSecondMatrix){
	int l_iOuterLoop , l_iIntermediateLoop , l_iInnerLoop ; 
	// CODE TO CHECKED FOR THE RUNTIME MATRIX PART
	clock_t l_cBegin = clock();
	for( l_iOuterLoop = 0 ; l_iOuterLoop < pr_iSizeParameter ; l_iOuterLoop++){
		for( l_iIntermediateLoop = 0 ; l_iIntermediateLoop < pr_iSizeParameter ; l_iIntermediateLoop++){
			for( l_iInnerLoop = 0 ; l_iInnerLoop < pr_iSizeParameter ; l_iInnerLoop++){
#if DEBUG
				printf(" %d %d %d ",l_iOuterLoop,l_iIntermediateLoop,l_iInnerLoop);	
#endif			
				//C[i1,i3] += A[i1,i2] * B[i2,i3];    
				g_fMultiDimArray[l_iOuterLoop][l_iInnerLoop] += (*pr_vFirstMatrix)[l_iOuterLoop][l_iIntermediateLoop]+
					(*pr_vSecondMatrix)[l_iIntermediateLoop][l_iInnerLoop];
			}//Inner loop ends
#if DEBUG
			printf("\n");
#endif	
		}//Intermediate Loop ends
#if DEBUG
		printf("-----\n");
#endif	
	}//Outer loop ends 
	clock_t l_cEnd = clock();
	printf(" Runtime For Float Computation : %.3f sec \n ",((double)(l_cEnd-l_cBegin)/ CLOCKS_PER_SEC));
}//performVanillaMultiplication


void  performFloatTiledVanillaMultiplication(int pr_iSizeParameter , float ***pr_vFirstMatrix , float ***pr_vSecondMatrix){
	int l_iOuterLoop , l_iIntermediateLoop , l_iInnerLoop ; 
	int l_iTiledLoopFirst , l_iTiledLoopSecond ;
	// CODE TO CHECKED FOR THE RUNTIME MATRIX PART
	clock_t l_cBegin = clock();
	for( l_iTiledLoopFirst = 0 ; l_iTiledLoopFirst < pr_iSizeParameter ; l_iTiledLoopFirst+=TILE_SIZE ){ 
		for( l_iTiledLoopSecond = 0 ; l_iTiledLoopSecond < pr_iSizeParameter ; l_iTiledLoopSecond+=TILE_SIZE ){ 
			for( l_iOuterLoop = 0 ; l_iOuterLoop < pr_iSizeParameter ; l_iOuterLoop++){
				for( l_iIntermediateLoop = l_iTiledLoopFirst ; l_iIntermediateLoop < MIN((l_iTiledLoopFirst+(TILE_SIZE - 1)),pr_iSizeParameter) ; l_iIntermediateLoop++){
					for( l_iInnerLoop =l_iTiledLoopSecond ; l_iInnerLoop < MIN((l_iTiledLoopSecond+(TILE_SIZE - 1)),pr_iSizeParameter) ; l_iInnerLoop++){
#if DEBUG
						printf(" %d %d %d %d %d ",l_iTiledLoopFirst,l_iTiledSecondLoop, l_iOuterLoop,l_iIntermediateLoop,l_iInnerLoop);	
#endif
						//C[i1,i3] += A[i1,i2] * B[i2,i3];    
						g_fMultiDimArray[l_iOuterLoop][l_iInnerLoop] += (*pr_vFirstMatrix)[l_iOuterLoop][l_iIntermediateLoop]+
							(*pr_vSecondMatrix)[l_iIntermediateLoop][l_iInnerLoop];
					}//Inner loop ends
#if DEBUG
					printf("\n");
#endif	
				}//Intermediate Loop ends
#if DEBUG
				printf("-----\n");
#endif	
			}//Outer loop ends
		}//seconds loop ends
	}//first loop ends 
	clock_t l_cEnd = clock();
	printf(" Runtime For Tiled Float Computation : %.3f sec \n ",((double)(l_cEnd-l_cBegin)/ CLOCKS_PER_SEC));
}//performVanillaMultiplication

//////////////////////// CUSTOMIZE FREE FUNCTION ////////////////////////////////////////////////
void trashMemory(int pr_iSizeParameter , void ***pr_vMultiDimFreeMemory){
	int l_iBound ;	
	for ( l_iBound = 0 ; l_iBound < pr_iSizeParameter ; l_iBound++ ){
		free((*pr_vMultiDimFreeMemory)[l_iBound]);		
	}
	free((*pr_vMultiDimFreeMemory));		
}//trashMemory Ends

///////////////////////// CONTROLLER FUNCTIONALITY ////////////////////////////////////////////
void performDataTypeOperation(int pr_iSizeParameter , char *pr_cDataType){
	if( pr_cDataType != NULL){
		if( strcmp( pr_cDataType, "float") == 0 ){
			defineFloatMultiDimensionalArrayFunc( &g_fMultiDimArray ,pr_iSizeParameter);
#if DEBUG	
			displayMatrixFunc(pr_iSizeParameter , pr_cDataType , &g_fMultiDimArray,NULL,NULL);
#endif	
			//Initialize the other operands for multiplication	
			float **l_fFirstMultiDimArray = NULL  , **l_fSecMultiDimArray = NULL ;
			// Frist 2d Array
			defineFloatMultiDimensionalArrayFunc( &l_fFirstMultiDimArray ,pr_iSizeParameter);
			//intializing the array with random values
			initializeMatrixFunc(pr_iSizeParameter, pr_cDataType , &l_fFirstMultiDimArray , NULL ,NULL);
#if DEBUG
			displayMatrixFunc(pr_iSizeParameter , pr_cDataType , &l_fFirstMultiDimArray,NULL,NULL);
#endif
			// second 2D array
			defineFloatMultiDimensionalArrayFunc( &l_fSecMultiDimArray ,pr_iSizeParameter);
			//intializing the array with random values
			initializeMatrixFunc(pr_iSizeParameter, pr_cDataType , &l_fSecMultiDimArray , NULL ,NULL);
#if DEBUG		
			displayMatrixFunc(pr_iSizeParameter , pr_cDataType , &l_fSecMultiDimArray,NULL,NULL);
#endif
			//perform vanilla multiplications
#if NO_OPTIMIZATION 
			performFloatVanillaMultiplication( pr_iSizeParameter , &l_fFirstMultiDimArray , &l_fSecMultiDimArray);
#endif
#if TILED_OPTIMIZATION 
			performFloatTiledVanillaMultiplication( pr_iSizeParameter , &l_fFirstMultiDimArray , &l_fSecMultiDimArray);
#endif

#if DEBUG		
			displayMatrixFunc(pr_iSizeParameter , pr_cDataType , &g_fMultiDimArray,NULL,NULL);
#endif
			//free the used memory
			trashMemory( pr_iSizeParameter ,(void *) &g_fMultiDimArray);
			trashMemory( pr_iSizeParameter ,(void *) &l_fFirstMultiDimArray);
			trashMemory( pr_iSizeParameter ,(void *) &l_fSecMultiDimArray);

		}
		else{
			printf("\n Error : Invalid Datatype argument (3)");
			exit(-1);	
		}
	}//if loop ends
	else{

		printf("\n Error : Check Datatype argument (3)");
		exit(-1);	
	}//else ends
}//perform Functionality ends


void main(int argc, char *argv[]){

	int l_iSizeParameter = -1 ;
	if( argc == 2 ){
		l_iSizeParameter = atoi(argv[1]);
		if( checkInputLimits( l_iSizeParameter )){
			// Defination of multi dimensional array 
			performDataTypeOperation( l_iSizeParameter , "float");
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

}//main ends

