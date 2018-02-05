//Author: 	Navdeep 79776167
//About	:	Vectorblox  matrix multiply
//Date	:	Nov 13 2017
/*Decription :
Input : An m*n matrix A and an n*p matrix B
Output: The product of AB
 */



//Header files
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "vbx.h"
#include "vbx_port.h"

//global const
#define MIN_SIZE 2
#define MAX_SIZE 4096
#define DEBUG 1
#define NO_OPTIMIZATION 0 
#define TILED_OPTIMIZATION 1 
#define TILE_SIZE 64 
#define MIN(x,y)((x)<(y)? x:y )
#define MAX(x,y)((x)>(y)? x:y )




//global vraibles
int *g_fMultiDimArray = NULL;

// main funtionality methods
int checkInputLimits(int pr_iSizeParameter){
	return ( MIN_SIZE <= pr_iSizeParameter && pr_iSizeParameter <= MAX_SIZE) ;
}//checkInputLimits ends

void defineFloatMultiDimensionalArrayFunc( int **pr_fMultiDimArray , int pr_iSizeParameter){
	// Defination of multi dimensional array 
	*pr_fMultiDimArray = (int *) calloc(  pr_iSizeParameter ,  sizeof(int));
	if(pr_fMultiDimArray == NULL){
		printf(" Error : Memory not allocated ");
		exit(0);	
	} 
}//defineFloatMultiDimensionalArray Ends

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
void initializeMatrixFunc(int pr_iSizeParameter , int **pr_fMultiDimArray ){
	int l_iOuterBound   ;
	for ( l_iOuterBound = 0 ; l_iOuterBound < (pr_iSizeParameter) ; l_iOuterBound++ ){
				(*pr_fMultiDimArray)[l_iOuterBound] = ( (int)rand()/(int) (RAND_MAX/10)  );
	}//outerBound ends
}//intializeMatrix ends
///////////////////////////// VANILLA MULTIPLICATION /////////////////////////////////////////////
void  performFloatTiledVanillaMultiplicationVector(int pr_iSizeParameter , int **pr_vFirstMatrix , int **pr_vSecondMatrix){
	int l_iOuterLoop , l_iIntermediateLoop , l_iColumnTiledLoop ; 
	int l_iTiledLoopFirst , l_iTiledLoopSecond , l_iTiledLoopThird ;
	vbx_timestamp_t l_vbxStart, l_vbxEnd;	
	int col, row;
	int *var = (int *)malloc( sizeof(int)) ; 
	clock_t l_cBegin = clock();
	unsigned l_uFrequency ;
#if VBX_SIMULATOR==1
        //initialize with 4 lanes,and 64kb of sp memory
        //word,half,byte fraction bits 16,15,4 respectively
        vbxsim_init( 4, 64, 256,6,5, 4 );
#endif
	// IT HAS TO GO INTO SEQUENCE OF THE 64* 1024 FOR TOTAL DMA OPERATIONS

	
	col = row = sqrt(pr_iSizeParameter);
	vbx_timestamp_start();
	l_uFrequency = vbx_timestamp_freq();

	vbx_word_t *va = vbx_sp_malloc( l_iNumSize);
	vbx_word_t *vb = vbx_sp_malloc( l_iNumSize);
	vbx_word_t *vr = vbx_sp_malloc( sizeof(int) );
	//vbx_dcache_flush( (pr_vFirstMatrix), l_iNumSize );
	l_vbxStart = vbx_timestamp();

	int l_iTileSize = row < TILE_SIZE ? row:TILE_SIZE; 
	int l_iNumSize = l_iTileSize * sizeof(int);
	vbx_set_vl( l_iTileSize );


     for (l_iTiledLoopFirst = 0; l_iTiledLoopFirst < row; l_iTiledLoopFirst += l_iTileSize) {
       // printf(" ==== ROW - [ %d ]  ==== " ,l_iTiledLoopFirst);  
	 for (l_iTiledLoopSecond = 0; l_iTiledLoopSecond < col; l_iTiledLoopSecond += l_iTileSize) {
              for (l_iTiledLoopThird = 0; l_iTiledLoopThird < row; l_iTiledLoopThird += l_iTileSize) {
       // printf(" l_iTiledLoopThird  - [ %d ] " ,l_iTiledLoopThird);  
 
                  for (l_iOuterLoop = l_iTiledLoopFirst; l_iOuterLoop < MIN( l_iTiledLoopFirst + l_iTileSize, row ); l_iOuterLoop++) {
      //  printf(" l_iOuterLoop  - [ %d ] " ,l_iOuterLoop);  
			vbx_dcache_flush(va,l_iNumSize);
			vbx_dma_to_vector( va,(*pr_vFirstMatrix)+( l_iOuterLoop * col + l_iTiledLoopThird ),  l_iNumSize);
                      for (l_iIntermediateLoop = l_iTiledLoopSecond; l_iIntermediateLoop < MIN( l_iTiledLoopSecond + l_iTileSize, col ); l_iIntermediateLoop++) {
      //  printf(" l_iIntermediateLoop  - [ %d ] " ,l_iIntermediateLoop);  
		//	printf(" %d",(*pr_vSecondMatrix)[(l_iTiledLoopThird * col)+l_iIntermediateLoop]); 

                          //for (l_iColumnTiledLoop = l_iTiledLoopThird; l_iColumnTiledLoop < MIN( l_iTiledLoopThird + l_iTileSize, row ); l_iColumnTiledLoop++) {
       // printf(" > %d " ,l_iColumnTiledLoop);  
			(*var) = 0 ;
				//vbx_sync();
                 //            printf("[ %d= %d * %d ] ",(x * col + y),(x * col + z),(y * col  + z));
                        //     res[ l_iOuterLoop * col + l_iIntermediateLoop ] +=  a[ l_iOuterLoop * col + l_iColumnTiledLoop ] * temp[ l_iIntermediateLoop * col  + l_iColumnTiledLoop  ];
			vbx_dcache_flush(vb,l_iNumSize);
				vbx_dma_to_vector( vb, (*pr_vSecondMatrix)+( l_iIntermediateLoop * col  + l_iTiledLoopThird ), l_iNumSize);
				vbx_acc( VVW, VMUL, vr, va, vb);
			vbx_dcache_flush(vr,sizeof(int));
			vbx_dma_to_host( var, vr, sizeof(int));
			vbx_sync();
			g_fMultiDimArray[ l_iOuterLoop * col + l_iIntermediateLoop ] += (*var);
#ifdef DEBUG
			printf("\n %d = [ %d * %d ] ---> %d ",(l_iOuterLoop * col + l_iTiledLoopThird  ),(*pr_vFirstMatrix)[(l_iOuterLoop * col + l_iIntermediateLoop)],(*pr_vSecondMatrix)[(l_iIntermediateLoop * col +l_iTiledLoopThird ) ],(*var));
#endif
		
		        }
                  }
              }
          }
      }


l_vbxEnd = vbx_timestamp();
vbx_sp_free();
#if DEBUG1	
vbxsim_get_stats();
vbxsim_print_stats();
vbxsim_print_stats_extended();
#endif
	printf("\n  %f Time Consumed - %f " ,(float) (l_vbxEnd-l_vbxStart), ((float)(l_vbxEnd-l_vbxStart)/l_uFrequency * 1000 ));
	clock_t l_cEnd = clock();
	printf(" Runtime For Tiled Float Computation : %f sec \n ",((double)(l_cEnd-l_cBegin)/ CLOCKS_PER_SEC));
}//performVanillaMultiplication

/////////////////////// PERFORM ORDER ACTIVITY ////////////////////////////
void changeMatrixOrder( int pr_iSizeParameter,  int **pr_vSecondMatrix){
	int l_iOuterLoop , l_iInnerLoop ;

	for( l_iOuterLoop = 0 ; l_iOuterLoop < pr_iSizeParameter ; l_iOuterLoop++){
		for( l_iInnerLoop = 0 ; l_iInnerLoop < pr_iSizeParameter ; l_iInnerLoop++){
			(*pr_vSecondMatrix)[l_iOuterLoop * pr_iSizeParameter + l_iInnerLoop ] = (*pr_vSecondMatrix)[ l_iOuterLoop+ l_iInnerLoop * pr_iSizeParameter];
		}//INerr loop ends
	}//outer loop ends
}//changeMatrixOrder


///////////////////////// CONTROLLER FUNCTIONALITY ////////////////////////////////////////////
void performDataTypeOperation(int pr_iSizeParameter  ){


		int l_iTotalData = (pr_iSizeParameter * pr_iSizeParameter);
		defineFloatMultiDimensionalArrayFunc( &g_fMultiDimArray ,l_iTotalData);
#if DEBUG	
		displayMatrixFunc(l_iTotalData  , &g_fMultiDimArray);
#endif	
		//Initialize the other operands for multiplication	
		int *l_fFirstMultiDimArray = NULL  , *l_fSecMultiDimArray = NULL ;
		// Frist 2d Array
		defineFloatMultiDimensionalArrayFunc( &l_fFirstMultiDimArray ,l_iTotalData);
		//intializing the array with random values
		initializeMatrixFunc(l_iTotalData , &l_fFirstMultiDimArray  );
#if DEBUG
		displayMatrixFunc(l_iTotalData  , &l_fFirstMultiDimArray);
#endif
		// second 2D array
		defineFloatMultiDimensionalArrayFunc( &l_fSecMultiDimArray ,l_iTotalData);
		//intializing the array with random values
		initializeMatrixFunc(l_iTotalData , &l_fSecMultiDimArray );
#if DEBUG		
	//	displayMatrixFunc(l_iTotalData  , &l_fSecMultiDimArray);
#endif
		changeMatrixOrder(pr_iSizeParameter , &l_fSecMultiDimArray);
#if DEBUG		
		displayMatrixFunc(l_iTotalData  , &l_fSecMultiDimArray);
#endif
		//perform vanilla multiplications
		performFloatTiledVanillaMultiplicationVector( l_iTotalData  , &l_fFirstMultiDimArray , &l_fSecMultiDimArray);

#if DEBUG		
		displayMatrixFunc(l_iTotalData  , &g_fMultiDimArray);
#endif
		//free the used memory
		free(g_fMultiDimArray);
		free(l_fFirstMultiDimArray);
		free(l_fSecMultiDimArray);

		g_fMultiDimArray = NULL;
		l_fFirstMultiDimArray = NULL;
		l_fSecMultiDimArray = NULL;
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

}//main ends

