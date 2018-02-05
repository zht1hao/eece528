//Author: 	Navdeep 79776167
//About	:	Vectorblox  matrix multiply
//Date	:	Nov 13 2017
/*Decription :
Input : An m*n matrix A and an n*p matrix B
Output: The product of AB
 */

//-Ofast -fforward-propagate -fthread-jumps:q

//Header files
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "vbx.h"
#include "vbx_port.h"

//global const
#define MIN_SIZE 2
#define MAX_SIZE 4096
#define DEBUG 0
#define DEBUG1 0
#define TILE_SIZE 2 
#define MIN(x,y)((x)<(y)? x:y )
#define MAX(x,y)((x)>(y)? x:y )




//global vraibles
int *g_fMultiDimArray = NULL;

// main funtionality methods
int checkInputLimits(int pr_iSizeParameter){
	return ( MIN_SIZE <= pr_iSizeParameter && pr_iSizeParameter <= MAX_SIZE) ;
}//checkInputLimits ends

void defineIntMultiDimensionalArrayFunc( int **pr_fMultiDimArray , int pr_iSizeParameter){
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
	int k, j, i ; 
	int jj, kk ;
	vbx_timestamp_t l_vbxStart, l_vbxEnd;	
	int col, row;
	int *var = (int *)malloc( sizeof(int)) ; 
	unsigned l_uFrequency ;
#if VBX_SIMULATOR==1
	//initialize with 4 lanes,and 64kb of sp memory
	//word,half,byte fraction bits 16,15,4 respectively
	vbxsim_init( 4, 64, 256,6,5, 4 );
#endif
	// IT HAS TO GO INTO SEQUENCE OF THE 64* 1024 FOR TOTAL DMA OPERATIONS
	col = row = sqrt(pr_iSizeParameter);
	int l_iTileSize = row < TILE_SIZE ? row:TILE_SIZE; 
	int l_iNumSize = l_iTileSize * sizeof(int);
	vbx_word_t *va = vbx_sp_malloc( l_iNumSize);
	vbx_word_t *vb = vbx_sp_malloc( l_iNumSize);
	vbx_word_t *vr = vbx_sp_malloc( sizeof(int) );
	vbx_set_vl( l_iTileSize );
	vbx_set_2D( row, row*sizeof(vbx_word_t), row*sizeof(vbx_word_t), row*sizeof(vbx_word_t) );

	vbx_timestamp_start();
	l_uFrequency = vbx_timestamp_freq();
	l_vbxStart = vbx_timestamp();
#if DEBUG
	printf("[%d,%d]- Tile- %d Size- %d\n",row,col,l_iTileSize,l_iNumSize);
#endif    
	printf("[%d,%d]- Tile- %d Size- %d\n",row,col,l_iTileSize,l_iNumSize);





 for (k = 0; k < col; k += l_iTileSize)
        {    
vbx_dma_to_vector_2D( va, (*pr_vFirstMatrix)+( k * col ), num_bytes, l_iTileSize/*noOfRows*/, /* row element */l_iTileSize*sizeof(vbx_word_t) ,/*row element*/ l_iTileSize* sizeof(int) );
printf("\n K = %d ",k); 
                for (j = 0; j < row; j += l_iTileSize)
                {  
vbx_dma_to_vector_2D( vb, (*pr_vFirstMatrix)+( j * col ), num_bytes, l_iTileSize/*noOFRows*/, /* row element */l_iTileSize*sizeof(vbx_word_t), /*row element*/ l_iTileSize*sizeof(int) );
  //              vbx_2D( VVW, VMUL,
    //                            vr,
      //                           va,
        //                        vb);
//vbx_dma_to_host_2D( C, vr, num_bytes, 4, 4*sizeof(vbx_word_t), 4*sizeof(vbx_word_t) );
//vbx_dma_to_host( C, vv, vector_len);
vbx_sync();

	}
}

#if 0
	for (k = 0; k < col; k += l_iTileSize)
	{    
printf("\n K = %d ",k); 
		for (j = 0; j < row; j += l_iTileSize)
		{       
          printf("\n J = %d ",j); 
			for (i = 0; i < row; ++i)
			{ 
//printf("\n ==>> i = %d ",i); 
				for (jj = j; jj < min(j + l_iTileSize, row); ++jj)
				{
		//		printf("\n");
//printf("\n JJ = %d ",jj); 
					for (kk = k; kk < min(k + l_iTileSize, col); ++kk)

					{
		
//printf("\t KK = %d ",kk); 
		//		printf(" %d =  [ %d * %d ]",(i*col+ jj),(*pr_vFirstMatrix)[i*col+kk],(*pr_vSecondMatrix)[kk+jj*col]);
						 g_fMultiDimArray[i*col+ jj] += (*pr_vFirstMatrix)[i*col+kk] * (*pr_vSecondMatrix)[kk+jj* col];
					}
				}
			}
		}
	}
#endif
#if 0 
	for (l_iTiledLoopFirst = 0; l_iTiledLoopFirst < row; l_iTiledLoopFirst += l_iTileSize) {
#if DEBUG
		printf("BLock Size = %d \n",l_iTiledLoopFirst);
#endif
		for (l_iTiledLoopSecond = 0; l_iTiledLoopSecond < col; l_iTiledLoopSecond += l_iTileSize) {
			for (l_iTiledLoopThird = 0; l_iTiledLoopThird < row; l_iTiledLoopThird += l_iTileSize) {
				printf(" l_iTiledLoopThird  - [ %d ] " ,l_iTiledLoopThird);  			
				for (l_iOuterLoop = l_iTiledLoopFirst; l_iOuterLoop < MIN( l_iTiledLoopFirst + l_iTileSize, row ); l_iOuterLoop++) {

					printf(" BLOCK [ %d ]  COLUMN OUTER - %d \n ", l_iOuterLoop, (*pr_vFirstMatrix)[( l_iOuterLoop * col + l_iTiledLoopThird )]);


					//			vbx_dcache_flush(va,l_iNumSize);
					//					vbx_dma_to_vector_2D( va,(*pr_vFirstMatrix)+( l_iOuterLoop * col + l_iTiledLoopThird ),  l_iNumSize,
					//						 row/*noOfRows*/, /* row element */row*sizeof(vbx_word_t) ,/*row element*/ row* sizeof(int) );

					//vbx_dma_to_vector( va,(*pr_vFirstMatrix)+( l_iOuterLoop * col + l_iTiledLoopThird ),  l_iNumSize);
					//vbx_sync();
					for (l_iIntermediateLoop = l_iTiledLoopSecond; l_iIntermediateLoop < MIN( l_iTiledLoopSecond + l_iTileSize, col ); l_iIntermediateLoop++) {
						(*var) = 0 ;

						printf(" BLOCK [ %d ] INNER - %d \n ", l_iIntermediateLoop, (*pr_vSecondMatrix)[( l_iIntermediateLoop * col + l_iTiledLoopThird )]);
				
						//			vbx_dcache_flush(vb,l_iNumSize);
						//						vbx_dma_to_vector_2D( va,(*pr_vSecondMatrix)+( l_iIntermediateLoop * col + l_iTiledLoopThird ),  l_iNumSize,
						//						 row/*noOfRows*/, /* row element */row*sizeof(vbx_word_t) ,/*row element*/ row* sizeof(int) );

						//vbx_dma_to_vector( vb, (*pr_vSecondMatrix)+( l_iIntermediateLoop * col  + l_iTiledLoopThird ), l_iNumSize);
						//						vbx_acc_2D( VVW, VMUL, vr, va, vb);
						//vbx_dcache_flush(vr,sizeof(int));
						//						vbx_dma_to_host( var, vr, sizeof(int));
						//						vbx_sync();

						//						g_fMultiDimArray[ l_iOuterLoop * col + l_iIntermediateLoop ] += (*var);
#if DEBUG1
						printf("\n %d = [ %d * %d ] ---> %d ",(l_iOuterLoop * col + l_iTiledLoopThird  ),(*pr_vFirstMatrix)[(l_iOuterLoop * col + l_iIntermediateLoop)],(*pr_vSecondMatrix)[(l_iIntermediateLoop * col +l_iTiledLoopThird ) ],(*var));
#endif

					}//l_iIntermediateLoop ends
				}//l_iOuterLoop ends
			}//l_iTiledLoopThird ends
		}//l_iTiledLoopSecond ends
	}//l_iTiledLoopFirst ends

#endif
	l_vbxEnd = vbx_timestamp();
	vbx_sp_free();

	#if DEBUG1	
	vbxsim_get_stats();
	vbxsim_print_stats();
	vbxsim_print_stats_extended();
	#endif
	printf("\n  %f Time Consumed - %f " ,(float) (l_vbxEnd-l_vbxStart), ((float)(l_vbxEnd-l_vbxStart)/l_uFrequency * 1000 ));
	clock_t l_cEnd = clock();
	//printf(" Runtime For Tiled Float Computation : %f sec \n ",((double)(l_cEnd-l_cBegin)/ CLOCKS_PER_SEC));
}//performVanillaMultiplication

/////////////////////// PERFORM ORDER ACTIVITY ////////////////////////////
void changeMatrixOrder( int pr_iSizeParameter,  int **pr_vSecondMatrix, int **pr_fSecMultiDimArrayTrans){
	int l_iOuterLoop , l_iInnerLoop ;

	for( l_iOuterLoop = 0 ; l_iOuterLoop < pr_iSizeParameter ; l_iOuterLoop++){
		for( l_iInnerLoop = 0 ; l_iInnerLoop < pr_iSizeParameter ; l_iInnerLoop++){
			(*pr_vSecondMatrix)[l_iOuterLoop * pr_iSizeParameter + l_iInnerLoop ] = (*pr_fSecMultiDimArrayTrans)[ l_iOuterLoop+ l_iInnerLoop * pr_iSizeParameter];
		}//INerr loop ends
	}//outer loop ends
}//changeMatrixOrder


///////////////////////// CONTROLLER FUNCTIONALITY ////////////////////////////////////////////
void performDataTypeOperation(int pr_iSizeParameter  ){


	int l_iTotalData = (pr_iSizeParameter * pr_iSizeParameter);
	defineIntMultiDimensionalArrayFunc( &g_fMultiDimArray ,l_iTotalData);
#if DEBUG1	
	displayMatrixFunc(l_iTotalData  , &g_fMultiDimArray);
#endif	
	//Initialize the other operands for multiplication	
	int *l_fFirstMultiDimArray = NULL  , *l_fSecMultiDimArray = NULL, *l_fSecMultiDimArrayTrans= NULL ;
	// Frist 2d Array
	defineIntMultiDimensionalArrayFunc( &l_fFirstMultiDimArray ,l_iTotalData);
	//intializing the array with random values
	initializeMatrixFunc(l_iTotalData , &l_fFirstMultiDimArray  );
#if DEBUG1
	displayMatrixFunc(l_iTotalData  , &l_fFirstMultiDimArray);
#endif
	displayMatrixFunc(l_iTotalData  , &l_fFirstMultiDimArray);
	// second 2D array
	defineIntMultiDimensionalArrayFunc( &l_fSecMultiDimArray ,l_iTotalData);
	defineIntMultiDimensionalArrayFunc( &l_fSecMultiDimArrayTrans ,l_iTotalData);
	//intializing the array with random values
	initializeMatrixFunc(l_iTotalData , &l_fSecMultiDimArray );
	memcpy(l_fSecMultiDimArrayTrans, l_fSecMultiDimArray, sizeof(int)*l_iTotalData); 
#if DEBUG1		
	displayMatrixFunc(l_iTotalData  , &l_fSecMultiDimArray);
#endif
//	displayMatrixFunc(l_iTotalData  , &l_fSecMultiDimArray);
	displayMatrixFunc(l_iTotalData  , &l_fSecMultiDimArrayTrans);
	changeMatrixOrder(pr_iSizeParameter , &l_fSecMultiDimArray, &l_fSecMultiDimArrayTrans);
#if DEBUG1		
	displayMatrixFunc(l_iTotalData  , &l_fSecMultiDimArray);
#endif
	displayMatrixFunc(l_iTotalData  , &l_fSecMultiDimArray);
	//perform vanilla multiplications
	performFloatTiledVanillaMultiplicationVector( l_iTotalData  , &l_fFirstMultiDimArray , &l_fSecMultiDimArray);

#if DEBUG1		
	displayMatrixFunc(l_iTotalData  , &g_fMultiDimArray);
#endif
	displayMatrixFunc(l_iTotalData  , &g_fMultiDimArray);
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
