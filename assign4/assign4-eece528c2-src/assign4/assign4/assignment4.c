//Author: 	Navdeep 79776167
//About	:	Vectorblox  matrix multiply
//Date	:	Nov 13 2017
/*Decription :
Input : An m*n matrix A and an n*p matrix B
Output: The product of AB
 */
//NORMAL RUN -	./assignmwnt4.elf matrix_size 
//-Ofast -fforward-propagate -fthread-jumps:q
//	TEST DATA - set the TEST_DATA flag
//RUN -	 ./assignment4.elf 
//Header files
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "vbx.h"
#include "vbx_port.h"

//global const
#define MIN_SIZE 64
#define MAX_SIZE 4096
#define DEBUG 0
#define DEBUG1 0
#define TILE_SIZE 32
#define TEST_DATA 0
#define TEST_ELEMENT 4
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

void checkAnswer( int  **prStrSharedData ){

        FILE *l_ptrFileOutput ;
        l_ptrFileOutput = fopen("result.txt","r+");
        if( l_ptrFileOutput == NULL ){
                printf("\n Error : In opening file ");
        }else{
                int* result_vector = (int*) malloc(sizeof(int)*TEST_ELEMENT);
                int row_sum;

                for (int j=0; j<TEST_ELEMENT; j++){
                        row_sum = 0;
                        for (int k=0; k<=TEST_ELEMENT; k++){
                                row_sum += (*prStrSharedData)[j*TEST_ELEMENT+k];
                        }
                        result_vector[j] = row_sum;
                        printf(" \n result_vector[ %d ]",result_vector[j]);
                }

                float sumOfSquares = 0;
                int entryOfResidual;
                int element = 0 ;
                for (int i=0; i<TEST_ELEMENT; i++){
                        int scn = fscanf(l_ptrFileOutput,"%d", &element) ;
                        if ( scn == -1)
				printf("\n Error in benchmark file");
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
void initializeMatrixFunc(int pr_iSizeParameter , int **pr_fMultiDimArray ){
	int l_iOuterBound   ;
	if(TEST_DATA){
		for ( l_iOuterBound = 0 ; l_iOuterBound < (pr_iSizeParameter) ; l_iOuterBound++ ){
			(*pr_fMultiDimArray)[l_iOuterBound] =  l_iOuterBound ;
		}//outerBound ends
		
	}else{
		for ( l_iOuterBound = 0 ; l_iOuterBound < (pr_iSizeParameter) ; l_iOuterBound++ ){
			(*pr_fMultiDimArray)[l_iOuterBound] = ( (int)rand()/(int) (RAND_MAX/10)  );
		}//outerBound ends
	}
}//intializeMatrix ends
///////////////////////////// VANILLA MULTIPLICATION /////////////////////////////////////////////
void  performFloatTiledVanillaMultiplicationVector(int pr_iSizeParameter , int **pr_vFirstMatrix , int **pr_vSecondMatrix){
	int l_iOuterLoop , l_iIntermediateLoop  ; 
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
	// MOdifying the tile size 
	// in case if the parameters are less than the tile 
	int l_iTileSize = row < TILE_SIZE ? row:TILE_SIZE; 
	int l_iNumSize = l_iTileSize * sizeof(int);
	// defining the parameter on the scratch memory 
	vbx_word_t *va = vbx_sp_malloc( l_iNumSize);
	vbx_word_t *vb = vbx_sp_malloc( l_iNumSize);
	vbx_word_t *vr = vbx_sp_malloc( sizeof(int) );
	// setting the numbe rof elements to the send onto scractch pad 
	vbx_set_vl( l_iTileSize );

	vbx_timestamp_start();
	l_uFrequency = vbx_timestamp_freq();
	l_vbxStart = vbx_timestamp();
#if DEBUG
	printf("[%d,%d]- Tile- %d Size- %d\n",row,col,l_iTileSize,l_iNumSize);
#endif    
	/*
		At each instance single row been trnasferesd with corresponding block on the scratch pad for matrix A
		After dat, the ssecond matrix will be fecth to get the repsective column elements of the row
		Its bsically , 1-N mapping for matrix A to B	
	*/
 
	for (l_iTiledLoopFirst = 0; l_iTiledLoopFirst < row; l_iTiledLoopFirst += l_iTileSize) {
		printf("BLock Size = %d \n",l_iTiledLoopFirst);
		for (l_iTiledLoopSecond = 0; l_iTiledLoopSecond < col; l_iTiledLoopSecond += l_iTileSize) {
			for (l_iTiledLoopThird = 0; l_iTiledLoopThird < row; l_iTiledLoopThird += l_iTileSize) {

				for (l_iOuterLoop = l_iTiledLoopFirst; l_iOuterLoop < MIN( l_iTiledLoopFirst + l_iTileSize, row ); l_iOuterLoop++) {
					//fecthing new blocked row of MAtrix A to scracth pad
					vbx_dcache_flush(va,l_iNumSize);
					vbx_dma_to_vector( va,(*pr_vFirstMatrix)+( l_iOuterLoop * col + l_iTiledLoopThird ),  l_iNumSize);
					for (l_iIntermediateLoop = l_iTiledLoopSecond; l_iIntermediateLoop < MIN( l_iTiledLoopSecond + l_iTileSize, col ); l_iIntermediateLoop++) {
						(*var) = 0 ;
						// fetching the new blocked row of matrix B to scratch pad
						vbx_dcache_flush(vb,l_iNumSize);
						vbx_dma_to_vector( vb, (*pr_vSecondMatrix)+( l_iIntermediateLoop * col  + l_iTiledLoopThird ), l_iNumSize);
						vbx_acc( VVW, VMUL, vr, va, vb);
						//getting the accumulated result out of the vectorblox
						vbx_dma_to_host( var, vr, sizeof(int));
						vbx_sync();
						// setting the result into the resultant matrix
						g_fMultiDimArray[ l_iOuterLoop * col + l_iIntermediateLoop ] += (*var);
#if DEBUG1
						printf("\n %d = [ %d * %d ] ---> %d ",(l_iOuterLoop * col + l_iTiledLoopThird  ),(*pr_vFirstMatrix)[(l_iOuterLoop * col + l_iIntermediateLoop)],(*pr_vSecondMatrix)[(l_iIntermediateLoop * col +l_iTiledLoopThird ) ],(*var));
#endif

					}//l_iIntermediateLoop ends
				}//l_iOuterLoop ends
			}//l_iTiledLoopThird ends
		}//l_iTiledLoopSecond ends
	}//l_iTiledLoopFirst ends

	l_vbxEnd = vbx_timestamp();
	vbx_sp_free();
	//getting the stats of the vectorblox
	vbxsim_get_stats();
	vbxsim_print_stats();
	vbxsim_print_stats_extended();
	printf("\n  %f Time Consumed - %f " ,(float) (l_vbxEnd-l_vbxStart), ((float)(l_vbxEnd-l_vbxStart)/l_uFrequency * 1000 ));
	clock_t l_cEnd = clock();
	printf(" Runtime For Tiled Float Computation : %f sec \n ",((double)(l_cEnd-l_cBegin)/ CLOCKS_PER_SEC));
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
	defineFloatMultiDimensionalArrayFunc( &g_fMultiDimArray ,l_iTotalData);
#if DEBUG1	
	displayMatrixFunc(l_iTotalData  , &g_fMultiDimArray);
#endif	
	//Initialize the other operands for multiplication	
	int *l_fFirstMultiDimArray = NULL  , *l_fSecMultiDimArray = NULL , *l_fSecMultiDimArrayTrans= NULL ;;
	// Frist 2d Array
	defineFloatMultiDimensionalArrayFunc( &l_fFirstMultiDimArray ,l_iTotalData);
	//intializing the array with random values
	initializeMatrixFunc(l_iTotalData , &l_fFirstMultiDimArray  );
#if DEBUG1
	displayMatrixFunc(l_iTotalData  , &l_fFirstMultiDimArray);
#endif
	// second 2D array
	defineFloatMultiDimensionalArrayFunc( &l_fSecMultiDimArray ,l_iTotalData);
	//intializing the array with random values
	initializeMatrixFunc(l_iTotalData , &l_fSecMultiDimArray );

	defineFloatMultiDimensionalArrayFunc( &l_fSecMultiDimArrayTrans ,l_iTotalData);
	memcpy(l_fSecMultiDimArrayTrans, l_fSecMultiDimArray, sizeof(int)*l_iTotalData);
#if DEBUG1		
	displayMatrixFunc(l_iTotalData  , &l_fSecMultiDimArray);
#endif
	changeMatrixOrder(pr_iSizeParameter , &l_fSecMultiDimArray,  &l_fSecMultiDimArrayTrans);
#if DEBUG1		
	displayMatrixFunc(l_iTotalData  , &l_fSecMultiDimArray);
#endif
	//perform vanilla multiplications
	performFloatTiledVanillaMultiplicationVector( l_iTotalData  , &l_fFirstMultiDimArray , &l_fSecMultiDimArray);

#if DEBUG1		
	displayMatrixFunc(l_iTotalData  , &g_fMultiDimArray);
#endif

	if( TEST_DATA )
		checkAnswer( &g_fMultiDimArray);

	//free the used memory
	free(g_fMultiDimArray);
	free(l_fFirstMultiDimArray);
	free(l_fSecMultiDimArray);
	free(l_fSecMultiDimArrayTrans);

	g_fMultiDimArray = NULL;
	l_fFirstMultiDimArray = NULL;
	l_fSecMultiDimArray = NULL;
	l_fSecMultiDimArrayTrans = NULL;
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

		if (TEST_DATA)
		{
			performDataTypeOperation( TEST_ELEMENT );
			printf("\n sucess ");

		}else{
			printf("\n Need Matrix size (Two argument expected i.e size & dataType)!");
			exit(-1);
		}	
	}	

}//main ends

