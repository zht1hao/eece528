#include <stdio.h>
#include "vbx.h"


int A[] = {1, 1, 1, 1,
		2, 2, 2,2,
		 3,3,3,3, 
		4,4,4,4 };
int B[] = {5,5,5,5, 
		6,6,6,6, 
		7,7,7,7, 
		8,8,8,8};
int C[] = {-1, -1, -1, -1,
		-1,-1,-1,-1,
		-1,-1,-1,-1,
		-1,-1,-1,-1};


int vbw_mtx_xp(vbx_half_t *v_dst, vbx_half_t *v_src, const int INROWS, const int INCOLS )
{
        vbx_set_vl( 1 );
        vbx_set_2D(  INCOLS, INROWS*sizeof(vbx_half_t),       sizeof(vbx_half_t), 0 );
        vbx_set_3D( INROWS,        sizeof(vbx_half_t), INCOLS*sizeof(vbx_half_t), 0 );
        vbx_3D( VVH, VMOV, v_dst, v_src, NULL);
        return 1;
}

 int vbw_mtx_mul( vbx_half_t *v_out, vbx_half_t *v_in1, const int rows1, const int cols1, vbx_half_t *v_in2, const int rows2, const int cols2 )
{       
        // TODO check for overlapping.
        
        const int OUTROWS = rows1;
        const int OUTCOLS = cols2;
        
        vbx_sp_push();
        
        int in2_bytes = cols2*rows2*sizeof(vbx_half_t);
        vbx_half_t *v_in2_trans = (vbx_half_t *)vbx_sp_malloc(in2_bytes);
        if (v_in2_trans==NULL){
                return  -1;
        }
        vbw_mtx_xp(v_in2_trans, v_in2, rows2, cols2);
        vbx_set_vl(cols1);
        vbx_set_2D(cols2, sizeof(vbx_half_t), 0, cols2*sizeof(vbx_half_t));
        vbx_set_3D(OUTROWS, OUTCOLS*sizeof(vbx_half_t), cols1*sizeof(vbx_half_t), 0);
        vbx_acc_3D(VVH,VMUL,v_out,v_in1,v_in2_trans);
        vbx_sp_pop();
        return 1;
}

int main(){
#if VBX_SIMULATOR==1
        //initialize with 4 lanes,and 64kb of sp memory
        //word,half,byte fraction bits 16,15,4 respectively
        vbxsim_init( 4, 64, 256,6,5, 4 );
#endif

	 int vector_len = 4 ;
    int num_bytes = vector_len * sizeof(vbx_word_t);

vbx_word_t *va = vbx_sp_malloc( num_bytes);
vbx_word_t *vb = vbx_sp_malloc( num_bytes);
vbx_word_t *vr = vbx_sp_malloc( num_bytes);
vbx_word_t *vv = vbx_sp_malloc( sizeof(vbx_word_t));



/*
vbw_mtx_mul(  vr,
 va, 4 , 4 , 
 vb, 4 , 4 );
*/

printf("\n ======== A ======== \n");
for ( int i = 0 ; i < 16; i ++)
{if ( i % 4 == 0 ) printf("\n");  
	printf(" %d ",A[i]);	
}

printf("\n ======== B ======== \n");
for ( int i = 0 ; i < 16; i ++)
{if ( i % 4 == 0 ) printf("\n");  
	printf(" %d ",B[i]);	
}
printf("\n SIZE = %d ",num_bytes);   

#if 1 
//for ( int i = 0 ; i < 4; i++)
{ 
vbx_set_vl( vector_len );
vbx_set_2D( 2, 2*sizeof(vbx_word_t), 2*sizeof(vbx_word_t), 2*sizeof(vbx_word_t) );
//vbx_dma_to_vector( va, A, 2*sizeof(vbx_word_t));
vbx_dma_to_vector_2D( va, A, num_bytes, 2/*noOfRows*/, /* row element */2*sizeof(vbx_word_t) ,/*row element*/ 2* sizeof(int) );
vbx_dma_to_vector_2D( vb, B, num_bytes, 2/*noOFRows*/, /* row element */2*sizeof(vbx_word_t), /*row element*/ 2*sizeof(int) );

/* vbx_set_2D( 2,sizeof(int),sizeof(int),sizeof(int));
*/
                vbx_2D( VVW, VMUL,
                                vr,
                                 va,
				vb);
vbx_dma_to_host_2D( C, vr, num_bytes, 2, 2*sizeof(vbx_word_t), 2*sizeof(vbx_word_t) );
//vbx_dma_to_host( C, vv, vector_len);
vbx_sync();

}
#endif
#if 0 
vbx_set_vl( vector_len );
 vbx_dma_to_vector( vb, B ,num_bytes  );
 vbx_dma_to_vector( va, A ,num_bytes  );
 vbx_acc( VVW, VMUL, vr, va, vb);

vbx_dma_to_host( C, vr, num_bytes);
vbx_sync();
#endif
printf(" \n ======== C ======== \n");
#if 0
for ( int i = 0 ; i < 4 ;  i++)
{
                vbx_set_2D( BLOCK_HEIGHT,
                                sizeof(uint16_t),
                                (BLOCK_WIDTH+SEARCH_WIDTH)*sizeof(uint8),
                                BLOCK_WIDTH*sizeof(uint8) );
                vbx_set_3D( SEARCH_WIDTH,
                                BLOCK_HEIGHT*sizeof(uint16_t),
                                sizeof(uint8_t),
                                0 );
                // 3D vector operation with reduction accumulate
printf(":===> %d",vr);	
}
#endif
for ( int i = 0 ; i < 16; i ++)
{if ( i % 4 == 0 ) printf("\n");  
	printf(" %d ",C[i]);	
}
#if 0
    /* step 1 */
    vbx_word_t * v_a = vbx_sp_malloc( num_bytes );
    vbx_word_t * v_b = vbx_sp_malloc( num_bytes );
    vbx_word_t * v_c = vbx_sp_malloc( num_bytes );
printf("\n Passed");
    /* step 2 */
    vbx_dma_to_vector( v_a, A, num_bytes );
    vbx_dma_to_vector( v_b, B, num_bytes );
printf("\n Passed");

    /* step 3 */
    vbx_set_vl( vector_len );
    vbx( VVW, VADD, v_c, v_a, v_b );
printf("\n Passed");

    /* step 4 */
    vbx_dma_to_host( C, v_c, num_bytes );
printf("\n Passed");

    /* step 5 */
    vbx_sp_free();
printf("\n Passed");

    vbx_sync();
printf("\n Passed");
    printf( "C[] = %d, %d, %d, %d\n",
    C[0], C[1], C[2], C[3] );

#endif
#if 0
int BLOCK_WIDTH = 2, BLOCK_HEIGHT = 2, SEARCH_HEIGHT = 2;
	for( y = 0; y < SEARCH_HEIGHT; y++ ) {
		// Set vector parameters to compute the SAD for each
		// row in a block across the search space width
		vbx_set_vl( BLOCK_WIDTH );
		vbx_set_2D( BLOCK_HEIGHT,
				sizeof(uint16_t),
				(BLOCK_WIDTH+SEARCH_WIDTH)*sizeof(uint8),
				BLOCK_WIDTH*sizeof(uint8) );
		vbx_set_3D( SEARCH_WIDTH,
				BLOCK_HEIGHT*sizeof(uint16_t),
				sizeof(uint8_t),
				0 );
		// 3D vector operation with reduction accumulate
		vbx_acc_2D( VVBHU, VABSDIFF,
				v_row_sad,
				v_img+y*(BLOCK_WIDTH+SEARCH_WIDTH),
				v_block );
		//Set vector parameters for accumulating the final SAD
		vbx_set_vl( BLOCK_HEIGHT/2 );
		vbx_set_2D( SEARCH_WIDTH,
				sizeof(uint8_t),
				BLOCK_HEIGHT*sizeof(uint16_t),
				BLOCK_HEIGHT*sizeof(uint16_t) );
		//2D vector operation to produces final SAD values
		vbx_acc_2D( VVHWU, VADD,
				v_result+y*SEARCH_WIDTH,
				v_row_sad,
				v_row_sad+(BLOCK_HEIGHT/2) );
		//Transfer the line back to the host
		vbx_dma_to_host( result+y*SEARCH_WIDTH,
				v_result+y*SEARCH_WIDTH,
				SEARCH_WIDTH*sizeof(output_type) );
	}
#endif


}
