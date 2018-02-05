#include <stdio.h>
#include "vbx.h"


const int num_elements=10;

int main()
{

#if VBX_SIMULATOR==1
	//initialize with 4 lanes,and 64kb of sp memory
	//word,half,byte fraction bits 16,15,4 respectively
	vbxsim_init( 4, 64, 256,6,5, 4 );
#endif

	//Allocate vectors in scratchpad
	vbx_word_t* a = vbx_sp_malloc( num_elements*sizeof(vbx_word_t) );
	vbx_word_t* b = vbx_sp_malloc( num_elements*sizeof(vbx_word_t) );
	vbx_word_t* c = vbx_sp_malloc( num_elements*sizeof(vbx_word_t) );

	//Set vector length, then compute 4*[1,2,3,...,10]
	vbx_set_vl(num_elements);
	vbx(SEW,VADD,a,1,0); //a = [1,2,3,...,10]
	vbx(SVW,VMOV,b,4,0); //b = [4,4,....,4]
	vbx(VVW,VMUL,c,a,b); //c = a * b

	//wait for all vector instructions to finish
	vbx_sync();

	//print out vector c
	printf( "%8d", c[0] );
	int i;
	for( i=1; i<num_elements; i++ ) {
		printf( ",%8d", c[i] );
	}
	printf( "\n" );

	return 0;
}
