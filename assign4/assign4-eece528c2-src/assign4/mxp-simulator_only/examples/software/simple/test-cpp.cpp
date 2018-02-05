#include <stdio.h>
#include "vbx.h"
#include "Vector.hpp"

const int num_elements=10;

int main()
{

#if VBX_SIMULATOR==1
	//initialize with 4 lanes,and 64kb of sp memory
	//word,half,byte fraction bits 16,15,4 respectively
	vbxsim_init( 4, 64, 256, 16, 15, 4 );
#endif

	//Allocating vectors in scratchpad
	VBX::Vector<vbx_word_t> a(num_elements);
	VBX::Vector<vbx_word_t> b(num_elements);
	VBX::Vector<vbx_word_t> c(num_elements);

	a = 1 + VBX::ENUM; //a = [1,2,3,...,10]
	b = 4;             //b = [4,4,....,4]
	c = a * b;

	//wait for all vector instructions to finish
	vbx_sync();

	//print out vector c
	c.printVec();
	printf("\n");

	return 0;
}
