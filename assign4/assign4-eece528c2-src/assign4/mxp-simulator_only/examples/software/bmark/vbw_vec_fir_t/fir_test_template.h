//define VBX_TEMPLATE_T first then include this file
// signed   VBX_BYTESIZE_DEF    VBX_HALFSIZE_DEF    VBX_WORDSIZE_DEF
// unsigned VBX_UBYTESIZE_DEF   VBX_UHALFSIZE_DEF   VBX_UWORDSIZE_DEF

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "vbx.h"
#include "vbx_port.h"
#include "vbx_common.h"
#include "vbx_test.h"

#include "scalar_vec_fir.h"
#include "vbw_vec_fir.h"

//Depends on VBX_TEMPLATE_T
#include "vbw_template_t.h"

//#define USE_TRANSPOSE
#define USE_2D

#define NTAPS     16
#define SAMP_SIZE 0x1000

#define STRINGIFY_INT(SOMETHING) #SOMETHING
#define STRINGIFY(SOMETHING) STRINGIFY_INT(SOMETHING)

#ifdef USE_TRANSPOSE
static double vbx_mm_t_test_vector_transpose( vbx_mm_t *vector_out, vbx_mm_t *sample, vbx_mm_t *coeffs, double scalar_time )
{
	int retval;
	vbx_timestamp_t time_start, time_stop;
	printf("\nExecuting MXP vector transpose FIR.... \n");fflush(stdout);

	vbx_timestamp_start();
	time_start = vbx_timestamp();
	retval = VBX_T(vbw_vec_fir_transpose_ext)( vector_out, sample, coeffs, SAMP_SIZE, NTAPS );
	time_stop = vbx_timestamp();

	printf("...done retval:%X\n", retval);fflush(stdout);
	return vbx_print_vector_time( time_start, time_stop, scalar_time );
}
#endif //ifdef USE_TRANSPOSE


static double vbx_mm_t_test_vector( vbx_mm_t *vector_out, vbx_mm_t *sample, vbx_mm_t *coeffs, double scalar_time )
{
	int retval;
	vbx_timestamp_t time_start, time_stop;
	printf("\nExecuting MXP vector FIR with Accum 2D....\n");fflush(stdout);

	vbx_timestamp_start();
	time_start = vbx_timestamp();
	retval = VBX_T(vbw_vec_fir_ext)( vector_out, sample, coeffs, SAMP_SIZE, NTAPS );
	time_stop = vbx_timestamp();

	printf("...done retval:%X\n", retval);fflush(stdout);
	return vbx_print_vector_time( time_start, time_stop, scalar_time );
}


static double vbx_mm_t_test_scalar( vbx_mm_t *scalar_out, vbx_mm_t *scalar_sample, vbx_mm_t *scalar_coeffs)
{
	vbx_timestamp_t time_start, time_stop;
	printf("\nExecuting scalar vector FIR...\n");fflush(stdout);

	vbx_timestamp_start();
	time_start = vbx_timestamp();
	VBX_T(scalar_vec_fir)( scalar_out, scalar_sample, scalar_coeffs, SAMP_SIZE, NTAPS );
	time_stop = vbx_timestamp();

	printf("...done\n");fflush(stdout);
	return vbx_print_scalar_time( time_start, time_stop );
}



int TEST_NAME()
{
	double scalar_time;
	int errors=0;

	vbx_mm_t *scalar_sample = malloc( (SAMP_SIZE+NTAPS)*sizeof(vbx_mm_t) );
	vbx_mm_t *scalar_coeffs = malloc(             NTAPS*sizeof(vbx_mm_t) );
	vbx_mm_t *scalar_out    = malloc(         SAMP_SIZE*sizeof(vbx_mm_t) );

	vbx_mm_t *sample     = vbx_shared_malloc( (SAMP_SIZE+NTAPS)*sizeof(vbx_mm_t) );
	vbx_mm_t *coeffs     = vbx_shared_malloc(             NTAPS*sizeof(vbx_mm_t) );
	vbx_mm_t *vector_out = vbx_shared_malloc(         SAMP_SIZE*sizeof(vbx_mm_t) );

	if((scalar_out == NULL) || (vector_out == NULL)){
		printf("\nMalloc failed\n");fflush(stdout);
		VBX_TEST_END(-1);
	}

	printf("\n\n**** " STRINGIFY(vbx_mm_t) " test ****\n");fflush(stdout);
	
	VBX_T(test_zero_array)( scalar_out, SAMP_SIZE );
	VBX_T(test_zero_array)( vector_out, SAMP_SIZE );

	VBX_T(test_init_array)( scalar_sample, SAMP_SIZE, 0xff );
	VBX_T(test_copy_array)( sample, scalar_sample, SAMP_SIZE );
	VBX_T(test_init_array)( scalar_coeffs, NTAPS, 1 );
	VBX_T(test_copy_array)( coeffs, scalar_coeffs, NTAPS );

	VBX_T(test_zero_array)( scalar_sample+SAMP_SIZE, NTAPS );
	VBX_T(test_zero_array)( sample+SAMP_SIZE, NTAPS );

	printf("\nSamples:\n");fflush(stdout);
	VBX_T(test_print_array)( scalar_sample, min(SAMP_SIZE,MAX_PRINT_LENGTH) );
	printf("\nCoefficients:\n");fflush(stdout);
	VBX_T(test_print_array)( scalar_coeffs, min(NTAPS,MAX_PRINT_LENGTH) );

	scalar_time = vbx_mm_t_test_scalar( scalar_out, scalar_sample, scalar_coeffs);
	VBX_T(test_print_array)( scalar_out,  min(SAMP_SIZE,MAX_PRINT_LENGTH) );


	#ifdef USE_TRANSPOSE
	vbx_mm_t_test_vector_transpose( vector_out, sample, coeffs, scalar_time );
	VBX_T(test_print_array)( vector_out,  min(SAMP_SIZE,MAX_PRINT_LENGTH) );
	errors += VBX_T(test_verify_array)( scalar_out, vector_out, SAMP_SIZE-NTAPS );
	#endif //USE_TRANSPOSE

	#ifdef USE_2D
	vbx_mm_t_test_vector( vector_out, sample, coeffs, scalar_time );
	VBX_T(test_print_array)( vector_out,  min(SAMP_SIZE,MAX_PRINT_LENGTH) );
	errors += VBX_T(test_verify_array)( scalar_out, vector_out, SAMP_SIZE-NTAPS );
	#endif //USE_2D

	free(scalar_sample);
	free(scalar_coeffs);
	free(scalar_out);

	vbx_shared_free(sample);
	vbx_shared_free(coeffs);
	vbx_shared_free(vector_out);

	return errors;
}
