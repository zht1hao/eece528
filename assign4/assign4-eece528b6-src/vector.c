#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <vbx.h>
#include <vbx_port.h>

#define n 4096


int main()
{
	//initialize with 4 lanes,and 64kb of sp memory
    //word,half,byte fraction bits 16,15,4 respectively
    vbxsim_init( 4, 64, 256, 16, 15, 4);
    
    int i,j,k;
    int *A =(int*)malloc(n*n*sizeof(int));
    int *B =(int*)malloc(n*n*sizeof(int));
    int *C =(int*)malloc(n*n*sizeof(int));
    int *D=(int*)malloc(n*n*sizeof(int));
    int *B_trans=(int*)malloc(n*n*sizeof(int));
    _Bool validate=0;

    srand(0);
    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            *(A+i*n+j)=rand()%4096;
            *(B+i*n+j)=rand()%4096;
            *(C+i*n+j)=0;
            *(D+i*n+j)=0;
        }
    }

    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            *(B_trans+j*n+i) = *(B+i*n+j);
        }
    }

    

    vbx_dcache_flush_all();

    vbx_word_t *A_vec=(vbx_word_t*) vbx_sp_malloc(n*sizeof(int));
    vbx_word_t *B_vec=(vbx_word_t*) vbx_sp_malloc(n*sizeof(int));
    vbx_word_t *C_vec=(vbx_word_t*) vbx_sp_malloc(n*sizeof(int));

    vbx_set_vl(n);

    for (i=0; i<n; i++)
    {
	vbx_dma_to_vector(A_vec, A+i*n, n*sizeof(int));

        for(j=0; j<n; j++)
        {
	        vbx_dma_to_vector(B_vec, B_trans+j*n, n*sizeof(int));
            vbx_acc(VVWS, VMUL, C_vec+j, A_vec, B_vec);
        }
        vbx_dma_to_host(C+i*n, C_vec, n*sizeof(int));
        vbx_sync();
    }


    
    if(validate)
    {
    	//test the result
	    for(i=0; i<n; i++)
    	{
        	for(k=0; k<n; k++)
	        {
            	for(j=0; j<n; j++)
            	{
                	*(D+i*n+j)+=*(A+i*n+k) * *(B+k*n+j);
            	}
        	}

    	}

    	for (i=0; i<n; i++)
    	{
        	for (j=0; j<n; j++)
        	{
            	if(*(D+i*n+j) != *(C+i*n+j))
            	{
                	printf("Error!\n");
                	break;
            	}
        	}
    	}
    }


    vbxsim_print_stats_extended();
    return 0;
}
