#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "vbx.h"
#include "vbx_port.h"

#define type int
#define N 4096


int main()
{

    #if VBX_SIMULATOR==1
	//initialize with 4 lanes,and 64kb of sp memory
	//word,half,byte fraction bits 16,15,4 respectively
	vbxsim_init( 4, 64, 256, 6, 5, 4);
    #endif

    type *A =(type*)malloc(N*N*sizeof(type));
    type *B =(type*)malloc(N*N*sizeof(type));
    type *C =(type*)malloc(N*N*sizeof(type));
    type *BT=(type*)malloc(N*N*sizeof(type));
    int i,j,k;
    int err;

    for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        {
            *(A+i*N+j)=rand()%1024;
            *(B+i*N+j)=rand()%1024;
            *(C+i*N+j)=0;
        }
    }

    for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        {
            *(BT+j*N+i) = *(B+i*N+j);
        }
    }


    vbx_dcache_flush_all();

    vbx_word_t*adata;
    vbx_word_t*bdata;
    vbx_word_t*cdata;

    adata=(vbx_word_t*) vbx_sp_malloc(N*sizeof(type));
    bdata=(vbx_word_t*) vbx_sp_malloc(N*sizeof(type));
    cdata=(vbx_word_t*) vbx_sp_malloc(N*sizeof(type));

    vbx_set_vl(N);

    for (i=0; i<N; i++)
    {
        vbx_dma_to_vector(bdata, BT+i*N, N*sizeof(type));
        for(j=0; j<N; j++)
        {
            vbx_dma_to_vector(adata, A+j*N, N*sizeof(type));
            vbx_acc(VVWS, VMUL, cdata+j, adata, bdata);
        }
        vbx_dma_to_host(C+i*N, cdata, N*sizeof(type));
        vbx_sync();
    }

    type *D=(type*)malloc(N*N*sizeof(type));
    type sum;
    for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        {
            sum=0;
            for(k=0; k<N; k++)
            {
                sum=sum+ *(A+i*N+k) * *(B+k*N+j);
            }
            *(D+i*N+j)=sum;
        }

    }

    for (j=0; j<N; j++)
    {
        for (i=0; i<N; i++)
        {
            if(*(D+j*N+i)!= *(C+i*N+j))
            {
                err=1;
            }
        }
    }

    vbxsim_print_stats_extended();

    if(err)
    {
        printf("Error!\n");
    }

    return 0;
}
