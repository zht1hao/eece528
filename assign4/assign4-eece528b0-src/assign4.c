#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <assert.h>
#include <vbx.h>
#include "vbx_port.h"
#undef VBX_TEMPLATE_T
#define VBX_TEMPLATE_T VBX_WORDSIZE_DEF
#define ms 4096


// the scratchpad can only hold 3 lines of matrix B ,and treat each element in A as a single scalar
int *a;
int b[ms][ms];
int c[ms][ms];
int i,j,k;
int b1[ms];
int b2[ms];
int b0[ms];

void print(int *a,int b[ms][ms],int c[ms][ms],int n);

int main()
{
    struct timeval start;//start time
    struct timeval finish;//finish time
    gettimeofday(&start,NULL);

    vbxsim_init(4,64,256,6,5,4);
    a = (int*)malloc(ms*ms*sizeof(int));
    //b = (int*)malloc(ms*ms*sizeof(int));
    //c = (int*)malloc(ms*ms*sizeof(int));
    srand(time(0));
    for(i=0;i<ms;i++)
    {
        for(j=0;j<ms;j++)
        {
            a[i*ms+j]=rand();
            b[i][j]=rand();
            c[i][j]=0;//initialize the output array
        }
    }

    printf("before calculation\n");
    print(a,b,c,4);

    vbx_dcache_flush_all();

    //Allocate vectors in scratchpad
    vbx_word_t* vb0=(vbx_word_t*)vbx_sp_malloc(sizeof(int)*ms);
    vbx_word_t* vb1=(vbx_word_t*)vbx_sp_malloc(sizeof(int)*ms);
    vbx_word_t* vb2=(vbx_word_t*)vbx_sp_malloc(sizeof(int)*ms);
    vbx_word_t* vc=(vbx_word_t*)vbx_sp_malloc(sizeof(int)*ms);

    vbx_set_vl(ms);
    for (k=0; k<ms; k++)
    {
        vbx_dma_to_vector(vc,c[k],ms*sizeof(int));

        for (i=0; i<ms; i+=3)
        {
            vbx_dma_to_vector(vb0,b[i],  ms*sizeof(int));
            vbx(SVW, VMULLO,vb0,a[k*ms+i],vb0);
            vbx(VVW,VADD,vc,vb0,vc);
            if(i+1<ms)
            {
                vbx_dma_to_vector(vb1,b[i+1],ms*sizeof(int));
                vbx(SVW, VMULLO,vb1,a[k*ms+i+1],vb1);
                vbx(VVW,VADD,vc,vb1,vc);
            }
            if(i+2<ms)
            {
                vbx_dma_to_vector(vb2,b[i+2],ms*sizeof(int));
                vbx(SVW, VMULLO,vb2,a[k*ms+i+2],vb2);
                vbx(VVW,VADD,vc,vb2,vc);
            }
        }
        vbx_dma_to_host(c[k],vc,ms*sizeof(int));
    }

    vbx_sync();

    vbxsim_print_stats_extended();

    printf("\nafter calculation\n");
    print(a,b,c,4);

    vbx_sp_free();
    //vbxsim_destroy();

    gettimeofday(&finish,NULL);
    double dtime=(double)(finish.tv_sec-start.tv_sec)+(double)(finish.tv_usec-start.tv_usec)/1000000;
    printf("*****\ngtime=%g\n*****\n",dtime);
}

void print(int *a,int b[ms][ms],int c[ms][ms],int n)
{
    printf("A\n");
    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
        {
            printf("%d\t",a[i*ms+j]);
        }
        printf("\n");
    }

    printf("\nB\n");
    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
        {
        printf("%d\t",b[i][j]);
        }
        printf("\n");
    }
    printf("\nC\n");
    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
        {
        printf("%d\t",c[i][j]);
        }
        printf("\n");
    }
}
