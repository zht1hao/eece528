#undef VBX_TEMPLATE_T
#define VBX_TEMPLATE_T VBX_WORDSIZE_DEF
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <vbx.h>
#define ms 4096

#include "vbx_port.h"
#include <assert.h>


// the scratchpad can only hold 3 lines of matrix B ,and treat each element in A as a single scalar
int *a;
int b[ms][ms];
int c[ms][ms];
int c1[ms][ms];
int diff[ms][ms];
int i,j,k;
int b1[ms];
int b2[ms];
int b0[ms];
int main()
{
    vbxsim_init( 4, 64, 256,6,5, 4 );
    a = (int*)malloc(ms*ms*sizeof(int));
   srand(time(0));
        for(i=0;i<ms;i++)
    {
        for(j=0;j<ms;j++)
        {
            a[i*ms+j]=rand()%10;
            b[i][j]=rand()%10;
            c[i][j]=0;//initialize the output array
            c1[i][j]=0;
            diff[i][j]=-1;
        }
    }

//printf("before calculation\n");

///////////////////////////////
/*
printf("A\n");
for (i=0; i<ms; i++){
printf("\n");
for (j=0; j<ms; j++)
{
printf("%d\t",a[i*ms+j]);
}
}

///////////////////////////////

printf("B\n");
for (i=0; i<ms; i++){
printf("\n");
for (j=0; j<ms; j++)
{
printf("%d\t",b[i][j]);
}
}
///////////////////////////////

printf("C\n");
for (i=0; i<ms; i++){
printf("\n");
for (j=0; j<ms; j++)
{
printf("%d\t",c[i][j]);
}
}

*/
vbx_dcache_flush_all();


//Allocate vectors in scratchpad

vbx_word_t* vb0=(vbx_word_t*)vbx_sp_malloc(sizeof(int)*ms);
    vbx_word_t* vb1=(vbx_word_t*)vbx_sp_malloc(sizeof(int)*ms);
    vbx_word_t* vb2=(vbx_word_t*)vbx_sp_malloc(sizeof(int)*ms);
    vbx_word_t* vc=(vbx_word_t*)vbx_sp_malloc(sizeof(int)*ms);

//measure time
     vbx_timestamp_t time_start, time_stop;

     vbx_timestamp_start();
     time_start = vbx_timestamp();
/*
for (i=0; i<ms; i++)
{
 b0[i]=b[0][i];
 b1[i]=b[1][i];
 b2[i]=b[2][i];
}
*/


for (k=0; k<ms; k++)
{
    vbx_dma_to_vector(vc,c[k],  ms*sizeof(int));

    for (i=0; i<ms; i+=3)
    {
    if(i<ms-2){
   vbx_dma_to_vector(vb0,b[i],  ms*sizeof(int));
   vbx_dma_to_vector(vb1,b[i+1],ms*sizeof(int));
   vbx_dma_to_vector(vb2,b[i+2],ms*sizeof(int));
    vbx_set_vl(ms);
    vbx(SVW, VMULLO,vb0,a[i+k*ms],vb0);
    vbx(SVW, VMULLO,vb1,a[i+1+k*ms],vb1);
    vbx(SVW, VMULLO,vb2,a[i+2+k*ms],vb2);

    vbx( VVW, VADD, vc, vb0, vc );
    vbx( VVW, VADD, vc, vb1, vc );
   vbx( VVW, VADD, vc, vb2, vc );
  }
else if ((i==ms-2)&&(ms%3==2)){
   vbx_dma_to_vector(vb0,b[i],  ms*sizeof(int));
   vbx_dma_to_vector(vb1,b[i+1],ms*sizeof(int));
   vbx_set_vl(ms);
    vbx(SVW, VMULLO,vb0,a[i+k*ms],vb0);
    vbx(SVW, VMULLO,vb1,a[i+1+k*ms],vb1);

    vbx( VVW, VADD, vc, vb0, vc );
    vbx( VVW, VADD, vc, vb1, vc );

    }
else if ((i==ms-1)&&(ms%3==1))
{

vbx_dma_to_vector(vb0,b[i],  ms*sizeof(int));

    vbx_set_vl(ms);
    vbx(SVW, VMULLO,vb0,a[i+k*ms],vb0);

    vbx( VVW, VADD, vc, vb0, vc );
}



}
  vbx_dma_to_host( c[k], vc, ms*sizeof(int));

}
vbx_sync();
vbxsim_print_stats_extended();

// printf("after calculation\n");
/*  for(i= 0; i<ms; i++)
 {
  printf("%d\t",c[0][i]);
  printf("%d\t",b[0][i] );
 }
*/

////////////////////////////////
/*
printf("A\n");
for (i=0; i<ms; i++){
printf("\n");
for (j=0; j<ms; j++)
{
printf("%d\t",a[i*ms+j]);
}
}

///////////////////////////////

printf("\nB\n");
for (i=0; i<ms; i++){
printf("\n");
for (j=0; j<ms; j++)
{
printf("%d\t",b[i][j]);
}
}
///////////////////////////////

printf("\nC\n");
for (i=0; i<ms; i++){
printf("\n");
for (j=0; j<ms; j++)
{
printf("%d\t",c[i][j]);
}
}
*/

    vbx_sp_free();
   //return VBW_SUCCESS;
    time_stop = vbx_timestamp();
double runtime = time_stop-time_start;
printf("\nruntime=%f\n", runtime);

/*---------------- verification of the matrix multiplication result----------*/
/*
int p,q,d;

for (p=0; p< ms; p++){
     for (q = 0; q < ms; q++){
            for(d =0; d < ms; d++){
        c1[p][d]+=a[p*ms+q]*b[q][d];
      }

  }
}


for (p=0; p< ms; p++){
     for (q = 0; q < ms; q++){
            diff[p][q]= c[p][q] - c1[p][q];
           printf("%d\t",diff[p][q]);
    }
    printf("\n");
 }

*/
}
