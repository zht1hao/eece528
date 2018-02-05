#include <stdio.h>
#include "vbx.h"

int main(){

const int length = 8;
int A[length] = {1,2,3,4,5,6,7,8};
int B[length] = {10,20,30,40,50,60,70,80};
int C[length] = {100,200,300,400,500,600,700,800};
int D[length];
vbx_dcache_flush_all();

const int data_len= length * sizeof(int);
vbx_word_t*va= (vbx_word_t*)vbx_sp_malloc( data_len);
vbx_word_t*vb= (vbx_word_t*)vbx_sp_malloc( data_len);
vbx_word_t*vc= (vbx_word_t*)vbx_sp_malloc( data_len);

vbx_dma_to_vector( va, A, data_len);
vbx_dma_to_vector( vb, B, data_len);
vbx_dma_to_vector( vc, C, data_len);

vbx_set_vl( length );
vbx( VVW, VADD, vb, va, vb);
vbx( VVW, VADD, vc, vb, vc);
vbx_dma_to_host( D, vc, data_len);

vbx_sync();
vbx_sp_free();


}
