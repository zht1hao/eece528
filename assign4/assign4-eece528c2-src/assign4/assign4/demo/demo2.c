#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <vbx.h>
//#include <vbx_asm_or_sim.h>

#define DIM 1024
#define DATA_TYPE int
#define MAX_NUM 5

void Malloc1d(DATA_TYPE **a) {
        (*a) = malloc(DIM * DIM * sizeof(DATA_TYPE));
};

void Init1d(DATA_TYPE **a) {
        for(int i=0; i<DIM*DIM; ++i)
                (*a)[i]=(DATA_TYPE)rand()/(DATA_TYPE)(RAND_MAX/MAX_NUM);
};

void PrintMatrix(DATA_TYPE *a, int xdim) {
        printf("Matrix: \n");
        int i,j;
        for (i=0; i<xdim; i++) {
                for (j=0; j<DIM; j++) {
                        printf("%d ", a[i*DIM+j]);
                };
                printf("\n");
        };
};


void swap_pointers(vbx_word_t *p0, vbx_word_t *p1) {
        vbx_word_t *temp;
        temp = p0;
        p0 = p1;
        p1 = temp;
}

int main() {
#if VBX_SIMULATOR==1
        vbxsim_init(256, 64, 256,6,5, 4 );
#endif

        DATA_TYPE *a, *b, *c;
        int i,j;
        //TRANSPOSE matrix B
        Malloc1d(&a); Malloc1d(&b); Malloc1d(&c);
        Init1d(&a); Init1d(&b);
        //PrintMatrix(a, DIM);
        //PrintMatrix(b, DIM);

        vbx_word_t *row_a0 = vbx_sp_malloc(DIM*sizeof(vbx_word_t));
        //vbx_word_t *row_a1 = vbx_sp_malloc(DIM*sizeof(vbx_word_t));
        vbx_word_t *row_b0 = vbx_sp_malloc(DIM*sizeof(vbx_word_t));
        //vbx_word_t *row_b1 = vbx_sp_malloc(DIM*sizeof(vbx_word_t));
        vbx_word_t *row_c = vbx_sp_malloc(DIM*sizeof(vbx_word_t));

        //vbx_dma_to_vector(row_a0, a, DIM*4);
        //vbx_dma_to_vector(row_b0, b, DIM*sizeof(vbx_word_t));

        vbx_set_vl(DIM);
        for(i=0; i<DIM; i++) {
                vbx_dma_to_vector(row_a0, &a[i*DIM], DIM*sizeof(vbx_word_t));
                for(j=0; j<DIM; j++) {
                        vbx_dma_to_vector(row_b0, &b[j*DIM], DIM*sizeof(vbx_word_t));

                        vbx_acc(VVW, VMUL, &row_c[j], row_a0, row_b0);
                        //printf("MARIA calculated c[i][j] = %d\n", row_c[j]);          
                //swap_pointers(row_b0, row_b1);

                };
                //swap_pointers(row_a0, row_a1);                
                vbx_sync();
                vbx_dma_to_host(&c[i*DIM], (vbx_ubyte_t*)row_c, DIM*sizeof(vbx_word_t));

        };
        vbx_sync();
        //PrintMatrix(c, DIM);
        //struct simulator_statistics results = vbxsim_get_stats();
        vbxsim_print_stats_extended();
        return 0;
}
