#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <vbx.h>
//#include <vbx_asm_or_sim.h>

//#define DIM 4096
#define DATA_TYPE int
#define MAX_NUM 5

void Malloc1d(DATA_TYPE **a);
void Init1d(DATA_TYPE **a);
void PrintMatrix(DATA_TYPE *a, int xdim);
void swap_pointers(vbx_word_t **p0, vbx_word_t **p1);
void MatrixMul(DATA_TYPE *a, DATA_TYPE *b, DATA_TYPE *c);

int main() {
#if VBX_SIMULATOR==1
        //initialize with 4 lanes,and 64kb of sp memory
        //word,half,byte fraction bits 16,15,4 respectively
        vbxsim_init(256, 64, 256,6,5, 4 );
#endif

	DATA_TYPE *a, *b, *c, *cref; 
	int i,j;
	//TRANSPOSE matrix B
	Malloc1d(&a); Malloc1d(&b); Malloc1d(&c); Malloc1d(&cref);	
	Init1d(&a); Init1d(&b);
#ifdef PRINT	
	PrintMatrix(a, DIM);
	PrintMatrix(b, DIM);
#endif	
	
	vbx_word_t *row_a0 = vbx_sp_malloc(DIM*sizeof(vbx_word_t));
	//vbx_word_t *row_a1 = vbx_sp_malloc(DIM*sizeof(vbx_word_t));
	vbx_word_t *row_b0 = vbx_sp_malloc(DIM*sizeof(vbx_word_t));
	vbx_word_t *row_b1 = vbx_sp_malloc(DIM*sizeof(vbx_word_t));
	vbx_word_t *row_c = vbx_sp_malloc(DIM*sizeof(vbx_word_t));
	
	//vbx_dma_to_vector(row_a0, a, DIM*sizeof(vbx_word_t));
	vbx_dma_to_vector(row_b0, b, DIM*sizeof(vbx_word_t));
	
	vbx_set_vl(DIM);
	for(i=0; i<DIM; i++) {
		//if (i<(DIM-1)) {
		//	vbx_dma_to_vector(row_a1, &a[(i+1)*DIM], DIM*sizeof(vbx_word_t));
		//}
		vbx_dma_to_vector(row_a0, &a[i*DIM], DIM*sizeof(vbx_word_t));
		for(j=0; j<DIM; j++) {
			if (j<(DIM-1)) {
				vbx_dma_to_vector(row_b1, &b[(j+1)*DIM], DIM*sizeof(vbx_word_t));
			}			
			vbx_acc(VVW, VMUL, &row_c[j], row_a0, row_b0);
			//printf("MARIA calculated c[i][j] = %d\n", row_c[j]);
			//vbx_sync();		
			swap_pointers(&row_b0, &row_b1);			
			//vbx_sync();
		}; 
		//vbx_sync();		
		vbx_dma_to_host(&c[i*DIM], (vbx_ubyte_t*)row_c, DIM*sizeof(vbx_word_t));	
		//swap_pointers(&row_a0, &row_a1);		
	};
	vbx_sync();
#ifdef PRINT	
	PrintMatrix(c, DIM);
#endif	
	vbxsim_print_stats_extended();
	//verify
	MatrixMul(a,b,cref);
#ifdef PRINT
	printf("Matrix verify: \n");
	PrintMatrix(cref, DIM);
#endif
	int rand_x = (int)rand()%DIM;
	int rand_y = (int)rand()%DIM;
	if (c[rand_x*DIM+rand_y] != cref[rand_x*DIM+rand_y]) {
		printf("c[%d][%d] (%d) != cref[%d][%d] (%d)", rand_x, rand_y, c[rand_x*DIM+rand_y], rand_x, rand_y, cref[rand_x*DIM+rand_y]);
		return 0;
	};	
	return 0;	
}; 

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

void swap_pointers(vbx_word_t **p0, vbx_word_t **p1) {
	vbx_word_t *temp;
	temp = (*p0);
	(*p0) = (*p1);
	(*p1) = temp;
}

void MatrixMul(DATA_TYPE *a, DATA_TYPE *b, DATA_TYPE *c) {
	int i,j,k;
	for (i=0; i<DIM; i++) {
		for (j=0; j<DIM; j++) {
			c[i*DIM+j] = 0;
			for (k=0; k<DIM; k++) {
				c[i*DIM+j] += a[i*DIM+k] * b[j*DIM+k];
			};
		};
	};
}
