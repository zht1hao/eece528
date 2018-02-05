#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <sys/time.h>
#include <sys/types.h>
#include <math.h>

//#define DIM 32
//#define DATA_TYPE float
//#define NTHREADS 16 //max 1024
#define MAX_THREADS_PER_BLOCK 1024
#define MAX_NUM 10

void Malloc2d(DATA_TYPE ***a, int dimx, int dimy) {
        *a = (DATA_TYPE**)malloc(dimx * sizeof(DATA_TYPE*));
        if (*a == NULL){
                printf("ERROR: out of memory\n");
        };
        for (int i=0; i<dimx; i++) {
                (*a)[i] = (DATA_TYPE*)malloc(dimy * sizeof(DATA_TYPE));
                if ((*a)[i] == NULL){
                        printf("ERROR: out of memory\n");
                };
        };
};

void Init1d(DATA_TYPE *a, DATA_TYPE *a_ref) {
        int i,j;
        for(i=0; i<DIM; ++i)
         	 for(j=0; j<DIM; ++j) {
                a[i*DIM+j] = (DATA_TYPE)rand()/(DATA_TYPE)(RAND_MAX/MAX_NUM);
                a_ref[i*DIM+j] = a[i*DIM+j];
            };
};

void PrintMatrix(DATA_TYPE *a, int xdim) {
	int i,j;
	for (i=0; i<xdim; i++) {
		for (j=0; j<DIM; j++) {
			printf("%f ", a[i*DIM+j]);
		};
		printf("\n");
	};
};

__device__  int ModifyRows(DATA_TYPE *a, int start_row, int end_row, int col, DATA_TYPE *b) {
#ifdef HIGH_VERB	
	printf("ModifyRows %d-%d with row %d\n", start_row, end_row, col);
#endif	
	int i,k;
	for (i=start_row; i<=end_row; i++) {
		if (a[i*DIM+col]==0) {
			continue;
		};
		if (b[col]==0) {
			printf("a[%d][%d]=0!!!", i, col);
			return 1;
		};
		__shared__ DATA_TYPE coeff;
		coeff = a[i*DIM+col]/b[col];
		for (k=col+1; k<DIM; k++) {
			a[i*DIM+k] -= b[k] * coeff;
		};
		a[i*DIM+col]=0;
#ifdef HIGH_VERB
		printf("row%d = row%d - row%d * %f\n",i, i, col, coeff);
		for (k=0; k<DIM; k++) {
			printf("%f ", a[i*DIM+k]);
		};
		printf("\n");
#endif	
	};
	
	return 0;
};

__global__ void ModifyMatrix(DATA_TYPE *a, int step) {
    int offset = MAX_THREADS_PER_BLOCK*blockIdx.x + step*threadIdx.x;
    int j,i;
    //__shared__ DATA_TYPE pivot_row[DIM];
	//printf("ModifyMatrix called offset=%d bx=%d tx=%d step=%d\n", offset, blockIdx.x, threadIdx.x, step);
	for (j=0; j<DIM; j++) {
		//for (i=offset; i<offset+step; i++) {
		//	pivot_row[i]=a[j*DIM+i];
		//}
		__syncthreads();
		if (j<(offset+step)) {
#ifdef HIGH_VERB		
			printf("task=%d Modifying all rows using pivot from row %d\n", offset, j); 
#endif	
			int start_row = j>=offset? j+1 : offset;
			if (ModifyRows(a, start_row, offset+step-1, j, &a[j*DIM])) {
				break;
			}	
		};	
		__syncthreads();
	};
	for (j=offset+1; j<offset+step; j++) {
#ifdef HIGH_VERB		
		printf("task=%d Modifying all rows using pivot from row %d\n", thread_id, j-1); 
#endif		
		if (ModifyRows(a, j, offset+step-1, j-1, &a[(j-1)*DIM+0])) {
			break;
		};
	__syncthreads();
	};
}

void sequential_gaussian(DATA_TYPE *a) {
	int i,j,k;
	for (j=1; j<DIM; j++) {

		int col = j-1;
		for (i=j; i<DIM; i++) {

			if (a[i*DIM+col]==0) {
				continue;
			};
			DATA_TYPE coeff = a[i*DIM+col]/a[col*DIM+col];
						
			for (k=col+1; k<DIM; k++) {
				a[i*DIM+k] -= a[col*DIM+k] * coeff;
			};
			a[i*DIM+col] = 0;
		};
	};
}

int main() {
	//int i;
	srand(time(NULL));
	DATA_TYPE *a = (DATA_TYPE*)malloc(DIM*DIM*sizeof(DATA_TYPE));
	DATA_TYPE *a_ref = (DATA_TYPE*)malloc(DIM*DIM*sizeof(DATA_TYPE));
	Init1d(a, a_ref);
#ifdef PRINT	
	printf("Matrix before: \n");
	PrintMatrix(a, DIM);
#endif	
	DATA_TYPE *a_gpu;
	cudaMalloc((void **)&a_gpu, DIM*DIM*sizeof(DATA_TYPE));
	cudaMemcpy(a_gpu, a, DIM*DIM*sizeof(DATA_TYPE), cudaMemcpyHostToDevice);
	struct timeval start, stop;
	dim3 dimGrid (1, 1, 1);
  	dim3 dimBlock(NTHREADS, 1, 1);
	printf("starting to measure runtime\n");
	gettimeofday(&start, 0);
	ModifyMatrix<<<dimGrid,dimBlock>>> (a_gpu, DIM/NTHREADS);
	cudaDeviceSynchronize();
	gettimeofday(&stop, 0);
	printf("Time = %.6f\n", (stop.tv_sec+stop.tv_usec*1e-6)-(start.tv_sec+start.tv_usec*1e-6));
	printf("copying from gpu to cpu\n");
	cudaMemcpy(a, a_gpu, DIM*DIM*sizeof(DATA_TYPE), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
	printf("freeing the gpu mem\n");
	cudaFree(a_gpu);
#ifdef PRINT
	printf("Matrix after: \n");
	PrintMatrix(a, DIM);
#endif
#ifdef VERIFY
	sequential_gaussian(a_ref);
#ifdef PRINT
	printf("Matrix verify: \n");
	PrintMatrix(a_ref, DIM);
#endif
	int rand_x = (int)rand()%DIM;
	int rand_y = (int)rand()%DIM;
	if (a[rand_x*DIM+rand_y] != a_ref[rand_x*DIM+rand_y]) {
		printf("a[%d][%d] (%f) != a_ref[%d][%d] (%f)", rand_x, rand_y, a[rand_x*DIM+rand_y], rand_x, rand_y, a_ref[rand_x*DIM+rand_y]);
		return 0;
	};
#endif
	return 0;
}


