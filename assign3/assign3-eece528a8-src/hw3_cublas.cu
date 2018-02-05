#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <sys/time.h>
#include <sys/types.h>
#include <math.h>
#include <cuda_runtime.h>
#include "cublas_v2.h"

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
	srand(time(NULL));
	DATA_TYPE *a = (DATA_TYPE*)malloc(DIM*DIM*sizeof(DATA_TYPE));
	DATA_TYPE *a_ref = (DATA_TYPE*)malloc(DIM*DIM*sizeof(DATA_TYPE));
	Init1d(a, a_ref);
#ifdef PRINT	
	printf("Matrix before: \n");
	PrintMatrix(a, DIM);
#endif	
	cudaError_t cudaStat ; // cudaMalloc status
	cublasStatus_t stat ; // CUBLAS functions status
	cublasHandle_t handle ; // CUBLAS context
	//DATA_TYPE *a_gpu;
	//cudaMalloc((void **)&a_gpu, DIM*DIM*sizeof(DATA_TYPE));
	//cudaMemcpy(a_gpu, a, DIM*DIM*sizeof(DATA_TYPE), cudaMemcpyHostToDevice);
	//cublasSetMatrix (DIM, DIM, sizeof(*a), a, DIM, a_gpu, DIM);
	DATA_TYPE *pivot_row_d;
	DATA_TYPE *curr_row_d;
	cudaMalloc((void **)&pivot_row_d, DIM*sizeof(DATA_TYPE));
	cudaMalloc((void **)&curr_row_d, DIM*DIM*sizeof(DATA_TYPE));
	stat = cublasCreate (&handle ); // initialize CUBLAS context
	int i,j,k;
	struct timeval start, stop;
	printf("starting to measure runtime\n");
	gettimeofday(&start, 0);
	for (j=1; j<DIM; j++) {
		int col = j-1;
		stat = cublasSetVector (DIM-col, sizeof(DATA_TYPE) ,&a[col*DIM+col] ,1 ,pivot_row_d ,1); // cp x- >d_x
		//pivot_row_d = &a_gpu[col*DIM+col];
		for (i=j; i<DIM; i++) {
			
			if (a[i*DIM+col]==0) {
				continue;
			};
			stat = cublasSetVector (DIM-col, sizeof(DATA_TYPE) ,&a[i*DIM+col] ,1 ,curr_row_d ,1); // cp x- >d_
			//curr_row_d = &a_gpu[i*DIM+col];
			DATA_TYPE coeff = -a[i*DIM+col]/a[col*DIM+col];
			stat=cublasSaxpy(handle, DIM-col, &coeff, pivot_row_d, 1, curr_row_d, 1);
			//cudaDeviceSynchronize();			
			stat = cublasGetVector(DIM-col, sizeof(DATA_TYPE), curr_row_d, 1, &a[i*DIM+col], 1); // cp d_y - >y
		};
	};
	//cudaDeviceSynchronize();
	gettimeofday(&stop, 0);
	printf("Time = %.6f\n", (stop.tv_sec+stop.tv_usec*1e-6)-(start.tv_sec+start.tv_usec*1e-6));
	printf("copying from gpu to cpu\n");
	//cudaMemcpy(a, a_gpu, DIM*DIM*sizeof(DATA_TYPE), cudaMemcpyDeviceToHost);
	//cudaDeviceSynchronize();
	//cublasGetMatrix (DIM, DIM, sizeof(*a), a_gpu, DIM, a, DIM);
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
	if ((int)a[rand_x*DIM+rand_y] != (int)a_ref[rand_x*DIM+rand_y]) {
		printf("a[%d][%d] (%f) != a_ref[%d][%d] (%f)", rand_x, rand_y, a[rand_x*DIM+rand_y], rand_x, rand_y, a_ref[rand_x*DIM+rand_y]);
		return 0;
	};
#endif
	return 0;
}


