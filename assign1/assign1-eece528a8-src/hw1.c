#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
//#define DIM 32
//#define DATA_TYPE float
#define MAX_NUM 10

void PrintMatrix(DATA_TYPE *a, int xdim) {
	int i,j;
	for (i=0; i<xdim; i++) {
		for (j=0; j<DIM; j++) {
			printf("%f ", a[i*DIM+j]);
		};
		printf("\n");
	};
};

void Malloc2d(DATA_TYPE **a, int dimx, int dimy) {
        *a = malloc(dimx * dimy * sizeof(DATA_TYPE));
};

void Init2d(DATA_TYPE *a) {
        int i,j;
        for(i=0; i<DIM; ++i)
                for(j=0; j<DIM; ++j)
                        a[i*DIM+j] = (DATA_TYPE)rand()/(DATA_TYPE)(RAND_MAX/MAX_NUM);
};

int ModifyRows(DATA_TYPE *a, int num_rows, int offset, int col, DATA_TYPE *b) {
	int i,k;	
	for (i=0; i<num_rows; i++) {
		if (a[i*DIM+col]==0) {
			continue;
		};
		if (b[col]==0) {
			printf("a[%d][%d]=0!!!", offset+i, col);
			return 1;
		};
		DATA_TYPE coeff = a[i*DIM+col]/b[col];
		a[i*DIM+col] = 0;
		for (k=col+1; k<DIM; k++) {
			a[i*DIM+k] -= b[k] * coeff;
		};
#ifdef HIGH_VERB
		printf("row%d = row%d - row%d * %f\n",offset+i, offset+i, col, coeff);
		for (k=0; k<DIM; k++) {
			printf("%f ", a[i*DIM+k]);
		};
		printf("\n");
#endif	
	};
	
	return 0;
};

void swap_rows(DATA_TYPE *a, int dest, int src) {
	int j;
	for (j=0; j<DIM; j++) {
		DATA_TYPE temp = a[dest*DIM+j];
		a[dest*DIM+j] = a[src*DIM+j];
		a[src*DIM+j]  = temp;
	};
};

int first_row_with_non_zero_coeff(DATA_TYPE *a, int col) {
	int result = col;
	for (int i=col+1; i<DIM; i++) {
		if (a[i*DIM+col] != 0) {
			result = i;
			break;
		};	
	};
	return result;
};

int main() {
	srand(time(NULL));
	int i,j,k;
	DATA_TYPE *a;
	Malloc2d(&a, DIM, DIM);
	Init2d(a);
#ifdef PRINT	
	printf("Matrix before: \n");
	PrintMatrix(a, DIM);
#endif
	struct timeval start, stop;
	gettimeofday(&start, 0);
	for (j=1; j<DIM; j++) {
		
		int col = j-1;
		if (ModifyRows(&a[j*DIM], DIM-j, j, col, &a[(col)*DIM])) {
			break;
		};
	};
	gettimeofday(&stop, 0);
	printf("Time = %.6f\n", (stop.tv_sec+stop.tv_usec*1e-6)-(start.tv_sec+start.tv_usec*1e-6));
#ifdef PRINT
	printf("Matrix after: \n");
	PrintMatrix(a, DIM);
#endif	
	return 0;
};


