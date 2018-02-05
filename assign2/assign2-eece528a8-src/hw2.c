#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include <stdbool.h>
#include <sys/time.h>
#include <sys/types.h>

//#define DIM 32
//#define DATA_TYPE float
//#define NTHREADS 16
#define MAX_NUM 10

struct ThreadFunctionInput {
	DATA_TYPE **a;
	int offset;
	int step;
	int thread_id;
};

pthread_cond_t row_done[DIM]; //DIM
pthread_mutex_t row_lock[DIM];
int row_done_state[DIM];

void Malloc2d(DATA_TYPE ***a, int dimx, int dimy) {
        *a = malloc(dimx * sizeof(DATA_TYPE*));
        if (*a == NULL){
                printf("ERROR: out of memory\n");
        };
        for (int i=0; i<dimx; i++) {
                (*a)[i] = malloc(dimy * sizeof(DATA_TYPE));
                if ((*a)[i] == NULL){
                        printf("ERROR: out of memory\n");
                };
        };
};

void Init2d(DATA_TYPE **a, DATA_TYPE **a_ref) {
	int i,j;
    for(i=0; i<DIM; ++i) {
        for(j=0; j<DIM; ++j) {
            a[i][j] = (DATA_TYPE)rand()/(DATA_TYPE)(RAND_MAX/MAX_NUM);
            a_ref[i][j] = a[i][j];
        };
    };
};

void PrintMatrix(DATA_TYPE **a, int xdim) {
	int i,j;
	for (i=0; i<xdim; i++) {
		for (j=0; j<DIM; j++) {
			printf("%f ", a[i][j]);
		};
		printf("\n");
	};
};

int ModifyRows(DATA_TYPE **a, int start_row, int end_row, int col, DATA_TYPE *b) {
#ifdef HIGH_VERB	
	printf("ModifyRows %d-%d with row %d\n", start_row, end_row, col);
#endif	
	int i,k;	
	for (i=start_row; i<=end_row; i++) {
		if (a[i][col]==0) {
			continue;
		};
		if (b[col]==0) {
			printf("a[%d][%d]=0!!!", i, col);
			return 1;
		};
		DATA_TYPE coeff = a[i][col]/b[col];
		a[i][col] = 0;
		for (k=col+1; k<DIM; k++) {
			a[i][k] -= b[k] * coeff;
		};
#ifdef HIGH_VERB
		printf("row%d = row%d - row%d * %f\n",i, i, col, a[i][col]/b[col]);
		for (k=0; k<DIM; k++) {
			printf("%f ", a[i][k]);
		};
		printf("\n");
#endif	
	};
	
	return 0;
};

void* ModifyMatrix(void  *thread_input_ext) {
	int i,j;
	struct ThreadFunctionInput *thread_input = (struct ThreadFunctionInput*)thread_input_ext;
	int thread_id = (*thread_input).thread_id;
	int offset = (*thread_input).offset;
	int step = (*thread_input).step;
	DATA_TYPE **a = (*thread_input).a;
	for (j=0; j<offset; j++) {
#ifdef HIGH_VERB		
		printf("thread %d waiting for row %d\n", thread_id, j);
#endif		
		pthread_mutex_lock(&row_lock[j]);
		while (row_done_state[j]==0) {
#ifdef HIGH_VERB			
			printf("thread %d cond_wait on row_lock for row %d\n", thread_id, j);
#endif			
			pthread_cond_wait(&row_done[j], &row_lock[j]);

		}
#ifdef HIGH_VERB		
		printf("thread %d got row ready for row %d, unlocking mutex and starting updating the matrix\n", thread_id, j);
#endif		
		pthread_mutex_unlock(&row_lock[j]);
		if (ModifyRows(a, offset, offset+step-1, j, &a[j][0])) {
			break;
		};
	};	
#ifdef HIGH_VERB	
	printf("%d broadcasting row %d\n", thread_id, offset);
#endif	
	pthread_mutex_lock(&row_lock[offset]);
	row_done_state[offset]=1;
	pthread_cond_broadcast(&row_done[offset]);
	pthread_mutex_unlock(&row_lock[offset]);
	for (j=offset+1; j<offset+step; j++) {
#ifdef HIGH_VERB		
		printf("task=%d Modifying all rows using pivot from row %d\n", thread_id, j-1); 
#endif		
		if (ModifyRows(a, j, offset+step-1, j-1, &a[j-1][0])) {
			break;
		};
#ifdef HIGH_VERB		
		printf("%d broadcasting row %d\n", thread_id, j);
#endif		
		pthread_mutex_lock(&row_lock[j]);
		row_done_state[j]=1;
		pthread_cond_broadcast(&row_done[j]);
		pthread_mutex_unlock(&row_lock[j]);
	};
};

void sequential_gaussian(DATA_TYPE **a) {
	int i,j,k;
	for (j=1; j<DIM; j++) {

		int col = j-1;
		for (i=j; i<DIM; i++) {

			if (a[i][col]==0) {
				continue;
			};
			DATA_TYPE coeff = a[i][col]/a[col][col];
						
			for (k=col+1; k<DIM; k++) {
				a[i][k] -= a[col][k] * coeff;
			};
			a[i][col] = 0;
		};
	};
}

int main() {
	int i;
	srand(time(NULL));
	DATA_TYPE **a;
	Malloc2d(&a, DIM, DIM);
	DATA_TYPE **a_ref;
	Malloc2d(&a_ref, DIM, DIM);
	Init2d(a, a_ref);
	for (i=0; i<DIM; i++){
        pthread_mutex_init(&row_lock[i], NULL);
        pthread_cond_init(&row_done[i], NULL);
        row_done_state[i] = 0;
    }
#ifdef PRINT	
	printf("Matrix before: \n");
	PrintMatrix(a, DIM);
#endif	
	pthread_t threads[NTHREADS];
	struct ThreadFunctionInput thread_inputs[NTHREADS];
	struct timeval start, stop;
	int step = DIM/NTHREADS;
	gettimeofday(&start, 0);
	for (i=0; i<(NTHREADS); i++) {
		thread_inputs[i].a = a;
		thread_inputs[i].offset = step*i;
		thread_inputs[i].step = step;
		thread_inputs[i].thread_id = i;
		if (i>0) {
			pthread_create(&threads[i], NULL, ModifyMatrix, (void*)&thread_inputs[i]);
		}
	}
	ModifyMatrix((void*)&thread_inputs[0]);
	for(i=1; i<NTHREADS; i++) {
        pthread_join(threads[i], NULL);
	}
	gettimeofday(&stop, 0);
	for (i=0; i<DIM; i++){
        pthread_mutex_destroy(&row_lock[i]);
        pthread_cond_destroy(&row_done[i]);
    }
#ifdef PRINT
	printf("Matrix after: \n");
	PrintMatrix(a, DIM);
#endif
	printf("Time = %.6f\n", (stop.tv_sec+stop.tv_usec*1e-6)-(start.tv_sec+start.tv_usec*1e-6));		
	//verify
	sequential_gaussian(a_ref);
#ifdef PRINT
	printf("Matrix verify: \n");
	PrintMatrix(a_ref, DIM);
#endif
	int rand_x = (int)rand()%DIM;
	int rand_y = (int)rand()%DIM;
	if (a[rand_x][rand_y] != a_ref[rand_x][rand_y]) {
		printf("a[%d][%d] (%f) != a_ref[%d][%d] (%f)", rand_x, rand_y, a[rand_x][rand_y], rand_x, rand_y, a_ref[rand_x][rand_y]);
		return 0;
	};
	return 0;
}