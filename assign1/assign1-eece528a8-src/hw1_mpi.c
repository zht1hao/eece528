#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <stdbool.h>
#include <sys/time.h>

//#define DIM 32
//#define DATA_TYPE float
#define MAX_NUM 10
#define DIV_ROUND_UP(a, b) (((a) / (b)) + (((a) % (b)) > 0 ? 1 : 0))

void Malloc2dint(int ***a, int dimx, int dimy) {
        *a = malloc(dimx * sizeof(int*));
        if (*a == NULL){
                printf("ERROR: out of memory\n");
        };
        for (int i=0; i<dimx; i++) {
                (*a)[i] = malloc(dimy * sizeof(int));
                if ((*a)[i] == NULL){
                        printf("ERROR: out of memory\n");
                };
        };
};

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

void swap_rows(DATA_TYPE ***a, int dest, int src) {
	int j;
	for (j=0; j<DIM; j++) {
		DATA_TYPE temp = (*a)[dest][j];
		(*a)[dest][j] = (*a)[src][j];
		(*a)[src][j]  = temp;
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

int first_row_with_non_zero_coeff(DATA_TYPE **a, int col) {
	int result = col;
	for (int i=col+1; i<DIM; i++) {
		if (a[i][col] != 0) {
			result = i;
			break;
		};	
	};
	return result;
};

int** calc_rows_distribution(int num_slaves) {
	int **result;
	Malloc2dint(&result, DIM, num_slaves);
	for (int j=0; j<DIM; j++) {
		int num_rows = DIM-(j+1);
		int res = num_rows % num_slaves;
		int *indices_to_add;
		if (res>0) {
			indices_to_add = malloc(res*sizeof(int));
		};
		for (int i=0; i<res; i++) {
			indices_to_add[i] = rand() % num_slaves;
		};
		for (int i=0; i<num_slaves; i++) {
			result[j][i] = num_rows/num_slaves;
		};
		for (int i=0; i<res; i++) {
			result[j][indices_to_add[i]] += 1;
		};
	};
#ifdef FULL_VERB
	printf("Generated rows distribution between slaves (sort of load balancer)\n");
	int i,j;
	for (i=0; i<DIM; i++) {
		for (j=0; j<num_slaves; j++) {
			printf("%d ", result[i][j]);
		};
		printf("\n");
	};
#endif
	return result;
};

MPI_Status status;

int ModifyRows(DATA_TYPE **a, int start_row, int end_row, int offset, int col, DATA_TYPE *b) {
	
	int i,k;	
	for (i=start_row; i<=end_row; i++) {
		if (a[i][col]==0) {
			continue;
		};
		if (b[col]==0) {
			printf("a[%d][%d]=0!!!", offset+i, col);
			return 1;
		};
		DATA_TYPE coeff = a[i][col]/b[col];
		a[i][col] = 0;
		for (k=col+1; k<DIM; k++) {
			a[i][k] -= b[k] * coeff;
		};
#ifdef HIGH_VERB
		printf("row%d = row%d - row%d * %f\n",offset+i, offset+i, col, a[i][col]/b[col]);
		for (k=0; k<DIM; k++) {
			printf("%f ", a[i][k]);
		};
		printf("\n");
#endif	
	};
	
	return 0;
};

void ModifyMatrix(DATA_TYPE **a, int offset, int step, int* pivot_owner, int num_tasks, int task_id) {
	int i,j;
	if (task_id>0) {
		for (j=0; j<offset; j++) {
			int col; //should be equal to j		
			DATA_TYPE b[DIM];
			//printf("task %d waiting to recv col=%d from %d\n", task_id, j, pivot_owner[j]);
			MPI_Recv(&col,    1, 		MPI_INT,   pivot_owner[j], 3, MPI_COMM_WORLD, &status);
			//printf("task %d received row %d from task %d\n", task_id, col, pivot_owner[j]);
			if (col!=j) {
				printf("wrong order!!!\n");
				return;
			};
			MPI_Recv(&b,      DIM, 		MPI_FLOAT, pivot_owner[j], 3, MPI_COMM_WORLD, &status); //row j
			if (ModifyRows(a, 0, step-1, offset, j, b)) {
				break;
			};
		};	
	};
	
	for (int s=task_id+1; s<num_tasks; s++) {	
		//printf("%d sending row %d to task %d\n", task_id, offset, s);
		MPI_Send(&offset, 	1, 		MPI_FLOAT, s, 3, MPI_COMM_WORLD);					
		MPI_Send(&a[0][0], 	DIM, 		MPI_FLOAT, s, 3, MPI_COMM_WORLD);
	};
	for (j=1; j<step; j++) {
		//printf("task=%d Modifying all rows using pivot from row %d\n", task_id, j-1); 
		int col_to_send = offset+j;
		if (ModifyRows(a, j, step-1, offset, offset+j-1, &a[j-1][0])) {
			break;
		};
		for (int s=task_id+1; s<num_tasks; s++) {
			//printf("task %d sending row %d to task %d\n", task_id, col_to_send, s);
			MPI_Send(&col_to_send, 		1, 		MPI_FLOAT, s, 3, MPI_COMM_WORLD);					
			MPI_Send(&a[j][0], 	DIM, 		MPI_FLOAT, s, 3, MPI_COMM_WORLD);
		};
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

int main(int argc, char *argv[]) {
	srand(time(NULL));
	int i, j, s, col, task_id, num_tasks, num_slaves, num_of_rows, step, offset;
	struct timeval start, stop;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &task_id);
	MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
	num_slaves = num_tasks - 1;
	//if (num_slaves == 0) {
	//	printf("No slaves, cant run mpi\n");
	//	return 1;
	//};
	DATA_TYPE **a;
	step = DIM/num_tasks;
	offset = step;
	int pivot_owner[DIM];
	for (i=0; i<num_tasks; i++) {
		for (j=0; j<step; j++) {
			pivot_owner[i*step+j] = i;
		};
	}; 
///////////////////////////////////////////MASTER//////////////////////////////////////
	if (task_id==0)  {//master 
		printf("Gaussian elimination calc. Num of processors: %d Step: %d \n", num_slaves, step);
		Malloc2d(&a, DIM, DIM);
		DATA_TYPE **a_ref;
		Malloc2d(&a_ref, DIM, DIM);
		Init2d(a, a_ref);
#ifdef PRINT	
		printf("Matrix before: \n");
		PrintMatrix(a, DIM);
#endif		
		for (s=1; s<num_tasks; s++) {
			
			MPI_Send(&offset, 	1, 		MPI_INT,   s, 1, MPI_COMM_WORLD);
			
			for (int k=0; k<step; k++) {
				MPI_Send(&a[offset+k][0], DIM, 	MPI_FLOAT, s, 1, MPI_COMM_WORLD);
			};	
			
			offset += step;
		};

		gettimeofday(&start, 0);

		ModifyMatrix(a, 0, step, pivot_owner, num_tasks, 0);

		offset = step;	
		for (s=1; s<num_tasks; s++) {
#ifdef HIGH_VERB 
			printf("waiting for results from slave %d column %d offset=%d size=%d \n", s, j, offset, step*DIM); 
#endif 
			for (int k=0; k<step; k++) {					
				MPI_Recv(&a[offset+k][0], DIM, 	MPI_FLOAT, s, 2, MPI_COMM_WORLD, &status);
			};
#ifdef HIGH_VERB
			printf("got results from slave %d column %d offset=%d size=%d \n", s, j, offset, step*DIM);
#endif 
#ifdef FULL_VERB
			printf ("matrix after stage %d \n", j);
			PrintMatrix(a, DIM);
#endif
			offset += step;
		};
		gettimeofday(&stop, 0);
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
	};
//////////////////////////////////////////SLAVES/////////////////////////////////////
	if (task_id>0) {
		//printf("task %d waiting for offset and data \n", task_id);	
		MPI_Recv(&offset,   1, 		MPI_INT,   0, 1, MPI_COMM_WORLD, &status);
		Malloc2d(&a, step, DIM);
		for (int k=0; k<step; k++) {
                	MPI_Recv(&a[k][0],      DIM, 	MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &status);
		};
#ifdef FULL_VERB 
		printf("task_id=%d received matrix to work on:\n", task_id);
		PrintMatrix(a, step);
#endif				
		ModifyMatrix(a, offset, step, pivot_owner, num_tasks, task_id);
#ifdef HIGH_VERB
		printf("Sending results to master. task_id=%d column=%d offset=%d size=%d\n", task_id, j, offset, step*DIM);			
#endif		
		for (int k=0; k<step; k++) {				
			MPI_Send(&a[k][0], DIM, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
		};
	};
	MPI_Finalize();
};

