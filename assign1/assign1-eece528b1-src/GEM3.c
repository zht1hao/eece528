#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>
#include<time.h>

 
#define send_data_tag 10000
#define return_data_tag 10001
#define check_normal 1
#define check_end 0
main(int argc, char **argv){

int ierr, my_id, num_procs,i,j,k;
float *array,div;
int Matrix_size,*spread;
clock_t begin,end;
double total_time;
ierr = MPI_Init(&argc, &argv);
ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

sscanf(argv[1],"%d",&Matrix_size);

MPI_Bcast(&Matrix_size,1,MPI_INT,0,MPI_COMM_WORLD);

array = (float*) malloc(Matrix_size * Matrix_size * sizeof(float));

spread = (int*)malloc(Matrix_size*sizeof(int));

if (my_id == 0){
	for(i=0; i<Matrix_size; i++){
		for(j=0; j<Matrix_size; j++){
			*(array + i*Matrix_size + j) = rand()/(float)(RAND_MAX)*1000;
		}
	}
	begin = clock();
printf("begin\n"); 

}
	MPI_Bcast(array,Matrix_size*Matrix_size,MPI_FLOAT,0,MPI_COMM_WORLD);

	for(i=0; i<Matrix_size;i++)
		*(spread +i) = i % num_procs; 
	MPI_Bcast(spread, Matrix_size, MPI_INT, 0, MPI_COMM_WORLD);

for(i=0; i<Matrix_size;i++){
	MPI_Bcast((array + i*Matrix_size + i),Matrix_size - i,MPI_FLOAT,*(spread+i),MPI_COMM_WORLD);
	for(j=i + 1; j<Matrix_size; j++){
		if(*(spread+j) == my_id){
			div = *(array + j*Matrix_size + i) / *(array + i*Matrix_size + i);
			for(k=i; k<Matrix_size; k++){
				*(array + j*Matrix_size + k) = *(array + j*Matrix_size + k) - div*(*(array + i*Matrix_size + k)); 
			}
		}
	}
	if(my_id == 0){
	for(j=0; j<i; j++)
		*(array + i*Matrix_size +j) = 0;
	}
}
	if(my_id == 0){
	printf("end\n");
		end = clock();
/*		for(j=0; j<Matrix_size; j++){
			for(i=0; i<Matrix_size; i++){
				printf("%*.*f",15,2,*(array+j*Matrix_size+i));
			}
			printf("\n");
		}*/
		total_time =(double)(end - begin)/CLOCKS_PER_SEC;
		printf("time = %lf\n", total_time);	
	}	

ierr = MPI_Finalize();
}
