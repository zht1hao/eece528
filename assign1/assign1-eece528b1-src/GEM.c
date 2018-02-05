#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>
#include<time.h>

 
#define send_data_tag 10000
#define return_data_tag 10001
#define check_normal 1
#define check_end 0
int Matrix_size;
main(int argc, char **argv){

MPI_Status status;
int ierr, my_id, num_procs,check_send ,i,j,k,prime_loop,inner_loop,column_loop;
float *array, div;
int size,counter,num_per_core_send,num_per_core,num_for_core,split_line;
clock_t bin,end;
double total_time;
ierr = MPI_Init(&argc, &argv);
ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
int n = 0;
sscanf(argv[1],"%d",&Matrix_size);


if(my_id == 0){

	size = Matrix_size;
	array = (float*)malloc(sizeof(float)*size*size);
	
	for (i = 0; i < size*size; i++)
		*(array+i) = rand()/(float)RAND_MAX*400000;

	bin = clock();
	for(prime_loop = 0; prime_loop < size; prime_loop++){
		for(inner_loop = prime_loop+1; inner_loop < size; inner_loop++){
			div = *(array + inner_loop * size + prime_loop)/ *(array + prime_loop * size+ prime_loop);
			num_per_core = (size - prime_loop) / num_procs ;
			split_line = size % num_procs;
			if(split_line == 0){
				counter = num_per_core;
			}else{
				counter = num_per_core+1;
			}
			for(i = prime_loop; i < counter; i++){
				*(array + inner_loop*size + i)=*(array + inner_loop*size + i) - *(array + prime_loop*size +i)*div;
			}
			for(i = 1; i < num_procs; i++){
				if (i < split_line){
					num_per_core_send = num_per_core+1;	
					}else{
					num_per_core_send = num_per_core;
				}
					check_send = check_normal;
					ierr = MPI_Send(&check_send,1,MPI_INT,i,send_data_tag,MPI_COMM_WORLD);
					ierr = MPI_Send(&div,1,MPI_FLOAT,i,send_data_tag,MPI_COMM_WORLD);
					ierr = MPI_Send(&num_per_core_send,1,MPI_INT,i,send_data_tag,MPI_COMM_WORLD);
					ierr = MPI_Send(array+prime_loop*size+counter,num_per_core_send,MPI_FLOAT,i,send_data_tag,MPI_COMM_WORLD);
					ierr = MPI_Send(array+inner_loop*size+counter,num_per_core_send,MPI_FLOAT,i,send_data_tag,MPI_COMM_WORLD);
					counter = counter + num_per_core_send;
					
			}
			if(split_line == 0){
				counter = num_per_core;
			}else{
				counter = num_per_core + 1;
			}
			for(i = 1; i < num_procs; i++){
				if(i < split_line){
					num_per_core_send = num_per_core + 1;
				}else{
					num_per_core_send = num_per_core;
				}
		
				ierr = MPI_Recv(array + inner_loop*size + counter, num_per_core_send, MPI_FLOAT, i, MPI_ANY_TAG,MPI_COMM_WORLD, &status);
				counter = counter + num_per_core_send;
			}	
		}

	}
	end = clock();
	total_time = (double)(end-bin)/(CLOCKS_PER_SEC);
	check_send = check_end;
	for(i=1;i<num_procs;i++)
		ierr = MPI_Send(&check_send, 1, MPI_INT, i, send_data_tag, MPI_COMM_WORLD);
	printf("total time: %lf\n",total_time);
}
else{
	float *line1,*line2,div_got;
	int check_got,i;
//      printf("%d",Matrix_size);
	line1 = (float*)malloc(sizeof(float)*Matrix_size);
	line2 = (float*)malloc(sizeof(float)*Matrix_size);
	while(1){
		ierr = MPI_Recv(&check_got,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		if(check_got == check_normal){
			ierr = MPI_Recv(&div_got,1,MPI_FLOAT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			ierr = MPI_Recv(&num_for_core,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			ierr = MPI_Recv(line1,num_for_core,MPI_FLOAT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			ierr = MPI_Recv(line2,num_for_core,MPI_FLOAT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			for(i = 0; i < num_for_core; i++){
				*(line2+i) = *(line2+i) - *(line1+i)*div_got;
			}
			ierr = MPI_Send(line2,num_for_core,MPI_INT,0,send_data_tag,MPI_COMM_WORLD);
		}else{
			break;
		}
	}	
}
ierr = MPI_Finalize();

}
