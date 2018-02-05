#include<stdio.h>
#include<stdlib.h>
#include<pthread.h>
#include<time.h>
pthread_mutex_t fl;
pthread_cond_t conditionf;

float *array;
int row_proced=0;
int size,num_procs;
void *calculate(void *ptr){
	int *data_get;
	int my_id,i,j,k;
	float div;
	my_id = (void*)ptr;
	for(i = 0; i < size-1; i++){
//		printf("id:%d i = %d\n",my_id,i);
		for(j = i+1; j < size; j++){
			if(my_id == j %(num_procs-1)){
//				printf("my_id = %d,j = %d\n",my_id,j);
				div = *(array+j*size+i)/ *(array+i*size+i);
				for(k = i; k < size; k++){
					if(k == i)
						*(array+j*size+k) = 0;	
					else
						*(array+j*size+k) -= *(array+i*size+k)*div;
					}
//				printf("id:%d  passed k, i = %d, j = %d\n",my_id,i,j);
				if(j == i+1){
					pthread_mutex_lock(&fl);
					row_proced= row_proced +1;
//					printf("id:%d added row_proced, now it's %d\n",my_id, row_proced);
					pthread_cond_broadcast(&conditionf);
					pthread_mutex_unlock(&fl);
				}
				
			}
		
		}
	//	printf("id:%d before if row_proced = %d\n",my_id, row_proced);
		if(row_proced<=i){
//			printf("id:%d inside if row proced = %d i = %d \n",my_id,row_proced, i);
			pthread_mutex_lock(&fl);
			while(row_proced <= i){
//				printf("id:%d inside while\n",my_id);
				pthread_cond_wait(&conditionf,&fl);
			}
			pthread_mutex_unlock(&fl);
		}
	}
}
void initialize_array(){
	int i,j;
	srand(time(NULL));
	for (i=0; i<size; i++)
		for(j=0; j<size; j++)
			*(array+i*size+j)=(float)rand()/RAND_MAX*1000;
}
void print_array(){
	int i,j;
	for(i=0; i<size; i++){
		for(j=0; j<size; j++)
			printf("%f ",*(array+i*size+j));
	printf("\n");
	}
}
void main(int argc, char **argv){
	pthread_t thread[32];
	struct timespec begin,end;
	clock_t begin2,end2;
	double time_total,time_total2;
	int i,j,k;
	int iret;
	float *div;
	div = (float*)malloc(sizeof(float));
	if (argc != 3){
		printf("enter a matrix size and threads number  pls\n like: ./GEp 4096 32\n");
		exit(0);	
	}else{
		sscanf(argv[1],"%d",&size);
		sscanf(argv[2],"%d",&num_procs);
	}
	if(num_procs > 32){
		printf("Max threads number is 32, pls enter a number less than 32\n");
		exit(0);	
	}
//	printf("%d \n%d \n",size,num_procs);
	array = (float*)malloc(size*size*sizeof(float));
	initialize_array();
//	print_array();
	begin2 = clock();
	clock_gettime(CLOCK_MONOTONIC,&begin);
	for (i = 0; i <num_procs; i++){
		iret = pthread_create(&thread[i],NULL,calculate,(void*)i);	
	}
	for(i = 0; i < num_procs; i++){
		pthread_join(thread[i],NULL);
	}
//	print_array();	
	clock_gettime(CLOCK_MONOTONIC,&end);
	end2 = clock();
	time_total = (end.tv_sec - begin.tv_sec);
	time_total += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;
	time_total2 = (double)(end2 - begin2)/(CLOCKS_PER_SEC * num_procs);
	printf("time by CLOCK_MONOTONIC = %lf\n time by clock() = %lf\n",time_total,time_total2);
	exit(0);
}


