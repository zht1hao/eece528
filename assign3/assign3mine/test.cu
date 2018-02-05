#include<stdio.h>

__global__
void Row(int size,int line1,float *array,float *array2){
	int i = blockIdx.x *blockDim.x + threadIdx.x;
	int line2,col,condition;
	float div;
	condition = (line1+1)*size;
	
	if (i + 1 > condition){
		if( i % size <= line1){
			*(array+i)=0;
		}
		else{
			line2 = i / size;
			col = i % size;
			div = *(array2 + size*line2 + line1) / *(array2 + size * line1 + line1);
			*(array+i) = *(array+i) - *(array + size * line1 + col) * div;
		}
	}
}

int main(int argc, char** argv){
	int size;
	float total_time;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	srand(time(NULL));
	if(argc != 2) {
		printf("pls use this execution with size of matrix\n");
		printf("for example:./GEN 4096\n");
		exit(0);
	}else{
		sscanf(argv[1],"%d",&size);
	}
	printf ("%d\n",size);
	float *d_array,*array,*d2_array;
	array = (float*)malloc(size*size*sizeof(float));
	
	cudaMalloc(&d_array, size*size*sizeof(float));
	cudaMalloc(&d2_array, size*size*sizeof(float));
//	initialize<<<(size+511)/512,512>>>(d_array);
//	cudaMemcpy(d_array, array, size*sizeof(float),cudaMemcpyDeviceToHost);
	for (int i = 0; i < size; i++){
		for (int j = 0; j < size; j++)
			*(array+i*size+j) = i+j+1;//(double)rand()/RAND_MAX*1000; 
	}
//	for (int i = 0; i < size; i++){
//		for (int j = 0; j < size; j++)
//			printf("%f ", *(array+i*size+j));
//		printf("\n");
//	}
	cudaMemcpy(d2_array,array,size*size*sizeof(float),cudaMemcpyHostToDevice);
	cudaMemcpy(d_array,array,size*size*sizeof(float),cudaMemcpyHostToDevice);
	cudaEventRecord(start);
	for (int i = 0; i < size-3; i++){
		Row<<<(size*size+511)/512,512>>>(size,i,d_array,d2_array);

	}
	cudaEventRecord(stop);
	cudaMemcpy(array,d_array, size*size*sizeof(float),cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();	
//	for (int i = 0; i < size; i++){
//		for (int j = 0; j < size; j++)
//			printf("%lf ", *(array+i*size+j));
//		printf("\n");
//	}
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&total_time,start,stop);
	total_time = (float) total_time / 1000;
	printf("total_time = %lf sec\n",total_time);
	cudaFree(d_array);
	cudaFree(d2_array);
	free(array);
}
