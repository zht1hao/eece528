#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<cuda.h>
__global__ void add(int *a, int *b, int *c){
	*c = *a + *b;
}
int main(){
	int a,b,c;
	int *d_a,*d_b,*d_c;	
	int size = sizeof(int);
	
	cudaMalloc((void **)&d_a, size);
	cudaMalloc((void **)&d_b, size);
	cudaMalloc((void **)&d_c, size);
	
	a = 2;
	b = 7;
	cudaMemcpy(d_a, &a, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, &b, size, cudaMemcpyHostToDevice);

	add<<<1,1>>>(d_a, d_b, d_c);
	
	cudaMemcpy(&c,d_c,size,cudaMemcpyDeviceToHost);
	
	printf("c = %d\n",c);	
	printf("a+b = %d\n",a+b);
	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_c);
	return 0;
}