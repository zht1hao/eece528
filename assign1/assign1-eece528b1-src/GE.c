#include<stdio.h>
#include<time.h>
#include<stdlib.h>

main(int argc, char *argv[]){
int i,j,k,size;
float *array,div,*array2;
clock_t begin,end;
double time_spent;
srand((unsigned)time(NULL));

if (argc!=2){
 	printf("you may entered nothing or thoo much numner, just enter one number pls.\n");
	exit(1);
}else{
	if(sscanf(argv[1],"%d",&size) == 0){
	printf("seems you got a wrong input type\n");
	exit(1);
	}

}

array = (float*)malloc(sizeof(float)*size*size);
for (i = 0; i<size*size; i++){
	*(array+i) = rand()/(float)(RAND_MAX)*4000000;
}

begin = clock();
for (i = 0;i < size;i++)
{
	for(j = i + 1;j < size;j++){
	
		div = *(array+j*size+i) / *(array+size*i+i);	
//		printf("%f\n",div);
	
		for(k = i; k < size;k++){
			*(array+j*size+k) = *(array+j*size+k) - div * *(array+i*size+k);
		}		
	}
}
end = clock();
time_spent = (double)(end - begin)/CLOCKS_PER_SEC;
printf("time spent: %lf \n",time_spent);
//for(i=0;i<size;i++){
//	for(j=0;j<size;j++){
//		printf("%f ",*(array+i*size+j));
//	}
//printf("\n");
}

