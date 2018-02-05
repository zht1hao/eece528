// EECE528 assignment0 attempt code
// last modify @ Sep 22ed 5:35 AM 
// Haotian Zhang @ Desktop 
// Double version
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#define Data_type double
//#define Allocating_ways 1;// 0 for int a[m*n]; 1 for int a[m][n] ;2 for int *a; 3 for int **a;
//calculate the runtime for the multiply
double Calculate_Run_Time(clock_t* beginP, clock_t* endP) {
	double time;
	time = (double)(*endP - *beginP) / CLOCKS_PER_SEC;
	return time;
}
//get a random float from 0 to 300000000 
//float range 2.22E-308 to 1.79E+308 
Data_type Randome_double() {
	return (1.12E+38) * ((Data_type)rand() / RAND_MAX);
}
//Initialize matrix 
//randomSwitch 0 for all 0, 1 for random
int Initialize_Matrix(Data_type* aP, int size, int randomSwitch) { 
	int i, j;
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++) {
			if (!randomSwitch)
				*(aP + i*size + j) = 0;
			else
				if(rand()%2)								 //negative numbers covered
					*(aP + i*size + j) = -Randome_double();
				else
					*(aP + i*size + j) = Randome_double();
			 //printf("%d ",*(aP + i*size + j));
		}
		//printf("\n");
	}
	return 0;
}
//multiply matrix
//a b for input, c for output
void Multiply_matrix(Data_type *aP, Data_type *bP, Data_type *cP, int size) {
	int i, j, k;
	Data_type *tmp;
	tmp = (Data_type*)malloc(size*size * sizeof(Data_type));
	for (i = 1; i < size; i++)
		for (j = 0; j < size; j++)
			*(tmp + i*size + j) = *(bP + j*size + i);

	printf("\n \n \n");
	printf("start multiply size:%d \n", size);
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			for (k = 0; k < size; k++) {
				*(cP + i*size + j) = *(cP + i*size + j) + *(aP + i*size + k) * *(bP + k*size + j);
			}
		}
		printf("still running! %d %%\r", i * 100 / size);
	}
	printf("\n");
}
//simple matrix output
//for check use
int Matrix_output(Data_type *aP, int size) {
	int i, j;
	for (i = 0; i <= size - 1; i++)
	{
		for (j = 0; j <= size - 1; j++)
			//printf("%d  ", *(aP + i*size + j));
		//printf("\n");
	}
	return 0;
}
//the main function
void main()
{
	printf("pls enter a number as the matrix size(max 99999):  ");
	int MatrixSize = 0;
	char s[5];
	fgets(s, 5, stdin);
	MatrixSize = atoi(s);
	clock_t begin, end; //clock vairables
	double totalTime = 0; // time output
	srand((unsigned)time(0));
	//define dynamic matrix for thows from chkstk.asm after 256*256 matrix
	int *a = NULL, *b = NULL, *c = NULL;
	a = (Data_type*)malloc(MatrixSize*MatrixSize * sizeof(Data_type));

	b = (Data_type*)malloc(MatrixSize*MatrixSize * sizeof(Data_type));

	c = (Data_type*)malloc(MatrixSize*MatrixSize * sizeof(Data_type));

	//initialize matrixs
	Initialize_Matrix(a, MatrixSize, 1);
	Initialize_Matrix(b, MatrixSize, 1);
	Initialize_Matrix(c, MatrixSize, 0);
	//the matrix multiply begin
	begin = clock(); //begin to count time
	Multiply_matrix(a, b, c, MatrixSize);
	//Matrix_output(c, Matrix_size);
	end = clock();//end the counter
	totalTime = Calculate_Run_Time(&begin, &end);
	printf("%lf", totalTime);
	free(a);
	free(b);
	free(c);
}

