#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<math.h>

int main()
{
	
	float **A,*B,*X;
	int i, j, k, size;
		
	printf("\nEnter the order of matrix: ");
	if(scanf("%d",&size)){};

	A= (float**) calloc(size, sizeof(float *));
	for(i=0;i<size;i++)
	{
		A[i]=(float *)calloc(size,sizeof(float));
	}
	B=calloc(size,sizeof(float));
	X=calloc(size,sizeof(float));

	for(i=0; i<size; i++)
	{
		for(j=0; j<size; j++)
		{
			A[i][j]=((float)rand()/RAND_MAX)*10;
			/*printf("%f ",A[i][j]);*/
		}
		/*printf("\n");*/
	}
	/* printf("B array is: ");*/

	for(i=0; i<size; i++)
	{
		B[i]=((float)rand()/RAND_MAX)*10;
		/*  printf("%f ",B[i]);*/
	}

	
/* **********************Forward elimination**********************  */

	clock_t start, end;
    start=clock();
	float multiplier;

	for (k = 0; k < size - 1; k++) 
	{
		for (i = k + 1; i < size; i++)
		{
			multiplier = A[i][k] / A[k][k];
			for (j = k; j < size; j++) 
			{
				A[i][j] -= A[k][j] * multiplier;
			}	
		B[i] -= B[k] * multiplier;
		}
	}

	end=clock();

/* **********************Backward Substitution**********************  */

	for (i = size - 1; i >= 0; i--) 
	{
		X[i] = B[i];
		for (j = size-1; j > i; j--) 
		{
			X[i] -= A[i][j] * X[j];
		}
		X[i] /= A[i][i];
	}

	/*printf("\nx is: size");
	for(i=1;i<=size;i++)
	printf("%f\t",X[i]);*/

	printf("\nExecution time: %f\ns",(double)(end-start)/CLOCKS_PER_SEC);
	return(0);
}
