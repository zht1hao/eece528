#include<stdlib.h>
#include<stdio.h>
#include<math.h>

#define TEST_ELEMENT 3
void setAnswer(float** A, float* b){

	FILE *l_ptrFileOutput ;
	l_ptrFileOutput = fopen("result.txt","w+");
	if( l_ptrFileOutput == NULL ){
		printf("\n Error : In opening file ");
	}else{
		float* result_vector = (float*) malloc(sizeof(float)*TEST_ELEMENT);
		float row_sum; 

		for (int j=0; j<TEST_ELEMENT; j++){
			row_sum = 0;
			for (int k=0; k<=TEST_ELEMENT; k++){
				row_sum += A[j][k];
			}
			result_vector[j] = row_sum;
			printf(" \n result_vector[ %f ]",result_vector[j]);
		}

		float sumOfSquares = 0;
		float entryOfResidual;
		for (int i=0; i<TEST_ELEMENT; i++){
			entryOfResidual = result_vector[i] - b[i];
			sumOfSquares += entryOfResidual*entryOfResidual;
			fprintf(l_ptrFileOutput,"%.20f\n", result_vector[i]);
		}
		sumOfSquares = sqrt(sumOfSquares);
		printf("\nThe L2-Norm of the result vector from Ax-b is: %.20f\n", sumOfSquares);
		//fprintf(l_ptrFileOutput,"%.20f\n",sumOfSquares);

		free(result_vector);
		fclose(l_ptrFileOutput);
	}
}

void checkAnswer(float** A ){

	FILE *l_ptrFileOutput ;
	l_ptrFileOutput = fopen("result.txt","r+");
	if( l_ptrFileOutput == NULL ){
		printf("\n Error : In opening file ");
	}else{
		float* result_vector = (float*) malloc(sizeof(float)*TEST_ELEMENT);
		float row_sum; 

		for (int j=0; j<TEST_ELEMENT; j++){
			row_sum = 0;
			for (int k=0; k<=TEST_ELEMENT; k++){
				row_sum += A[j][k];
			}
			result_vector[j] = row_sum;
			printf(" \n result_vector[ %f ]",result_vector[j]);
		}

		float sumOfSquares = 0;
		float entryOfResidual;
		float element = 0 ;
		for (int i=0; i<TEST_ELEMENT; i++){
			fscanf(l_ptrFileOutput,"%f", &element) ;
			entryOfResidual = result_vector[i] - element;
			sumOfSquares += entryOfResidual*entryOfResidual;
		}
		float squares = 0 ;	
		//fscanf(l_ptrFileOutput,"%f", &squares) ;
		sumOfSquares = sqrt(sumOfSquares);
		printf("squares === %.20f \n", squares);	

		if( sumOfSquares == 0 )
		 printf("\n SUCESSS :: The L2-Norm of the result vector from Ax-b is: %.5f\n", sumOfSquares);
		else
		 printf("\n FAIL :: The Bench Mark for L2-Norm  %.20f & result is %.5f \n", squares, sumOfSquares);
		
		free(result_vector);
		fclose(l_ptrFileOutput);
	}
}

int main()
{
	float **a,*temp,app,sum,mult;
	int i,j,k,n,p;
	//take no of terms
	// printf("Enter n : ");scanf("%d",&n);
	n =  3;

	//memory allocation
	a = (float**)malloc(n*sizeof(float*));
	for(i = 0; i < n; i++)
		a[i] = (float*)malloc(n*sizeof(float));
	temp = (float*)malloc(n*sizeof(float));
	//take input from user
	printf("Intialize the elements of augmended matrix rowwise\n");
	a[0][0] = 10;
	a[0][1] = -7;
	a[0][2] = 3;
	a[1][0] = -6;
	a[1][1] = 8;
	a[1][2] = 4;
	a[2][0] = 2;
	a[2][1] = 6;
	a[2][2] = 9;

for(i=0;i<n;i++){
                for(j=0;j<n;j++)
                        printf("%.2f\t",a[i][j]);
                printf("\n");
        }
	//generation of scalar matrix
	for(i=0;i<(n);i++){

		//calculating triangular matrix
		for(j=i+1;j<n;j++){
			mult = a[j][i]/a[i][i];
			for(k=0;k<n;k++)
				a[j][k] -= mult*a[i][k];
		}
	}
	//for calculating value of z,y,x via backward substitution method
	for(i=n-1;i>=0;i--)
	{
		sum = 0;
		for(j=i+1;j<n;j++)
			sum += a[i][j]*temp[j];
		temp[i] = (a[i][n]-sum)/a[i][i];
	}
	printf("****The matrix is : ***\n");
	for(i=0;i<n;i++){
		for(j=0;j<n;j++)
			printf("%.2f\t",a[i][j]);
		printf("\n");
	}
	//display solution
	printf("-------------The solution is ----------\n");
	for(i=0;i<n;i++)
		printf("X[%d] = %.2f\n",i+1,temp[i]);
	setAnswer(a, temp);
checkAnswer(a);
	//free allocated memory
	for(i = 0; i < n; i++)
		free(a[i]);
	free(a);
	free(temp);
	return 0;
}
