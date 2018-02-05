#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <time.h>

int test = 0; /* Keep this value 1 in order to conduct a sanity check of the correctness of the parallel version, zero otherwise*/
float **A, *B, *X;/* A * B = X
                        where A is the coefficient matrix, B is the constant matrix and X is the variables matrix */
float **Aseq, *Bseq, *Xseq; /*These matrices are only allocated if the test condition is true*/
int N; /*size of the square matrix*/


/*Function to check the sanity of the parallel and sequential versions*/
void sequentialsanitycheck()
{
        int i,j,k;
        float xel;
        /*Forward elimination*/
        for (k = 0; k < N - 1; k++)
        {
                for (i = k + 1; i < N; i++)
                {
                        xel = Aseq[i][k] / Aseq[k][k];
                        for (j = k; j < N; j++)
                        {
                                Aseq[i][j] -= Aseq[k][j] * xel;
                        }
                Bseq[i] -= Bseq[k] * xel;
                }
        }
        /*Backsubstitution*/
        for (i = N - 1; i >= 0; i--)
        {
                Xseq[i] = Bseq[i];
                for (j = N-1; j > i; j--)
                {
                        Xseq[i] -= Aseq[i][j] * Xseq[j];
                }
                Xseq[i] /= Aseq[i][i];
        }
        /*Check if X and Xseq are eual*/
        for(i=0;i<N;i++){
 if(X[i]!=Xseq[i]){
                        break;  /* break out of the loop whenever the values are uneuqal*/
                }
        }
        if(i==N){ /*If did not break out of the loop i.e. all values Of X and Xseq match*/
                printf("Both the parallel and sequential versions of X matrices have same values.\n");
        }
        else{
                printf("The parallel and sequential versions of X matrices have different values.\n");
        }
}

void backsubstitution()
{
        int i,j;
         for (i = N - 1; i >= 0; i--) {
                X[i] = B[i];
                for (j = N-1; j > i; j--) {
                        X[i] -= A[i][j] * X[j];
                }
                X[i] /= A[i][i];
        }

        if(N<6){
                for(i = 0; i < N; i++){
                        printf("\nX[%d]= %f\n",i,X[i]);
                }
        }

}
int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
    int i, j, k, rank, numProcessors;
    float xel;
    clock_t start, end;
    double eliminationTime;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcessors);
        start=(clock_t)0;
    if (!rank)
 {
        printf("Enter the order of the matrix: \n");
        if(scanf("%d", &N)){};
                if(N>=6){printf("\nSize too large, won't be able to display the matrices. If you want to see a sample result, please enter a size less than or equal to 5\n");}
                printf("Processing Gaussian elimination...\n");
        }
                MPI_Bcast (&N,1,MPI_INT,0,MPI_COMM_WORLD); /*Broadcast the value of N to all the procesors*/
                A= (float**) calloc(N, sizeof(float *));        /*Allocate the matrices A, B and X in the memory */
                for(i=0;i<N;i++)
                {
                        A[i]=(float *)calloc(N,sizeof(float));
                }
                B=calloc(N,sizeof(float));
                X=calloc(N,sizeof(float));
                if(test==1){ /*If sanity check is required, allocate memory for Aseq, Bseq and Xseq matrices*/
                        Aseq= (float**) calloc(N, sizeof(float *));
                        for(i=0;i<N;i++)
                        {
                                Aseq[i]=(float *)calloc(N,sizeof(float));
                        }
                        Bseq=calloc(N,sizeof(float));
                        Xseq=calloc(N,sizeof(float));
                }
        if (!rank){
                if(N<6){printf("A matrix: \n");}
                for(i=0;i<N;i++){
                        for(j=0;j<N;j++){
                                A[i][j]=((float)rand()/RAND_MAX)*10;
                                if(N<6){printf("%f \t", A[i][j]);}
                                if(test){Aseq[i][j]=A[i][j];}
                        }
                        if(N<6){printf("\n");}
                        B[i]=((float)rand()/RAND_MAX)*10;
                        if(test){Bseq[i]=B[i];}

                }
                if(N<6){printf("\nB matrix:\n");}
                for(i=0;i<N;i++){
                        if(N<6){printf("%f\t",B[i]);}
                }
 if(N<6){printf("\n");}
                start = clock();
/*A[0][0]=8.401877;
A[0][1]=3.943829;
A[0][2]=7.830992;
A[1][0]=7.984400;
A[1][1]=9.116474;
A[1][2]=1.975514;
A[2][0]=3.352227;
A[2][1]=7.682296;
A[2][2]=2.777747;
B[0]=5.539700;
B[1]=4.773971;
B[2]=6.288709;
*/

    }
        /* Broadcast the A and B matrices to all the processors*/
        for(i=0;i<N;i++){
                MPI_Bcast (&A[i][0],N,MPI_FLOAT,0,MPI_COMM_WORLD);
        }
        MPI_Bcast (B,N,MPI_FLOAT,0,MPI_COMM_WORLD);

    for(i=0;i<N;i++)
    {
        MPI_Bcast (&(A[i][i]),N-i,MPI_FLOAT,i % numProcessors,MPI_COMM_WORLD);
        MPI_Bcast (&B[i] ,1,MPI_FLOAT,i % numProcessors,MPI_COMM_WORLD);
        for(j=i+1; j<N; j++){
                        /*Parallel execution by multiple processors*/
            if(j % numProcessors == rank){
                xel = A[j][i] / A[i][i];
                for(k=i; k<N; k++)
                {
                    A[j][k]-= xel*A[i][k];
                }
                B[j] -= xel* B[i];
            }
        }
        if(!rank){memset(&A[i][0],0,i*sizeof(float)); }
    }
    if(!rank)
 {
        end= clock();
        eliminationTime=(double)(end-start)/CLOCKS_PER_SEC;
        backsubstitution();
        printf("Forward elimination time=%fs\n", eliminationTime);
        if(test){  /*conduct a sanity test if test is set to 1*/
                printf("Performing sequential sanity check......\n");
                sequentialsanitycheck();
        }
    }
    MPI_Finalize();
    return 0;

}
