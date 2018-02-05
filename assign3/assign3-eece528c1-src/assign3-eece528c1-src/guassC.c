#include <stdlib.h>
#include <stdio.h>
#include <time.h>

int main(int argc, char** argv)
{

    //Declarations
    srand ( time(NULL) );
    int i,j,k,n;
    float c,sum=0.0,temp=0.0;
    float **A = NULL;
    float *x = NULL;
    clock_t t, t1, t2;

    //Command line arguments
    n = atoi(argv[1]);

    A = malloc(sizeof(float*)*n);
    x = (float*)malloc(sizeof(float)*n);

    for(i=0;i<n;i++)
        A[i] = malloc(sizeof(float)*n);

    for(i=0; i<n; i++)
    {
        for(j=0; j<(n+1); j++)
        {
            A[i][j] = ((float)rand()/(float)(RAND_MAX)) * 5.0;
        }
    }

    t1 = clock();
    for(j=0; j<n; j++) /* loop for the generation of upper triangular matrix*/
    {
        for(i=0; i<n; i++)
        {
            if(i>j)
            {
                c=A[i][j]/A[j][j];

                for(k=0; k<n+1; k++)
                {
                        A[i][k]=A[i][k]-c*A[j][k];
                }
            }
        }
    }
    t2 = clock();
    printf("\nThe time taken is %f seconds\n",(float)(t2-t1)/CLOCKS_PER_SEC);
    x[n-1]=A[n-1][n]/A[n-1][n-1];

    /* this loop is for backward substitution*/
    for(i=n-2; i>=0; i--)
    {
        sum=0;
        for(j=i; j<n; j++)
        {
                sum+= A[i][j]*x[j];
        }
        x[i]=(A[i][n]-sum)/A[i][i];
    }
    return(0);
}
