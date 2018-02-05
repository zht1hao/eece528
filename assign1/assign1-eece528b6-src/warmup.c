#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define n 4096


float A[n][n];
float temp;
int i,j,k;


void mat_init()
{
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            A[i][j]=rand();
        }
    }
}

void swap(float a,float b)
{
    temp=a;
    a=b;
    b=temp;
}

void row_swap(float A[n][n])
{
    int max;
    for(i=0;i<n-1;i++)
    {
        max=i;
        for(j=i+1;j<n;j++)
        {
            if(abs(A[i][i]) < abs(A[j][i]))
            {
                max=j;
            }
        }
        for(k=0;k<n;k++)
        {
            swap(A[i][k],A[max][k]);
        }
    }
}

void gaus_elim(float A[n][n])
{
    for(i=0;i<n-1;i++)
    {
        for(j=i+1;j<n;j++)
        {
            temp=A[j][i]/A[i][i];
            for(k=i+1;k<n;k++)
            {
                A[j][k]-=temp*A[i][k];
            }
            A[j][i]=0;
        }
    }
}

void print_mat(float A[n][n])
{
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            printf("%6.0f ",A[i][j]);
        }
        printf("\n");
    }
}

int main()
{
    //initialize
    srand(0);
    mat_init(A);

    clock_t start= clock();

    //row swapping
    row_swap(A);

    //elimination
    gaus_elim(A);

    clock_t end = clock();

    //print_mat(A);
    float runtime = (float)(end-start)/CLOCKS_PER_SEC;
    printf("runtime: %f\n",runtime);
    return 0;
}

