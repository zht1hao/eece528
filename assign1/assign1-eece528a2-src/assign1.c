/*
Vanilla Matrix gaussian elimination in C
EECE 528
Guangyu Zhang
17550138
*/

#include <time.h>
#include <stdio.h>
#include <stdlib.h>

//#define size 1024
#define type_matrix float
#define randbig 1000


int main()
{
    int size=0;
    printf("Please input an integer value for matrix size: \n");
    scanf("%d", &size);
    type_matrix *A;
    A=(type_matrix*)malloc(size*size*sizeof(type_matrix));
    int i,j,k;
    float c;
    printf("Program Start! \n");
    //generate A
    for(i=0;i<size;i++)
    {
        for(j=0;j<size;j++)
        {
            *(A+i*size+j)=((type_matrix)rand()/(type_matrix)(RAND_MAX)) * randbig;
        }
    }
/*
    for (j=0;j<size;j++)
    {
        for (i=0;i<size;i++)
        {
           printf("%*.*f",10,1,*(A+j*size+i));
        }
        printf("\n");
    }
    */
    //calculate C
    clock_t start_calculation = clock();
    for(i=0; i<size-1; i++)
    {
        for(j=i+1; j<size; j++)
        {
            c=*(A+j*size+i) / *(A+i*size+i);
            //printf("%f \n", c);
            for(k=i; k<size; k++)
            {
                *(A+j*size+k)=*(A+j*size+k)-c*(*(A+i*size+k));
            }
        }
    }
    clock_t end_calcualtion= clock();
    for (j=0;j<size;j++)
    {
       /* for (i=0;i<size;i++)
        {
           printf("%*.*f",10,1,*(A+j*size+i));
        }
        printf("\n");*/
    }
    double time_passed=(double)(end_calcualtion-start_calculation)/CLOCKS_PER_SEC;
    printf("Gaussian Elimination time=%f seconds. \n", time_passed);
    printf("End of the program! \n");
    return 0;
}
