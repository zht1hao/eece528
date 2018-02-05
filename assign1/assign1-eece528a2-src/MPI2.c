/*
MPI Matrix gaussian elimination in C
EECE 528
Guangyu Zhang
17550138
*/


#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define type_matrix float
#define randbig 1000.0

int main(int argc, char** argv) {

    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
    int size=4;
    int i,j,k;
    type_matrix c;
    type_matrix *A;
    int *map;
    clock_t start_calculation;
    clock_t end_calcualtion;
    // Find out rank, size
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (world_rank==0)
    {
        printf("Please input an integer value for matrix size: \n");
        scanf("%d", &size);

    }
    MPI_Bcast (&size,1,MPI_INT,0,MPI_COMM_WORLD);
    A=(type_matrix*)malloc(size*size*sizeof(type_matrix));
   // c=(type_matrix*)malloc(size*sizeof(type_matrix));
    map=(int*)malloc(size*sizeof(int));
   // MPI_Barrier(MPI_COMM_WORLD);

    if (world_rank==0)
    {
            for(i=0;i<size;i++)
            {
                for(j=0;j<size;j++)
                {
                    *(A+i*size+j)=((type_matrix)rand()/(type_matrix)(RAND_MAX)) * randbig;
                }
            }
    }
    MPI_Bcast (A,size*size,MPI_FLOAT,0,MPI_COMM_WORLD);
    //MPI_Barrier(MPI_COMM_WORLD);

    for(i=0; i<size; i++)
    {
        *(map+i)= i % world_size;
    }
    MPI_Bcast (map,size,MPI_INT,0,MPI_COMM_WORLD);

    if (world_rank==0)
    {
        start_calculation = clock();
    }

    for(i=0;i<size;i++)
    {
        MPI_Bcast ((A+i*size+i),size-i,MPI_FLOAT,*(map+i),MPI_COMM_WORLD);
        for(j=i+1; j<size; j++)
        {
            if(*(map+j)==world_rank)
            {
                c=*(A+j*size+i) / *(A+i*size+i);
                //printf("%f \n", c);
                for(k=i; k<size; k++)
                {
                    *(A+j*size+k)=*(A+j*size+k)-c*(*(A+i*size+k));
                }
            }
        }
        if(world_rank==0)
        {
	    for(j=0;j<i;j++)
		*(A+i*size+j)=0;
        }
    }



    if(world_rank==0)
    {
        end_calcualtion= clock();
       /* for (j=0;j<size;j++)
        {
            for (i=0;i<size;i++)
            {
               printf("%*.*f",15,2,*(A+j*size+i));
            }
            printf("\n");
        }*/
        double time_passed=(double)(end_calcualtion-start_calculation)/CLOCKS_PER_SEC;
        printf("Gaussian Elimination time=%f seconds. \n", time_passed);
        printf("End of the program! \n");
    }
    MPI_Finalize();
    return 0;

}
