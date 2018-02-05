#include<stdio.h>
#include<stdlib.h>
#include <time.h>
#include<mpi.h>
#define n 4096

typedef struct{
    float value;
    int id;
} value_id;


float A[n][n];
float separate_A[n][n];
float global_max[n];
int pivoting[n];

int main(int argc, char *argv[])
{
    float maxvalue,tmp;
    double start,end,runtime;
    int i,j,k,myid,numprocs,maxindex;
    value_id me,max;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

    int numtasks=n/numprocs;
    for (i=0; i<numtasks; i++)
    {
        pivoting[i]=0;
    }

    //initiate
    if (myid==0)
    {
        srand(0);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                A[i][j] = rand();
            }
        }
    }

    MPI_Bcast (&A[0][0],n*n,MPI_FLOAT,0,MPI_COMM_WORLD);

    k=0;
    for(i=0;i<n;i++)
    {
        if(myid==i%numprocs)
        {
            for (j=0; j<n;j++)
            {
                separate_A[k][j]=A[i][j];
            }
            k++;
        }
    }


    start=MPI_Wtime();

    for(i=0;i<n;i++)
    {
        maxvalue=0;
        for(j=0;j<numtasks;j++)
        {
            if (pivoting[i]==0&&abs(separate_A[j][i])>maxvalue)
            {
                maxvalue = abs(separate_A[j][i]);
                maxindex = j;
            }
        }

        if(pivoting[maxindex]==0)
        {
            me.value=separate_A[maxindex][i];
        }else
        {
            me.value=0;
        }
        me.id=myid;

        MPI_Allreduce(&me,&max,1,MPI_FLOAT_INT,MPI_MAXLOC,MPI_COMM_WORLD);

        if(myid==max.id)
        {
            for(j=0;j<n;j++)
            {
                global_max[j]=separate_A[maxindex][j];
            }
            pivoting[maxindex]=1;
        }

        MPI_Bcast(global_max,n,MPI_FLOAT,max.id,MPI_COMM_WORLD);

        for(j=0;j<numtasks;j++)
        {
            if(pivoting[j]==0)
            {
                tmp=separate_A[j][i]/global_max[i];
                for(k=i+1;k<n;k++)
                {
                    separate_A[j][k]-=global_max[k]*tmp;
                }
            }
        }
    }

    end=MPI_Wtime();
    runtime=end-start;

    if(myid==0)
    {
    	printf("runtime: %f\n",runtime);
    }
    
    
    MPI_Finalize();
    return 0;
}

