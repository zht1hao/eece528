#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <pthread.h>
#define n 4096
#define thread_num 32

float A[n][n];
int flag[n];
int done[n];
int thread_id[thread_num];
pthread_t threads[thread_num];

void mat_init(float A[n][n])
{
    int i,j;
    for (i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            A[i][j] = rand();
        }
    }
}
void *gaus_elim(void *thread_id)
{
    int i,j,k;
    int *myid=(int*)thread_id;
    float temp;
    for(i=0;i<n;i++)
    {
        while(done[i]==0)
        {
            usleep(10);
        }
        for(j=(*myid);j<n;j+=thread_num)
        {
            if(j==i)
            {
                flag[j]=1;
            }
            if(flag[j]==0)
            {
                temp=A[j][i]/A[i][i];
                for(k=i;k<n;k++)
                {
                    A[j][k]-=temp*A[i][k];
                }
                A[j][i]=0;
            }
            if(j==i+1 && i<n-1)
            {
                done[i+1]=1;
            }
        }
    }
}



int main(int argc, char *argv[])
{
    int i,j,k;
    for(i=0;i<thread_num;i++)
    {
        thread_id[i]=i;
    }
    for(i=0;i<n;i++)
    {
        flag[i]=0;
        done[i]=0;
    }
    flag[0]=1;
    done[0]=1;

    //initialize
    srand(0);
    mat_init(A);
    

    struct timeval starttime,endtime;
    gettimeofday(&starttime,NULL);

    for (i=0;i<thread_num;i++)
    {
        pthread_create(&threads[i],NULL,gaus_elim,(void*)&thread_id[i]);
    }

    for (i=0;i<thread_num;i++)
    {
        pthread_join(threads[i],NULL);
    }

    gettimeofday(&endtime,NULL);
    double runtime=(double)(endtime.tv_sec-starttime.tv_sec)+(double)(endtime.tv_usec-starttime.tv_usec)/1000000;
    printf("runtime=%f\n",runtime);

    return 0;
}