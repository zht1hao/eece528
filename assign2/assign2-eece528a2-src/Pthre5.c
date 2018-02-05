/*
Pthread Matrix gaussian elimination in C
EECE 528
Guangyu Zhang
17550138
*/

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/time.h>

//#define size 1024
#define type_matrix float
#define randbig 1000
#define total_thread 32
#define test 0
#define MIN(a,b) (((a)<(b))?(a):(b))

//global variables
type_matrix *A;
int size;
int thread1,thread0;
int *map1;
int *map2;
//int next_barrier1;
pthread_barrier_t mybarrier;
pthread_mutex_t eli_lock;
pthread_cond_t eli_cond;
//pthread_mutex_t ini_lock;
//pthread_cond_t ini_cond;

void *forward_eli(void *threadarg){
	int *thread_count = (int*)threadarg;
        int i,j,k;
        float c;
        *map1=1;
    for(i=0; i<size-1; i++)
    {
        if(*(map1+i)==0)
        {
            pthread_mutex_lock(&eli_lock);
            while(*(map1+i)==0)
            {
                pthread_cond_wait(&eli_cond, &eli_lock);
            }
            pthread_mutex_unlock(&eli_lock);
        }

        for (j=i+1; j<size; j++)
        {
            if(*(map2+j)==*thread_count)
            {
                // printf("line %d \n",j);
                c=*(A+j*size+i) / *(A+i*size+i);
                for (k=i; k<size; k++)
                {
                    // printf("line %d \n",k);
                *(A+j*size+k)=*(A+j*size+k)-c*(*(A+i*size+k));
                }
                if(j==i+1)
                {
                    pthread_mutex_lock(&eli_lock);
                    pthread_cond_broadcast(&eli_cond);
                    *(map1+j)=1;
                    pthread_mutex_unlock(&eli_lock);

                }
            }
        }





    if (test)
        {
            printf("thread %d \n",*thread_count);
        }
        //pthread_barrier_wait(&mybarrier);
    }
}

int main()
{
    pthread_t threads[total_thread];
    pthread_cond_init(&eli_cond, NULL);
    pthread_mutex_init(&eli_lock, NULL);
    //pthread_cond_init(&ini_cond, NULL);
    //pthread_mutex_init(&ini_lock, NULL);
    //pthread_barrier_init(&mybarrier, NULL, total_thread);
    struct timeval start_calculation,end_calcualtion;
    int j,i;
    int thread_count=0;
    printf("Please input an integer value for matrix size: \n");
    scanf("%d", &size);
    A=(type_matrix*)malloc(size*size*sizeof(type_matrix));
    printf("Program Start! \n");
    //generate A
    for(i=0;i<size;i++)
    {
        for(j=0;j<size;j++)
        {
            *(A+i*size+j)=((type_matrix)rand()/(type_matrix)(RAND_MAX)) * randbig;
        }
    }

   /* for (j=0;j<size;j++)
    {
        for (i=0;i<size;i++)
        {
           printf("%*.*f",10,1,*(A+j*size+i));
        }
        printf("\n");
    }*/

     //calculate C
     int *map = (int*)malloc(total_thread*sizeof(int));
     map1 = (int*)malloc(size*sizeof(int));
     map2 = (int*)malloc(size*sizeof(int));
	for(j = 0; j < total_thread; j++)
    {
    	map[j]= j;
    }
    for(j = 0; j < size; j++)
    {
    	*(map1+j) = 0;
    }
    for(j = 0; j < size; j++)
    {
    	*(map2+j) = j%total_thread;
    }
 //   next_barrier=1;
    gettimeofday(&start_calculation, NULL);

    for (thread_count = 0; thread_count < total_thread; thread_count++)
    {
        pthread_create(&threads[thread_count], NULL, forward_eli, (void*)&map[thread_count]);
    }

    for (thread_count = 0; thread_count < total_thread; thread_count++)
    {
        pthread_join(threads[thread_count], NULL);
    }

    //pthread_barrier_destroy(&mybarrier);
    pthread_mutex_destroy(&eli_lock);
    pthread_cond_destroy(&eli_cond);
    //pthread_mutex_destroy(&ini_lock);
    //pthread_cond_destroy(&ini_cond);

    gettimeofday(&end_calcualtion, NULL);
    if (test)
    {
        for (j=0;j<size;j++)
        {
            for (i=0;i<size;i++)
            {
            printf("%*.*f",10,1,*(A+j*size+i));
            }
            printf("\n");
        }
    }
    double time_passed=((end_calcualtion.tv_sec  - start_calculation.tv_sec) * 1000000u + end_calcualtion.tv_usec - start_calculation.tv_usec) / 1.e6;
    printf("Gaussian Elimination time=%f seconds. \n", time_passed);
    printf("End of the program! \n");
    return 0;
}
