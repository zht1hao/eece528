#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/time.h>

int test = 1; /* Keep this value 1 in order to conduct a sanity check of the correctness of the parallel version, zero otherwise*/
int N;  /*size of the square matrix*/
float **A, *B, *X;/* A * B = X
                        where A is the coefficient matrix, B is the constant matrix and X is the variables matrix */
float **Aseq, *Bseq, *Xseq; /*These matrices are only allocated if the test condition is true*/

int *isProcessed,*thread_num, numOfThreads;
pthread_mutex_t mut;
pthread_cond_t condition;
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
        /*Check if X and Xseq are equal*/
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
                //printf("\nX[%d]= %f\n",i,X[i]);
                }
        X[i] /= A[i][i];
        }
        if(N<6){
                        printf("\nX matrix: ");
                        for(i = 0; i < N; i++){
                                        printf("%f ",X[i]);
                        }
                        printf("\n");
        }
}
void waitToProcess(int index)
{
        if(*(isProcessed+index)==0){
                pthread_mutex_lock(&mut);
                while(*(isProcessed+index)==0){
                        pthread_cond_wait(&condition, &mut);
                }
                pthread_mutex_unlock(&mut);
        }
}



void *threadedElimination(void *tdparam){
        int i,j,k;
        //struct threadData *td = (struct threadData*)tdparam;
        int *threadID = (int *)tdparam;
        float xel;
        for(i = 0; i < N - 1 ; i++){
                for (j = i + 1; j < N; j++){
  if( j % numOfThreads == *threadID)
                        {
                                xel = A[j][i] / A[i][i];
                                for (k = i; k < N; k++){
                                        A[j][k] -= xel * A[i][k];
                                }
                                B[j] -= B[i] * xel;
                                if(j==i+1){
                                        pthread_mutex_lock(&mut);
                                        pthread_cond_broadcast(&condition);
                                        *(isProcessed+i)=1;
                                        pthread_mutex_unlock(&mut);
                                }
                        }
                }
				waitToProcess(i);/*wait until the row i has been processed*/
				//pthread_exit(0);
        }
        return NULL;
}

int main()
{
        pthread_cond_init(&condition, NULL);
        pthread_mutex_init(&mut, NULL);
        double eliminationTime;
        struct timeval start,end;
        int j,i;
        int t=0;
        printf("\n***************************************Enter the order of matrix: ***********************************************");
        if(scanf("%d",&N)){};
        printf("Enter the number of threads in the program: ");
        if(scanf("%d", &numOfThreads)){};
        if(N>=6){printf("\nSize too large, won't be able to display the matrices. If you want to see a sample result, please enter a size less than or equal to 5\n");}
 printf("Processing Gaussian elimination...\n");
        A= (float**) malloc(N * sizeof(float *));
        for(i=0;i<N;i++)
        {
                A[i]=(float *)malloc(N * sizeof(float));
        }
        B = malloc(N * sizeof(float));
        X = malloc(N * sizeof(float));
        if(test==1){ /*If sanity check is required, allocate memory for Aseq, Bseq and Xseq matrices*/
            Aseq= (float**) calloc(N, sizeof(float *));
            for(i=0;i<N;i++)
            {
                    Aseq[i]=(float *)calloc(N,sizeof(float));
            }
            Bseq=calloc(N,sizeof(float));
            Xseq=calloc(N,sizeof(float));
        }
        thread_num = malloc(numOfThreads * sizeof(float));

        for(i = 0; i < N; i++){
                for(j = 0; j < N; j++){
                        A[i][j]=((float)rand()/RAND_MAX) * 10;
                        if(test){Aseq[i][j]=A[i][j];}
                }
                B[i]=((float)rand()/RAND_MAX) * 10;
                if(test){Bseq[i]=B[i];}
                X[i]=0.0f;
        }

/*A[0][0]=8.401877;
A[0][1]=3.943829;
A[0][2]=7.830992;
A[1][0]=8.984400;
A[1][1]=9.116474;
A[1][2]=1.975514;
A[2][0]=5.352227;
A[2][1]=7.682296;
A[2][2]=2.777747;
B[0]=5.539700;
B[1]=4.773971;
B[2]=6.288709;*/

        if(N<6){
                printf("A matrix entered: \n");
                for (i=0; i<N; i++){
                        for (j=0; j<N; j++){
                                printf("%f ",A[i][j]);
                        }
        printf("\n");
                }
        }
        isProcessed = (int*)calloc(N , sizeof(int));
        for(i=0;i<numOfThreads;i++){ thread_num[i]=i;}
        gettimeofday(&start, NULL);
        pthread_t *thread;
        thread = malloc(numOfThreads * sizeof (pthread_t));
        //struct threadData *td = malloc (N*sizeof (struct threadData));
        for (t = 0; t < numOfThreads; t++){
                pthread_create(&thread[t], NULL, threadedElimination, (void *)&thread_num[t]);
        }

        for (t = 0; t < numOfThreads; t++){
                pthread_join(thread[t], NULL);
        }

        pthread_mutex_destroy(&mut);
        pthread_cond_destroy(&condition);

        gettimeofday(&end, NULL);
        if (N<6)
        {
                printf("A matrix upper triangle \n");
                for (i = 0; i < N; i++){
                        for (j = 0; j < N; j++){
                                printf("%f ", A[i][j]);
                        }
                        printf("\n");
                }
                printf("\nB matrix: ");
                for (i = 0; i < N; i++){
                        printf("%f ", B[i]);
 }
                printf("\n");
        }
        eliminationTime=((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
        printf("Forward elimination execution time=%f seconds. \n", eliminationTime);
        backsubstitution();
        if(test){  /*conduct a sanity test if test is set to 1*/
                printf("Performing sequential sanity check......\n");
                sequentialsanitycheck();
        }
        return 0;
}
