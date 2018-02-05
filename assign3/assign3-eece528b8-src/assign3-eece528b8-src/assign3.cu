#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "cuda.h"
#include <string.h>

#define N 4096
#define threadNum 1024
#define blockNum N/threadNum

float *a;

void Forward();
void Initialize_input(float *ary);
void PrintMat(float *ary);
void checkCUDAError(const char *msg);

unsigned int totalKernelTime = 0;

int main()
{
    a = (float *)malloc(N * N * sizeof(float));
    Initialize_input(a);
    Forward();

    printf("runtime = %f s\n",totalKernelTime * 1e-6);
    //PrintMat(a);
}

__global__ void elimination(float *a, int t)
{
    __shared__ float multiple;
    if (threadIdx.x == 0)
    {
        multiple = a[blockIdx.y * N  + t] / a[t * N + t];
    }
    __syncthreads();
    int xidx =  blockIdx.x * blockDim.x + threadIdx.x;

    if (blockIdx.y > t && t == xidx)
    {
        a[N * blockIdx.y + xidx] = 0;
    }
    if (blockIdx.y > t && xidx > t)
    {
        a[N * blockIdx.y +xidx] -= multiple * a[N * t + xidx];
    }
}

void Forward()
{
    int t;
    float *cuda_a;

    cudaMalloc((void **) &cuda_a, N * N * sizeof(float));

    struct timeval time_start, time_end;
    gettimeofday(&time_start, NULL);

    cudaMemcpy(cuda_a, a, N * N * sizeof(float), cudaMemcpyHostToDevice);

    dim3 dimBlock(threadNum, 1, 1);
    dim3 dimGrid(blockNum, N, 1);

    for (t = 0; t < N; t++)
    {
        elimination<<<dimGrid,dimBlock>>>(cuda_a, t);
        cudaThreadSynchronize();
    }

    cudaMemcpy(a, cuda_a, N * N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaFree(cuda_a);

    gettimeofday(&time_end, NULL);
    totalKernelTime = (time_end.tv_sec * 1000000 + time_end.tv_usec) - (time_start.tv_sec * 1000000 + time_start.tv_usec);
}

void Initialize_input(float *ary)
{
    int i, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N ; j++)
        {
            ary[i * N + j] = rand() % 32768;
            while (ary[i * N + j] == 0)
                ary[i * N + j] = rand() % 32768;
        }
    }
}

void PrintMat(float *ary)
{
    int i, j;
    for (i = N - 5; i < N; i++)
    {
        for (j = N - 5; j < N ; j++)
        {
            printf("%8.2f ", ary[i * N + j]);
        }
        printf("\n");
    }
    printf("\n");
}

void checkCUDAError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err)
    {
        fprintf(stderr, "Cuda error: %s: %s.\n", msg,
                                  cudaGetErrorString( err) );
        exit(EXIT_FAILURE);
    }
}
