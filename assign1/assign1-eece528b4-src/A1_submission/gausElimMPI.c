#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <assert.h>

#include <x86intrin.h>

// the following define should be passed as a compile-time parameter
// #define DIM 2048

// quick max macro funcion
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

int main(int argc, char *argv[])
{
  // mpi setup stuff
  MPI_Init(&argc, &argv);

  int rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  // saved myself some trouble here. Let's assume the following is true
  assert(DIM % nprocs == 0);

  // split up the work among all threads
  int myRows = DIM / nprocs;
  int startRow = rank * myRows;
  int endRow = startRow + myRows;

  // timing stuff
  clock_t start, diff;
  int elapsed_time_mm;

  // 
  float** a;
  // allocate the row pointers into the memory
  (a) = (float **)malloc(DIM*sizeof(float*));

  // set up the pointers into the contiguous memory
  for (int i=0; i<DIM; i++)
  {
    (a)[i] = (float*)malloc(DIM*sizeof(float));// &(p[i*DIM]);
  }
  // et voila, now I can use a[i][j] and a is guarantted to be contiguos in memory

  float f;

  // only master thread randomizes the matrix
  if (rank == 0)
  {
    // randomize
    srand(time(NULL));  

    for (int i=0; i < DIM; i++)
    {
      for (int j=0; j < DIM; j++)
      {
         a[i][j] = (float)rand();
      }
    }
  }

  // print the matrix for debugging
  // before any operations are done onthe matrix
  /*
  for (int i = 0; i < DIM; i++)
  {
    for (int j = 0; j < DIM; j++)
    {
      printf("%i ", (int)a[i][j]);
    }
    printf("\n");
  }

  printf("\n\n\n"); 
  */

  // master handles the clock:
  if (rank == 0)
  {
    start = clock();
  }

  // broadcast the corresponding chunk to each of the matrices
  // master does the sending
  if (rank == 0)
  {
    // send chunk to processes
    for (int i = 1; i < nprocs; i++)
    {
      MPI_Send(&a[i*myRows][0], DIM*myRows, MPI_FLOAT, i, i, MPI_COMM_WORLD);
    }
  }
  else
  {
    // each process will receive it's chunk from master
    MPI_Recv(&a[rank*myRows][0], DIM*myRows, MPI_FLOAT, 0, rank, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }
  // make sure all processes are up to speed
  // MPI_Barrier(MPI_COMM_WORLD); we dont need this because broadcast is a barrier 

  // matrix has been chunked up and sent. Now do the computation

  for(int j=0; j<DIM; j++)
  { 
    // broadcast the current row
    // need trick to see which process needs to broadcast: (int)(j / myRows)
    MPI_Bcast(&a[j][0], DIM, MPI_FLOAT, (int)(j / myRows), MPI_COMM_WORLD);
    
    // we need the max function here to determine where to start the loop
    for(int i = max(j+1, startRow); i < startRow + myRows; i++)
    {
      f = a[i][j]/a[j][j];
      // assert(f); // this assert may be triggered on a bad initialization, however generally it was useful for debugging
      
      // create a vector with f 
      __m128 f_scalar = _mm_set1_ps(-f);
      
      for(int k=j; k<DIM; k++)
      {
        // if k % 4 we can vecorize the rest of the loop
        // we need this check because SSE float vectors are 16-byte aligned
	if (k % 4 == 0)
        {
          int kk;
          for(kk=k; kk<DIM; kk+=4)
          {
            // load 4 elements starting at a[i][k] 
            __m128 one = _mm_load_ps(&a[i][k]);
            
            // load 4 elements starting at a[j][k]
            __m128 two =  _mm_load_ps(&a[j][k]);
            // one - f * two (scalar times vector)
            __m128 add = _mm_add_ps(one, _mm_mul_ps(two, f_scalar));

            // store the result at a[i][k]
            _mm_store_ps(&a[i][k], add);
            k+=4;
          }
	  //break;
        } 
        else
        {
          // do the normal way 1x1 element
          a[i][k] = a[i][k] - f*a[j][k];
        }
      }
    }
  }

  // print the matrix for debugging and equivalence check
  //
  // 
  //
  // enable the following lines of code, but redirect the output to file. 
  // Otherwise the command line will get flooded
  // Save the output file as the reference output.
  //
  // Trick to avoid float rounding errors: cast the output to int, only take what is in front of the comma
  // this may be problematic for small amtrices as high precision matters
  // but for large matrices this is irrelevant and avoids error due to rounding
  /*
  for (int i = 0; i < DIM; i++)
  {
    for (int j = 0; j < DIM; j++)
    {
      printf("%i ", (int)a[i][j]);
    }
    printf("\n");
  }
  */
  
  // again, let master handle the clock
  if (rank == 0)
  {
    // end the clock here, we only measure forward substitution
    diff = clock() - start;
    // print out the time tally
    elapsed_time_mm = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time tally (%i - float - %i) : %d milliseconds \n", DIM, nprocs, elapsed_time_mm);
  }

  // done
  MPI_Finalize(); 

  return 0;
} 
