
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

#include <x86intrin.h>

// the following define should be passed by the makefile
// may be useful for debugging
// #define DIM 1024

int main()
{
  srand(time(NULL));  

  // timing stuff
  clock_t start, diff;
  int elapsed_time_mm;

  // allocate the arrays on the heap
  // for ease of use, allocate them such that we can write a[i][j]
  // NOTE: a is of dimension DIM X DIM
  float** a;
  a = (float**)malloc(DIM * sizeof(float*));
  
  for (int i = 0; i < DIM; i++)
  {
    a[i] = (float*)malloc(DIM * sizeof(float));
  }

  // f helper
  float f = 0.f;

  // fill a and d with random numbers
  for (int i=0; i < DIM; i++)
  {
    for (int j=0; j < DIM; j++)
    {
      a[i][j] = (float)rand();
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
  
  //
  // forward iteration 
  //

  // start the clock
  start = clock();

  // create the upper triangular matrix
  for(int j=0; j<DIM; j++)
  {
    for(int i=j+1; i<DIM; i++)
    {
      f = a[i][j]/a[j][j];
      //assert(f);
      
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
        } 
        else 
        {
          // do the normal way 1x1 element
          a[i][k] = a[i][k] - f*a[j][k];
        }
      }
    }
  }

  diff = clock() - start;
  
  
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

  // end the clock here, we only measure forward substitution
  diff = clock() - start;

  // print out the time tally
  elapsed_time_mm = diff * 1000 / CLOCKS_PER_SEC;
  printf("Time tally (%i - float) : %d milliseconds \n", DIM, elapsed_time_mm);

  // free the matrix
  for (int i = 0; i < DIM; i++)
  {
    free(a[i]);
  }
  free(a);

  // done!
  return 0;
} 
