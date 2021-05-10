# include <iostream>
# include <cstdio>
# include <cstdlib>
# include "papi.h"
#include <cmath>
#include "basic_operations.h"

int main(int argc, char **argv)
{
    int N = std::atoi(argv[1]); // Tamaño de la matriz
    int BLOCK = std::atoi(argv[2]); // Tamaño del bloque

    // Matrix declaration : Modeled as 1D array
    matriz A(N*N);
    matriz AT(N*N);
    // Declare as pointers and ask for memory to use the heap
    // initialize matrices
    for (int ii = 0; ii < N; ++ii) {
      for (int jj = 0; jj < N; ++jj) {
	A[ii*N + jj] = ii + jj + 0.99;
	AT[ii*N + jj] = 0.0;
      }
    }
    // PAPI vars
    float real_time, proc_time,mflops;
    long long flpops;
    float ireal_time, iproc_time, imflops;
    long long iflpops;
    int retval;
    // PERFOMANCE MEASURE
    // start PAPI counters
    if((retval=PAPI_flops_rate(PAPI_FP_OPS,&ireal_time,&iproc_time,&iflpops,&imflops)) < PAPI_OK)
      {
	printf("Could not initialise PAPI_flops \n");
	printf("Your platform may not support floating point operation event.\n");
	printf("retval: %d\n", retval);
	exit(1);
      }
    transpose_blocking(A, AT, N, BLOCK);
    if((retval=PAPI_flops_rate(PAPI_FP_OPS,&real_time, &proc_time, &flpops, &mflops))<PAPI_OK)
      {
	printf("retval: %d\n", retval);
	exit(1);
      }
    printf("Real_time: %f Proc_time: %f Total flpops: %lld MFLOPS: %f\n",
	   real_time, proc_time,flpops,mflops);
    // Do something here, like computing the average of the resulting matrix, to avoid the optimizer deleting the code
    printf("%.15e\n", AT[0]);

    return 0;
  }
