# include <iostream>
# include <cstdio>
# include <cstdlib>
# include <eigen3/Eigen/Dense>
# include "papi.h"
//# include <Eigen/Dense>

int code_to_be_measured(const Eigen::MatrixXd & M, Eigen::MatrixXd & MT);


int main(int argc, char **argv)
{
  const int N = std::atoi(argv[1]);
  // Matrix declaration : Modeled as 1D array
  // Declare as pointers and ask for memory to use the heap
  // Matrix declaration
  Eigen::MatrixXd A(N, N), AT(N, N);
  
  // initialize matrices
  for (int ii = 0; ii < N; ++ii) {
    for (int jj = 0; jj < N; ++jj) {
      A(ii*N + jj) = ii/N + jj/N + 0.99/N; //evitar overflow
      AT(ii*N + jj) = 0.0;
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
  
  code_to_be_measured(A, AT);
  if((retval=PAPI_flops_rate(PAPI_FP_OPS,&real_time, &proc_time, &flpops, &mflops))<PAPI_OK)
    {
      printf("retval: %d\n", retval);
      exit(1);
    }
  printf("Real_time: %f Proc_time: %f Total flpops: %lld MFLOPS: %f\n",
	 real_time, proc_time,flpops,mflops);
  // Do something here, like computing the average of the resulting matrix, to avoid the optimizer deleting the code
  std::cout <<AT.sum()<<"\n";
  return 0;
}

int code_to_be_measured(const Eigen::MatrixXd & M, Eigen::MatrixXd & MT)
{
  MT = 2.3456*M.transpose().eval();
  return 0;
}
