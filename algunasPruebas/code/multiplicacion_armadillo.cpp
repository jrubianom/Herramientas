# include <iostream>
# include <cstdio>
# include <cstdlib>
# include <armadillo>
# include "papi.h"

int code_to_be_measured(const arma::dmat & M1,const arma::dmat & M2, arma::dmat & P);

int main(int argc, char **argv){
    const int N = std::atoi(argv[1]);

// Matrix declaration : Modeled as 1D array
// Declare as pointers and ask for memory to use the heap
// Matrix declaration

    arma::dmat A(N, N), B(N, N),C(N, N);

// initialize matrices, armadillo handles column-major order

    for (int jj = 0; jj < N; ++jj) {
        for (int ii = 0; ii < N; ++ii) {
            A(ii + jj*N) = ii/N + jj/N + 0.99/N;
            B(ii + jj*N) = ii/N + jj/N - 0.99/N;
            C(ii + jj*N) = 0.0;
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

    code_to_be_measured(A, B, C);

    if((retval=PAPI_flops_rate(PAPI_FP_OPS,&real_time, &proc_time, &flpops, &mflops))<PAPI_OK)
    {
        printf("retval: %d\n", retval);
        exit(1);
    }

    printf("Real_time: %f Proc_time: %f Total flpops: %lld MFLOPS: %f\n",
           real_time, proc_time,flpops,mflops);

// Do something here, like computing the average of the resulting matrix, to avoid the optimizer deleting the code
    printf("%.15e\n", accu(C));
    return 0;
}

int code_to_be_measured(const arma::dmat & M1, const arma::dmat & M2, arma::dmat & P)
{
    P = M1*M2;
    return 0;
}
