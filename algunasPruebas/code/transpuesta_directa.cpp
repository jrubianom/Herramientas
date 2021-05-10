#include <iostream>
#include <cstdlib>
#include <vector>
#include "basic_operations.h"

int main(int argc, char **argv)
  {
    const int N = std::atoi(argv[1]);
    // Matrix declaration : Modeled as 1D array
    matriz A(N*N);
    matriz AT(N*N);
    // initialize matrices
    for (int ii = 0; ii < N; ++ii) {
      for (int jj = 0; jj < N; ++jj) {
        A[ii*N + jj] = ii + jj + 0.99;

      }
    }
    // Block size


    transpose_direct(A, AT, N);
    print_matriz(A,N,N);
    print_matriz(AT,N,N);


    return 0;
  }
