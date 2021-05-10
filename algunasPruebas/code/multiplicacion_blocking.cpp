#include <iostream>
#include <cstdlib>
#include <vector>
#include "basic_operations.h"

int main(int argc, char **argv)
{

     const int N = std::atoi(argv[1]); //Tamaño de la matriz

     const int BLOCK = std::atoi(argv[2]); //Tamaño del bloque

     
    // Matrix declaration : Modeled as 1D array
     matriz A(N*N);
     matriz B(N*N);
     matriz C(N*N);
    // Declare as pointers and ask for memory to use the heap
    // initialize matrices
    for (int ii = 0; ii < N; ++ii) {
      for (int jj = 0; jj < N; ++jj) {
      A[ii*N + jj] = ii*jj;
      B[ii*N + jj] = jj+2;
      C[ii*N + jj] = 0.0;
      }
    }
   // A = {2,0,1,3,0,0,5,1,1};
   // B = {1,0,1,1,2,1,1,1,0};
    multiply_blocking(A, B, C, N, BLOCK);
    print_matriz(C,N,N);

    return 0;
  }
