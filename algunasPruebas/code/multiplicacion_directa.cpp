#include <iostream>
#include <cstdlib>
#include <vector>
#include "basic_operations.h"

int main(int argc, char **argv)
  {
    const int N = std::atoi(argv[1]);
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

        for (int ii=0; ii < N; ++ii){
        for (int jj=0; jj < N; ++jj){
            std::cout << B[N*ii + jj] << '\t';
        }
        std::cout << '\n';
    }
    std::cout<<'\n';

    // multiply_direct(A, B, C, N);
    // print_matriz(A,N,N);
    // std::cout<<'\n';
    // print_matriz(B,N,N);
    // std::cout<<'\n';
    // print_matriz(C,N,N);
    // std::cout<<'\n';


    return 0;
  }
