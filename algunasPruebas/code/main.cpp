#include <iostream>
#include <vector>
#include "basic_operations.h"

int main(){
    int m=4, p=5, n=6;
    matriz A(m*p);
    matriz B(p*n);

    int count = 0;
    for (int ii=0; ii<m; ++ii){
        for (int jj=0; jj<p; ++jj){
            A[ii*p + jj] = count;
            count++;
        }
    }
    count = 0;
    for (int ii=0; ii<p; ++ii){
        for (int jj=0; jj<n; ++jj){
            B[ii*n + jj] = count;
            count++;
        }
    }
    print_matriz(A, m, p);
    std::cout << '\n';
    print_matriz(B, p, n);
    std::cout << '\n';

    matriz E = multiply_direct(A, B, m, p, n);
    print_matriz(E, m, n);
    std::cout << '\n';

    matriz F = copy(A, m, p, 1, 3, 2, 4);
    print_matriz(F, 3, 3);

    return 0;
}
