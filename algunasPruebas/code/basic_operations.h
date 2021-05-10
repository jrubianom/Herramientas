#ifndef _BASIC_OPERATIONS_
#define _BASIC_OPERATIONS_

#include <iostream>
#include <vector>

typedef std::vector<double> matriz;

void print_matriz(matriz &A, int m, int n);
matriz add(matriz &A, matriz &B, int m, int n);
matriz substract(matriz &A, matriz &B, int m, int n);
void multiply_direct(matriz &A, matriz &B, matriz &C, int m);
void multiply_blocking(matriz &A, matriz &B, matriz &C, int m, int s);
matriz copy(matriz &origin, int m, int n, int row1, int row2, int col1, int col2);
void transpose_direct(matriz &A, matriz &B, int m);
void transpose_blocking(matriz &A, matriz &AT, int n, int BLOCK);
#endif
