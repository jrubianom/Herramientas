#include <iostream>
#include <vector>
#include "basic_operations.h"

void print_matriz(matriz &A, int m, int n){
    for (int ii=0; ii < m; ++ii){
        for (int jj=0; jj < n; ++jj){
            std::cout << A[n*ii + jj] << '\t';
        }
        std::cout << '\n';
    }
}
matriz add(matriz &A, matriz &B, int m, int n){
    matriz temp(m*n);
    for (int ii=0; ii < m; ++ii){
        for (int jj=0; jj < n; ++jj){
            temp[n*ii + jj] = (A[n*ii + jj] + B[n*ii + jj]);
        }
    }
    return temp;
}

matriz substract(matriz &A, matriz &B, int m, int n){
    matriz temp(m*n);
    for (int ii=0; ii < m; ++ii){
        for (int jj=0; jj < n; ++jj){
            temp[n*ii + jj] = (A[n*ii + jj] - B[n*ii + jj]);
        }
    }
    return temp;
}

matriz copy(matriz &origin, int m, int n, int row1, int row2, int col1, int col2){
    /*
      Copia el bloque delimitado por las filas row1 y row2, y las columnas col1 y col2.
      m y n indican las dimensiones de la matriz original origin (m*n)
    */
    int rows = row2-row1+1, cols = col2-col1+1;

    matriz temp(rows*cols);

    for (int ii=0; ii < rows; ++ii){
        for (int jj=0; jj < cols; ++jj){
            temp[cols*ii + jj] = origin[(row1+ii)*n + col1 + jj];
        }
    }
    return temp;
}

void multiply_direct(matriz &A, matriz &B, matriz &C, int m){
    /*
      Multiplica las matrices A y B de tamaños m*p y p*n, respectivamente.
    */
    for (int ii=0; ii<m; ii++){
        for (int kk=0; kk<m; kk++){
            //double sum = 0;
            for (int jj=0; jj<m; jj++){
                C[m*ii + kk] += 1.0*A[m*ii + jj]*B[m*jj + kk];
            }

        }
    }
}

void transpose_direct(matriz &A, matriz &B,  int m)
{

  // procesarla : transponerla
  for(int ii=0; ii < m; ++ii){
    for(int jj=0; jj < m; ++jj){
      B[jj*m +ii] = 1.0*A[ii*m + jj];
    }
  }
}

void transpose_blocking(matriz &A, matriz &AT, int n, int BLOCK){
    for (int i = 0; i < n; i+=BLOCK) {
        for (int j = 0; j < n; j+=BLOCK) {
            for (int ib = 0; ib < BLOCK && i + ib < n; ++ib) {
                for (int jb = 0; jb < BLOCK && j + jb < n; ++jb) {
                    AT[(j+jb) * n + (i+ib)] = 1.0*A[(i+ib) * n + (j+jb)];
                }
            }
        }
    }
}



void multiply_blocking(matriz &A, matriz &B, matriz &C, int m, int s){

//ESTA MÉTODO DE MULTIPLICACIÓN POR BLOQUES SOLO ESTÁ VERIFICADO PARA MATRICES CUADRADA
  
  double temp;

  //multiplicacion por bloques
  for(int jj=0;jj<m;jj+= s){
    for(int kk=0;kk<m;kk+= s){
      for(int i = 0;i<m;i++){
        for(int j = jj; j<std::min(m, jj+s); j++){
          temp = 0;
          for(int k = kk; k<std::min(m, kk+s); k++){
            temp += A[m*i+k]*B[m*k+j];
          }
          C[m*i+j] += temp;
        }
      }
    }
  }
}
