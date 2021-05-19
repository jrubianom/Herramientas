#include<eigen3/Eigen/Dense>
#include<iostream>
#include<cstdlib>
#include "omp.h"

void print_matrix(double* A, int n, int m);
void init_row_major(double* A,Eigen::MatrixXd &X,int seed);
void init_column_major(double* A,Eigen::MatrixXd &X,int seed);
void Multiply_serial(double* A,double* B, double* C,int n, int m, int p,int jmin, int jmax);
void Multiply_parallel(double* A,double* B, double* C, int n, int m,int p,int nthrds);
void make_copy(double* A, double* A_local,int n,int m);
void make_copy_cmo(double* B, double* B_local,int m,int jmin, int jmax);
void fill_matrix(double* C, double* C_local,int n, int jmin, int jmax);

int main(int argc,char** argv){

    int nthrds = std::atoi(argv[2]);
    int N = std::atoi(argv[1]);
    int n = N,m=N,p=N;
    //int n = 2, m = 4, p = 3;
    double t1,t2;
    double A[n*m], B[m*p], C[n*p]={0.0};
    Eigen::MatrixXd X(n,m),Y(m,p),Z(n,p);
    init_row_major(A,X,2);
    init_column_major(B,Y,3);
    Z = X*Y;
    t1 = omp_get_wtime();
    Multiply_parallel(A,B,C,n,m,p,nthrds);
    t2 = omp_get_wtime();
//    print_matrix(C,n,p);
//    std::cout<<Z<<std::endl;
    std::cout<<"Time:\t"<<t2-t1<<std::endl;
//    init_row_major(A,X,2);

}

void print_matrix(double A[], int n, int m){

    std::cout<<std::endl;
    for(int i = 0;i<n;i++){
        for(int j = 0;j<m;j++){
            std::cout<<"\t"<<A[j*n+i]<<"\t";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

void init_row_major(double* A,Eigen::MatrixXd &X,int seed){
    srand(seed);
    int n = X.rows(); int m = X.cols();
    X = Eigen::MatrixXd::Random(n,m);
    for(int i = 0; i<n; i++){
        for(int j=0; j<m;j++){
            A[i*m+j] = X(i,j);
        }
    }
}

void init_column_major(double* A,Eigen::MatrixXd &X,int seed){
    srand(seed);
    int n = X.rows(); int m = X.cols();
    X = Eigen::MatrixXd::Random(n,m);
    for(int j = 0; j<m; j++){
        for(int i=0; i<n;i++){
            A[j*n+i] = X(i,j);
        }
    }
}

void Multiply_serial(double* A,double* B, double* C,int n, int m,int p,int jmin, int jmax){
    for(int col = jmin;col<jmax;col++){
        for(int i=0;i<n;i++){
            for(int j = 0; j<m;j++){
                C[(col-jmin)*n+i] += (A[i*m+j]*B[col*m+j]);
            }
        }
    }
}

void make_copy(double* A, double* A_local,int n, int m){
    for(int i = 0; i<n; i++){
        for(int j=0; j<m;j++){
            A_local[i*m+j] = A[i*m+j];
        }
    }
}

void make_copy_cmo(double* B, double* B_local,int m, int jmin, int jmax){
    for(int j = jmin; j<jmax;j++){
        for(int i=0; i<m;i++){
            B_local[(j-jmin)*m+i] = B[j*m+i];
        }
    }
}

void fill_matrix(double* C, double* C_local,int n, int jmin, int jmax){
    for(int j = jmin; j < jmax ; j++){
        for(int i = 0; i < n ; i++){
            C[j*n+i] = C_local[(j-jmin)*n+i];
        }
    }
}

void Multiply_parallel(double* A,double* B, double* C, int n, int m,int p,int nthrds){

    int N = p/nthrds;
    int r = p % nthrds;
#pragma omp parallel num_threads(nthrds)
    {
    int ID = omp_get_thread_num();
    int index_0 = N*ID; int index_f = N*(ID+1);
    double C_local[n*N]={0.0};

Multiply_serial(A,B,C_local,n,m,N,index_0,index_f);
#pragma omp critical
{
fill_matrix(C,C_local,n,index_0,index_f);
}
  }
}
