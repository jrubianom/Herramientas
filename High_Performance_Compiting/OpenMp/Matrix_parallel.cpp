#include<eigen3/Eigen/Dense>
#include<iostream>
#include<cstdlib>
#include "omp.h"

void print_matrix(double* A, int n, int m);
void init_row_major(double* A,Eigen::MatrixXd &X,int seed);
void init_column_major(double* A,Eigen::MatrixXd &X,int seed);
void Mv(double* A,double* B, double* C,int n, int m, int col);
void multiplication(double* A,double* B, double* C, int n, int m,int p,int nthrds);

int main(int argc,char** argv){

    int nthrds = std::atoi(argv[1]);
    int n = 2,m=3,p=4;
    double A[n*m], B[m*p], C[n*p]={0.0};
    Eigen::MatrixXd X(n,m),Y(m,p),Z(n,p);
    init_row_major(A,X,2);
    init_column_major(B,Y,3);
    Z = X*Y;
    multiplication(A,B,C,n,m,p,nthrds);
    print_matrix(C,n,p);
    std::cout<<Z<<std::endl;

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

void Mv(double* A,double* B, double* C,int n, int m,int col){
    for(int i=0;i<n;i++){
        for(int j = 0; j<m;j++){
            C[col*n+i] += (A[i*m+j]*B[col*m+j]);
        }
    }
}

void multiplication(double* A,double* B, double* C, int n, int m,int p,int nthrds){

    int N = p/nthrds;
    int r = p % nthrds;
#pragma omp parallel num_threads(nthrds)
    {
    int ID = omp_get_thread_num();
    int index_0 = N*ID; int index_f = N*(ID+1);
    for (int col=index_0; col<index_f;col++){
        Mv(A,B,C,n,m,col);
    }
    }
}
