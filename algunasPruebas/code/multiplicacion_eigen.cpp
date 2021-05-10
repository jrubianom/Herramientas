# include <iostream>
# include <cstdio>
# include <cstdlib>
# include <eigen3/Eigen/Dense>
//# include <Eigen/Dense>

int code_to_be_measured(const Eigen::MatrixXd & A, const Eigen::MatrixXd & B, Eigen::MatrixXd & C);


int main(int argc, char **argv)
{
  const int N = std::atoi(argv[1]);
  // Matrix declaration : Modeled as 1D array
  // Declare as pointers and ask for memory to use the heap
  // Matrix declaration
  Eigen::MatrixXd A(N, N), B(N, N), C(N, N);

  // initialize matrices
  for (int ii = 0; ii < N; ++ii) {
    for (int jj = 0; jj < N; ++jj) {
      A(ii*N + jj) = ii*jj;
      B(ii*N + jj) = jj+2;
      C(ii*N + jj) = 0.0;
    }
  }
//  A << 2, 0, 1,
//    3, 0, 0,
//         5, 1, 1;
//  B << 1, 0, 1,
//    1, 2, 1,
//    1, 1, 0;

  code_to_be_measured(A, B, C);
  std::cout<<A<<'\n'<<B<<'\n'<<C<<std::endl;
  return 0;
}

int code_to_be_measured(const Eigen::MatrixXd & A, const Eigen::MatrixXd & B, Eigen::MatrixXd & C)
{
  C = A*B;
  return 0;
}
