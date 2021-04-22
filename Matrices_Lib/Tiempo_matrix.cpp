#include <iostream>
#include<cstdlib>
#include <eigen3/Eigen/Dense>
#include<chrono>

int create_arrays(int N);

int main(int argc, char** argv)
{
    using namespace std::chrono;
    int dimensions[] = {4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048};
    int N,len = sizeof(dimensions)/sizeof(dimensions[0]);
    const int SEED = 3;
    srand(SEED); //control seed

    std::cout<<"tamaño de la matriz\t"
             <<"Solve takes: [s]\t"
             <<"Eigen vectors takes: [s]\t"
             <<"Accuracy of solution\t"
             <<"The second component of an eigenvector is:"
             <<std::endl;

    for(int i = 0;i<len;i++){
        N = dimensions[i];
        Eigen::MatrixXd A = Eigen::MatrixXd::Random(N,N);
        Eigen::VectorXd b = Eigen::VectorXd::Random(N);
        //tik
        steady_clock::time_point t1 = steady_clock::now();

        Eigen::VectorXd sol = A.colPivHouseholderQr().solve(b);
        //tok
        steady_clock::time_point t2 = steady_clock::now();

        Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(A);
        //tik
        steady_clock::time_point t3 = steady_clock::now();

        Eigen::VectorXcd v = eigensolver.eigenvectors().col(0);
        duration<double> time_sol = duration_cast<duration<double>>(t2 - t1);
        duration<double> time_eigen = duration_cast<duration<double>>(t3 - t2);
        std::cout <<N<<"\t"
                  <<time_sol.count()<<"\t"
                  <<time_eigen.count()<<'\t'
                  <<(A*sol-b).norm()<<'\t'
                  <<v(1)<<std::endl;
    }

    return 0;
}
