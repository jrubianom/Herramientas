#include<iostream>
#include<cstdio>
#include "omp.h"

double Integral(int steps);
double integral(int steps);

int main(int argc, char** argv){
    int NTHRDS = std::atoi(argv[1]);
    std::cout.precision(15); std::cout.setf(std::ios::scientific);
    const int steps = 100000000;
    omp_set_num_threads(NTHRDS);
    double result = Integral(steps);
double dif = result-(10*10*10.0/3);
std::cout <<dif<<std::endl;

}

double integral(int steps){
    double sum = 0.0;
   #pragma omp parallel
    {
    double x,dx = 10.0/steps;
    #pragma omp for reduction(+:sum)
    for(int i = 0; i<steps; i++){
        x = i*dx;
        sum += dx*x*x;
    }
    }
    return sum;
}

double Integral(int steps){
    double sum = 0.0,dx = 10.0/steps;
    int i;
   #pragma omp parallel for reduction(+:sum)
    for(i = 0; i<steps; i++){
        sum += dx*i*i*dx*dx;
    }
    return sum;
}
