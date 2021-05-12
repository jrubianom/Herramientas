#include <iostream>
#include<cstdio>
#include "omp.h"

static long num_steps = 1000000;
double step;
int main(int argc,char** argv)
{
    int NUM_THREADS = std::atoi(argv[1]);
    double pi=0.0,t2,t1;//sum[NUM_THREADS]
    int nthreads;
    step = 1.0/(double) num_steps;
    omp_set_num_threads(NUM_THREADS);
    t1 = omp_get_wtime();
#pragma omp parallel
{
    double x,sum=0.0; int ID,nthrds;
    ID = omp_get_thread_num();
    nthrds = omp_get_num_threads();
    if(ID==0) nthreads = nthrds;
    for(int i=ID;i<num_steps;i = i+nthrds){
        x = (i+0.5)*step;
        sum += 4.0/(1.0+x*x);
    }
#pragma omp critical
    pi += sum*step;
}

t2 = omp_get_wtime();

std::cout<<"PI:\t"<<pi<<"\tTime:\t"<<t2-t1<<std::endl;
return 0;

}
