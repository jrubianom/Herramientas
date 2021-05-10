# include <iostream>
# include <cstdio>
# include <cstdlib>
# include "papi.h"
#include <cmath>
#include <numeric>
#include <fstream>
#include "basic_operations.h"

int main(int argc, char **argv)
{
    std::cout.precision(15);
    std::cout.setf(std::ios::scientific);

    std::ofstream fout("datos_transpuesta_blocking.txt");

    const int N_ITER = 10; // Numero de veces que se mide el desempeño de la funcion

    const int N = std::atoi(argv[1]); // Tamaño de la matriz

    float data_mflops[N_ITER], data_time[N_ITER]; // arreglos para los datos brutos de mflops y tiempo

    // Matrix declaration : Modeled as 1D array
    matriz A(N*N);
    matriz AT(N*N);

    // inicializar matrices
    for (int ii = 0; ii < N; ++ii) {
      for (int jj = 0; jj < N; ++jj) {
	A[ii*N + jj] = ii + jj + 0.99;
	AT[ii*N + jj] = 0.0;
      }
    }
    // PAPI vars
    float real_time, proc_time,mflops;
    long long flpops;
    float ireal_time, iproc_time, imflops;
    long long iflpops;
    int retval;
    // PERFOMANCE MEASURE
    // start PAPI counters
    if((retval=PAPI_flops_rate(PAPI_FP_OPS,&ireal_time,&iproc_time,&iflpops,&imflops)) < PAPI_OK)
      {
	printf("Could not initialise PAPI_flops \n");
	printf("Your platform may not support floating point operation event.\n");
	printf("retval: %d\n", retval);
	exit(1);
      }

    // BLOCK: Tamaño del bloque
    for (int BLOCK = 1; BLOCK <= N && BLOCK <= 4096; BLOCK = 2*BLOCK){

        for (int ii=0; ii<N_ITER; ++ii){
	    // Reinicia contadores
            if((retval=PAPI_flops_rate(PAPI_FP_OPS,&real_time, &proc_time, &flpops, &mflops))<PAPI_OK)
            {
                printf("retval: %d\n", retval);
                exit(1);
            }

	    // Codigo que se va a medir
            transpose_blocking(A, AT, N, BLOCK);

	    // Retorna valor de los contadores
            if((retval=PAPI_flops_rate(PAPI_FP_OPS,&real_time, &proc_time, &flpops, &mflops))<PAPI_OK)
            {
                printf("retval: %d\n", retval);
                exit(1);
            }

	    // Almacena en arreglos el valor de los contadores para luego procesar
	    data_time[ii] = real_time;
	    data_mflops[ii] = mflops;

        /*
	 printf("Real_time: %f Proc_time: %f flpops: %lld MFLOPS: %f\n",
               real_time, proc_time,flpops,mflops);
	*/

	    // Do something here, like computing the average of the resulting matrix, to avoid the optimizer deleting the code
	    printf("%.15e\n\n", AT[0]);

            }

	// Se calcula la media
        float mean_time = std::accumulate(data_time, data_time + N_ITER, 0.0)/(N_ITER);
        float mean_mflops = std::accumulate(data_mflops, data_mflops + N_ITER, 0.0)/(N_ITER);

        fout << BLOCK << '\t' <<  mean_time << '\t' << mean_mflops << '\n';

    }

    fout.close();

    return 0;
  }
