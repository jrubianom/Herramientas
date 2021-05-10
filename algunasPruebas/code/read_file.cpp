//Este programa lee de un archivo datos.txt los tiempo y mflops medidos, y calcula los promedios y desviaciones estandar

#include <iostream>
#include <cmath>
#include <fstream> //liberia para manejo de archivos
#include <string>
#include <stdio.h>
#include <gsl/gsl_statistics.h> //gsl para hallar promedio y desviacion estandar

int get_number_of_lines(std::string fname); //Obtener numero de lineas en los .txt
double mean (double data[], int N); //hallar promedio de los valores de un vector
double standar_deviation_m (double data[], int N, double avrg); // hallar varianza de los valores de un vector, habiendo hallado el promedio


int main(int argc, char **argv){

    std::string file_name = "datos_mflops.txt";
    std::ifstream infile(file_name);

    int N = get_number_of_lines(file_name);

    //arreglos para los datos, real time, process time y mflops
    double real_t [N];
    //double proc_t [N];
    double mflops [N];

    int row = 0;

    // Se llenan los arreglos por filas. Primera columna real time,
    // segunda columna process time, y tercera columna mflops
    while(!infile.eof()){

        infile >> real_t[row];
	// infile >> proc_t[row];
        infile >> mflops[row];
        row += 1;
    }

    double average_time = mean (real_t,N);
    double average_mflops = mean (mflops,N);
    double sd_time = standar_deviation_m (real_t, N, average_time);
    double sd_mflops = standar_deviation_m (mflops, N, average_mflops);

    std::cout << "Tiempo promedio y desv. est." << "\n";
    std::cout << average_time<< "\t"<<sd_time<<"\n \n";
    std::cout << "Mflops promedio y desv. est." << "\n";
    std::cout << average_mflops<< "\t \t"<<sd_mflops<<"\n";

    return 0;
}

int get_number_of_lines(std::string fname){
    std::ifstream file(fname);
    int n=0;
    std::string row;
    while ( std::getline(file, row) )
        ++n;
    return n;
}

double mean (double data[], int N){
  double valor = gsl_stats_mean(data, 1, N);
  return valor;
}

double standar_deviation_m (double data[], int N, double avrg){
  double variance = gsl_stats_variance_m(data, 1, N, avrg);
  double valor = sqrt(variance);
  return valor;
}

