#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <fstream>

const double k_b =  1.380649*10e-23; //cte Boltzmann
double f_mb(double v, double m, double T);

int main(auto Particles){
  double sum_velocities = 0;
  double sum_kinetic_energy = 0;
  double avg_temperature = 0;
  double nparticles = 10;
  double single_temperature = 0;
  double avg_kinetic_energy = 0;
  double sum_mass = 0;
  double avg_mass;
  
  for (auto & val: Particles){
    sum_velocities += val.velocity;
    sum_kinetic_energy += 0.5*val.mass*val.velocity*val.velocity;
    sum_mass += val.mass;
  }
  avg_kinetic_energy = sum_kinetic_energy/nparticles;
  avg_temperature = 2*avg_kinetic_energy/(3*k_b);
  avg_mass = sum_mass/nparticles;

  //Para graficar. No es la velocidad max/min de las particulas
  double vmin = 0;
  double vmax = Particles.vel.max();
  double N = 20;
  double delta = (vmax-vmin)/N;
  double v;  
  std::ofstream fout ("theoric_values.txt"); 
  for (int ii = 0; ii<N; ii++){
    v = vmin+ii*delta;
    fout << v << "\t"
	 << f_mb(v, avg_mass, avg_temperature)<< "\n";//Save data in file theoric_values.txt
  }
  return 0;
}

double f_mb(double v, double m, double T){
  double valor = 8*pow(m/(2*M_PI*k_b*T),1.5) * exp(-m*v*v/(2*k_b*T));
  return valor;
}
