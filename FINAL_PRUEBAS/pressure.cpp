#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include "vector_operations.h"
#include "Motor.h"

double presion_z (std::vector<Particula> P);
double presion_y (std::vector<Particula> P);
double presion_x (std::vector<Particula> P);
double presion_general (std::vector<Particula> P);

int main(int argc, char**argv)
{
  presion_x(particulas);
  presion_y(particulas);
  presion_z(particulas);
  presion_general;
}


double presion_z (std::vector<Particula> P)
{
  double N = P.size();
  double G = 0.0; //Here will be made the sum of the masses*vz^2
  for(int ii = 0; ii < N; ii++){
    Pu = P[ii]
      G += P[ii].m*P[ii].vel[2]*P[ii].vel[2];
  }
  double p = G/(3.0*V);
}


double presion_y (std::vector<Particula> P)
{
  double N = P.size();
  double G = 0.0; //Here will be made the sum of the masses*vy^2
  for(int ii = 0; ii < N; ii++){
    Pu = P[ii]
      G += P[ii].m*P[ii].vel[1]*P[ii].vel[1];
  }
  double p = G/(3.0*V);
}

double presion_x (std::vector<Particula> P)
{
  double N = P.size();
  double G = 0.0; //Here will be made the sum of the masses*vx^2
  for(int ii = 0; ii < N; ii++){
    Pu = P[ii]
      G += P[ii].m*P[ii].vel[0]*P[ii].vel[0];
  }
  double p = G/(3.0*V);
}

double presion_general (std::vector<Particula> P)
{
  double N = P.size();
  double G = 0.0; //Here will be made the sum of the masses*velocity^2
  for(int ii = 0; ii < N; ii++){
    G += P[ii].m*inn_prod(P[ii].vel, P[ii].vel);
  }
  double p = G/(3.0*V);
}
