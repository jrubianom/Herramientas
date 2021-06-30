#include <random>
#include <iostream>
#include "Motor.h"
#include "vector_operations.h"

void initial_conditions(int seed_0,int N_particulas, std::vector<Particula> &Ps, vector &dbox){
  //inputs
  double a = dbox[0],b = dbox[1], c = dbox[2]; // dimensiones de la caja
  double max_speed = 0.1;
  double min_radius = 0.001;
  double max_radius = 0.02;
  double min_mass = 0.000001;
  double max_mass = 0.001;
  
  int seed = seed_0; 
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> random_posx(0+max_radius, a-max_radius);
  std::uniform_real_distribution<double> random_posy(0+max_radius, b-max_radius);
  std::uniform_real_distribution<double> random_posz(0+max_radius, c-max_radius);
  std::uniform_real_distribution<double> random_vel(0, max_speed);
  std::uniform_real_distribution<double> random_phi(0,2*M_PI);
  std::uniform_real_distribution<double> random_theta(-0.5*M_PI, 0.5*M_PI);
  std::uniform_real_distribution<double> random_radius(min_radius,max_radius);
  std::uniform_real_distribution<double> random_mass(min_mass, max_mass);

  double vx,vy,vz,m,r,speed,phi,theta;
  vector vel0(3);
  vector pos0(3);
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  for(int ii = 0; ii < N_particulas; ii++){

    x = 

    
    pos0 = {random_posx(gen),random_posy(gen),random_posz(gen)};
    speed = random_vel(gen);
    phi = random_phi(gen);
    theta = random_theta(gen);
    vx = speed*std::sin(theta)*std::cos(phi);
    vy = speed*std::sin(phi)*std::sin(theta);
    vz = speed*std::cos(theta);
    r =  random_radius(gen);
    m =  random_mass(gen);
    vel0 = {vx,vy,vz};
    Particula P(pos0,vel0,m,r);
    Ps[ii] = P;
  }
}
