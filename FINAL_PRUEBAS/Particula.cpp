#include "Particula.h"
#include "vector_operations.h"

void Particula::init(vector pos0, vector vel0,double m, double r){
  pos = pos0;
  vel = vel0;
  masa = m;
  radio = r;
}

double Particula::Energia(){
  double speed2 = inn_prod(vel,vel);
  double E = 0.5*masa*speed2;
  return E;
}

double Particula::Rapidez(){
  double speed = Magnitude(vel);
  return speed;
}


void Particula::Evolucionar(double dt){  //Evolucionar actualiza X --> X = X_0 + V*dt
  for (int ii = 0; ii<3; ii++){
    pos[ii] += vel[ii]*dt;
  }
}