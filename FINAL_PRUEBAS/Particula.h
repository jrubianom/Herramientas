#ifndef _PARTICULA_H
#define _PARTICULA_H

#include "vector_operations.h"

struct Particula
{
  vector pos {0.0, 0.0, 0.0};  //posicion
  vector vel{0.0, 0.0, 0.0};  //velocidad
  double masa = 0;
  double radio = 0;
  double Rapidez();
  double Energia();
  void Evolucionar(double dt);
  void init(vector pos0, vector vel0,double m, double r);
};

#endif