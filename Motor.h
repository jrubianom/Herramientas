#ifndef _MOTOR_H
#define _MOTOR_H

#include"vector_operations.h"
#include<vector>

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

void initial_conditionsp(int seed_0,int N_particulas,
                        std::vector<Particula> &Ps, vector dbox,
                        double m, double r, double speed0, bool flag_z,
                         bool flag_speed, int thid, int nth);

void initial_conditions(int seed_0,int N_particulas,
                        std::vector<Particula> &Ps, vector dbox,
                        double m, double r, double speed0,
                        bool flag_z, bool flag_speed);



double Tiempo_colision_pared(Particula P,vector dimensiones,double tiempo_actual);
void Colision_pared(Particula &P,vector dimensiones); //Actualiza la velocidad de P al chocar con la pared
double Tiempo_colision_particulas(Particula & p1, Particula & p2,double t_actual);
void Colision_particulas(Particula & p1, Particula & p2);

#endif
