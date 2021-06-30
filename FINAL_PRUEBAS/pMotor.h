#ifndef _PMOTOR_H
#define _PMOTOR_H

#include "Particula.h"
#include "vector_operations.h"
#include <vector>

void initial_conditions(int seed_0,int N_particulas,
                        std::vector<Particula> &Ps, vector dbox,
                        double m, double r, double speed0, bool flag_z, bool flag_speed);;

double Tiempo_colision_pared(Particula P,vector dimensiones,double tiempo_actual);
void Colision_pared(Particula &P,vector dimensiones); //Actualiza la velocidad de P al chocar con la pared
double Tiempo_colision_particulas(Particula & p1, Particula & p2,double t_actual);
void Colision_particulas(Particula & p1, Particula & p2);

void colision(std::vector<Particula> &Ps,int id1,int id2, vector dbox);

void Evolve_system(double dt,int N,std::vector<Particula> &Ps);

#endif
