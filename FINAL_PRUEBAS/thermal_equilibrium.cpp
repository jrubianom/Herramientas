#include <iostream>
#include <cmath>
#include <vector>
#includo <queue>
#include "vector_operations.h"
#include "Motor.h"

struct evento
{
  Particula *p1, *p2;
  double tiempo;
  bool flag;
  evento(Particula & pa1, Particula & pa2, double tiempo_actual, bool fg_wall);
  void update_time(vector &dbox, double t_actual){
    if(flag == true) tiempo = Tiempo_colision_pared(*p1,dbox,t_actual);
    else tiempo = Tiempo_colision_particulas(*p1,*p2,t_actual);
  }
};

evento::evento(Particula & pa1, Particula & pa2, double tiempo_actual, bool fg_wall){
  p1=&pa1;
  p2=&pa2;
  tiempo=tiempo_actual;
  flag=fg_wall;
}

//void in_lista();
void in_eventos();
void Evolve_system(double dt,int N,std::vector<Particula*> particulas);
void colision(Event &E,vector dbox);

int main (int argc, char **argv)
{
  int Dim;
  int N_particulas;
  int N_eventos;
  double t_actual=0.0;
  double dt;
  evento E_actual;
  std::vector<double> Box(Dim);
  std::vector<Particula> prtla(N_particulas);
  std::vector<Particula*> ptr(N_particulas); 
  std::vector<evento> eventos(N_eventos);// Aqui deberia ir la lista de prioridad
  Queue lista;
  
  initial_conditions(seed, N_particulas, prtla, Box);
  in_eventos(eventos, prtla, lista, N_particulas, Box, t_actual);

  for (int ii=0; ii< Ncol; ii++){
    
    E_actual = lista.top(); lista.pop();

    dt = E_actual.time - t_actual;
    
    evolucionar_sistema(dt, N_particulas, ptr)
    //update velocities of particles in this event
    colision(E_actual,Box);
    //update priority-queue
    t_actual = E_actual.time;
    
    actualizar_evento(E_actual, t_actual, eventos);
    
    actualizar_lista();
  }
  
  return 0;
}


void in_eventos(std::vector<eventos> &Es, std::vector<Particula> & particulas, Queue Lista, int N_par, vector &Box, double t_actual )
{
  for(int i=0;i < N_par;i++){
        evento E(*particulas[i],*particulas[i],true,t_actual);
        Es.push_back(E);
        Lista.push(E);
        //std::cout<<i+1<<"-"<<i+1<<"\t"<<E.time<<std::endl;
        for(int j = i+1; j < N_par; j++){
	  Event E(*particulas[i],*particulas[j],false,t_actual);
	  E.update_time(Box, t_actual);
	  if(E.tiempo != -1){
	    Es.push_back(E);
	    Lista.push(E);
	  }
	  //std::cout<<i+1<<"-"<<j+1<<"\t"<<E.time<<std::endl;
        }
  }
}

void evolucionar_sistema(double dt,int N,std::vector<Particula*> particulas){
    for(int i = 0; i < N; i++){
        particulas[i]->Evolucionar(dt);
    }
}

void colision(Event &E,vector dbox){
    if(E.flag == true) Colision_pared(*E.p1,dbox);
    else Colision_particulas(*E.p1,*E.p2);
}


  //void in_lista(std::vector<Particula> pp, /*std::lista de prioridad*/ list, double at)
  /*{
  
  for (int ii=0; ii<pp.size(), ii++){
    for(int jj=ii+1; jj<pp.size(), jj++){
      evento par_col;
      if(Tiempo_colision_particulas(pp[ii], pp[jj],at)>0.0){
	par_col(p[ii], p[jj], Tiempo_colision_particulas(pp[ii], pp[jj],at))
	  list.push(par_col)
	  }
    }
  }
  }*/
