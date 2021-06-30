#include "Particula.h"
#include "Motor.h"
#include "Event.h"
#include "vector_operations.h"
#include <algorithm>

void Event::update_time(std::vector<Particula> &Ps, double t_actual, vector dbox){
    if(id1 == id2) time = Tiempo_colision_pared(Ps[id1],dbox,t_actual);
    else time = Tiempo_colision_particulas(Ps[id1],Ps[id2],t_actual);
}

bool CompareByTime(const Event &a,const Event &b){
     /*
    Determina si el evento a ocurre antes que el evento b.
    Siguiendo la convención de Motor::Tiempo_choque_pared y
    Motor::Tiempo_choque_particula, toma el tiempo de choque -1
    como más grande a cualquiera.
    */
    bool order = true;
    if(a.time == -1.0){
         order = false;
    } else if(b.time != -1){
        order = (a.time) < (b.time);
     }
     return order;
}

void init_Events(std::vector<Event> &Events,int N_particulas,
                 std::vector<Particula> &Ps, vector dbox){
    /*
    Inicializa la lista de eventos Events.
    Events es un arreglo unidimensional que modela una matriz
    triangular superior. Las entradas de la diagonal (i,i)
    representan un choque con una pared; las otras entradas (i,j),
    el choque entre las particulas i y j.
    Note que solo existe la mitad superior de la "matriz"
    */
    int index;
    for(int i=0;i < N_particulas;i++){
        for(int j = i; j < N_particulas; j++){
            index = i*N_particulas + j - (i*(i+1))/2;
            Event E;
            E.id1 = i; E.id2 = j;
            E.update_time(Ps,0,dbox);
            Events[index] = E;
        }
    }
}

vector get_next_time(std::vector<Event> Events){
    /*
    Extrae el evento más próximo de la lista de eventos.
    Retona un vector con las id de las particulas que chocan
    y el tiempo en que chocan.
    */
    std::vector<Event> E = Events;
    std::sort(E.begin(), E.end(),CompareByTime);
    double id1 = (double)E[0].id1;
    double id2 = (double)E[0].id2;
    double time = E[0].time;
    vector data = {id1,id2,time};
    return data;
}

void actualizar_eventos(int id1,int id2, std::vector<Event> &Events,
                         std::vector<Particula> &Ps,int N_particulas, 
                         vector dbox)
{ 
    /*
    Actualiza la lista de eventos, después de que ocurre un choque.
    id1 e id2 indican qué particulas chocaron.
    */
  int N = 2,ii = id1;// N == 2 -> chocaron dos partículas
  int index = id1*N_particulas+id2-(id1*(id1+1))/2;
  if(id1 == id2) N = 1; // N == 1 -> chocó una partícula con una pared
  double t_actual = Events[index].time;
	
  for(int i = 0; i < N; i++){
    if(i == 1) ii = id2;
    for(int jj = 0; jj < ii; jj++){
	index = jj*N_particulas+ii-(jj*(jj+1))/2;
	Events[index].update_time(Ps,t_actual, dbox);
      }
    for(int jj = ii; jj < N_particulas; jj++){
	index = ii*N_particulas+ jj - (ii*(ii+1))/2;
	Events[index].update_time(Ps,t_actual, dbox);
      }
    }
}