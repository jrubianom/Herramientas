#include <iostream>
#include <algorithm>
#include <queue>
#include <deque>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <cmath>
#include "Motor.h"
#include "Event.h"
#include "vector_operations.h"
#include <omp.h>

const vector dbox = {1,1,1};



void parallel_init_particulas_y_eventos(int SEED,int N_particulas,
                                        std::vector<Particula> &vec_particulas,
                                        std::vector<Event> &Events, vector dbox,
                                        double m, double r,double vel_inicial,
                                        bool flagz, bool flag_s);

void init_Eventsp(std::vector<Event> &Events,int N_particulas,
                  std::vector<Particula> &Ps, int thid, int nth);

void init_Events(std::vector<Event> &Events,int N_particulas,
                 std::vector<Particula> &Ps);

void Evolve_system(double dt,int N,std::vector<Particula> &Ps);

void colision(std::vector<Particula> &Ps,int id1,int id2);

vector get_next_time(std::vector<Event> Events);
vector parallel_get_next_time(std::vector<Event> Events);

void actualizar_eventos(int id1,int id2, std::vector<Event> &Events,
                        std::vector<Particula> &Ps,int N_particulas);

void actualizar_eventosp(int id1,int id2, std::vector<Event> &Events,
                        std::vector<Particula> &Ps,int N_particulas);

double get_energy(std::vector<Particula> Ps,int N_particulas);

double avg_speed(std::vector<Particula> Ps,int N_particulas);
double avg_speedp(std::vector<Particula> Ps,int N_particulas);


void push_histogram_data(std::vector<Particula> &Ps, int N_particulas);

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv){

    int N_particulas = std::atoi(argv[1]),N_col = std::atoi(argv[2]);
    int N_eventos = (N_particulas*(N_particulas+1))/2;

    double m = 0.01, r = 0.01, t_actual = 0.0;

    int SEED = 3;

    double vel_inicial = 1; // Velocidad inicial con que se inicializan las particulas

    std::vector<Particula> vec_particulas(N_particulas);
    std::vector <Event> Events(N_eventos);


    //inicializacion paralela de los vectores particulas y eventos
    //(primero particulas, luego eventos)

    initial_conditions(3,N_particulas,vec_particulas,dbox, m, r, 1, 0, 0);

#pragma omp parallel
    {

    int thid = omp_get_thread_num();
    int nth = omp_get_num_threads();
    init_Eventsp(Events,N_particulas,vec_particulas,thid,nth);

    }

    /*

    parallel_init_particulas_y_eventos(SEED,N_particulas,vec_particulas,
                                       Events,dbox, m, r,
                                     vel_inicial, 0, 0);

    */
    vector current_event_info;

    double dt = 0, Next_time;
    int id1,id2;

    for(int i=0; i < N_col;i++){
      
        //Obtiene el vector {id1,id2,tiempo_proximo}
        //current_event_info = get_next_time(Events);
        current_event_info = parallel_get_next_time(Events);

        id1 = (int) current_event_info[0]; id2 = (int)current_event_info[1];
        Next_time = current_event_info[2];

        dt = Next_time - t_actual;

        //Evoluciona hasta el momento del choque
        Evolve_system(dt,N_particulas,vec_particulas);


        //Ejecuta el choque y actualiza las velocidades del las partículas involucradas
        colision(vec_particulas,id1,id2);


        //Actualiza tiempo de eventos
        t_actual = Next_time;
        actualizar_eventosp(id1,id2,Events,vec_particulas,N_particulas);

    }
    push_histogram_data(vec_particulas, N_particulas);

    std::cout << "La velocidad cuadratica promedio final es " <<
        avg_speedp(vec_particulas,N_particulas) << std::endl;
}


void init_Events(std::vector<Event> &Events,int N_particulas,
                 std::vector<Particula> &Ps){
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




void init_Eventsp(std::vector<Event> &Events,int N_particulas,
                 std::vector<Particula> &Ps,int thid,int nth){
    /*
    Inicializa la lista de eventos Events.
    Events es un arreglo unidimensional que modela una matriz
    triangular superior. Las entradas de la diagonal (i,i)
    representan un choque con una pared; las otras entradas (i,j),
    el choque entre las particulas i y j.
    Note que solo existe la mitad superior de la "matriz"
    */
    int index;
    int i0 = thid*N_particulas/nth; // indice inicial local para el thread thid
    int N_loc = (thid+1)*N_particulas/nth; // indice final local
    for(int i=i0;i < N_loc;i++){
        for(int j = i; j < N_particulas; j++){
            index = i*N_particulas + j - (i*(i+1))/2;
            Event E;
            E.id1 = i; E.id2 = j;
            E.update_time(Ps,0,dbox);
            Events[index] = E;
        }
    }
}


void Evolve_system(double dt,int N,std::vector<Particula> &Ps){
    /*
    Actualiza las posiciones de las particulas al tiempo t_actual + dt.
    */
#pragma omp parallel
    {
        int thid = omp_get_thread_num();
        int nth = omp_get_num_threads();
        int i0 = thid*N/nth;
        int N_loc = (thid+1)*N/nth;
        for(int i = i0; i < N_loc; i++){
            Ps[i].Evolucionar(dt);
        }
    }

}



void colision(std::vector<Particula> &Ps,int id1,int id2){
    /*
    Ejecuta una colisión.
    id1 == id2 -> Partícula - Pared
    else -> Partícula - Partícula
    */
        if(id1 == id2) Colision_pared(Ps[id1],dbox);
        else Colision_particulas(Ps[id1],Ps[id2]);
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


vector parallel_get_next_time(std::vector<Event> Events){

    std::vector<Event> E = Events;
    std::sort(E.begin(), E.end(),CompareByTime);
    double id1 = (double)E[0].id1;
    double id2 = (double)E[0].id2;
    double time = E[0].time;
    vector data = {id1,id2,time};
    return data;
}



void parallel_init_particulas_y_eventos(int SEED,int N_particulas,
                                        std::vector<Particula> &vec_particulas,
                                        std::vector<Event> &Events, vector dbox,
                                        double m, double r,double vel_inicial,
                                        bool flagz, bool flag_s){

#pragma omp parallel shared(SEED,N_particulas,dbox, m, r,vel_inicial, flagz, flag_s)
    {
        int thid = omp_get_thread_num();
        int nth = omp_get_num_threads();
        initial_conditionsp(SEED,N_particulas,vec_particulas,dbox, m, r,
                            vel_inicial, flagz, flag_s,thid,nth);
#pragma omp barrier

        init_Eventsp(Events,N_particulas,vec_particulas,thid,nth);

    }
}




void actualizar_eventos(int id1,int id2, std::vector<Event> &Events,
                         std::vector<Particula> &Ps,int N_particulas)
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
	Events[index].update_time(Ps,t_actual,dbox);
      }
    for(int jj = ii; jj < N_particulas; jj++){
	index = ii*N_particulas+ jj - (ii*(ii+1))/2;
	Events[index].update_time(Ps,t_actual,dbox);
      }
    }
}


double get_energy(std::vector<Particula> Ps,int N_particulas){
    // Energía total del sistema
    double E_tot = 0;
    for(int i = 0;i < N_particulas; i++){
        E_tot += Ps[i].Energia();
    }
    return E_tot;
}


double avg_speed(std::vector<Particula> Ps,int N_particulas){
    // Energía total del sistema
    double avg = 0;
    for(int i = 0;i < N_particulas; i++){
        avg += Ps[i].Rapidez()*Ps[i].Rapidez();
    }
    return avg/N_particulas;
}


double avg_speedp(std::vector<Particula> Ps,int N_particulas){
    // Energía total del sistema
    double avg = 0;
#pragma omp parallel for reduction(+:avg)
        for(int i = 0;i < N_particulas; i++){
            avg += Ps[i].Rapidez()*Ps[i].Rapidez();
        }

    return avg/N_particulas;
}

void push_histogram_data(std::vector<Particula> &vec_particulas, int N_particulas){
  std::ofstream vfout ("velocities_data.txt");
  for (int ii = 0; ii<N_particulas; ii++){
    vfout << vec_particulas[ii].Rapidez()<< "\n";
  }
}





////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

/*

void actualizar_eventosp(int id1,int id2, std::vector<Event> &Events,
                         std::vector<Particula> &Ps,int N_particulas)
{

  int N = 2;// N == 2 -> chocaron dos partículas
  int index = id1*N_particulas+id2-(id1*(id1+1))/2;
  if(id1 == id2) N = 1; // N == 1 -> chocó una partícula con una pared
  int ii;
  double t_actual = Events[index].time;

#pragma omp parallel shared(N,id1,id2,dbox,t_actual,N_particulas) private(ii)
  {
      ii = id1;
      int thid = omp_get_thread_num();
      int nth = omp_get_num_threads();
      int i0 = thid*N_particulas/nth;
      int N_loc = (thid+1)*N_particulas/nth;
      int Thidii = int( (nth*ii)/N_particulas ); // el thid que contiene al indice ii

      for(int i = 0; i < N; i++){
          if(i == 1) {
              ii = id2;
              Thidii = int( (nth*ii)/N_particulas );
          }
          if(thid != Thidii){
              if(thid < ii){
                  for(int jj = i0; jj < N_loc; jj++){
                      index = jj*N_particulas+ii-(jj*(jj+1))/2;
                      Events[index].update_time(Ps,t_actual, dbox);
                  }
              }else{
                  for(int jj = i0; jj < N_loc; jj++){
                      index = ii*N_particulas+ jj - (ii*(ii+1))/2;
                      Events[index].update_time(Ps,t_actual, dbox);
                  }
              }
          }else{
              for(int jj = i0; jj < ii; jj++){
                  index = jj*N_particulas+ii-(jj*(jj+1))/2;
                  Events[index].update_time(Ps,t_actual, dbox);
              }
              for(int jj = ii; jj < N_loc; jj++){
                  index = ii*N_particulas+ jj - (ii*(ii+1))/2;
                  Events[index].update_time(Ps,t_actual, dbox);
              }
          }

      }
  }

}

*/

void actualizar_eventosp(int id1,int id2, std::vector<Event> &Events,
                         std::vector<Particula> &Ps,int N_particulas)
{




  int N = 2,ii = id1;// N == 2 -> chocaron dos partículas
  int index = id1*N_particulas+id2-(id1*(id1+1))/2;
  if(id1 == id2) N = 1; // N == 1 -> chocó una partícula con una pared
  double t_actual = Events[index].time;

  for(int i = 0; i < N; i++){
      if(i == 1) ii = id2;
      for(int jj = 0; jj < ii; jj++){
          index = jj*N_particulas+ii-(jj*(jj+1))/2;
          Events[index].update_time(Ps,t_actual,dbox);
      }
  }

#pragma omp parallel shared(N,id1,id2,dbox,t_actual,N_particulas) private(ii)
  {
      int ii = id1;
      int thid = omp_get_thread_num();
      int nth = omp_get_num_threads();
      int i0 = thid*N_particulas/nth;
      int N_loc = (thid+1)*N_particulas/nth;
      int Thidii = int( (nth*ii)/N_particulas ); // el thid que contiene al indice ii

      for(int i = 0; i < N; i++){
          if(i == 1) {
              ii = id2;
              Thidii = int( (nth*ii)/N_particulas );
          }
          if(thid > Thidii){
              for(int jj = i0; jj < N_loc; jj++){
                  index = ii*N_particulas+ jj - (ii*(ii+1))/2;
                  Events[index].update_time(Ps,t_actual, dbox);
              }
          } else if(thid == Thidii){

              for(int jj = ii; jj < N_loc; jj++){
                  index = ii*N_particulas+ jj - (ii*(ii+1))/2;
                  Events[index].update_time(Ps,t_actual, dbox);
              }
          }

      }

  }
}
