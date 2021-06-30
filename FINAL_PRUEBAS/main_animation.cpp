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
#include "vector_operations.h"

const vector dbox = {10,10,10};

struct Event{
    /*
    Contiene el tiempo de choque entre las partículas id1 e id2.
    Si id1 == id2, el evento representa un choque con una pared.

    -----------------------------------------------------------

    Métodos:

    update_time: Actualiza el tiempo de choque entre las partículas.
    */
    double time;
    int id1,id2;
    void update_time(std::vector<Particula> &Ps, double t_actual){
        if(id1 == id2) time = Tiempo_colision_pared(Ps[id1],dbox,t_actual);
        else time = Tiempo_colision_particulas(Ps[id1],Ps[id2],t_actual);
    }

};

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
                 std::vector<Particula> &Ps);

void Evolve_system(double dt,int N,std::vector<Particula> &Ps);

void colision(std::vector<Particula> &Ps,int id1,int id2);

vector get_next_time(std::vector<Event> Events);


void actualizar_eventos(int id1,int id2, std::vector<Event> &Events,
                        std::vector<Particula> &Ps,int N_particulas);

double get_energy(std::vector<Particula> Ps,int N_particulas);


void push_histogram_data(std::vector<Particula> &Ps, int N_particulas);

void push_animation_data(std::vector<Particula> &vec_particulas, int N_particulas, std::ofstream& xfout, std::ofstream& yfout);

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv){
    int N_particulas = std::atoi(argv[1]),N_col = std::atoi(argv[2]);
    int N_eventos = (N_particulas*(N_particulas+1))/2;
    vector dbox = {10,10,10};
    long double m = 0.01, r = 0.1149, t_actual = 0.0;

    std::vector<Particula> vec_particulas(N_particulas);
    initial_conditions(3,N_particulas,vec_particulas,dbox, m, r, 0.1, 0, 0); //El ultimo en 0 da aleatoriedad

    std::vector <Event> Events(N_eventos),Eventos2(3);

    init_Events(Events,N_particulas,vec_particulas);

    vector current_event_info;

    long double dt = 0, Next_time;
    int id1,id2;
	
    //Inicializar lo necesario para guardar los datos de la animacion
    const long double dt_animation = 0.08;
    long double time_aux;
    int n_anim_steps = 0;;
    std::ofstream xfout ("xanimation_data.txt");
    std::ofstream yfout ("yanimation_data.txt");

    for(int ii=0; ii < N_col;ii++){
      
        //Obtiene el vector {id1,id2,tiempo_proximo}
        current_event_info = get_next_time(Events);

        id1 = (int) current_event_info[0]; id2 = (int)current_event_info[1];
        Next_time = current_event_info[2];

        dt = Next_time - t_actual;
	n_anim_steps = int (dt/dt_animation);
   
	
	//Hallar datos entre choque y choque para la animacion
	for (int kk = 0; kk < n_anim_steps; kk++){
	  time_aux = dt/n_anim_steps;
	  push_animation_data(vec_particulas, N_particulas, xfout, yfout);
	  Evolve_system(time_aux,N_particulas,vec_particulas); 	  	  
	}

	//Evolucionar el sistema un pequeño delta para mitigar los errores de redondeo
	Evolve_system(dt - (time_aux*(n_anim_steps)),N_particulas,vec_particulas);
	
	
        //Ejecuta el choque y actualiza las velocidades del las partículas involucradas
        colision(vec_particulas,id1,id2);


	//Actualiza tiempo de eventos
        t_actual = Next_time;
        actualizar_eventos(id1,id2,Events,vec_particulas,N_particulas);
	
    }
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
            E.update_time(Ps,0);
            Events[index] = E;
        }
    }
}


void Evolve_system(double dt,int N,std::vector<Particula> &Ps){
    /*
    Actualiza las posiciones de las particulas al tiempo t_actual + dt.
    */
    for(int i = 0; i < N; i++){
        Ps[i].Evolucionar(dt);
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
	Events[index].update_time(Ps,t_actual);
      }
    for(int jj = ii; jj < N_particulas; jj++){
	index = ii*N_particulas+ jj - (ii*(ii+1))/2;
	Events[index].update_time(Ps,t_actual);
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

void push_histogram_data(std::vector<Particula> &vec_particulas, int N_particulas){
  std::ofstream vfout ("velocities_data.txt");
  for (int ii = 0; ii<N_particulas; ii++){
    vfout << vec_particulas[ii].Rapidez()<< "\n";
  }
}

void push_animation_data(std::vector<Particula> &Ps, int N_particulas, std::ofstream& xfout, std::ofstream& yfout){
  for (int ii = 0; ii<N_particulas; ii++){
    xfout << Ps[ii].pos[0]<< "\t";
    yfout << Ps[ii].pos[1]<< "\t";
  }
  xfout << "\n";
  yfout << "\n";
}
