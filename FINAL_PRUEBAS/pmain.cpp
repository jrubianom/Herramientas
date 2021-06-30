#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <cmath>
#include "Particula.h"
#include "Motor.h"
#include "Event.h"
#include "vector_operations.h"

const vector dbox = {3,1,1};

double get_energy(std::vector<Particula> Ps,int N_particulas);

double avg_speed(std::vector<Particula> Ps,int N_particulas);

void push_histogram_data(std::vector<Particula> &Ps, int N_particulas);

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv){

    std::cout.precision(15);
    std::cout.setf(std::ios::scientific);
    
    int N_particulas = std::atoi(argv[1]),N_col = std::atoi(argv[2]);
    int N_eventos = (N_particulas*(N_particulas+1))/2;

    double m = 0.01, r = 0.01, t_actual = 0.0;

    int SEED = 3;

    double vel_inicial = 10; // Velocidad inicial con que se inicializan las particulas

    std::vector<Particula> vec_particulas(N_particulas);

    initial_conditions(SEED,N_particulas,vec_particulas,dbox, m, r, vel_inicial, 0, 1);

    std::vector <Event> Events(N_eventos);

    init_Events(Events,N_particulas,vec_particulas, dbox);

    vector current_event_info;

    double dt = 0, Next_time;
    int id1,id2;

    for(int i=0; i < N_col;i++){
      
        //Obtiene el vector {id1,id2,tiempo_proximo}
        current_event_info = get_next_time(Events);

        id1 = (int) current_event_info[0]; id2 = (int)current_event_info[1];
        Next_time = current_event_info[2];

        dt = Next_time - t_actual;

        //Evoluciona hasta el momento del choque
        Evolve_system(dt,N_particulas,vec_particulas);


        //Ejecuta el choque y actualiza las velocidades del las partículas involucradas
        colision(vec_particulas,id1,id2,dbox);


        //Actualiza tiempo de eventos
        t_actual = Next_time;
        actualizar_eventos(id1,id2,Events,vec_particulas,N_particulas,dbox);

    }
    push_histogram_data(vec_particulas, N_particulas);

    std::cout << "Velocidad final: " << avg_speed(vec_particulas, N_particulas) << '\n';
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

void push_histogram_data(std::vector<Particula> &vec_particulas, int N_particulas){
  std::ofstream vfout ("velocities_data.txt");
  for (int ii = 0; ii<N_particulas; ii++){
    vfout << vec_particulas[ii].Rapidez()<< "\n";
  }
}
