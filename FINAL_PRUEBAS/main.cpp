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

const vector dbox = {1,1,1};

double get_energy(std::vector<Particula> Ps,int N_particulas);

double avg_speed(std::vector<Particula> Ps,int N_particulas);

void push_histogram_data(std::vector<Particula> &Ps, int N_particulas, int count);

void print_params(double avg_s2, double m);

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

    double vel_inicial = std::atof(argv[3]); // Velocidad inicial con que se inicializan las particulas

    std::vector<Particula> vec_particulas(N_particulas);

    initial_conditions(SEED,N_particulas,vec_particulas,dbox, m, r, vel_inicial, 1, 1);

    std::vector <Event> Events(N_eventos);

    init_Events(Events,N_particulas,vec_particulas, dbox);

    vector current_event_info;

    double dt = 0, Next_time;
    int id1,id2;

    int count = 0;

    int i = 0;

    for(i=0; i < N_col;i++){
      
        if (i % 100 == 0){
            push_histogram_data(vec_particulas, N_particulas, i/100);
        }

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

    push_histogram_data(vec_particulas, N_particulas, i/100);

    std::cout << "Velocidad final: " << avg_speed(vec_particulas, N_particulas) << '\n';
    print_params(avg_speed(vec_particulas, N_particulas),m);
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

void print_params(double avg_s2, double m){
    double mu = 2*std::sqrt(2*avg_s2/(M_PI*3));
    double moda = std::sqrt(2*avg_s2/3);
    double T = m*avg_s2/3;
    std::cout << "La velocidad promedio final es " << mu << std::endl;
    std::cout << "La moda (el pico) es en : " << moda << std::endl;
    std::cout << "La Temperatura de equilibrio por K_b es : " << T << std::endl;
}



void push_histogram_data(std::vector<Particula> &vec_particulas, int N_particulas, int count){
  std::string name = "velocities_data" + std::to_string(count) + "-" + std::to_string(N_particulas) + ".txt";
  std::ofstream vfout (name);
  for (int ii = 0; ii<N_particulas; ii++){
    vfout << vec_particulas[ii].Rapidez()<< "\n";
  }
}
