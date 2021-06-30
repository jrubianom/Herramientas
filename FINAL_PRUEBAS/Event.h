#ifndef _EVENT_H
#define _EVENT_H

#include "Particula.h"
#include "Motor.h"
#include "vector_operations.h"

struct Event
{
     /*
    Contiene el tiempo de choque entre las partículas id1 e id2.
    Si id1 == id2, el evento representa un choque con una pared.

    -----------------------------------------------------------

    Métodos:

    update_time: Actualiza el tiempo de choque entre las partículas.
    */
    double time;
    int id1,id2;
    void update_time(std::vector<Particula> &Ps, double t_actual, vector dbox);
};

bool CompareByTime(const Event &a,const Event &b);

void init_Events(std::vector<Event> &Events,int N_particulas,
                 std::vector<Particula> &Ps, vector dbox);

vector get_next_time(std::vector<Event> Events);

void actualizar_eventos(int id1,int id2, std::vector<Event> &Events,
                        std::vector<Particula> &Ps,int N_particulas,
                        vector dbox);

#endif
