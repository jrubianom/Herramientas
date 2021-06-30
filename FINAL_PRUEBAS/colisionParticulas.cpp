
#include <iostream>
#include <vector>
#include "vector_operations.h"
#include <cmath> // ¿Podemos poner este include en vector_operations?

double Tiempo_colision_particulas(Particula & p1, Particula & p2, double t_actual);

double Tiempo_colision_particulas(Particula & p1, Particula & p2, double t_actual){
	/*
	Determina si p1 y p2 choca, y retorna el momento en que lo hacen.
	Hay dos condiciones.
	Primero, que el producto punto de la velocidad relativa y la
	posición relativa sea menor que cero; esto es, que se estén acercando.
	Segundo, que el discriminante de la ecuación |pos1 - pos2|^2 <= (r1 + r2)^2 sea
	mayor que cero; esto es, que el tiempo de colisión sea real.
	Si no colisionan retorna -1.
 	*/
	double R = p1.radio + p2.radio;
    vector pos1 = p1.pos, pos2 = p2.pos; // Posiciones iniciales
    vector v1 = p1.vel, v2 = p2.vel; // Velocidades iniciales

	vector vr(3,0); subs_vec(v2, v1, vr); // Velocidad de 2 con respecto a 1
	vector posr(3,0); subs_vec(pos2, pos1, posr); // Posición de 2 con respescto a 1

	double vr_posr = inn_prod(vr, posr); // Producto punto de la velocidad y posición relativa

   	double mag_vr = Magnitude(vr); // Magnitud de vr
   	double mag_posr = Magnitude(posr); // Magnitud de posr

   	double Dis = vr_posr*vr_posr - mag_vr*mag_vr*(mag_posr*mag_posr - R*R); // Si chocan, Dis es pequeño por desigualda de Cauchy - Schwarz

	if (Dis > 0 && vr_posr < 0) { // No hace falta considerar la precisión
		return t_actual + ( -vr_posr - std::sqrt(Dis) ) / (mag_vr*mag_vr); // La raíz negativa corresponde al primer contacto entre las esferas. Note que vr_posr es negativa
	} else {
		return -1;
	}

}

