#include <iostream>
#include <random>
#include <cmath>
#include "Motor.h"
#include "vector_operations.h"
#include <vector>

const double Tol = 1e-10; //Tolerancia precision distancias
const double Tolv = 1e-15; // Tolerancia precision velocidades
const int D = 3;

////////////////////////////////////////////////////////////////
///Inicializacion

void initial_conditions(int seed_0,int N_particulas,
                        std::vector<Particula> &Ps, vector dbox, double m, double r, double speed0, bool flag_z, bool flag_speed){
  //inputs
  double a = dbox[0],b = dbox[1], c = dbox[2]; // dimensiones de la caja
  double max_speed = speed0;
  int nz = 0;
  int ny = 4;
  int nx = 5;
  int NZ = 10; //in case you want to distribute the positions also in z, NZ tells the program how many slices you want to have in z 
  if(flag_z == false)  ny = (N_particulas/nx)+1; 
  else{ ny = ((N_particulas/NZ)/(nx))+1; nz = (N_particulas/(nx*ny))+1;}
  
  int seed = seed_0;
  std::mt19937 gen(seed);
  /*std::uniform_real_distribution<double> random_posx(0+r, a-r);
  std::uniform_real_distribution<double> random_posy(0+r, b-r);
  std::uniform_real_distribution<double> random_posz(0+r, c-r);*/
  std::uniform_real_distribution<double> random_vel(0, max_speed);
  std::uniform_real_distribution<double> random_phi(0,2*M_PI);
  std::uniform_real_distribution<double> random_theta(0, M_PI);

  double x,y,z,vx,vy,vz,speed,phi,theta;
  vector vel0(3);
  vector pos0(3);
  int ix = 0;
  int iy = 0;
  int iz = 0;

  /*flag_speed == 0 -> gives random velocities with the random library taking the entry of the function as the maximum speed of the random distribution
    flag speed == 1 -> gives the constant velocity given in the entries of the function
  */
  
  if(flag_speed == 0){
    if(flag_z == false){ //if we do not want to distribute the particles in z
      for(int ii = 0; ii < N_particulas; ii++){
	/*the particles will be distributed in the various axis in the priorty order
	  x and then y */
	ix = (ii % nx)+1;
	iy = (ii/nx)+1;
	x = (a/(nx+2))*ix;
	y = (b/(ny+2))*iy;
	z = r;
	pos0 = {x,y,z};
	speed = random_vel(gen);
	phi = random_phi(gen);
	theta = random_theta(gen);
	vx = speed*std::cos(phi);
	vy = speed*std::sin(phi);
	vz = 0;
	vel0 = {vx,vy,vz};
	Particula P;
	P.init(pos0,vel0,m,r);
	Ps[ii] = P;
      }
    }
    
    if(flag_z == true){ //if we want to distribute the particles also in z
      for(int ii = 0; ii < N_particulas; ii++){
	/*the particles will be distributed in the various axis in the priorty order
	  x, then y, and then z*/
	ix = (ii % nx)+1;
	iy = (ii/nx % ny)+1;
	iz = (ii/(nx*ny))+1;
	x = (a/(nx+2))*ix;
	y = (b/(ny+2))*iy;
	z = (c/(nz+2))*iz;
	pos0 = {x,y,z};
	speed = random_vel(gen);
	phi = random_phi(gen);
	theta = random_theta(gen);
	vx = speed*std::sin(theta)*std::cos(phi);
	vy = speed*std::sin(phi)*std::sin(theta);
	vz = speed*std::cos(theta);
	vel0 = {vx,vy,vz};
	Particula P;
	P.init(pos0,vel0,m,r);
	Ps[ii] = P;
      }
    }
  }

  if(flag_speed == 1){
    if(flag_z == false){ //if we do not want to distribute the particles in z
      for(int ii = 0; ii < N_particulas; ii++){
	/*the particles will be distributed in the various axis in the priorty order
	  x and then y */
	ix = (ii % nx)+1;
	iy = (ii/nx)+1;
	x = (a/(nx+2))*ix;
	y = (b/(ny+2))*iy;
	z = r;
	pos0 = {x,y,z};
	speed = speed0;
	phi = random_phi(gen);
	theta = random_theta(gen);
	vx = speed*std::cos(phi);
	vy = speed*std::sin(phi);
	vz = 0;
	vel0 = {vx,vy,vz};
	Particula P;
	P.init(pos0,vel0,m,r);
	Ps[ii] = P;
      }
    }

    if(flag_z == true){ //if we want to distribute the particles also in z
      for(int ii = 0; ii < N_particulas; ii++){
	/*the particles will be distributed in the various axis in the priorty order
	  x, then y, and then z*/
	ix = (ii % nx)+1;
	iy = (ii/nx % ny)+1;
	iz = (ii/(nx*ny))+1;
	x = (a/(nx+2))*ix;
	y = (b/(ny+2))*iy;
	z = (c/(nz+2))*iz;
	pos0 = {x,y,z};
	speed = speed0;
	phi = random_phi(gen);
	theta = random_theta(gen);
	vx = speed*std::sin(theta)*std::cos(phi);
	vy = speed*std::sin(phi)*std::sin(theta);
	vz = speed*std::cos(theta);
	vel0 = {vx,vy,vz};
	Particula P;
	P.init(pos0,vel0,m,r);
	Ps[ii] = P;
      }
    }
  }

  
}


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

/*Se calcula el tiempo con la primera pared con la que chocaria (analizando la proyección), por ejemplo, viendo la velocidad en x
 y que tanto tarda la posicion en x en llegar a alguno de los dos planos x=0 o x =a, para esto ultimo se mira el signo de la velocidad en x.
Se hace lo mismo con Y y Z, y retorna el menor de estos tiempos. En caso de que la partciula esté quieta, retorna -1, pues nunca chocará con ninguna pared.
*/

double Tiempo_colision_pared(Particula P,vector dimensiones,double tiempo_actual){
  double tmin,t;
  double v_d; //velocidad en direccion d
  double r = P.radio;
  int flag = 0;
  for(int d = 0; d < D; d++){
    v_d = P.vel[d];
    if(std::abs(v_d) > Tolv){
      if(v_d < 0){
        t = (r-P.pos[d])/P.vel[d];
      } else if(v_d > 0){
        t = (dimensiones[d]-r-P.pos[d])/P.vel[d];
      }
      if(flag == 0){
        tmin = t;
        flag = 1;
      }
      else if(t<=tmin){
        tmin = t;
      }

    }
  }
  tmin += tiempo_actual;
  if(P.Rapidez()<Tolv){
    tmin = -1;
  }
  return tmin;

}


void Colision_pared(Particula &P,vector dimensiones){
  double a,b,c;
  a = dimensiones[0]; b = dimensiones[1]; c = dimensiones[2];
  double r = P.radio;
   double dis[6] = {r,a-r,r,b-r,r,c-r}; /*distancias en las que el centro de una particula de radio r está al tocar los planos {x=0,x=a,y=0,y=b,z=0,z=c}
   dir es signo(convencion usual)*1 de la dirección de la normal del plano actual (en este caso x = 0);
   Por ejemplo, si n^ = k^, entonces dir = 1, si n^ = -k^, entonces dir = -1
                                        */
   int dir = 1;
   for(int plano = 0; plano < 6; plano++){
       /*si x1 = x, x2= y, x3=z. La condicion es : si posicion de P está tocando el plano Y SU VELOCIDAD en la direccion normal del plano es hacia el plano,
      entonces ha chocado y se actualiza su velocidad
       */
     if(std::abs(P.pos[plano/2]-dis[plano])<=Tol && (P.vel[plano/2]*dir)<0){
           P.vel[plano/2] *=-1.0 ;
       }
       dir *= -1;
   }
}


double Tiempo_colision_particulas(Particula & p1, Particula & p2,double t_actual){
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

	if (Dis >= 0 && vr_posr < 0) { // No hace falta considerar la precisión
		return t_actual + ( -vr_posr - std::sqrt(Dis) ) / (mag_vr*mag_vr); // La raíz negativa corresponde al primer contacto entre las esferas. Note que vr_posr es negativa
	} else {
		return -1;
	}

}


void Colision_particulas(Particula & p1, Particula & p2){

    vector dP(3,0);
    double m1 = p1.masa,m2 = p2.masa;
    vector pr(3,0); // p1.pos-p2.pos
    vector vr(3,0); // p1.vel - p2.vel
    subs_vec(p1.pos,p2.pos,pr);
    Scalar_vec(pr,1.0/Magnitude(pr));
    subs_vec(p1.vel,p2.vel,vr);
    double q1 = -2*m2/(m1+m2)*inn_prod(vr,pr);
    Scalar_vec(pr,q1);
    sum_vec(pr,p1.vel,p1.vel);
    Scalar_vec(pr,-m1/m2);
    sum_vec(pr,p2.vel,p2.vel);

}

void colision(std::vector<Particula> &Ps,int id1,int id2, vector dbox){
    /*
    Ejecuta una colisión.
    id1 == id2 -> Partícula - Pared
    else -> Partícula - Partícula
    */
        if(id1 == id2) Colision_pared(Ps[id1],dbox);
        else Colision_particulas(Ps[id1],Ps[id2]);
}

void Evolve_system(double dt,int N,std::vector<Particula> &Ps){
    /*
    Actualiza las posiciones de las particulas al tiempo t_actual + dt.
    */
    for(int i = 0; i < N; i++){
        Ps[i].Evolucionar(dt);
    }
}
