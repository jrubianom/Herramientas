#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch2/catch.hpp"
//#include "catch.hpp"
#include "Motor.h"
#include <cmath>

vector Test_colision_pared(vector dbox,vector x_0,
                           vector v_0,double r,int N){
    Particula P; P.init(x_0,v_0,0,r);
    double dt = 0;
    for(int i = 0; i<= N; i++){
        P.Evolucionar(dt);
        Colision_pared(P,dbox);
        dt = Tiempo_colision_pared(P,dbox,0);
    }

    return P.pos;

}


double Test_tiempo_colision_particulas(vector pos1, vector pos2, vector vel1,
                                       vector vel2, double r, double t_act)
{
	Particula P1; P1.init(pos1,vel1,0,r);
	Particula P2; P2.init(pos2,vel2,0,r);

	return Tiempo_colision_particulas(P1, P2, t_act);
}

//Test 1, partcula va de esquina a esquina del cubo

double r;
vector dbox,x_0,v_0,x_f;
int N;
double d; //variable util, en ocasiones será una dimension de la caja menos el radio

TEST_CASE( "Particula-pared, esquina-esquina" ) {

    r = 1.5;
    dbox = {30,30,30};
    x_0 = {r,r,30-r};
    v_0 = {4,4,-4};
    N = 16;
    d = dbox[0]-r;
    x_f = Test_colision_pared(dbox,x_0,v_0,r,N);
    REQUIRE(x_f[0] == r);
    REQUIRE(x_f[1] == r);
    REQUIRE(x_f[2] == d);
}


//Test 2, Particula va en direecion y, golpeando los planos y=0, y= b

TEST_CASE( "Particula-pared,velocidad solo en direccion -y" ) {

    r = 4.5;
    dbox = {10.2,13.82,15};
    x_0 = {4.6,4.9,9.0};
    v_0 = {0,-1.1,0};
    N = 19;
    d = dbox[1]-r;
    x_f = Test_colision_pared(dbox,x_0,v_0,r,N);
    REQUIRE(x_f[0] == x_0[0]);
    REQUIRE(x_f[1] == r);
    REQUIRE(x_f[2] == x_0[2]);
    N = 32;
    x_f = Test_colision_pared(dbox,x_0,v_0,r,N);
    REQUIRE(x_f[0] == x_0[0]);
    REQUIRE(x_f[1] == d);
    REQUIRE(x_f[2] == x_0[2]);
}

//Test 3, Particula va en direecion x, golpeando los planos x=0, x= a

TEST_CASE( "Particula-pared,velocidad solo en direccion -x" ) {

    r = 4.5;
    dbox = {10.2,13.82,15};
    x_0 = {4.6,4.9,9.0};
    v_0 = {-1.1,0,0};
    N = 13;
    d = dbox[0]-r;
    x_f = Test_colision_pared(dbox,x_0,v_0,r,N);
    REQUIRE(x_f[0] == r);
    REQUIRE(x_f[1] == x_0[1]);
    REQUIRE(x_f[2] == x_0[2]);
    N = 42;
    x_f = Test_colision_pared(dbox,x_0,v_0,r,N);
    REQUIRE(x_f[0] == d);
    REQUIRE(x_f[1] == x_0[1]);
    REQUIRE(x_f[2] == x_0[2]);
}


//Test 4, Particula va en direecion z, golpeando los planos z=0, z= c

TEST_CASE( "Particula-pared,velocidad solo en direccion +z" ) {

    r = 4.5;
    dbox = {10.2,13.82,15};
    x_0 = {4.6,4.9,9.0};
    v_0 = {0,0,1.1};
    N = 50;
    d = dbox[2]-r;
    x_f = Test_colision_pared(dbox,x_0,v_0,r,N);
    REQUIRE(x_f[0] == x_0[0]);
    REQUIRE(x_f[1] == x_0[1]);
    REQUIRE(x_f[2] == r);
    N = 71;
    x_f = Test_colision_pared(dbox,x_0,v_0,r,N);
    REQUIRE(x_f[0] == x_0[0]);
    REQUIRE(x_f[1] == x_0[1]);
    REQUIRE(x_f[2] == d);
}



//Test 5, Partcula con centro en el plano z = r,velocidad solo en
// X-Y, rebota en posiciones centrales

TEST_CASE("P-Pared, velocidad paralela plano XY,trayetoria cuadrada ") {
    r = 2.3;
    dbox = {18.0,18.0,15};
    double med = dbox[0]/2; //mitad de la longitud en x ( o y, es lo mismo)
    x_0 = {r,med,3.0};
    v_0 = {1.2,1.2,0};
    d = 18.0-r;
    N = 20;
    x_f = Test_colision_pared(dbox,x_0,v_0,r,N);
    REQUIRE(x_f[0] == r);
    REQUIRE(x_f[1] == med);
    REQUIRE(x_f[2] == x_0[2]);
    N = 21;
    x_f = Test_colision_pared(dbox,x_0,v_0,r,N);
    REQUIRE(x_f[0] == med);
    REQUIRE(x_f[1] == d);
    REQUIRE(x_f[2] == x_0[2]);
    N = 22;
    x_f = Test_colision_pared(dbox,x_0,v_0,r,N);
    REQUIRE(x_f[0] == d);
    REQUIRE(x_f[1] == med);
    REQUIRE(x_f[2] == x_0[2]);
    N = 23;
    x_f = Test_colision_pared(dbox,x_0,v_0,r,N);
    REQUIRE(x_f[0] == med);
    REQUIRE(x_f[1] == r);
    REQUIRE(x_f[2] == x_0[2]);
}

//Test 6, particula con velocidad 0

TEST_CASE() {
    r = 0.006;
    dbox = {1,0.5,0.9};
    x_0 = {0.3,0.02,0.03};
    v_0 = {0.0,0.0,0.0};
    N = 19;
    x_f = Test_colision_pared(dbox,x_0,v_0,r,N);
    REQUIRE(x_f[0] == x_0[0]);
    REQUIRE(x_f[1] == x_0[1]);
    REQUIRE(x_f[2] == x_0[2]);
}

//Test 7, tiempo de colisión para particulas

TEST_CASE("Particula-Particula, velocidad solo en x") {
	double r = 1;
	vector pos1 = {0,0,0}, pos2 = {10,0,0};
	vector vel1 = {1,0,0}, vel2 = {-1,0,0};
	REQUIRE(std::fabs(Test_tiempo_colision_particulas(pos1, pos2, vel1, vel2, r, 0) - 4) < 1e-3);
}

TEST_CASE("Particula-Particula, velocidad solo en y") {
	double r = 1;
	vector pos1 = {0,0,0}, pos2 = {0,10,0};
	vector vel1 = {0,1,0}, vel2 = {0,-1,0};
	REQUIRE(std::fabs(Test_tiempo_colision_particulas(pos1, pos2, vel1, vel2, r, 0) - 4) < 1e-3);
}


TEST_CASE("Particula-Particula, velocidad solo en z") {
	double r = 1;
	vector pos1 = {0,0,0}, pos2 = {0,0,10};
	vector vel1 = {0,0,1}, vel2 = {0,0,-1};
	REQUIRE(std::fabs(Test_tiempo_colision_particulas(pos1, pos2, vel1, vel2, r, 0) - 4) < 1e-3);
}

// Velocidades no paralelas

TEST_CASE("Particula-Particula, vectores velocidad oblicuos") {
    double r = 1;
    vector pos1 = {0,0,0}, pos2 = {10,0,0};
    vector vel1 = {1,-1,0}, vel2 = {-1,-1,0};
    REQUIRE(std::fabs(Test_tiempo_colision_particulas(pos1, pos2, vel1, vel2, r, 0) - 4) < 1e-3);
}

TEST_CASE("Particula-Particula, vectores velocidad perpendiculares") {
    double r = 1;
    vector pos1 = {0,4,0}, pos2 = {4,0,0};
    vector vel1 = {1,0,0}, vel2 = {0,0.5,0};
    REQUIRE(std::fabs(Test_tiempo_colision_particulas(pos1, pos2, vel1, vel2, r, 0) - 4) < 1e-3);
}
