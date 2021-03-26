#include <cmath>
#include <iostream>

double sin_pade(double x);
double sin_aux(double u);

int main(void) {
  std::cout.precision(16);              // 16 decimal places
  std::cout.setf(std::ios::scientific); // use scientific notation

  double xmin =0.01, xmax = 2*M_PI,delta = 0.01;
  double Dif,exact;
  for(double x=xmin;x<=xmax;x+=delta){
    exact = std::sin(x);
    Dif = std::abs((sin_pade(x)-exact))/exact;
    std::cout<<x<<"\t"<<Dif<<"\n";
  }

  return 0;
}
  double sin_pade(double x) {
      double xmod,senmod,sen,residuo,xabs,u;
      xabs = std::abs(x);
        xmod = std::fmod(xabs,2*M_PI); // residuo con valor absoluto de x, se usará que sen es impar
        if(M_PI<=xmod && xmod<=3*M_PI/2){
          xmod = 3*M_PI-xmod;
        }
        if(3*M_PI/2<=xmod && xmod<=2*M_PI){
          xmod = xmod-2*M_PI;
          if(std::abs(xmod)<=M_PI/6){
            sen = sin_aux(xmod);
          }
          if(std::abs(xmod) >M_PI/6){
            u = xmod/3;
            sen = (3-4*sin_aux(u)*sin_aux(u))*sin_aux(u);
          }
        }
        if(M_PI/2<=xmod && xmod<M_PI){
          xmod = M_PI-xmod;
        }
        if(0<=xmod && xmod<=M_PI/2){
          if(xmod<=M_PI/6){
            sen = sin_aux(xmod);
          }
          if(xmod>M_PI/6){
            u = xmod/3;
            sen = (3-4*sin_aux(u)*sin_aux(u))*sin_aux(u);
          }
        }
        if(std::abs(xmod)<1e-8){
          sen = xmod;
        }
        if(x<0){
          sen = -sen;
        }
      return sen;

    }

double sin_aux(double u) {
  double nom,den,s;
  nom = 1-(29593.0/207636.0)*u*u+(34911.0/7613320.0)*std::pow(u,4)-(479249.0/11511339840.0)*std::pow(u,6);
  den = 1+(1671.0/69212.0)*u*u+(97.0/351384.0)*std::pow(u,4)+(2623.0/1644477120.0)*std::pow(u,6);
  return u*nom/den;
}
