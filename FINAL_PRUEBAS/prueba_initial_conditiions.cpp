#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include "Motor.h"
#include "vector_operations.h"

int main(int argc, char**argv){

  int NP = 17;
  vector BOX(3);
  BOX = {1.0, 1.0, 1.0};
  int seed = 5;
  double M = 0.01;
  double R = 0.5;
  
  std::vector<Particula> p(NP);

  initial_conditions(seed, NP, p, BOX, M, R, 5, 0, 1);
  
  for(int ii = 0; ii < NP; ii++){
    std::cout << p[ii].pos[0] << "\t" << p[ii].pos[1] << "\t" << p[ii].pos[2] << "\n";
  }
  
  return 0;
}
