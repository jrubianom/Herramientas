#include <iostream>

void print_array(const double data[], const int & size);

int main()
{
  const int NX = 2, NY = 3, NZ = 4;
  double *x, y[NY], z[NZ];
  x = new double [NX];
  int ii, jj, kk;
  //print_array(x, NX);
  //print_array(y, NY);
 // print_array(z, NZ);
  std::cout << std::endl;

  for (ii = 0; ii < NX; ++ii) {
    x[ii] = ii;
  }
  print_array(x, NX);
 // print_array(y, NY);
 //print_array(z, NZ);
  std::cout << std::endl;

  for (jj = 0; jj < NY; ++jj) {
//    x[jj] = ii;  bound limit
    y[jj] = jj; 
  }
  print_array(x, NX);
  print_array(y, NY);
//  print_array(z, NZ);
  std::cout << std::endl;

  for (kk = 0; kk < NZ; ++kk) {
  //x[kk] = ii;
  // y[kk] = jj;
    z[kk] = kk;
  }
  print_array(x, NX); 
  print_array(y, NY); 
  print_array(z, NZ);
  std::cout << std::endl;

  return EXIT_SUCCESS;
}

void print_array(const double data[], const int & size)
{
  std::cout << std::endl;
  for (int ii = 0; ii < size; ++ii) {
    std::cout << data[ii] << "  " ;
  }
}
