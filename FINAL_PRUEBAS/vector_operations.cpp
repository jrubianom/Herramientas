#include <iostream>
#include <algorithm>
#include <vector>
#include <functional>
#include <numeric>
#include <cmath>
#include "vector_operations.h"


// Retorna a+b

void sum_vec(vector & a, vector & b, vector & result)
{
  std::transform(a.begin(), a.end(), b.begin(), result.begin() , std::plus<double>());
}

// Retorna a-b

void subs_vec(vector & a, vector & b, vector & result)
{
  std::transform(a.begin(), a.end(), b.begin(), result.begin() , std::minus<double>());

}

// Retorna a°b (producto  punto)

double inn_prod(vector & a, vector & b)
{
  double init = 0.0;
  return  std::inner_product(a.begin(), a.end(), b.begin(), init); 
}

// Retorna la magnitud del vector del argumento |a|

double Magnitude(vector & a)
{
  return std::sqrt(inn_prod(a,a));
}

//producto por escalar

void Scalar_vec(vector &a, double k){
  std::transform(a.begin(), a.end(),a.begin(),
                 [&k](auto& c){return c*k;});
}

//print imprime un vector en la consola

void print(vector &data){
    int n = data.size();
    for(int i = 0;i<n;i++){
        std::cout<<"\t"<< data[i];
    }
    std::cout<<std::endl;
}
