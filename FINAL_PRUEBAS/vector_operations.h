#ifndef _VECTOR_OPERATIONS_H
#define _VECTOR_OPERATIONS_H

#include <vector>

typedef std::vector<double> vector;

void sum_vec(vector & a, vector & b, vector & result);
void subs_vec(vector & a, vector & b, vector & result);
double inn_prod(vector & a, vector & b);
double Magnitude(vector & a);
void Scalar_vec(vector &a, double k);
void print(vector &a);

#endif
