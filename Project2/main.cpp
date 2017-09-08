#include <iostream>
#include <cstdio>
#include <cmath>
#include <armadillo>
#include <fstream>

using namespace std;
using namespace arma;

double func (double x_i)
{
  return 100*exp(-10*x_i);
}


int main (int argc, char *argv[])
{
    mat A = randu<mat>(5,5);
    mat B = randu<mat>(5,5);

    cout << A*B << endl;

    return 0 ;
}
