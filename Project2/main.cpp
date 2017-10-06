#include <iostream>
#include <cstdio>
#include <cmath>
#include <armadillo>
#include <fstream>
#include <boost/timer.hpp>

using namespace std;
using namespace arma;

double func_f (double x_i);
mat LU_solver(mat f_tilde,mat A);
mat specific_solver(mat b, mat f_tilde,int n);
mat general_solver(mat a,mat b,mat c, mat f_tilde,int n);
int write_to_file(string datafil,mat f_tilde,int n,double h);

int main (int argc, char *argv[])
{
    int n = atoi(argv[1]); double x_start = 0.0; double x_slutt = 1.0;
    double h = (x_slutt-x_start)/(n+1); //x_i = i*h

    mat a = ones<vec>(n);a[0] = 0;
    mat b = (ones<vec>(n))*(-2);
    mat c = ones<vec>(n);c[n-1] = 0;
    mat f_tilde = zeros<vec>(n);

    for (int i = 0; i<n;i++){
        f_tilde[i] = -(pow(h,2))*func_f((i+1)*h);
    }

    //armadillo solves everything
    mat A = zeros<mat>(n,n);
    for (int i=0; i<n; i++) {
        A(i,i)=b[i];
        if (i!=n-1) A(i,i+1) = c[i];
        if (i!=n-1) A(i+1,i) = a[i+1];
    }
    vec ar = solve(A,f_tilde);
    //solving done

    boost::timer t_test;
    t_test.elapsed();

    boost::timer t;
    mat general_answer = general_solver(a,b,c,f_tilde,n);
    cout << "Time usage general: " << t.elapsed() << endl;
    write_to_file("General.dat",general_answer,n,h);

    boost::timer t_;
    mat specific_answer = specific_solver(b,f_tilde,n);
    cout << "Time usage specific: " << t_.elapsed() << endl;
    write_to_file("Specific.dat",specific_answer,n,h);

    boost::timer t__;
    mat LU_answer = LU_solver(f_tilde,A);
    cout << "Time usage LU: " << t__.elapsed() << endl;
    write_to_file("LU.dat",LU_answer,n,h);
    return 0 ;
}

mat general_solver(mat a,mat b,mat c, mat f_tilde,int n)
{
    //forward substitution
    for (int i = 1; i<n;i++){
        double alpha = a[i]/b[i-1];
        a[i] = a[i] - alpha*b[i-1];
        b[i] = b[i] - alpha*c[i-1];
        f_tilde[i] = f_tilde[i] - alpha*f_tilde[i-1];
    }

    //backward substitution
    for (int i = n-2; i>-1;i--){
        double c_ledd = c[i]/b[i+1];
        c[i] = c[i] - c_ledd*b[i+1];
        f_tilde[i] = f_tilde[i] - c_ledd*f_tilde[i+1];
    }

    //scaling
    for (int i = 0; i<n;i++){
        f_tilde[i] = f_tilde[i]/b[i];
        b[i] = b[i]/b[i];
    }

    return f_tilde;
}

mat specific_solver(mat b, mat f_tilde,int n)
{
    //forward substitution
    for (int i = 1; i<n;i++){
        double alpha = 1/b[i-1];
        b[i] = b[i] - alpha;
        f_tilde[i] = f_tilde[i] - alpha*f_tilde[i-1];
    }
    //backward substitution
    for (int i = n-2; i>-1;i--){
        f_tilde[i] = f_tilde[i] - (1./b[i+1])*f_tilde[i+1];
    }
    //scaling
    for (int i = 0; i<n;i++){
        f_tilde[i] = f_tilde[i]/b[i];
    }
    return f_tilde;
}

mat LU_solver(mat f_tilde, mat A){
    mat L,U;
    lu(L,U,A);
    mat Y = solve(L,f_tilde);
    mat V = solve(U,Y);
    return V;
}

int write_to_file(string datafil,mat f_tilde,int n,double h)
{
    ofstream outFile;
    outFile.open(datafil, std::ios::out);
    if (! outFile.is_open()) {
        cout << "Problem opening file." << endl;
        exit(1);
    }
    for (int i = 0; i < n; i++) {
        outFile << (i+1)*h << " " <<  f_tilde[i] << endl;
    }
    outFile.close();
    return 0;
}

double func_f (double x_i)
{
  return 100*exp(-10*x_i);
}
