#include <iostream>
#include <cstdio>
#include <cmath>
#include <armadillo>
#include <fstream>
#include <stdlib.h>
#include <boost/timer.hpp>

using namespace std;
using namespace arma;


bool test_jacobi();
vec jacobi(mat a, int n);
void jacobi_solver(mat &a,int n,int k,int l);
double biggest_func(mat a,int n, int &k, int &l);
int write_to_file(string datafil,mat f_tilde,int n,double h);





int main (int argc, char *argv[])
{

    if(!test_jacobi) {
        cout << "error" << endl;
        exit(0);
    }

    double rho_end = 5;
    //int n = atoi(argv[1]);
    int n = 60;
    double h = (rho_end)/n;


    //non interact
    mat V = zeros<mat>(n,n);
    for (int i=0; i<n; i++) {
        V(i,i)= (2/(h*h)) +(h*i*h*i);
        if (i < n-1){
          V(i+1,i) = -1/(h*h);
        }
        if (i > 0){
          V(i-1,i) = -1/(h*h);
        }
    }

    /*interact
    double w_r = 0.1;
    mat V = zeros<mat>(n,n);
    for (int i=0; i<n; i++) {
        V(i,i)= (2/(h*h)) +w_r*w_r*(h*i*h*i) + (1/(h*i));
        if (i < n-1){
          V(i+1,i) = -1/(h*h);
        }
        if (i > 0){
          V(i-1,i) = -1/(h*h);
        }
    }
    */
    cout << jacobi(V,n) << endl;

    //cout << sort(eig_sym(V)) << endl;
    // 3, 7, 11
    // omega1 non interact

    boost::timer t_test;
    t_test.elapsed();
}

vec jacobi(mat a, int n)
{
    int k;
    int l;
    double biggest = biggest_func(a,n,k,l);
    while (abs(biggest)>0.001){
        jacobi_solver(a,n,k,l);
        biggest = biggest_func(a,n,k,l);
    }
    vec eigs = sort(a.diag());
    vec smallEigs(3);
    smallEigs(0) = eigs(0);
    smallEigs(1) = eigs(1);
    smallEigs(2) = eigs(2);
    return smallEigs;
}



void jacobi_solver(mat &a,int n,int k,int l){
    double tau = (a(l,l) - a(k,k))/(2*a(k,l));
    double t_;
    if ( tau > 0 ) {
        t_ = 1.0/ (tau + sqrt(1.0 + tau*tau));
    } else {
        t_ = -1.0/ (-tau + sqrt(1.0 + tau*tau));
    }
    //double t_ = -tau + sqrt(1+(tau*tau));
    double c_ = 1/sqrt(1 + (t_*t_));//cos(theta)
    double s_ = c_*t_;//sin(theta)
    double a_kk = a(k,k);
    double a_ll = a(l,l);
    double a_kl = a(k,l);

    a(k,k) = (a_kk*c_*c_) - (2*a_kl*c_*s_) + (a_ll*s_*s_);
    a(l,l) = (a_ll*c_*c_) + (2*a_kl*c_*s_) + (a_kk*s_*s_);
    a(k,l) = ((a_kk - a_ll)*s_*c_) + (a_kl*(c_*c_ - s_*s_));
    a(l,k) = a(k,l);
    for(int i = 0;i<n;i++){
        if (i !=k && i!= l){
            double a_il = a(i,l);
            double a_ik = a(i,k);
            a(i,k) = a_ik*(c_) - a_il*(s_);
            a(i,l) = a_il*(c_) + a_ik*(s_);
            a(k,i) = a(i,k);
            a(l,i) = a(i,l);
        }
    }
}

double biggest_func(mat a,int n, int &k, int &l){
    double biggest = 0;
    for(int i = 0;i<n;i++){
        for(int j=0;j<n;j++){
            if(i!=j){
                if(abs(a(i,j)) > abs(biggest)){
                    biggest = a(i,j);
                    k = i; l = j;
                }
            }
        }
    }
    return biggest;
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


bool test_jacobi(){
    //mat a;
    //a = jacobi();
    return true;
}

