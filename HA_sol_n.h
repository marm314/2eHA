#ifndef _HA_SOL_N_H_
#define _HA_SOL_N_H_
#include<iostream>
#include<cstring>
#include<fstream>
#include<iomanip>
#include<stdio.h>
#include<stdlib.h>
#include"Numbers.h"
#include"Mathematical_Functions.h"
#include"legendre_quadrature.h"
#include"sphere_lebedev_rule.h"
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

using namespace std;
///////////////
//Begin class//
///////////////
class HA_SOL_n
{
private:
 int n,order_rad,order_ang;
 double l,m,nx,ny,nz;
 double e_double_prime,e_prime,sqrt_omega_r;
 string quadrature_name;
 void filling_coef();
 void calculate_norm_constant();
 double u_r_over_r(double);
 double a_n(double,int,double);
 void R_and_r(double [3],double [3],double [3],double [3]);
 void grid_avail(int &);
 void create_grid();
 void clear_grid();
public:
 HA_SOL_n();
 HA_SOL_n(int,int,int,int,int,int,int,string);
 double Ntotal();
 double Intracule(double [3]);
 double Extracule(double [3]);
 double PSY12(double [3],double [3]);
 void set_order_ang(int);
 double Rho_r(double [3]);
 double Int_dm1(double [3],double [3]);
 double Int_density();
~HA_SOL_n();
 double E,omega,N;
 double *an;
};
#endif // _HA_SOL_N_H_
