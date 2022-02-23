/////////////////////////////////////////////////////////////////////////////////////////////////////
// 2eHAan                                                                                          //
// ------                                                                                          //
// Analytic solutions of 2e- HA. (for both states, i.e., singlet and triplet for the lowest omega  //
// values).                                                                                        //
// 1) Finds the value of omega for which the power series is finite.                               //
// 2) Builds the coefficients for the finite power series solution.                                //
// 3) Finds the normalization constant for the wave function and builds the wave function.         //
// 4) Calculates several properties from analytic solutions.                                       //
// writen by: M.Sc. Mauricio Rodríguez M. for support and comments: marm314@gmail.com              //
/////////////////////////////////////////////////////////////////////////////////////////////////////
#include<iostream>
#include<string>
#include<locale>
#include<sstream>
#include<fstream>
#include<iomanip>
#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include"cuba.h"
#include"Numbers.h"
#include"Mathematical_Functions.h"
#include"Input_commands.h"
#include"HA_sol_n.h"
#include"D_read_calc_rho.h"
////////////////////////////////////////////////////////////////////////////////////////////////////
#define LAST 4
#define SEED 0
#define NVEC1 1
#define PCORE 10000
#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
#define NNEW 1000
#define FLATNESS 25.
#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0
#define KEY 0
#define SPIN NULL
////////////////////////////////////////////////////////////////////////////////////////////////////
using namespace std;
////////////////////////////////////////////////////////////////////////////////////////////////////
// Global Variables needed for cuba intervals and FCHK_WFN DIs                                    //
///////////////////////////////////////////////////////////////////////////////////////////////////7
int quant_3dho[3],quant_3dho2[3];
double R_extra,r_intra,r_min,r_max,theta_min,theta_max,phi_min,phi_max;
double rr2_min,rr2_max,thetat2_min,thetat2_max,phip2_min,phip2_max;
READ_FCHK_WFN *Read_fchk_wfn;
////////////////////////////////////////////////////////////////////////////////////////////////////
int Int_rho(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata);
int Int_rho2(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata);
int Int_rho3(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata);
int Int_2part_rho(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata);
int Int_2part_rho_fchk_wfn(const int *ndim, const double xx[],const int *ncomp,
double ff[],void *userdata);
int Int_intra(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata);
int Int_extra(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata);
int Int_vee(const int *ndim, const double xx[],const int *ncomp, double ff[],
void *userdata);
int Int_coef_fourier_wfn(const int *ndim, const double xx[],const int *ncomp,
double ff[],void *userdata);
double Int_intracule_sph(HA_SOL_n &ha_sol_n);
double Int_extracule_sph(HA_SOL_n &ha_sol_n);
double Int_Vee(HA_SOL_n &ha_sol_n);
double Int_density(HA_SOL_n &ha_sol_n);
double Int_density2(HA_SOL_n &ha_sol_n,double R_min,double R_max,double Theta_min,
double Theta_max,double Phi_min,double Phi_max);
double Int_density_fchk_wfn(double R_min,double R_max,double Theta_min,
double Theta_max,double Phi_min,double Phi_max);
double Int_dm2(HA_SOL_n &ha_sol_n,double R_min,double R_max,double Theta_min,double Theta_max,
double Phi_min,double Phi_max,double R2_min,double R2_max,double Theta2_min,double Theta2_max,
double Phi2_min,double Phi2_max);
double Int_dm2_fchk_wfn(HA_SOL_n &ha_sol_n,double R_min,double R_max,double Theta_min,
double Theta_max,double Phi_min,double Phi_max,double R2_min,double R2_max,
double Theta2_min,double Theta2_max,double Phi2_min,double Phi2_max);
double Int_wfn_apx(HA_SOL_n &ha_sol_n, int &ncores,int &fail);
void Quant_fill(int **Quant,int styp);
////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
 int i;
 double density,eval,r[3],INTRACULE[3]={ZERO};
 double rho_1,rho_2,rho_12,integral,intracule_scanned,Cnorm;
 ofstream results;
 if(argc==2)
 {
  string name_in(argv[1]);
  string name;
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  // Read input information                                                                        //
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  Input_commands IN(name_in);
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  // Set results file name                                                                         //
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  name="2eHA_";
  ostringstream convert;
  convert<<IN.n;
  name=name+convert.str();
  name=name+"_";
  ostringstream convert1;
  convert1<<IN.l;
  name=name+convert1.str();
  name=name+"_";
  ostringstream convert2;
  convert2<<IN.m;
  name=name+convert2.str()+".txt";
  results.open((name).c_str());
  HA_SOL_n ha_sol_n(IN.n,IN.l,IN.m,IN.nx,IN.ny,IN.nz,IN.order_rad,name);
  if(IN.order_ang!=0 || (IN.DI_an || IN.density || IN.dm1 || IN.scan_density)){ha_sol_n.set_order_ang(IN.order_ang);}
  results<<"##########################################################################"<<endl;
  results<<"##########################################################################"<<endl;
  results<<"# 2e- HA Program for Analytic Evaluations in C++                         #"<<endl;
  results<<"# written by M.Sc. Mauricio Rodríguez M email: marm3.14@gmail.com        #"<<endl;
  results<<"##########################################################################"<<endl;
  results<<"##########################################################################"<<endl;
  results<<endl;
  results<<endl;
  results<<"#*****************************************************************************#";
  results<<endl;
  results<<"# Copyright (C) 2014 M.Sc. Mauricio A. Rodríguez Mayorga                      #";
  results<<endl;
  results<<"# Ph.D. student at University of the Basque Country (Donostia) EHU/UPV        #";
  results<<endl;
  results<<"# for support and comments send an email to: marm314@gmail.com                #";
  results<<endl;
  results<<"#*****************************************************************************#";
  results<<endl;
  results<<"#  This program is free software: you can redistribute it and/or modify       #";
  results<<endl;
  results<<"#  it under the terms of the GNU General Public License as published by       #";
  results<<endl;
  results<<"#  the Free Software Foundation, either version 3 of the License, or          #";
  results<<endl;
  results<<"#  (at your option) any later version.                                        #";
  results<<endl;
  results<<"#  This program is distributed in the hope that it will be useful,            #";
  results<<endl;
  results<<"#  but WITHOUT ANY WARRANTY; without even the implied warranty of             #";
  results<<endl;
  results<<"#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #";
  results<<endl;
  results<<"#  GNU General Public License for more details.                               #";
  results<<endl;
  results<<"#  You should have received a copy of the GNU General Public License          #";
  results<<endl;
  results<<"#  along with this program.  If not, see <http://www.gnu.org/licenses/>.      #";
  results<<endl;
  results<<"#*****************************************************************************#";
  results<<endl;
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  // Evaluate the intracule                                                                        //
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  if(IN.intracule)
  {
   //Evaluates intracule density and properties
   results<<endl;
   results<<setprecision(5)<<fixed<<"x_intra y_intra z_intra |Phi(r12)|^2: "<<endl;
   results<<IN.r_intra[0]<<" "<<IN.r_intra[1]<<" "<<IN.r_intra[2]<<" ";
   eval=ha_sol_n.Intracule(IN.r_intra);
   results<<setprecision(10)<<fixed<<scientific<<eval*eval<<endl;
   r_intra=norm3D(IN.r_intra);
   results<<setprecision(5)<<fixed<<"|r_12| = u: "<<r_intra<<endl;
   results<<"(u = sqrt[ x_intra^2 + y_intra^2 + z_intra^2 ])"<<endl;
   eval=Int_intracule_sph(ha_sol_n);
   results<<setprecision(10)<<fixed<<scientific<<"h(u): "<<eval<<endl;
   results<<"(h(u) = Int |Phi(u sin(theta) cos(phi), u sin(theta) sin(phi), u cos(theta))|^2 sin(theta) dtheta dphi)"<<endl;
   results<<setprecision(10)<<fixed<<scientific<<"h(u) u^2 : "<<eval*r_intra*r_intra<<endl;
   eval=Int_Vee(ha_sol_n);
   results<<setprecision(10)<<fixed<<scientific<<"Vee: "<<eval<<endl;
   results<<"(Vee = Int (1/u) |Phi(u sin(theta) cos(phi), u sin(theta) sin(phi), u cos(theta))|^2 sin(theta) u^2 du dtheta dphi)"<<endl;
   results<<endl;
   results<<"##########################################################################"<<endl;
  }
  if(IN.scan_intracule)
  {
   //Scan the intracule density and properties
   results<<endl;
   results<<"# x_intra    y_intra     z_intra      u=|r_12|      |Phi(r12)|^2      h(u)             h(u) u^2 : "<<endl;
   integral=ZERO;
   for(intracule_scanned=ZERO;intracule_scanned<IN.last_intra;intracule_scanned=intracule_scanned+IN.step_intra)
   {
    if(IN.axis_intra=='X' || IN.axis_intra=='x')
    {
     INTRACULE[0]=intracule_scanned;
    }
    else if(IN.axis_intra=='Y' || IN.axis_intra=='y')
    {
     INTRACULE[1]=intracule_scanned;
    }
    else
    {
     INTRACULE[2]=intracule_scanned;
    }
    r_intra=norm3D(INTRACULE);
    results<<setprecision(5)<<fixed<<scientific<<INTRACULE[0]<<" "<<INTRACULE[1]<<" "<<INTRACULE[2]<<" "<<r_intra<<" ";
    eval=ha_sol_n.Intracule(INTRACULE);
    results<<setprecision(10)<<fixed<<scientific<<eval*eval<<" ";
    eval=Int_intracule_sph(ha_sol_n);
    eval=eval;
    results<<eval<<"  ";
    results<<eval*r_intra*r_intra<<" ";
    integral=integral+eval*r_intra*r_intra*IN.step_intra;
    results<<endl;
   }
   eval=Int_Vee(ha_sol_n);
   results<<"# "<<setprecision(10)<<fixed<<scientific<<endl;
   results<<"# Vee = Int (1/u) |Phi(u sin(theta) cos(phi), u sin(theta) sin(phi), u cos(theta))|^2 sin(theta) u^2 du dtheta dphi "<<eval<<endl;
   results<<"# Int h(u) u^2 du : "<<setprecision(5)<<fixed<<integral<<endl;
   results<<"##########################################################################"<<endl;
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  // Evaluate the extracule                                                                        //
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  if(IN.extracule)
  {
   //Evaluates extracule function (Not density)
   eval=ha_sol_n.Extracule(IN.r_extra);
   results<<endl;
   results<<setprecision(5)<<fixed<<"X_extra Y_extra Z_extra |Xi(R12)|^2: "<<endl;
   results<<IN.r_extra[0]<<" "<<IN.r_extra[1]<<" "<<IN.r_extra[2]<<" ";
   results<<setprecision(10)<<fixed<<scientific<<eval*eval<<endl;
   R_extra=norm3D(IN.r_extra);
   results<<setprecision(5)<<fixed<<"|R_12| = U: "<<R_extra<<endl;
   results<<"(U = sqrt[ X_extra^2 + Y_extra^2 + Z_extra^2 ])"<<endl;
   eval=Int_extracule_sph(ha_sol_n);
   eval=eval;
   results<<setprecision(10)<<fixed<<scientific<<"H(U): "<<eval<<endl;
   results<<"(H(U) = Int |Xi(U sin(theta) cos(phi), U sin(theta) sin(phi), U cos(theta))|^2 sin(Theta) dtheta dphi)"<<endl;
   results<<setprecision(10)<<fixed<<scientific<<"H(U) U^2 : "<<eval*R_extra*R_extra<<endl;
   results<<endl;
   results<<"##########################################################################"<<endl;
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  // Evaluate the wave function                                                                    //
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  if(IN.wave_function)
  {
   eval=ha_sol_n.PSY12(IN.r1_wf,IN.r2_wf);
   results<<endl;
   results<<setprecision(5)<<fixed<<"PSY(r1,r2): "<<endl;
   results<<IN.r1_wf[0]<<" "<<IN.r1_wf[1]<<" "<<IN.r1_wf[2]<<" "<<IN.r2_wf[0]<<" "<<IN.r2_wf[1]<<" "<<IN.r2_wf[2]<<" ";
   //Sqrt(TWO) To make a density that integrates to TWO
   results<<setprecision(10)<<fixed<<scientific<<sqrt(TWO)*eval<<endl;
   results<<endl;
   results<<"##########################################################################"<<endl;
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  // Evaluate the density                                                                          //
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  if(IN.density)
  {
   density=ha_sol_n.Rho_r(IN.r_dens);
   results<<endl;
   results<<setprecision(5)<<fixed<<"Density at: "<<IN.r_dens[0]<<"\t"<<IN.r_dens[1]<<"\t"<<IN.r_dens[2]<<"\t\t";
   //The number TWO is here to make the density integrate to TWO.
   results<<setprecision(10)<<fixed<<scientific<<TWO*density<<endl;
   results<<endl;
   results<<"##########################################################################"<<endl;
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  // Scan options                                                                                  //
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  if(IN.scan_density)
  {
   r[0]=ZERO;r[1]=ZERO,r[2]=ZERO;
   density=ha_sol_n.Rho_r(r);
   results<<"# Scanned density along "<<IN.axis<<"-axis"<<endl;
   results<<"# x \t\t y \t\t z \t\t Density"<<endl;
   results<<setprecision(5)<<fixed<<r[0]<<"\t\t"<<r[1]<<"\t\t"<<r[2]<<"\t\t";
   results<<setprecision(10)<<fixed<<scientific<<TWO*density<<endl;
   if(IN.axis=='x' ||IN.axis=='X')
   {i=0;}
   else if(IN.axis=='y' ||IN.axis=='Y')
   {i=1;}
   else
   {i=2;}
   for(r[i]=IN.step;r[i]<IN.last+IN.step;r[i]=r[i]+IN.step)
   {
    density=ha_sol_n.Rho_r(r);
    results<<setprecision(5)<<fixed<<r[0]<<"\t\t"<<r[1]<<"\t\t"<<r[2]<<"\t\t";
    results<<setprecision(10)<<fixed<<scientific<<TWO*density<<endl;
   }
   results<<endl;
   results<<"##########################################################################"<<endl;
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  // Integration options                                                                           //
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  if(IN.integrate_density)
  {
   results<<endl;
   eval=Int_density(ha_sol_n); //Result must be 1 and we make it 2
   //The number TWO is here to make the density integrate to TWO.
   results<<setprecision(10)<<fixed<<scientific<<"Integration of the density: "<<TWO*eval<<endl;
   results<<endl;
   results<<"##########################################################################"<<endl;
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  // DMN options                                                                                   //
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  if(IN.dmn)
  {
   results<<endl;
   if(IN.dm1)
   {
    eval=ha_sol_n.Int_dm1(IN.r1_dm1,IN.r1_prime_dm1);
    results<<"DM1"<<endl;
    results<<"          r1                        r1'          DM1 RHO(r1,r1'): "<<endl;
    //TWO to integrate the density to TWO
    results<<setprecision(5)<<fixed<<IN.r1_dm1[0]<<" "<<IN.r1_dm1[1]<<" "<<IN.r1_dm1[2]<<" "<<IN.r1_prime_dm1[0]<<" "<<IN.r1_prime_dm1[1]<<" "<<IN.r1_prime_dm1[2]<<" ";
    results<<setprecision(10)<<fixed<<scientific<<TWO*eval<<endl;
    results<<endl;
    results<<"##########################################################################"<<endl;
   }
   if(IN.dm2)
   {
    results<<"DM2"<<endl;
    results<<"          r1                        r1': "<<endl;
    results<<setprecision(5)<<fixed<<IN.r1_dm2[0]<<" "<<IN.r1_dm2[1]<<" "<<IN.r1_dm2[2]<<" "<<IN.r1_prime_dm2[0]<<" "<<IN.r1_prime_dm2[1]<<" "<<IN.r1_prime_dm2[2]<<endl;
    results<<"          r2                        r2': "<<endl;
    results<<setprecision(5)<<fixed<<IN.r2_dm2[0]<<" "<<IN.r2_dm2[1]<<" "<<IN.r2_dm2[2]<<" "<<IN.r2_prime_dm2[0]<<" "<<IN.r2_prime_dm2[1]<<" "<<IN.r2_prime_dm2[2]<<endl;
    results<<"Value of the DM2 at (r1,r2,r1',r2'): ";
    //TWO to make the density integrate to TWO.
    results<<setprecision(10)<<fixed<<scientific<<TWO*ha_sol_n.PSY12(IN.r1_dm2,IN.r2_dm2)*ha_sol_n.PSY12(IN.r1_prime_dm2,IN.r2_prime_dm2)<<endl;
    results<<endl;
    results<<"##########################################################################"<<endl;
   }
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  // Exact DI d(A,B) (ESI)                                                                         //
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  if(IN.DI_an)
  {
   results<<endl;
   results<<"##########################################################################"<<endl;
   results<<endl;
   results<<"DI calculated taking "<<IN.n<<" terms"<<endl;
   convert_angles(IN.theta1_min_an,IN.theta1_max_an,IN.phi1_min_an,IN.phi1_max_an);
   convert_angles(IN.theta2_min_an,IN.theta2_max_an,IN.phi2_min_an,IN.phi2_max_an);
   //d(A,B)=2 [int_A rho(1) d1 int_B rho(2) d2 - int_A int_B DM2(2,1) d2d1]
   //The M. Taut wavefunction is noralized to 1 and that is why we have some
   //numbers "2" for rho_1, rho_2  and rho_12.
   rho_1=TWO*Int_density2(ha_sol_n,IN.r1_min_an,IN.r1_max_an,IN.theta1_min_an,IN.theta1_max_an,IN.phi1_min_an,IN.phi1_max_an);
   results<<setprecision(10)<<fixed<<scientific<<" N(A)  : "<<rho_1<<endl;
   rho_2=TWO*Int_density2(ha_sol_n,IN.r2_min_an,IN.r2_max_an,IN.theta2_min_an,IN.theta2_max_an,IN.phi2_min_an,IN.phi2_max_an);
   results<<setprecision(10)<<fixed<<scientific<<" N(B)  : "<<rho_2<<endl;
   rho_12=TWO*Int_dm2(ha_sol_n,IN.r1_min_an,IN.r1_max_an,IN.theta1_min_an,IN.theta1_max_an,IN.phi1_min_an,IN.phi1_max_an,
        IN.r2_min_an,IN.r2_max_an,IN.theta2_min_an,IN.theta2_max_an,IN.phi2_min_an,IN.phi2_max_an);
   results<<setprecision(10)<<fixed<<scientific<<" N(A,B): "<<rho_12<<endl;
   eval=TWO*((rho_1*rho_2)-rho_12);
   results<<setprecision(10)<<fixed<<scientific<<" DI    : "<<eval<<endl;
   results<<endl;
   results<<"##########################################################################"<<endl;
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // DI using an FCHK or WFN file and some analytic DM2 and Rho(r) [generating a pair corr. function g(1,2) //
  // DM2(1,2)=Rho(1)*Rho(2)*g(1,2)] and using an FCHK or WFN to evaluate the densities [Rho(1) and Rho(2)]. //                                                                     //
  // DI=d(A,B)=2 [int_A rho(1) d1 int_B rho(2) d2 - int_A int_B DM2(2,1) d2d1]                              //
  //   = 2 [int_A rho(1) d1 int_B rho(2) d2 - int_A int_B Rho(1)*Rho(2)*g(1,2) d2d1]                        //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if(IN.DI_fchk_wfn)
  {
   Read_fchk_wfn=new READ_FCHK_WFN(IN.name_fchk_wfn,IN.name_log,IN.WFN,IN.log_file,IN.cas,IN.mult_in);
   results<<endl;
   results<<"##########################################################################"<<endl;
   results<<endl;
   results<<"DI calculated taking "<<IN.n<<" terms and "<<IN.name_fchk_wfn<<" file"<<endl;
   convert_angles(IN.theta1_min,IN.theta1_max,IN.phi1_min,IN.phi1_max);
   convert_angles(IN.theta2_min,IN.theta2_max,IN.phi2_min,IN.phi2_max);
   rho_1=Int_density_fchk_wfn(IN.r1_min,IN.r1_max,IN.theta1_min,IN.theta1_max,IN.phi1_min,IN.phi1_max);
   results<<setprecision(10)<<fixed<<scientific<<" N(A)  : "<<rho_1<<endl;
   rho_2=Int_density_fchk_wfn(IN.r2_min,IN.r2_max,IN.theta2_min,IN.theta2_max,IN.phi2_min,IN.phi2_max);
   results<<setprecision(10)<<fixed<<scientific<<" N(B)  : "<<rho_2<<endl;
   Cnorm=ONE;
   if(IN.cnormdi)
   {
    Cnorm=Int_dm2_fchk_wfn(ha_sol_n,0,1e99,0,PI,0,TWO*PI,0,1e99,0,PI,0,TWO*PI);
    Cnorm=ONE/sqrt(Cnorm);
    results<<" Cnorm : "<<Cnorm<<endl;
   }
   rho_12=((double)Read_fchk_wfn[0].nelectrons*((double)Read_fchk_wfn[0].nelectrons-ONE))
         *pow(Cnorm,TWO)*Int_dm2_fchk_wfn(ha_sol_n,
          IN.r1_min,IN.r1_max,IN.theta1_min,IN.theta1_max,IN.phi1_min,IN.phi1_max,
          IN.r2_min,IN.r2_max,IN.theta2_min,IN.theta2_max,IN.phi2_min,IN.phi2_max);
   results<<setprecision(10)<<fixed<<scientific<<" N(A,B): "<<rho_12<<endl;
   eval=TWO*((rho_1*rho_2)-rho_12);
   results<<setprecision(10)<<fixed<<scientific<<" DI    : "<<eval<<endl;
   results<<endl;
   results<<"##########################################################################"<<endl;
   Read_fchk_wfn[0].READ_FCHK_WFN::~READ_FCHK_WFN();
  }

  if(IN.apx_wfn)
  {
   int j,k,l,fail=1000,counter=0,terms,terms2;
   int **n_xyz,**n_xyz2;
   results<<endl;
   results<<"##########################################################################"<<endl;
   results<<endl;
   for(i=0;i<IN.nterms_3dho;i++)
   {
    terms=((i+1)*(i+2))/2;
    n_xyz=new int*[terms];
    for(j=0;j<terms;j++)
    {
     n_xyz[j]=new int[3];
    }
    Quant_fill(n_xyz,i);
    for(j=0;j<terms;j++)
    {
     quant_3dho2[0]=n_xyz[j][0];
     quant_3dho2[1]=n_xyz[j][1];
     quant_3dho2[2]=n_xyz[j][2];
     for(k=0;k<IN.nterms_3dho;k++)
     {
      terms2=((k+1)*(k+2))/2;
      n_xyz2=new int*[terms2];
      for(l=0;l<terms2;l++)
      {
       n_xyz2[l]=new int[3];
      }
      Quant_fill(n_xyz2,k);
      for(l=0;l<terms2;l++)
      {
       quant_3dho[0]=n_xyz2[l][0];
       quant_3dho[1]=n_xyz2[l][1];
       quant_3dho[2]=n_xyz2[l][2];
       counter++;
       integral=ZERO;
       fail=1000;
       results<<endl;
       integral=Int_wfn_apx(ha_sol_n,IN.ncores,fail);
       results<<setprecision(0)<<fixed;
       results<<quant_3dho[0]<<" "<<quant_3dho[1]<<" "<<quant_3dho[2]<<" \t ";
       results<<quant_3dho2[0]<<" "<<quant_3dho2[1]<<" "<<quant_3dho2[2]<<" \t ";
       if(integral>=ZERO){results<<" ";}
       results<<setprecision(10)<<fixed<<scientific<<integral<<endl;
      }
      for(l=0;l<terms2;l++)
      {
       delete[] n_xyz2[l];n_xyz2[l]=NULL;
      }
      delete[] n_xyz2; n_xyz2=NULL;
     }
    }
    for(j=0;j<terms;j++)
    {
     delete[] n_xyz[j];n_xyz[j]=NULL;
    }
    delete[] n_xyz;n_xyz=NULL;
   }
   results<<endl;
   results<<"# The number of evaluated terms is: "<<counter<<endl;
   results<<endl;
   results<<"##########################################################################"<<endl;
   results<<endl;
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  // Write the energy, the value of omega and the normalization constant                           //
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  cout<<"See "<<name<<" for results."<<endl;
  results<<endl;
  results<<"##########################################################################"<<endl;
  results<<endl;
  results<<setprecision(10)<<fixed<<scientific;
  results<<"Energy: "<<ha_sol_n.E<<" Omega: "<<ha_sol_n.omega<<endl;
  results<<"Coefficients of Pn(r12):"<<endl;
  for(i=0;i<IN.n;i++)
  {
   results<<"a_"<<setw(2)<<i<<": "<<setprecision(7)<<fixed<<scientific<<pow(sqrt(HALF*ha_sol_n.omega),(double)i)*ha_sol_n.an[i]*(pow(sqrt(HALF*ha_sol_n.omega),(double)IN.l+ONE)/
          (sqrt(HALF*ha_sol_n.omega)+pow(TEN,-TEN)))<<endl;
  }
  results<<setprecision(10)<<fixed<<scientific;
  results<<"Normalization constant: "<<ha_sol_n.N<<endl;
  results<<"[Normalization for the intracule function, Phi(r12), to make it integrate"<<endl;
  results<<" to ONE. The extracule, Xi(R12), integrates to ONE ]"<<endl;
  results<<"Total Normalization constant: "<<ha_sol_n.Ntotal()<<endl;
  results<<"[Normalization for the whole function, to make it integrate"<<endl;
  results<<" to ONE. All constants comming from the spherical harmonics of the intracule (Co from Ylm = Co P_l|m|(cos(tetha)) x function(m,phi)) "<<endl;
  results<<" and from the 3D-Gaussian of the extracule are included. ]"<<endl;
  results<<endl;
  results<<"##########################################################################"<<endl;
  results<<endl;
  results<<"          .--.                  Train yourself to let go"<<endl;
  results<<"::\\`--._,'.::.`._.--'/::        of everything you"<<endl;
  results<<"::::.  ` __::__ '  .::::        fear to lose."<<endl;
  results<<"::::::-:.`'..`'.:-::::::"<<endl;
  results<<"::::::::\\ `--' /::::::::                    -Yoda"<<endl;
  results<<endl;
  results<<"   Normal termination of 2e- HA program"<<endl;
  results<<endl;
  results<<"##########################################################################"<<endl;
  results.close();
 }
 else
 {cout<<"No input file given!"<<endl;}
 return 0;
}
///////////////////////////////////////////////
// Functions for evaluation of quantities    //
///////////////////////////////////////////////
////////////////////////////////////////////////////////
// Integrate the spherical averaged intracule         //
////////////////////////////////////////////////////////
double Int_intracule_sph(HA_SOL_n &ha_sol_n)
{
 int verbose,neval,fail,nregions;// This variable is not needed for Vegas!
 const int MINEVAL=0,MAXEVAL=1000000,NDIM=2,NCOMP=1;
 double integral[NCOMP], error[NCOMP], prob[NCOMP];
 const double EPSREL=1.0e-4, EPSABS=1.0e-5;
 void *SOLUTION=NULL;
 SOLUTION=&ha_sol_n;
 verbose=0;
 Cuhre(NDIM, NCOMP, Int_intra, SOLUTION, NVEC1,EPSREL, EPSABS, verbose | LAST,MINEVAL,
 MAXEVAL, KEY,STATEFILE,SPIN,&nregions, &neval, &fail, integral, error, prob);
 return integral[0];
}
////////////////////////////////////////////////////////
// Integrate the spherical averaged extracule         //
////////////////////////////////////////////////////////
double Int_extracule_sph(HA_SOL_n &ha_sol_n)
{
 int verbose,neval,fail,nregions;// This variable is not needed for Vegas!
 const int MINEVAL=0,MAXEVAL=1000000,NDIM=2,NCOMP=1;
 double integral[NCOMP], error[NCOMP], prob[NCOMP];
 const double EPSREL=1.0e-4, EPSABS=1.0e-5;
 void *SOLUTION=NULL;
 SOLUTION=&ha_sol_n;
 verbose=0;
 Cuhre(NDIM, NCOMP, Int_extra, SOLUTION, NVEC1,EPSREL, EPSABS, verbose | LAST,MINEVAL,
 MAXEVAL, KEY,STATEFILE,SPIN,&nregions, &neval, &fail, integral, error, prob);
 return integral[0];
}
////////////////////////////////////////////////////////
// Integrate the spherical averaged intracule         //
////////////////////////////////////////////////////////
double Int_Vee(HA_SOL_n &ha_sol_n)
{
 int verbose,neval,fail,nregions;// This variable is not needed for Vegas!
 const int MINEVAL=0,MAXEVAL=1000000,NDIM=3,NCOMP=1;
 double integral[NCOMP], error[NCOMP], prob[NCOMP];
 const double EPSREL=1.0e-4, EPSABS=1.0e-5;
 void *SOLUTION=NULL;
 SOLUTION=&ha_sol_n;
 verbose=0;
 Cuhre(NDIM, NCOMP, Int_vee, SOLUTION, NVEC1,EPSREL, EPSABS, verbose | LAST,MINEVAL,
 MAXEVAL, KEY,STATEFILE,SPIN,&nregions, &neval, &fail, integral, error, prob);
 return integral[0];
}
/////////////////////////////////////////////////////////////////////////
// Integrate the density  must get 1 since normalization is 1          //
/////////////////////////////////////////////////////////////////////////
double Int_density(HA_SOL_n &ha_sol_n)
{
 int verbose,neval,fail;//nregions; This variable is not needed for Vegas!
 const int MINEVAL=0,MAXEVAL=10000000,NDIM=6,NCOMP=1;
 double integral[NCOMP], error[NCOMP], prob[NCOMP];
 const double EPSREL=1.0e-4, EPSABS=1.0e-5;
 void *SOLUTION=NULL;
 verbose=0;
 SOLUTION=&ha_sol_n;
 Vegas(NDIM, NCOMP, Int_rho, SOLUTION, NVEC1,EPSREL, EPSABS, verbose, SEED,MINEVAL,
 MAXEVAL, NSTART, NINCREASE, NBATCH,GRIDNO, STATEFILE,SPIN,&neval, &fail, integral,
 error, prob);
 return integral[0];
}
/////////////////////////////////////////////////////////////////////////
// Integrate the density using regions for DI                          //
/////////////////////////////////////////////////////////////////////////
double Int_density2(HA_SOL_n &ha_sol_n,double R_min,double R_max,double Theta_min,
double Theta_max,double Phi_min,double Phi_max)
{
 int verbose,neval,fail;//nregions; This variable is not needed for Vegas!
 const int MINEVAL=0,MAXEVAL=10000000,NDIM=6,NCOMP=1;
 double integral[NCOMP], error[NCOMP], prob[NCOMP];
 const double EPSREL=1.0e-4, EPSABS=1.0e-5;
 void *SOLUTION=NULL;
 verbose=0;
 SOLUTION=&ha_sol_n;
 r_min=R_min;r_max=R_max;
 theta_min=Theta_min;theta_max=Theta_max;
 phi_min=Phi_min;phi_max=Phi_max;
 Vegas(NDIM, NCOMP, Int_rho2, SOLUTION, NVEC1,EPSREL, EPSABS, verbose, SEED,MINEVAL,
 MAXEVAL, NSTART, NINCREASE, NBATCH,GRIDNO, STATEFILE,SPIN,&neval, &fail, integral,
 error, prob);
 return integral[0];
}
/////////////////////////////////////////////////////////////////////////
// Integrate the density using regions for FCHK or WFN DI              //
/////////////////////////////////////////////////////////////////////////
double Int_density_fchk_wfn(double R_min,double R_max,double Theta_min,
double Theta_max,double Phi_min,double Phi_max)
{
 int verbose,neval,fail,nregions;// This variable is not needed for Vegas!
 const int MINEVAL=0,MAXEVAL=1000000,NDIM=3,NCOMP=1;
 double integral[NCOMP], error[NCOMP], prob[NCOMP];
 const double EPSREL=1.0e-4, EPSABS=1.0e-5;
 void *SOLUTION=NULL;
 verbose=0;
 r_min=R_min;r_max=R_max;
 theta_min=Theta_min;theta_max=Theta_max;
 phi_min=Phi_min;phi_max=Phi_max;
 Cuhre(NDIM, NCOMP, Int_rho3, SOLUTION, NVEC1,EPSREL, EPSABS, verbose | LAST,MINEVAL,
 MAXEVAL, KEY,STATEFILE,SPIN,&nregions, &neval, &fail, integral, error, prob);
 return integral[0];
}
/////////////////////////////////////////////////////////////////////////
// Integrate the DM2 using regions for DI                              //
/////////////////////////////////////////////////////////////////////////
double Int_dm2(HA_SOL_n &ha_sol_n,double R_min,double R_max,double Theta_min,
double Theta_max,double Phi_min,double Phi_max,double R2_min,double R2_max,
double Theta2_min,double Theta2_max,double Phi2_min,double Phi2_max)
{
 int verbose,neval,fail;//nregions; This variable is not needed for Vegas!
 const int MINEVAL=0,MAXEVAL=10000000,NDIM=6,NCOMP=1;
 double integral[NCOMP], error[NCOMP], prob[NCOMP];
 const double EPSREL=1.0e-4, EPSABS=1.0e-5;
 void *SOLUTION=NULL;
 verbose=0;
 SOLUTION=&ha_sol_n;
 r_min=R_min;r_max=R_max;
 theta_min=Theta_min;theta_max=Theta_max;
 phi_min=Phi_min;phi_max=Phi_max;
 rr2_min=R2_min;rr2_max=R2_max;
 thetat2_min=Theta2_min;thetat2_max=Theta2_max;
 phip2_min=Phi2_min;phip2_max=Phi2_max;
 Vegas(NDIM, NCOMP, Int_2part_rho, SOLUTION, NVEC1,EPSREL, EPSABS, verbose, SEED,MINEVAL,
 MAXEVAL, NSTART, NINCREASE, NBATCH,GRIDNO, STATEFILE,SPIN,&neval, &fail, integral,
 error, prob);
 return integral[0];
}
/////////////////////////////////////////////////////////////////////////
// Integrate the DM2 using regions for DI FCHK or FCI                  //
/////////////////////////////////////////////////////////////////////////
double Int_dm2_fchk_wfn(HA_SOL_n &ha_sol_n,double R_min,double R_max,double Theta_min,
double Theta_max,double Phi_min,double Phi_max,double R2_min,double R2_max,
double Theta2_min,double Theta2_max,double Phi2_min,double Phi2_max)
{
 int verbose,neval,fail;//,nregions;//nregions; This variable is not needed for Vegas!
 const int MINEVAL=0,MAXEVAL=1000000,NDIM=6,NCOMP=1;
 double integral[NCOMP], error[NCOMP], prob[NCOMP];
 const double EPSREL=1.0e-4, EPSABS=1.0e-5;
 void *SOLUTION=NULL;
 verbose=0;
 SOLUTION=&ha_sol_n;
 r_min=R_min;r_max=R_max;
 theta_min=Theta_min;theta_max=Theta_max;
 phi_min=Phi_min;phi_max=Phi_max;
 rr2_min=R2_min;rr2_max=R2_max;
 thetat2_min=Theta2_min;thetat2_max=Theta2_max;
 phip2_min=Phi2_min;phip2_max=Phi2_max;
 Vegas(NDIM, NCOMP, Int_2part_rho_fchk_wfn, SOLUTION, NVEC1,EPSREL, EPSABS, verbose, SEED,MINEVAL,
 MAXEVAL, NSTART, NINCREASE, NBATCH,GRIDNO, STATEFILE,SPIN,&neval, &fail, integral,
 error, prob);
 integral[0]=FOUR*PI*FOUR*PI*integral[0];
 return integral[0];
}
/////////////////////////////////////////////////////////////////////////
// Integrate the wave function expansion to obtain Fourier coeficients //
/////////////////////////////////////////////////////////////////////////
double Int_wfn_apx(HA_SOL_n &ha_sol_n, int &ncores,int &fail)
{
 int verbose,neval;//,nregions; //This variable is not needed for Vegas!
 const int MINEVAL=0,MAXEVAL=10000000,NDIM=6,NCOMP=1;
 double integral[NCOMP], error[NCOMP], prob[NCOMP];
 const double EPSREL=1.0e-4, EPSABS=1.0e-5;
 void *SOLUTION=NULL;
 verbose=0;
 SOLUTION=&ha_sol_n;
 cubacores(ncores,PCORE);
 Vegas(NDIM, NCOMP, Int_coef_fourier_wfn, SOLUTION, NVEC1,EPSREL, EPSABS, verbose, SEED,MINEVAL,
 MAXEVAL, NSTART, NINCREASE, NBATCH,GRIDNO, STATEFILE,SPIN,&neval, &fail, integral,
 error, prob);
 return integral[0];
}
///////////////////////////////////////////////
// Functions for calculating integrals       //
///////////////////////////////////////////////
/////////////////////////////////////////////////////////////
// Function integrable for the intracule PHI(r12) => h(u)  //
/////////////////////////////////////////////////////////////
int Int_intra(const int *ndim, const double xx[],const int *ncomp, double ff[],void *userdata)
{
 double AUX[3],teta1,phi1,res;
 HA_SOL_n *ha_sol_n;
 ha_sol_n=((HA_SOL_n*)userdata);
 teta1=acos(ONE-TWO*xx[0]);
 phi1=TWO*PI*xx[1];
 AUX[0]=r_intra*sin(teta1)*cos(phi1);
 AUX[1]=r_intra*sin(teta1)*sin(phi1);
 AUX[2]=r_intra*cos(teta1);
 res=ha_sol_n[0].Intracule(AUX);
 ff[0]=FOUR*PI*res*res;
 return 0;
}
/////////////////////////////////////////////////////////////
// Function integrable for the extracule Xi(R12) => H(U)  //
/////////////////////////////////////////////////////////////
int Int_extra(const int *ndim, const double xx[],const int *ncomp, double ff[],void *userdata)
{
 double AUX[3],teta1,phi1,res;
 HA_SOL_n *ha_sol_n;
 ha_sol_n=((HA_SOL_n*)userdata);
 teta1=acos(ONE-TWO*xx[0]);
 phi1=TWO*PI*xx[1];
 AUX[0]=R_extra*sin(teta1)*cos(phi1);
 AUX[1]=R_extra*sin(teta1)*sin(phi1);
 AUX[2]=R_extra*cos(teta1);
 res=ha_sol_n[0].Extracule(AUX);
 ff[0]=FOUR*PI*res*res;
 return 0;
}
///////////////////////////////////////////////////////////
// Function integrable for the intracule PHI(r12) => Vee //
///////////////////////////////////////////////////////////
int Int_vee(const int *ndim, const double xx[],const int *ncomp, double ff[],void *userdata)
{
 double AUX[3],r1,teta1,phi1,res;
 HA_SOL_n *ha_sol_n;
 ha_sol_n=((HA_SOL_n*)userdata);
 r1=xx[0]/(ONE-xx[0]);
 teta1=acos(ONE-TWO*xx[1]);
 phi1=TWO*PI*xx[2];
 AUX[0]=r1*sin(teta1)*cos(phi1);
 AUX[1]=r1*sin(teta1)*sin(phi1);
 AUX[2]=r1*cos(teta1);
 res=ha_sol_n[0].Intracule(AUX);
 // ff[0]=FOUR*PI*xx[0]*xx[0]*(ONE/pow(ONE-xx[0],FOUR))*res*res/r1;
 ff[0]=FOUR*PI*xx[0]*(ONE/pow(ONE-xx[0],THREE))*res*res;
 return 0;
}
/////////////////////////////////////////////////////////////////////////
// Function integrable of the one electron density   Rho=>Nelect       //
/////////////////////////////////////////////////////////////////////////
int Int_rho(const int *ndim, const double xx[],const int *ncomp, double ff[],void *userdata)
{
 double AUX[3],r1,teta1,phi1,AUX2[3],r2,teta2,phi2,res;
 HA_SOL_n *ha_sol_n;
 ha_sol_n=((HA_SOL_n*)userdata);
 r1=xx[0]/(ONE-xx[0]);
 teta1=acos(ONE-TWO*xx[1]);
 phi1=TWO*PI*xx[2];
 AUX[0]=r1*sin(teta1)*cos(phi1);
 AUX[1]=r1*sin(teta1)*sin(phi1);
 AUX[2]=r1*cos(teta1);
 r2=xx[3]/(ONE-xx[3]);
 teta2=acos(ONE-TWO*xx[4]);
 phi2=TWO*PI*xx[5];
 AUX2[0]=r2*sin(teta2)*cos(phi2);
 AUX2[1]=r2*sin(teta2)*sin(phi2);
 AUX2[2]=r2*cos(teta2);
 if(AUX[0]!=AUX2[0] && AUX[1]!=AUX2[1] && AUX[2]!=AUX2[2])
 {res=ha_sol_n[0].PSY12(AUX,AUX2);}
 else
 {
  AUX2[0]=AUX2[0]+pow(TEN,-EIGHT);
  AUX2[1]=AUX2[1]+pow(TEN,-EIGHT);
  AUX2[2]=AUX2[2]+pow(TEN,-EIGHT);
  res=ha_sol_n[0].PSY12(AUX,AUX2);
 }
 ff[0]=pow(FOUR*PI*res*xx[0]*xx[3]*pow(ONE/(ONE-xx[0]),TWO)*pow(ONE/(ONE-xx[3]),TWO),TWO);
 return 0;
}
/////////////////////////////////////////////////////////////////////////
// Function integrable of the one electron density  Rho=>Nelect for DI //
/////////////////////////////////////////////////////////////////////////
int Int_rho2(const int *ndim, const double xx[],const int *ncomp, double ff[],void *userdata)
{
 double AUX[3],r1,teta1,phi1,AUX2[3],r2,teta2,phi2,res;
 HA_SOL_n *ha_sol_n;
 ha_sol_n=((HA_SOL_n*)userdata);
 r1=xx[0]/(ONE-xx[0]);
 teta1=acos(ONE-TWO*xx[1]);
 phi1=TWO*PI*xx[2];
 AUX[0]=r1*sin(teta1)*cos(phi1);
 AUX[1]=r1*sin(teta1)*sin(phi1);
 AUX[2]=r1*cos(teta1);
 r2=xx[3]/(ONE-xx[3]);
 teta2=acos(ONE-TWO*xx[4]);
 phi2=TWO*PI*xx[5];
 AUX2[0]=r2*sin(teta2)*cos(phi2);
 AUX2[1]=r2*sin(teta2)*sin(phi2);
 AUX2[2]=r2*cos(teta2);
 if(AUX[0]!=AUX2[0] && AUX[1]!=AUX2[1] && AUX[2]!=AUX2[2])
 {res=ha_sol_n[0].PSY12(AUX,AUX2);}
 else
 {
  AUX2[0]=AUX2[0]+pow(TEN,-EIGHT);
  AUX2[1]=AUX2[1]+pow(TEN,-EIGHT);
  AUX2[2]=AUX2[2]+pow(TEN,-EIGHT);
  res=ha_sol_n[0].PSY12(AUX,AUX2);
 }
 if(r_min<=r1 && r1<=r_max)
 {
  if(theta_min<=teta1 && teta1<=theta_max)
  {
   if(phi_min<=phi1 && phi1<=phi_max)
   {
    ff[0]=pow(FOUR*PI*res*xx[0]*xx[3]*pow(ONE/(ONE-xx[0]),TWO)*pow(ONE/(ONE-xx[3]),TWO),TWO);
   }
   else
   {
    ff[0]=ZERO;
   }
  }
  else
  {
   ff[0]=ZERO;
  }
 }
 else
 {
  ff[0]=ZERO;
 }
 return 0;
}
/////////////////////////////////////////////////////////////////////////
// Function integrable of the one electron density  Rho=>Nelect for DI //
//                     using a FCHK o WFN file                         //
/////////////////////////////////////////////////////////////////////////
int Int_rho3(const int *ndim, const double xx[],const int *ncomp, double ff[],void *userdata)
{
 double AUX[3],r1,teta1,phi1,res;
 r1=xx[0]/(ONE-xx[0]);
 teta1=acos(ONE-TWO*xx[1]);
 phi1=TWO*PI*xx[2];
 AUX[0]=r1*sin(teta1)*cos(phi1);
 AUX[1]=r1*sin(teta1)*sin(phi1);
 AUX[2]=r1*cos(teta1);
 if(r_min<=r1 && r1<=r_max)
 {
  if(theta_min<=teta1 && teta1<=theta_max)
  {
   if(phi_min<=phi1 && phi1<=phi_max)
   {
    Read_fchk_wfn[0].rho_eval(AUX,res);
    ff[0]=FOUR*PI*res*pow(xx[0]*pow(ONE/(ONE-xx[0]),TWO),TWO);
   }
   else
   {
    ff[0]=ZERO;
   }
  }
  else
  {
   ff[0]=ZERO;
  }
 }
 else
 {
  ff[0]=ZERO;
 }
 return 0;
}
/////////////////////////////////////////////////////////////////////////
// Function integrable of the DM2 in regions for DI calculation        //
/////////////////////////////////////////////////////////////////////////
int Int_2part_rho(const int *ndim, const double xx[],const int *ncomp, double ff[],void *userdata)
{
 double AUX[3],r1,teta1,phi1,AUX2[3],r2,teta2,phi2,res;
 HA_SOL_n *ha_sol_n;
 ha_sol_n=((HA_SOL_n*)userdata);
 r1=xx[0]/(ONE-xx[0]);
 teta1=acos(ONE-TWO*xx[1]);
 phi1=TWO*PI*xx[2];
 AUX[0]=r1*sin(teta1)*cos(phi1);
 AUX[1]=r1*sin(teta1)*sin(phi1);
 AUX[2]=r1*cos(teta1);
 r2=xx[3]/(ONE-xx[3]);
 teta2=acos(ONE-TWO*xx[4]);
 phi2=TWO*PI*xx[5];
 AUX2[0]=r2*sin(teta2)*cos(phi2);
 AUX2[1]=r2*sin(teta2)*sin(phi2);
 AUX2[2]=r2*cos(teta2);
 if(r_min<=r1 && r1<=r_max)
 {
  if(theta_min<=teta1 && teta1<=theta_max)
  {
   if(phi_min<=phi1 && phi1<=phi_max)
   {
    if(rr2_min<=r2 && r2<=rr2_max)
    {
     if(thetat2_min<=teta2 && teta2<=thetat2_max)
     {
      if(phip2_min<=phi2 && phi2<=phip2_max)
      {
       if(AUX[0]!=AUX2[0] && AUX[1]!=AUX2[1] && AUX[2]!=AUX2[2])
       {res=ha_sol_n[0].PSY12(AUX,AUX2);}
       else
       {
        AUX2[0]=AUX2[0]+pow(TEN,-EIGHT);
        AUX2[1]=AUX2[1]+pow(TEN,-EIGHT);
        AUX2[2]=AUX2[2]+pow(TEN,-EIGHT);
        res=ha_sol_n[0].PSY12(AUX,AUX2);
       }
       ff[0]=pow(FOUR*PI*res*xx[0]*xx[3]*pow(ONE/(ONE-xx[0]),TWO)*pow(ONE/(ONE-xx[3]),TWO),TWO);
      }
      else
      {
       ff[0]=ZERO;
      }
     }
     else
     {
      ff[0]=ZERO;
     }
    }
    else
    {
     ff[0]=ZERO;
    }
   }
   else
   {
    ff[0]=ZERO;
   }
  }
  else
  {
   ff[0]=ZERO;
  }
 }
 else
 {
  ff[0]=ZERO;
 }
 return 0;
}

/////////////////////////////////////////////////////////////////////////
// Function integrable of the DM2 in regions for DI calculation        //
/////////////////////////////////////////////////////////////////////////
int Int_2part_rho_fchk_wfn(const int *ndim, const double xx[],const int *ncomp, double ff[],void *userdata)
{
 double AUX[3],r1,teta1,phi1,AUX2[3],r2,teta2,phi2,res,rho1,rho2,rho1_fchk_wfn,rho2_fchk_wfn;
 HA_SOL_n *ha_sol_n;
 ha_sol_n=((HA_SOL_n*)userdata);
 r1=xx[0]/(ONE-xx[0]);
 teta1=acos(ONE-TWO*xx[1]);
 phi1=TWO*PI*xx[2];
 r2=xx[3]/(ONE-xx[3]);
 teta2=acos(ONE-TWO*xx[4]);
 phi2=TWO*PI*xx[5];
 if(r_min<=r1 && r1<=r_max)
 {
  if(theta_min<=teta1 && teta1<=theta_max)
  {
   if(phi_min<=phi1 && phi1<=phi_max)
   {
    if(rr2_min<=r2 && r2<=rr2_max)
    {
     if(thetat2_min<=teta2 && teta2<=thetat2_max)
     {
      if(phip2_min<=phi2 && phi2<=phip2_max)
      {
       AUX[0]=r1*sin(teta1)*cos(phi1);
       AUX[1]=r1*sin(teta1)*sin(phi1);
       AUX[2]=r1*cos(teta1);
       AUX2[0]=r2*sin(teta2)*cos(phi2);
       AUX2[1]=r2*sin(teta2)*sin(phi2);
       AUX2[2]=r2*cos(teta2);
       if(AUX[0]!=AUX2[0] && AUX[1]!=AUX2[1] && AUX[2]!=AUX2[2])
       {res=ha_sol_n[0].PSY12(AUX,AUX2);}
       else
       {
        AUX2[0]=AUX2[0]+pow(TEN,-EIGHT);
        AUX2[1]=AUX2[1]+pow(TEN,-EIGHT);
        AUX2[2]=AUX2[2]+pow(TEN,-EIGHT);
        res=ha_sol_n[0].PSY12(AUX,AUX2);
       }
       Read_fchk_wfn[0].rho_eval(AUX,rho1);
       Read_fchk_wfn[0].rho_eval(AUX2,rho2);
       rho1=rho1/(double)Read_fchk_wfn[0].nelectrons;
       rho2=rho2/(double)Read_fchk_wfn[0].nelectrons;
       rho1_fchk_wfn=ha_sol_n[0].Rho_r(AUX);
       rho2_fchk_wfn=ha_sol_n[0].Rho_r(AUX2);
       if(rho1<pow(TEN,-TEN)){rho1=pow(TEN,-TEN);}
       if(rho1_fchk_wfn<pow(TEN,-TEN)){rho1_fchk_wfn=pow(TEN,-TEN);}
       if(rho2<pow(TEN,-TEN)){rho2=pow(TEN,-TEN);}
       if(rho2_fchk_wfn<pow(TEN,-TEN)){rho2_fchk_wfn=pow(TEN,-TEN);}
       //int rho1 d1 = 2 while int ha_sol_n.Rho_r(AUX) = 1 that is why we have to correct the FCHK or WFN density
       res=pow(res,TWO)*(rho1*rho2/(rho1_fchk_wfn*rho2_fchk_wfn));
       ff[0]=res*pow(xx[0]*xx[3]*pow(ONE/(ONE-xx[0]),TWO)*pow(ONE/(ONE-xx[3]),TWO),TWO);
      }
      else
      {
       ff[0]=ZERO;
      }
     }
     else
     {
      ff[0]=ZERO;
     }
    }
    else
    {
     ff[0]=ZERO;
    }
   }
   else
   {
    ff[0]=ZERO;
   }
  }
  else
  {
   ff[0]=ZERO;
  }
 }
 else
 {
  ff[0]=ZERO;
 }
 return 0;
}

////////////////////////////////////////////////////////////////////////
// Function integrable for the Fourier coeficients of the DM1 apx     //
////////////////////////////////////////////////////////////////////////
int Int_coef_fourier_wfn(const int *ndim, const double xx[],const int *ncomp, double ff[],void *userdata)
{
 double AUX[3],r1,teta1,phi1,AUX2[3],r2,teta2,phi2,res;
 HA_SOL_n *ha_sol_n;
 ha_sol_n=((HA_SOL_n*)userdata);
 r1=xx[0]/(ONE-xx[0]);
 teta1=acos(ONE-TWO*xx[1]);
 phi1=TWO*PI*xx[2];
 AUX[0]=r1*sin(teta1)*cos(phi1);
 AUX[1]=r1*sin(teta1)*sin(phi1);
 AUX[2]=r1*cos(teta1);
 r2=xx[3]/(ONE-xx[3]);
 teta2=acos(ONE-TWO*xx[4]);
 phi2=TWO*PI*xx[5];
 AUX2[0]=r2*sin(teta2)*cos(phi2);
 AUX2[1]=r2*sin(teta2)*sin(phi2);
 AUX2[2]=r2*cos(teta2);
 if(AUX[0]==AUX2[0] && AUX[1]==AUX2[1] && AUX[2]==AUX2[2])
 {
  AUX2[0]=AUX2[0]+pow(TEN,-EIGHT);
  AUX2[1]=AUX2[1]+pow(TEN,-EIGHT);
  AUX2[2]=AUX2[2]+pow(TEN,-EIGHT);
 }
 res=pow(TWO,HALF)*ha_sol_n[0].PSY12(AUX,AUX2)
    *Q_HARM_OSC_3D(quant_3dho2[0],quant_3dho2[1],quant_3dho2[2],ha_sol_n[0].omega,AUX)
    *Q_HARM_OSC_3D(quant_3dho[0],quant_3dho[1],quant_3dho[2],ha_sol_n[0].omega,AUX2);
 ff[0]=res*pow(FOUR*PI*xx[0]*xx[3]*pow(ONE/(ONE-xx[0]),TWO)*pow(ONE/(ONE-xx[3]),TWO),TWO);
 return 0;
}

//Generate n, l and m
void Quant_fill(int **Quant,int styp)
{
 int i,j,k,qfil;
 qfil=0;
 for(i=styp;i>-1;i=i-1)
 {
  for(j=i;j>-1;j=j-1)
  {
   for(k=j;k>-1;k=k-1)
   {
    if((i+j+k)==styp)
    {
     Quant[qfil][0]=i;
     Quant[qfil][1]=j;
     Quant[qfil][2]=k;
     qfil++;
     if(j!=k)
     {
      Quant[qfil][0]=i;
      Quant[qfil][1]=k;
      Quant[qfil][2]=j;
      qfil++;
     }
     if(i!=j)
     {
      Quant[qfil][0]=j;
      Quant[qfil][1]=i;
      Quant[qfil][2]=k;
      qfil++;
     }
     if((i!=j) && (j!=k))
     {
      Quant[qfil][0]=j;
      Quant[qfil][1]=k;
      Quant[qfil][2]=i;
      qfil++;
      Quant[qfil][0]=k;
      Quant[qfil][1]=i;
      Quant[qfil][2]=j;
      qfil++;
     }
     if(k!=i)
     {
      Quant[qfil][0]=k;
      Quant[qfil][1]=j;
      Quant[qfil][2]=i;
      qfil++;
     }
    }
   }
  }
 }
}

