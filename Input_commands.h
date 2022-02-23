#ifndef _INPUT_COMMANDS_H_
#define _INPUT_COMMANDS_H_
#include<iostream>
#include<algorithm>
#include<fstream>
#include<stdlib.h>
#include<string>
#include<stdio.h>
#include<iomanip>
#include"String_ops.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

using namespace std;
///////////////
//Begin class//
///////////////
class Input_commands
{
 public:
 Input_commands();
 Input_commands(string);
 int order_rad,order_ang,nterms_3dho;
 int n,l,m,nx,ny,nz,ncores;
 double step,last;
 double r_intra[3],r_extra[3],r1_wf[3],r2_wf[3],r_dens[3],r1_dm1[3],r1_prime_dm1[3];
 double r1_dm2[3],r1_prime_dm2[3],r2_dm2[3],r2_prime_dm2[3];
 double r1_min,r2_min,r1_max,r2_max,phi1_min,phi2_min,phi1_max,phi2_max,theta1_min,theta2_min,theta1_max;
 double theta2_max,last_intra,step_intra;
 double r1_min_an,r2_min_an,r1_max_an,r2_max_an,phi1_min_an,phi2_min_an,phi1_max_an,phi2_max_an,theta1_min_an,theta2_min_an;
 double theta1_max_an,theta2_max_an;
 bool radial_grid,intracule,scan_intracule,extracule,wave_function,density,scan_density,integrate_density,dmn,dm1,dm2,DI_an,DI_fchk_wfn;
 bool apx_wfn,cnormdi;
 char axis,axis_intra;
 string name_fchk_wfn,name_log;
 bool WFN,log_file,cas;
 int mult_in;
};
#endif // _INPUT_COMMANDS_H_
