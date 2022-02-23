#include"HA_sol_n.h"
   ///////////////////////////
   //Declare Class Functions//
   ///////////////////////////

/////////////////////////////////////////
/////////////////////////////////////////
// Creation and anihilation of objects //
// ----------------------------------- //
/////////////////////////////////////////
/////////////////////////////////////////

HA_SOL_n::HA_SOL_n(){cout<<"Not allowed default constructor"<<endl;}
//Build Density object from fchk or wfn (read all required variables)
HA_SOL_n::HA_SOL_n(int n_in,int l_in, int m_in, int nx_in, int ny_in, int nz_in,int order_rad_in,string name_in)
{
 if(n_in<2){n=2;}
 else{n=n_in;}
 an=new double[n];
 l=(double)l_in;m=(double)m_in;
 nx=(double)nx_in;ny=(double)ny_in;nz=(double)nz_in;
 order_rad=order_rad_in;
 quadrature_name=name_in.substr(0,(name_in.length()-4));
 ///////////////////////////////////////////////////////////////////////////////////////////////////
 // Calculate omega, the energy and the normalization constant                                    //
 ///////////////////////////////////////////////////////////////////////////////////////////////////
 e_double_prime=TWO*(l+(double)n)+ONE;
 filling_coef();
 e_prime=e_double_prime*omega/FOUR;
 E=TWO*e_prime+omega*((nx+HALF)+(ny+HALF)+(nz+HALF));
 create_grid(); //Legendre quadrature always created
 calculate_norm_constant();
}
//Delete object
HA_SOL_n::~HA_SOL_n()
{
 clear_grid(); //Legendre quadrature always deleted
 delete[] an;
 an=NULL;
}

///////////////////////////////////////////////
// Functions for building the 2eHA solution  //
///////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
// Find omega and the coefficients of the power series                 //
/////////////////////////////////////////////////////////////////////////
void HA_SOL_n::filling_coef()
{
 int i,counter=0;
 double error_abs,h=pow(TEN,-NINE);
 sqrt_omega_r=pow(TEN,-SEVEN);
 do
 {
  counter++;
  sqrt_omega_r=sqrt_omega_r-TWO*a_n(sqrt_omega_r,n,e_double_prime)*h/(a_n(sqrt_omega_r+h,n,e_double_prime)
  -a_n(sqrt_omega_r-h,n,e_double_prime));
  error_abs=abs(a_n(sqrt_omega_r,n,e_double_prime));
 }while(error_abs>pow(TEN,-TEN) && counter<10000);
 omega=TWO*pow(sqrt_omega_r,TWO);
 for(i=0;i<n;i++)
 {an[i]=a_n(sqrt_omega_r,i,e_double_prime);}
}
/////////////////////////////////////////////////////////////////////////
// Build the coefficients of the power series (see Taut 1993)          //
/////////////////////////////////////////////////////////////////////////
double HA_SOL_n::a_n(double sqrt_omega_r, int v, double e_double_prime)
{
 int v_1,v_2;
 double coef;
 if(v==0)
 {
  coef=ONE;
 }
 else if(v==1)
 {
  coef=ONE/(TWO*(l+ONE)*sqrt_omega_r);
 }
 else
 {
  v_1=v-1;
  v_2=v-2;
  coef=(ONE/(((double)v)*(((double)v)+TWO*l+ONE)))*((a_n(sqrt_omega_r,v_1,e_double_prime)/sqrt_omega_r)
  +a_n(sqrt_omega_r,v_2,e_double_prime)*(TWO*(l+((double)v))-ONE-e_double_prime));
 }
 return coef;
}
/////////////////////////////////////////////////////////////////////////
// Function u(r)/r (see Taut 1993)                                     //
/////////////////////////////////////////////////////////////////////////
double HA_SOL_n::u_r_over_r(double r)
{
 int i;
 double rho,u_rho;
 rho=r*sqrt_omega_r;
 u_rho=ZERO;
 for(i=0;i<n;i++)
 {u_rho=u_rho+an[i]*pow(rho,(double)i);}
 return (u_rho*pow(rho,l+ONE)*exp(-HALF*pow(rho,TWO)))/(rho+pow(TEN,-TEN));
}
/////////////////////////////////////////////////////////////////////////
// Normalization constant                                              //
/////////////////////////////////////////////////////////////////////////
void HA_SOL_n::calculate_norm_constant()
{
 int i;
 double r,integral=ZERO;
 double *coord_r,*w_radial;
 ifstream wr_read,r_read;
 coord_r=new double[order_rad];
 w_radial=new double[order_rad];
 //Grid for r
 wr_read.open((quadrature_name+"_w.txt").c_str());
 r_read.open((quadrature_name+"_x.txt").c_str());
 for(i=0;i<order_rad;i++)
 {
  r_read>>coord_r[i];
  wr_read>>w_radial[i];
 }
 wr_read.close();
 r_read.close();
 for(i=0;i<order_rad;i++)
 {
  //r is for entire space (0 to inf)
  r=coord_r[i]/(ONE-coord_r[i]);
  integral=integral+w_radial[i]*pow(r*u_r_over_r(r),TWO)*pow(ONE/(ONE-coord_r[i]),TWO);
 }
 delete[] w_radial;w_radial=NULL;
 delete[] coord_r;coord_r=NULL;
 N=ONE/pow(integral,HALF);
}
/////////////////////////////////////////////////////////////////////////
// Center of mass and intracule                                        //
/////////////////////////////////////////////////////////////////////////
void HA_SOL_n::R_and_r(double r1[3],double r2[3],double r[3],double R[3])
{
 int i;
 for(i=0;i<3;i++)
 {
  R[i]=(r1[i]+r2[i])/TWO;
  r[i]=(r2[i]-r1[i]);
 }
}

//////////////////////////////////////////////////////
// Functions for evaluations in the 2eHA solutions  //
//////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////
// Evaluate the intracule                                                                        //
///////////////////////////////////////////////////////////////////////////////////////////////////
double HA_SOL_n::Intracule(double r12[3])
{
 double eval,rnorm,rteta,rphi;
 rnorm=norm3D(r12);
 rteta=acos(r12[2]/(rnorm+pow(TEN,-TEN)));
 rphi=atan(r12[1]/(r12[0]+pow(TEN,-TEN)));
 eval=N*u_r_over_r(rnorm)*Y_lm(l,m,rteta,rphi);
 return eval;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
// Evaluate the extracule                                                                        //
///////////////////////////////////////////////////////////////////////////////////////////////////
double HA_SOL_n::Extracule(double R12[3])
{
 double eval,phi_nx,phi_ny,phi_nz;
 phi_nx=pow(pow(TWO,nx)*factorial((int)nx),-HALF)*pow(TWO*omega/PI,ONE/FOUR);
 phi_nx=phi_nx*H_n(nx,R12[0]*pow(TWO*omega,HALF))*exp(-omega*R12[0]*R12[0]);
 phi_ny=pow(pow(TWO,ny)*factorial((int)ny),-HALF)*pow(TWO*omega/PI,ONE/FOUR);
 phi_ny=phi_ny*H_n(ny,R12[1]*pow(TWO*omega,HALF))*exp(-omega*R12[1]*R12[1]);
 phi_nz=pow(pow(TWO,nz)*factorial((int)nz),-HALF)*pow(TWO*omega/PI,ONE/FOUR);
 phi_nz=phi_nz*H_n(nz,R12[2]*pow(TWO*omega,HALF))*exp(-omega*R12[2]*R12[2]);
 eval=phi_nx*phi_ny*phi_nz;
 return eval;
}
/////////////////////////////////////////////////////////////////////////
// Wave function                                                       //
/////////////////////////////////////////////////////////////////////////
double HA_SOL_n::PSY12(double r1[3],double r2[3])
{
 double eval,rnorm,rteta,rphi,phi_nx,phi_ny,phi_nz,r[3],R[3];
 R_and_r(r1,r2,r,R);
 phi_nx=pow(pow(TWO,nx)*factorial((int)nx),-HALF)*pow(TWO*omega/PI,ONE/FOUR);
 phi_nx=phi_nx*H_n(nx,R[0]*pow(TWO*omega,HALF))*exp(-omega*R[0]*R[0]);
 phi_ny=pow(pow(TWO,ny)*factorial((int)ny),-HALF)*pow(TWO*omega/PI,ONE/FOUR);
 phi_ny=phi_ny*H_n(ny,R[1]*pow(TWO*omega,HALF))*exp(-omega*R[1]*R[1]);
 phi_nz=pow(pow(TWO,nz)*factorial((int)nz),-HALF)*pow(TWO*omega/PI,ONE/FOUR);
 phi_nz=phi_nz*H_n(nz,R[2]*pow(TWO*omega,HALF))*exp(-omega*R[2]*R[2]);
 rnorm=norm3D(r);
 rteta=acos(r[2]/(rnorm+pow(TEN,-TEN)));
 rphi=atan(r[1]/(r[0]+pow(TEN,-TEN)));
 eval=N*u_r_over_r(rnorm)*Y_lm(l,m,rteta,rphi)*phi_nx*phi_ny*phi_nz;
 return eval;
}
/////////////////////////////////////////////////////////////////////////
// Total Nomralization constant for the Wave function                  //
/////////////////////////////////////////////////////////////////////////
double HA_SOL_n::Ntotal()
{
 double cte;
 cte=N*pow(pow(TWO,nx)*factorial((int)nx),-HALF)*pow(TWO*omega/PI,ONE/FOUR)*
          pow(pow(TWO,nx)*factorial((int)nx),-HALF)*pow(TWO*omega/PI,ONE/FOUR)*
          pow(pow(TWO,nx)*factorial((int)nx),-HALF)*pow(TWO*omega/PI,ONE/FOUR)*
          pow((TWO*l+ONE)*(double)(factorial((int)l-(int)abs(m)))/(FOUR*PI*(double)(factorial((int)l+(int)abs(m)))),HALF);
 if(m!=0){cte=cte*pow(TWO,HALF);}
 return cte;
}
///////////////////////////////////////////////////////////////////////
// Functions for evaluating quantities (Density, DM1, DM2, etc)      //
///////////////////////////////////////////////////////////////////////
void HA_SOL_n::set_order_ang(int n)
{
 order_ang=n;
 grid_avail(order_ang);
 cout<<"Angular grid chosen containning: "<<n<<" points"<<endl;
 cout<<"Angular grid available: "<<order_ang<<" points"<<endl;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
// Evaluate the density  Psi(1,2) => rho(r)                                                      //
///////////////////////////////////////////////////////////////////////////////////////////////////
double HA_SOL_n::Rho_r(double r1[3])
{
 int i,j;
 double r,r2[3],Integral=ZERO;
 double *coord_r,*w_radial,*x,*y,*z,*w_theta_phi,*theta,*phi,*temp,*temp1;
 ifstream wr_read,r_read;
 coord_r=new double[order_rad];
 w_radial=new double[order_rad];
 //Grid for r
 wr_read.open((quadrature_name+"_w.txt").c_str());
 r_read.open((quadrature_name+"_x.txt").c_str());
 for(i=0;i<order_rad;i++)
 {
  r_read>>coord_r[i];
  wr_read>>w_radial[i];
 }
 wr_read.close();
 r_read.close();
 //Grid for theta, phi
 x=new double[order_ang];
 y=new double[order_ang];
 z=new double[order_ang];
 w_theta_phi=new double[order_ang];
 theta=new double[order_ang];
 phi=new double[order_ang];
 ld_by_order(order_ang,x,y,z,w_theta_phi);
 temp=new double[1];
 temp1=new double[1];
 for(i=0;i<order_ang;i++)
 {
  xyz_to_tp(x[i],y[i],z[i],temp,temp1);
  phi[i]=temp[0];
  theta[i]=temp1[0];
 }
 delete[] temp; temp=NULL;
 delete[] temp1; temp1=NULL;
 //Loop over radial points
 for(i=0;i<order_rad;i++)
 {
  //r is for entire space (0 to inf)
  r=coord_r[i]/(ONE-coord_r[i]);
  for(j=0;j<order_ang;j++)
  {
   r2[0]=r*x[j];
   r2[1]=r*y[j];
   r2[2]=r*z[j];
   Integral=Integral+w_theta_phi[j]*w_radial[i]*pow(r*PSY12(r1,r2),TWO)*pow(ONE/(ONE-coord_r[i]),TWO);
  }
 }
 delete[] w_radial;w_radial=NULL;
 delete[] coord_r;coord_r=NULL;
 delete[] x;x=NULL;
 delete[] y;y=NULL;
 delete[] z;z=NULL;
 delete[] w_theta_phi; w_theta_phi=NULL;
 delete[] theta; theta=NULL;
 delete[] phi; phi=NULL;
 return FOUR*PI*Integral;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
// Evaluate the DM1  Psi(1,2) => DM1(r1,r1')                                                     //
///////////////////////////////////////////////////////////////////////////////////////////////////
double HA_SOL_n::Int_dm1(double r1[3],double r1_prime[3])
{
 int i,j;
 double r,r2[3],Integral=ZERO;
 double *coord_r,*w_radial,*x,*y,*z,*w_theta_phi,*theta,*phi,*temp,*temp1;
 ifstream wr_read,r_read;
 coord_r=new double[order_rad];
 w_radial=new double[order_rad];
 //Grid for r
 wr_read.open((quadrature_name+"_w.txt").c_str());
 r_read.open((quadrature_name+"_x.txt").c_str());
 for(i=0;i<order_rad;i++)
 {
  r_read>>coord_r[i];
  wr_read>>w_radial[i];
 }
 wr_read.close();
 r_read.close();
 //Grid for theta, phi
 x=new double[order_ang];
 y=new double[order_ang];
 z=new double[order_ang];
 w_theta_phi=new double[order_ang];
 theta=new double[order_ang];
 phi=new double[order_ang];
 ld_by_order(order_ang,x,y,z,w_theta_phi);
 temp=new double[1];
 temp1=new double[1];
 for(i=0;i<order_ang;i++)
 {
  xyz_to_tp(x[i],y[i],z[i],temp,temp1);
  phi[i]=temp[0];
  theta[i]=temp1[0];
 }
 delete[] temp; temp=NULL;
 delete[] temp1; temp1=NULL;
 //Loop over radial points
 for(i=0;i<order_rad;i++)
 {
  //r is for entire space (0 to inf)
  r=coord_r[i]/(ONE-coord_r[i]);
  for(j=0;j<order_ang;j++)
  {
   r2[0]=r*x[j];
   r2[1]=r*y[j];
   r2[2]=r*z[j];
   Integral=Integral+w_theta_phi[j]*w_radial[i]*r*r*PSY12(r1,r2)*PSY12(r1_prime,r2)*pow(ONE/(ONE-coord_r[i]),TWO);
  }
 }
 delete[] w_radial;w_radial=NULL;
 delete[] coord_r;coord_r=NULL;
 delete[] x;x=NULL;
 delete[] y;y=NULL;
 delete[] z;z=NULL;
 delete[] w_theta_phi; w_theta_phi=NULL;
 delete[] theta; theta=NULL;
 delete[] phi; phi=NULL;
 return FOUR*PI*Integral;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
// Integrate the wave function Psi(1,2) => Nelectrons                                            //
///////////////////////////////////////////////////////////////////////////////////////////////////
double HA_SOL_n::Int_density()
{
 int i,j,k,l;
 double real_1,real_2,r1[3],r2[3],Integral=ZERO;
 double *coord_r,*w_radial,*x,*y,*z,*w_theta_phi,*theta,*phi,*temp,*temp1;
 ifstream wr_read,r_read;
 coord_r=new double[order_rad];
 w_radial=new double[order_rad];
 //Grid for r
 wr_read.open((quadrature_name+"_w.txt").c_str());
 r_read.open((quadrature_name+"_x.txt").c_str());
 for(i=0;i<order_rad;i++)
 {
  r_read>>coord_r[i];
  wr_read>>w_radial[i];
 }
 wr_read.close();
 r_read.close();
 //Grid for theta, phi
 x=new double[order_ang];
 y=new double[order_ang];
 z=new double[order_ang];
 w_theta_phi=new double[order_ang];
 theta=new double[order_ang];
 phi=new double[order_ang];
 ld_by_order(order_ang,x,y,z,w_theta_phi);
 temp=new double[1];
 temp1=new double[1];
 for(i=0;i<order_ang;i++)
 {
  xyz_to_tp(x[i],y[i],z[i],temp,temp1);
  phi[i]=temp[0];
  theta[i]=temp1[0];
 }
 delete[] temp; temp=NULL;
 delete[] temp1; temp1=NULL;
 //Loop over radial points of the electron 1
 for(i=0;i<order_rad;i++)
 {
  //real_1 is for entire space for electron 1 (0 to inf)
  real_1=coord_r[i]/(ONE-coord_r[i]);
  for(j=0;j<order_ang;j++)
  {
   r1[0]=real_1*x[j];
   r1[1]=real_1*y[j];
   r1[2]=real_1*z[j];
   //Loop over radial points of the electron 1
   for(k=0;k<order_rad;k++)
   {
    //real_1 is for entire space for electron 1 (0 to inf)
    real_2=coord_r[k]/(ONE-coord_r[k]);
    for(l=0;l<order_ang;l++)
    {
     r2[0]=real_2*x[l];
     r2[1]=real_2*y[l];
     r2[2]=real_2*z[l];
     if(r2[0]==r1[0] && r2[1]==r1[1] && r2[2]==r1[2])
     {
      r2[0]=r2[0]+pow(TEN,-EIGHT);
      r2[1]=r2[1]+pow(TEN,-EIGHT);
      r2[2]=r2[2]+pow(TEN,-EIGHT);
     }
     Integral=Integral+w_theta_phi[l]*w_radial[k]*w_theta_phi[j]*w_radial[i]*pow(real_2*real_1*PSY12(r1,r2),TWO)*pow(ONE/(ONE-coord_r[i]),TWO)*pow(ONE/(ONE-coord_r[k]),TWO);
    }
   }
  }
 }
 delete[] w_radial;w_radial=NULL;
 delete[] coord_r;coord_r=NULL;
 delete[] x;x=NULL;
 delete[] y;y=NULL;
 delete[] z;z=NULL;
 delete[] w_theta_phi; w_theta_phi=NULL;
 delete[] theta; theta=NULL;
 delete[] phi; phi=NULL;
 return FOUR*FOUR*PI*PI*Integral;
}

//////////////////////////////////////////////
// Functions for Legendre quadrature        //
//////////////////////////////////////////////

void HA_SOL_n::create_grid()
{
 double r_inf,r_sup;
 r_inf=ZERO;
 r_sup=ONE;
 legendre_quadrature(quadrature_name,order_rad,r_inf,r_sup);
}
void HA_SOL_n::clear_grid()
{
 system(("rm "+quadrature_name+"_w.txt").c_str());
 system(("rm "+quadrature_name+"_x.txt").c_str());
 system(("rm "+quadrature_name+"_r.txt").c_str());
}
//////////////////////////////////////////////
// Functions for Sphere Lebedev GRIDNO      //
//////////////////////////////////////////////

void HA_SOL_n::grid_avail(int & Order)
{
  if(Order<=6)
  {
   Order =    6;
  }
  else if(Order<=14)
  {
   Order =   14;
  }
  else if(Order<=26)
  {
   Order =   26;
  }
  else if(Order<=38)
  {
   Order =   38;
  }
  else if(Order<=50)
  {
   Order =   50;
  }
  else if(Order<=74)
  {
   Order =   74;
  }
  else if(Order<=86)
  {
   Order =   86;
  }
  else if(Order<=110)
  {
   Order =  110;
  }
  else if(Order<=146)
  {
   Order =  146;
  }
  else if(Order<=170)
  {
   Order =  170;
  }
  else if(Order<=194)
  {
   Order =  194;
  }
  else if(Order<=230)
  {
   Order =  230;
  }
  else if(Order<=266)
  {
   Order =  266;
  }
  else if(Order<=302)
  {
   Order =  302;
  }
  else if(Order<=350)
  {
   Order =  350;
  }
  else if(Order<=434)
  {
   Order =  434;
  }
  else if(Order<=590)
  {
   Order =  590;
  }
  else if(Order<=770)
  {
   Order =  770;
  }
  else if(Order<=974)
  {
   Order =  974;
  }
  else if(Order<=1202)
  {
   Order = 1202;
  }
  else if(Order<=1454)
  {
   Order = 1454;
  }
  else if(Order<=1730)
  {
   Order = 1730;
  }
  else if(Order<=2030)
  {
   Order = 2030;
  }
  else if(Order<=2354)
  {
   Order = 2354;
  }
  else if(Order<=2702)
  {
   Order = 2702;
  }
  else if(Order<=3074)
  {
   Order = 3074;
  }
  else if(Order<=3470)
  {
   Order = 3470;
  }
  else if(Order<=3890)
  {
   Order = 3890;
  }
  else if(Order<=4334)
  {
   Order = 4334;
  }
  else if(Order<=4802)
  {
   Order = 4802;
  }
  else if(Order<=5294)
  {
   Order = 5294;
  }
  else
  {
   Order = 5810;
  }
}
