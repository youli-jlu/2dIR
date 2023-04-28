#include "response2d.h"
#include <numeric>
#include <fstream>
#include <ctime>
#include <string>
#include <complex>
#include <cmath>
//This functionn will calculate rephasing and non-rephasing response function
//There's no centain way to add population relaxation. I will use formula in ref 1.

using namespace std;

//Youli : use class method as response calculation
//reference:1: DOI: 10.1063/1.1633549 ; 2: DOI: 10.1063/1.1961472;
//3: 10.1021/jp907648y

//youli: some intial method for response2d class
response2d::response2d(para parameterin,transition transin1,transition transin2):
  parameter(parameterin),trans1(transin1),trans2(transin2) 
{
  tmax=parameter.tmax;
  ntotal=pow(tmax,parameter.dimension);
  dimension=parameter.dimension;
  gamma1=1.0/(trans1.T1/parameter.dt);
  gamma2=1.0/(trans2.T1/parameter.dt);
  R_r=new complex<double>[ntotal];
  R_nr=new complex<double>[ntotal];
}
response2d::~response2d() {};
void response2d::initial() {};

complex<double> & response2d::r(int t1,int t3)
{
  int index=t1*tmax+t3;
  return R_r[index];
}
complex<double> & response2d::nr(int t1,int t3)
{
  int index=t1*tmax+t3;
  return R_nr[index];
}
void response2d::printr(string filename)
{
  checkname check(filename);
  string fname=check.ifreplace();
  ofstream outfile(fname,ios::out);
  outfile << "re"  <<  "  "  << "im" << endl; 
  for (int i=0;i<tmax;i++)
  {
    for (int j=0;j<tmax;j++)
    {
      outfile <<i<< "  " << j<< " " << r(i,j).real()  <<  "  "  << r(i,j).imag() << endl; 
    }
  }
  outfile.close();
}
void response2d::printnr(string filename)
{
  checkname check(filename);
  string fname=check.ifreplace();
  ofstream outfile(fname,ios::out);
  outfile << "re"  <<  "  "  << "im" << endl; 
  for (int i=0;i<tmax;i++)
  {
    for (int j=0;j<tmax;j++)
    {
      outfile<<i<< "  " << j<< " "  << nr(i,j).real()  <<  "  "  << nr(i,j).imag() << endl; 
    }
  }
}
void response2d::zero_padding1d(complex<double> *&pR,int np)
{
  complex<double> *pR_new=nullptr;
  int tmax_new=tmax+np;
  int ntotal_new=pow(tmax_new,dimension);
  //!!!
  cout <<tmax_new<<endl;
  cout <<ntotal_new<<endl;
  pR_new= new complex<double>[ntotal_new];
  int index,index_new;
  for (int i=0;i<tmax;i++)
  {
    for (int j=0;j<tmax;j++)
    {
      index_new=j+tmax_new*i;
      index=j+tmax*i;
      pR_new[index_new]=pR[index];
    }
  }
  for (int i=tmax;i<tmax+np;i++)
  {
    for (int j=tmax;j<tmax+np;j++)
    {
      index=j+tmax_new*i;
      pR_new[index]=0.0+0.0*1i;
    }
  }
  delete [] pR;
  pR=pR_new;
}
void response2d::zero_padding()
{
  if (parameter.nzp>0){
    zero_padding1d(R_r,parameter.nzp);
    zero_padding1d(R_nr,parameter.nzp);
    tmax+=parameter.nzp;
    ntotal=pow(tmax,dimension);
  }
}
void response2d::trans_mean()
{
  trans1.fmean_cal();
  trans1.fmean_removal();
  trans2.fmean_cal();
  trans2.fmean_removal();
  return;
}
void response2d::initial_value()
{
  int i=0;
  for (i=0;i<tmax;i++)
  {
    r (0,i)/=2.0;
    nr(0,i)/=2.0;
  }
  for (i=1;i<tmax;i++)
  {
    r (i,0)/=2.0;
    nr(i,0)/=2.0;
  } 
}



//condon2d method
condon2d::condon2d(para parameterin,transition transin1,transition transin2)
  :response2d(parameterin,transin1,transin2) {} ;

double condon2d::gamma_TA(int t1,int t3) { 
  return exp( -(gamma1+gamma2)*t3/2.0 -gamma1*parameter.t2 - gamma1*t1/2.0 );}
double condon2d::gamma_SE(int t1,int t3) { 
  return exp( -(gamma1+gamma2)*t3/2.0 -gamma1*parameter.t2 - gamma1*t1/2.0 );}


void condon2d::calculation()
{
  //condon approximation,calculate mu10**4 and mu10**2*mu21**2
  double td10,td21;
  td10=trans1.td[0]*trans1.td[0]*trans1.td[0]*trans1.td[0];
  td21=trans1.td[0]*trans2.td[0]*trans1.td[0]*trans2.td[0];
  trans_mean();
  double domega=trans2.freq_mean-trans1.freq_mean;
  clock_t time[4];
  time[0] = clock();

  //*** calculate response function
  //*** based on ref1 and ref2
  double f=0.0,g10=0.0,g21=0.0;
  int tk=0;
  double delta;
  for(int k=0;k<parameter.np;k++)
  {
    tk=k*parameter.tgap;
    f=0.0;
    for(int t1=0;t1<tmax;t1++)
    {
      g10=0.0;
      g21=0.0;
      for(int t3=0;t3<tmax;t3++)
      {
        delta=domega*t3;
        r(t1,t3) += (2.0*td10*exp(double( f-g10)*1i)*gamma_SE(t1,t3) - td21*exp(double( f-g21+delta)*1i)*gamma_TA(t1,t3))/double(parameter.np);
        nr(t1,t3)+= (2.0*td10*exp(double(-f-g10)*1i)*gamma_SE(t1,t3) - td21*exp(double(-f-g21+delta)*1i)*gamma_TA(t1,t3))/double(parameter.np);
        g10+=trans1.freq[tk+t1+parameter.t2+t3];
        g21+=trans2.freq[tk+t1+parameter.t2+t3];
      }
      f+=trans1.freq[tk+t1];
    }
  }


  time[1]=clock();
  cout << "total time for response function:" 
    << double(time[1]-time[0]) / double(CLOCKS_PER_SEC) << "seconds"<< endl;

}

//non_condon2d method
non_condon2d::non_condon2d(para parameterin,transition transin1,transition transin2)
  :response2d(parameterin,transin1,transin2) {} ;
double non_condon2d::td10(int t1,int t3,int tk) { 
  return trans1.td[tk]*trans1.td[tk+t1]*trans1.td[tk+t1+parameter.t2]*trans1.td[tk+t1+parameter.t2+t3];}
double non_condon2d::td21(int t1,int t3,int tk) { 
  return trans1.td[tk]*trans1.td[tk+t1]*trans2.td[tk+t1+parameter.t2]*trans2.td[tk+t1+parameter.t2+t3];}

double non_condon2d::gamma_TA(int t1,int t3) { 
  return exp( -(gamma1+gamma2)*t3/2.0 -gamma1*parameter.t2 - gamma1*t1/2.0 );}
double non_condon2d::gamma_SE(int t1,int t3) { 
  return exp( -(gamma1+gamma2)*t3/2.0 -gamma1*parameter.t2 - gamma1*t1/2.0 );}

//non-condon approximation
void non_condon2d::calculation()
{
  trans_mean();
  double domega=trans2.freq_mean-trans1.freq_mean;

  clock_t time[4];
  time[0] = clock();

  //*** calculate response function
  //*** based on ref1 and ref2
  double f=0.0,g10=0.0,g21=0.0;
  int tk=0;
  double delta;
  for(int k=0;k<parameter.np;k++)
  {
    tk=k*parameter.tgap;
    f=0.0;
    for(int t1=0;t1<tmax;t1++)
    {
      g10=0.0;
      g21=0.0;
      for(int t3=0;t3<tmax;t3++)
      {
        delta=domega*t3;
        r(t1,t3) += (2.0*td10(t1,t3,tk)*exp(double( f-g10)*1i)*gamma_SE(t1,t3) 
            - td21(t1,t3,tk)*exp(double( f-g21+delta)*1i)*gamma_TA(t1,t3))/double(parameter.np);
        nr(t1,t3)+= (2.0*td10(t1,t3,tk)*exp(double(-f-g10)*1i)*gamma_SE(t1,t3) 
            - td21(t1,t3,tk)*exp(double(-f-g21+delta)*1i)*gamma_TA(t1,t3))/double(parameter.np);
        g10+=trans1.freq[tk+t1+parameter.t2+t3];
        g21+=trans2.freq[tk+t1+parameter.t2+t3];
      }
      f+=trans1.freq[tk+t1];
    }
  }

  time[1]=clock();
  cout << "total time for response function:" 
    << double(time[1]-time[0]) / double(CLOCKS_PER_SEC) << "seconds"<< endl;
}

//cummulant2d method
cummulant2d::cummulant2d(para parameterin,transition transin1,transition transin2)
  :response2d(parameterin,transin1,transin2) {} ;

double cummulant2d::gamma_TA(int t1,int t3) { 
  return exp( -(gamma1+gamma2)*t3/2.0 -gamma1*parameter.t2 - gamma1*t1/2.0 );}
double cummulant2d::gamma_SE(int t1,int t3) { 
  return exp( -(gamma1+gamma2)*t3/2.0 -gamma1*parameter.t2 - gamma1*t1/2.0 );}

double* cummulant2d::correlation_function(string filename,transition transin1,transition transin2)
{
  int tk;
  double *g;
  g=new double [(int) transin1.freq.size()];
  double c[(int) transin1.freq.size()]={0};
  //calculate frequency fluctuation correlation function (FFCF)
  int gtmax=2*parameter.tmax+parameter.t2;
  for (int t=0;t<gtmax;t++)
  {
    c[t]=0.0;
    for (int k=0;k<parameter.np;k++)
    {
      tk=k*parameter.tgap;
      c[t]+=transin1.freq[t+tk]*transin2.freq[tk];
    }
    c[t]/=parameter.np;
  }

  for (int t=0;t<gtmax;t++)
  {
    g[t]=0.0;
    for (int tau1=0;tau1<t;tau1++)
    {
      for (int tau2=0;tau2<tau1;tau2++)
      {
        g[t]+=c[tau2];
       // cout << t << "  " << tau1 << " " << tau2 << "  "   << g[t] << endl;
      }
    }
  }
  cout << (int) transin1.freq.size() << endl;
  ofstream check(filename,ios::out);
  //for (int t=0;t<(int) trans2.freq.size();t++) check << t << "  " <<  g[t] << endl;
  for (int t=0;t<gtmax;t++) check << t << "  " <<  c[t] << endl;
  check.close();
   return g;
}

void cummulant2d::calculation()
{
  //condon approximation,calculate mu**2
  double td10,td21;
  td10=trans1.td[0]*trans1.td[0]*trans1.td[0]*trans1.td[0];
  td21=trans1.td[0]*trans2.td[0]*trans1.td[0]*trans2.td[0];
  trans_mean();

  //calculate response function

  clock_t time[4];
  time[0] = clock();
  string filename="correlation11";
  double *g11=correlation_function(filename,trans1,trans1);
  filename="correlation12";
  double *g12=correlation_function(filename,trans1,trans2);
  filename="correlation22";
  double *g22=correlation_function(filename,trans2,trans2);
  double domega=trans2.freq_mean-trans1.freq_mean;


  double delta;
  double G1,G2,G3,G4;
  int t2=parameter.t2;
  for (int t1=0;t1<tmax;t1++)
  {
    for (int t3=0;t3<tmax;t3++)
    {
      //gamma_TA = 1.0;
      //gamma_SE = 1.0;
      delta=domega*t3;
      G1= -g11[t1] +g11[t2] -g11[t3] -g11[t1+t2] -g11[t2+t3] +g11[t1+t2+t3];
      G2= -g11[t1] +g12[t2] -g22[t3] -g12[t1+t2] -g12[t2+t3] +g12[t1+t2+t3];
      G3= -g11[t1] -g11[t2] -g11[t3] +g11[t1+t2] +g11[t2+t3] -g11[t1+t2+t3];
      G4= -g11[t1] -g12[t2] -g22[t3] +g12[t1+t2] +g12[t2+t3] -g12[t1+t2+t3];
      r(t1,t3) =(2.0*td10*exp(G1)*gamma_SE(t1,t3) - td21*exp(G2)*exp(1i*delta)*gamma_TA(t1,t3));
      nr(t1,t3)=(2.0*td10*exp(G3)*gamma_SE(t1,t3) - td21*exp(G4)*exp(1i*delta)*gamma_TA(t1,t3));
    }
  }
  time[1]=clock();
  cout << "total time for response function:"
    << double(time[1]-time[0]) / double(CLOCKS_PER_SEC) << "seconds"<< endl;
}






