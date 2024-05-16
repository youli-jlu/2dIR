#include "response1d.h"
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
//reference:1: doi: 10.1063/1.1514652

//youli: some intial method for response1d class
response1d::response1d(para parameterin,transition transin1):
  parameter(parameterin),trans1(transin1) 
{
  tmax=parameter.tmax;
  ntotal=pow(tmax,parameter.dimension);
  dimension=parameter.dimension;
  gamma1=1.0/(trans1.T1/parameter.dt);
  R_r=new complex<double>[ntotal];
}
response1d::~response1d() {};
void response1d::initial() {};

complex<double> & response1d::r(int t1)
{
  int index=t1;
  return R_r[index];
}

double response1d::td10(int t1,int tk) { 
  return trans1.tdx[tk]*trans1.tdx[tk+t1] 
    +trans1.tdy[tk]*trans1.tdy[tk+t1]
    +trans1.tdz[tk]*trans1.tdz[tk+t1]
    ;}

//void response1d::printr(string filename)
//{
//  checkname check(filename);
//  string fname=check.ifreplace();
//  ofstream outfile(fname,ios::out);
//  outfile << "t" << "  "  << "re"  <<  "  "  << "im" << endl; 
//  for (int i=0;i<tmax;i++)
//  {
//    outfile <<i*parameter.dt<< "  "  << r(i).real()  <<  "  "  << r(i).imag() << endl; 
//  }
//  outfile.close();
//  //return fname;
//}
void response1d::printr()
{
  ofstream outfile(parameter.responsefile,ios::out);
  outfile << "t" << "  "  << "re"  <<  "  "  << "im" << endl; 
  for (int i=0;i<tmax;i++)
  {
    outfile <<i*parameter.dt<< "  "  << r(i).real()  <<  "  "  << r(i).imag() << endl; 
  }
  cout << "saving response function in   " + parameter.responsefile << endl;
  outfile.close();
  //return fname;
}

void response1d::zero_padding1d(complex<double> *&pR,int np)
{
  complex<double> *pR_new=nullptr;
  int tmax_new=tmax+np;
  int ntotal_new=pow(tmax_new,dimension);

  pR_new= new complex<double>[ntotal_new];
  for (int i=0;i<tmax;i++)
  {
      pR_new[i]=pR[i];
  }
  for (int i=tmax;i<tmax+np;i++)
  {
      pR_new[i]=0.0+0.0*1i;
  }
  delete [] pR;
  pR=pR_new;
}
void response1d::zero_padding()
{
  if (parameter.nzp>0){
    zero_padding1d(R_r,parameter.nzp);
    tmax+=parameter.nzp;
    ntotal=pow(tmax,dimension);
  }
}
void response1d::trans_mean()
{
  trans1.fmean_cal();
  trans1.fmean_removal();
  return;
}
void response1d::initial_value()
{
  r(0)/=2.0;}
  



//condon1d method
condon1d::condon1d(para parameterin,transition transin1)
  :response1d(parameterin,transin1) {} ;

double condon1d::gamma(int t1) { return exp(-gamma1*t1);}

void condon1d::calculation()
{
  //condon approximation,calculate mu10**4 and mu10**2*mu21**2
  double td10_c;
  td10_c=td10(0,0);
  trans_mean();

  clock_t time[4];
  time[0] = clock();

  //*** calculate response function
  //*** based on ref1 and ref2
  double f=0.0;
  int tk=0;
  for(int k=0;k<parameter.np;k++)
  {
    tk=k*parameter.tgap;
    f=0.0;
    for(int t1=0;t1<tmax;t1++)
    {
      r(t1) += (td10_c*exp(f*1i)*gamma(t1))/double(parameter.np);
      f+=trans1.freq[tk+t1];
    }
  }


  time[1]=clock();
  cout << "total time for response function:" 
    << double(time[1]-time[0]) / double(CLOCKS_PER_SEC) << "seconds"<< endl;
}

//non_condon1d method
non_condon1d::non_condon1d(para parameterin,transition transin1)
  :response1d(parameterin,transin1) {} ;

double non_condon1d::gamma(int t1) { return exp(-gamma1*t1);}

//non-condon 
void non_condon1d::calculation()
{
  trans_mean();

  clock_t time[4];
  time[0] = clock();

  //*** calculate response function
  //*** based on ref1 and ref2
  double f=0.0;
  int tk=0;
  for(int k=0;k<parameter.np;k++)
  {
    tk=k*parameter.tgap;
    f=0.0;
    for(int t1=0;t1<tmax;t1++)
    {
      r(t1) += (td10(t1,tk)*exp(f*1i)*gamma(t1))/double(parameter.np);
      f+=trans1.freq[tk+t1];
    }
  }

  time[1]=clock();
  cout << "total time for response function:" 
    << double(time[1]-time[0]) / double(CLOCKS_PER_SEC) << "seconds"<< endl;
}

//cummulant1d method
cummulant1d::cummulant1d(para parameterin,transition transin1)
  :response1d(parameterin,transin1) {} ;

double cummulant1d::gamma(int t1) { return exp(-gamma1*t1);}
double* cummulant1d::correlation_function(string filename,transition transin1,transition transin2)
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
      }
    }
  }
  cout << (int) transin1.freq.size() << endl;
  ofstream check(filename,ios::out);
  for (int t=0;t<gtmax;t++) check << t << "  " <<  c[t] << endl;
  check.close();
   return g;
}

void cummulant1d::calculation()
{
  //condon approximation,calculate mu**2
  double td10_c;
  td10_c=td10(0,0);
  trans_mean();

  //calculate response function

  clock_t time[4];
  time[0] = clock();
  string filename="correlation11";
  double *g11=correlation_function(filename,trans1,trans1);


  for (int t1=0;t1<tmax;t1++)
  {
      r(t1) =td10_c*exp(-g11[t1])*gamma(t1);
      //cout << gamma(t1) << endl;
  }
  time[1]=clock();
  cout << "total time for response function:"
    << double(time[1]-time[0]) / double(CLOCKS_PER_SEC) << "seconds"<< endl;
}

