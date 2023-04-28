#include "dft.h"
#include <iomanip>

using namespace std;

//using FFTW to transform response function to frequency space.
// ref 3 (Zenni's book) provide basic concept about FFTW.


dft::dft(para parameterin):parameter(parameterin)
{
  //contain zero point
  if (parameter.nzp>0) tmax=parameter.tmax+parameter.nzp;
  else tmax=parameter.tmax;
  cout << " frequency grid = " << tmax << endl;
  dimension=parameter.dimension;
  ntotal=pow(tmax,dimension);
}
dft::~dft() {};


//*** 2d-DFT part
dft2d::dft2d(para parameterin,complex<double>* rin,complex<double>* nrin):
  dft(parameterin),r(rin),nr(nrin) {};

int dft2d::index(int i,int j)
{
  return j+i*tmax;
}


double* dft2d::time2freq(complex<double>* inc)
{
  int nt=tmax;
  int ntot=ntotal;
  complex<double>* in2=nullptr;
  in2 = new complex<double>[ntot];
  in2 = inc;
  cout << "response time" << nt << endl;
  cout << "ntot" << ntot << endl;
  //ofstream checkfile("check",ios::app);
  //print rephasing

  fftw_plan plan;
  fftw_complex *in,*out;
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntot);
  out= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntot);

  in=reinterpret_cast<fftw_complex*>(in2);
  plan = fftw_plan_dft_2d(nt,nt,in,out,FFTW_FORWARD,FFTW_ESTIMATE);

  fftw_execute(plan);

  //print rephasing freq
  //checkfile << "re  " << endl;
  //for (int i=0;i<nt;i++){
  // for (int j=0;j<nt;j++){
  //   checkfile << index(i,j) << " "<< in[index(i,j)][0] << "  "  << in[index(i,j)][1] << endl;
  //   checkfile << index(i,j) << " "<< out[index(i,j)][0] << "  "  << out[index(i,j)][1] << endl;
  // }}

  double* rw=nullptr;
  rw=new double[ntot];
  //reindex output and save only real part
  //for (int i=0;i<ntot;i++) checkfile << i << " " << out[i][0] << endl;
  int tw1,tw2;
  int j=0;
  for (int w1=0;w1<nt;w1++) 
  {
    tw1=w1+nt/2;
    if (tw1>=nt) tw1-=nt;
    for (int w2=0;w2<nt;w2++)
    {
      tw2=w2+nt/2;
      if (tw2>=nt) tw2-=nt;
      j=tw2+nt*tw1;
      rw[index(w1,w2)]=double(out[j][0]);
    }
  }
  fftw_destroy_plan(plan);
  fftw_free(in); 
  fftw_free(out);
  //checkfile.close();
  return rw;
}
  

void dft2d::spectra()
{

  double const pi=M_PI;
  clock_t time[2];
  time[0]=clock();
  double *rw=nullptr;
  double *nrw=nullptr;
  cout << "dft1" << endl;
  rw=time2freq(r);

  cout << "dft2" << endl;
  nrw=time2freq(nr);

  time[1]=clock();

  ofstream spectra(parameter.spectrafile,ios::out);
  spectra << "total DFI time:" << double(time[1]-time[0]) / double(CLOCKS_PER_SEC) << "seconds" << endl;
  cout << "total DFI time:" << double(time[1]-time[0]) / double(CLOCKS_PER_SEC) << "seconds" << endl;

  //remember our dt is in ps unit. 
  //change corresponding frequency from ps-1 to cm-1
  parameter.dw = 2.0*pi/(tmax*parameter.dt)*parameter.ps2cm;
  spectra << " frequency resolution:" << parameter.dw << endl;
  spectra << "frequency grid number:" << tmax << endl;
  cout  << " frequency resolution:" << parameter.dw << endl;
  cout  << "frequency grid number:" << tmax << endl;

  int nmax,nmin;
  nmax=int(parameter.omega_max/parameter.dw+tmax/2);
  nmin=int(parameter.omega_min/parameter.dw+tmax/2);
  if (nmax>tmax) nmax=tmax;
  if (nmin<0) nmin=0;

  cout << "print spectra" << endl;
  spectra << "w1    w3     rephasing   non-rephasing  total:" << tmax << endl;
  int indexr,indexnr;
  for (int w1=nmin;w1<nmax;w1++) 
  {
    for (int w3=nmin;w3<nmax;w3++) 
    {
      spectra<< setprecision(10)  << parameter.momega1+ parameter.dw*(w1-tmax/2) << "  " ;
      spectra<< setprecision(10)  << parameter.momega1+ parameter.dw*(w3-tmax/2) << "  " ;
      indexr= w3+(tmax-w1-1)*tmax;
      indexnr= w3+w1*tmax;
      spectra << rw[indexr]  << "  " << endl; 
      spectra << nrw[indexnr]  << "  " << endl; 
      spectra << rw[indexr] + nrw[indexnr]  << "  " << endl; 
    }
  }
  spectra.close();

}



//*** 1d-DFT part
dft1d::dft1d(para parameterin,complex<double>* rin):
  dft(parameterin),r(rin) {};
double* dft1d::time2freq(complex<double>* inc)
{
  int nt=tmax;
  int ntot=ntotal;
  complex<double>* in2=nullptr;
  in2 = new complex<double>[ntot];
  in2 = inc;
  cout << "response time" << nt << endl;
  cout << "ntot" << ntot << endl;
  //ofstream checkfile("checkdft",ios::app);
  //print rephasing

  fftw_plan plan;
  fftw_complex *in,*out;
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntot);
  out= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntot);

  in=reinterpret_cast<fftw_complex*>(in2);
  plan = fftw_plan_dft_1d(nt,in,out,FFTW_FORWARD,FFTW_ESTIMATE);

  fftw_execute(plan);

  //print rephasing freq
  //checkfile << "t inre  inim  outre outim" << endl;
  //for (int i=0;i<nt;i++){
  //   checkfile << i << " "<< in[i][0] << "  "  << in[i][1] << 
  //     " "<< out[i][0] << "  "  << out[i][1] << endl;
  //}

  double* rw=nullptr;
  rw=new double[ntot];
  //reindex output and save only real part
  int tw1;
  for (int w1=0;w1<nt;w1++) 
  {
    tw1=w1+nt/2;
    if (tw1>=nt) tw1-=nt;
    rw[w1]=double(out[tw1][0]);
  }
  fftw_destroy_plan(plan);
  fftw_free(in); 
  fftw_free(out);
  //checkfile.close();
  return rw;
}
  

void dft1d::spectra()
{

  double const pi=M_PI;
  clock_t time[2];
  time[0]=clock();
  double *rw=nullptr;
  cout << "dft1" << endl;
  rw=time2freq(r);

  time[1]=clock();

  ofstream spectra(parameter.spectrafile,ios::out);
  spectra << "total DFI time:" << double(time[1]-time[0]) / double(CLOCKS_PER_SEC) << "seconds" << endl;

  //remember our dt is in ps unit. 
  //change corresponding frequency from ps-1 to cm-1
  parameter.dw = 2.0*pi/(tmax*parameter.dt)*parameter.ps2cm;
  spectra << " frequency resolution:" << parameter.dw << endl;
  spectra << "frequency grid number:" << tmax << endl;
  cout  << " frequency resolution:" << parameter.dw << endl;
  cout  << "frequency grid number:" << tmax << endl;

  int nmax,nmin;
  nmax=int(parameter.omega_max/parameter.dw+tmax/2);
  nmin=int(parameter.omega_min/parameter.dw+tmax/2);
  if (nmax>tmax) nmax=tmax;
  if (nmin<0) nmin=0;

  cout << "spectra" << endl;
  for (int w1=nmin;w1<nmax;w1++) 
  {
      spectra<< setprecision(10)  << parameter.momega1+ parameter.dw*(w1-tmax/2) << "  " ;
      spectra << rw[w1] << endl;
  }
  spectra.close();

}

