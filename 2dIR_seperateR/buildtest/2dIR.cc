#include "2dIR.h"

using namespace std;

const double pi=3.141592653589793;

void spectra2d::calculation()
{
  //transition trans1,trans2;
  //read td and freq
  cout << "***reading frequency and transition***" << endl;
  readfreq();
  readtd(parameter.tdfile1,trans1);
  readtd(parameter.tdfile2,trans2);

  parameter.tunit(parameter.rtgap); //change unit to ps and also contain rtgap
  trans1.T1=parameter.trans1T1;
  trans2.T1=parameter.trans2T1;
  cout << "***basic parameters***" << endl;
  cout << "trans1 information:" << trans1.freq_max() << "  " << trans1.freq_min() << "  "<< trans1.fmean_cal() << endl;
  cout << "trans2 information:" << trans2.freq_max() << "  " << trans2.freq_min() << "  "<< trans2.fmean_cal() << endl;
  cout << "delta omega:" << trans2.freq_mean-trans1.freq_mean<< endl;
  parameter.momega1=trans1.fmean_cal();
  parameter.momega2=trans2.fmean_cal();


  //calculate ensemble parameter.
  parameter.get_np(parameter);
  cout << "ensemble max time:" << parameter.tmax << endl;
  //!!! End of reading part

  //!!! unit transformation
  //change time unit if needed
  //chang energy(frequency and time) to time-normaliztion unit
  for (int i=0;i<parameter.n;i++) {
    trans1.freq[i]=trans1.freq[i]/parameter.ps2cm*parameter.dt;
    trans2.freq[i]=trans2.freq[i]/parameter.ps2cm*parameter.dt;
  }

  //!!! End of unit transformation


  //!!! Response function calculation
  cout << "***begin response function calculation***"<< endl;
  response2d *R=nullptr;
  // use different derivative class based on response flag
  if (parameter.responseflag=="condon") {
    R= new condon2d(parameter,trans1,trans2); std::cout << "Doing 2d condon calculation "<< std::endl;}
  else if (parameter.responseflag=="non-condon") {
    R= new non_condon2d(parameter,trans1,trans2); std::cout << " Doing 2d non-condon calculation "<< std::endl;}
  else if (parameter.responseflag=="cumulant") {
    R= new cumulant2d(parameter,trans1,trans2); std::cout << " Doing 2d cumulant calculation "<< std::endl;}
  else {
    std::cerr << " wrong responseflag!!" << std::endl;
    exit(0);}

  R->initial();
  R->calculation();
  cout << "End of response calculation"<< endl;

  R->printr();
  //!!! End of response_function calculation

  //!!! Fourier transformation
  cout << "zero_padding numbers:" << parameter.nzp << endl;
  R->zero_padding();
  R->initial_value();

  cout << "***begin Fast Fourier Transformation***"<< endl;
  dft2d dft1(parameter,R->R1_r,R->R1_nr);
  dft1.spectra(parameter.spectrafile1);
  dft2d dft2(parameter,R->R2_r,R->R2_nr);
  dft2.spectra(parameter.spectrafile2);


  //!!! End of Fourier transformation
  delete R;
  return ;

}
 
void spectra2d::readfreq()
{
  ifstream freqfile(parameter.freqfile,ios::in); //freq-file
  if (!freqfile) {
    cerr << "open freqfile fail\n" ;
    cerr << "file:"<<parameter.freqfile<<"doesn't exist\n" << endl;
    exit(0);}
  parameter.n=0;                  //total time step
  double freq_t[2]; //tmple freq(t)
  string oneline; 
  getline(freqfile,oneline); // remove fisrt line
  while(freqfile.peek() != EOF ){
    parameter.n++;
    freqfile >> freq_t[0] >> freq_t[1] ;
    trans1.add_freq(freq_t[0]);
    trans2.add_freq(freq_t[1]);
    getline(freqfile,oneline); // remove other character
    for (int i=0;i<parameter.rtgap-1;i++) getline(freqfile,oneline); // skip oneline
  }
  freqfile.close();
}

void spectra2d::readtd(string filename,transition &transi)
{
  ifstream tdfile(filename,ios::in);
  string oneline;
  int nd=0;
  double td[3];
  if (!tdfile) {
    cerr << "open tdfile fail\n" ;
    cerr << "file:"<<parameter.tdfile1 <<"doesn't exist\n" << endl;
    exit(0);}
  getline(tdfile,oneline); // remove fisrt line
  if (parameter.responseflag=="non-condon")
  {
    while(tdfile.peek() != EOF ){
      nd++;
      tdfile >> td[0] >> td[1] >> td[2]  ;
      transi.tdx.push_back(td[0]); transi.tdy.push_back(td[1]); transi.tdz.push_back(td[2]);
      getline(tdfile,oneline); // remove other character
      for (int i=0;i<parameter.rtgap-1;i++) getline(tdfile,oneline); // skip oneline
    }
    if (nd<parameter.n) {
      cerr << "too less transition dipole\n" <<endl ;
      exit(0);
    }
  }
  else
  {
  //read only one transisition dipole
    tdfile >> td[0] >> td[1] >> td[2]  ;
    transi.tdx.push_back(td[0]); transi.tdy.push_back(td[1]); transi.tdz.push_back(td[2]);
  }
  tdfile.close();
}

