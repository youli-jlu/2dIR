#include "1dIR.h"

using namespace std;


void spectra1d::calculation()
{
  //transition trans1,trans2;
  const double pi=3.141592653589793;

  //read td and freq
  cout << "***reading frequency and transition***" << endl;
  readtd();
  readfreq();
  parameter.tunit(parameter.rtgap); //change unit to ps and also contain rtgap


  cout << "***basic parameters***" << endl;
  trans1.T1=parameter.trans1T1;
  cout << "trans information:" << trans1.freq_max() << "  " << trans1.freq_min() << "  "<< trans1.fmean_cal() << endl;
  parameter.momega1=trans1.fmean_cal();


  //calculate ensemble parameter.
  parameter.get_np(parameter);
  cout << "ensemble numbers:"<< parameter.np << endl;
  cout << "ensemble max time:" << parameter.tmax << endl;
  //!!! End of reading part

  //!!! unit transformation
  //parameter.tunit(parameter.rtgap*2.419e-5); //change unit to ps 
  //chang energy(frequency and time) to time-normaliztion unit
  for (int i=0;i<parameter.n;i++) {
    trans1.freq[i]=trans1.freq[i]/parameter.ps2cm*parameter.dt;
  }

  //!!! End of unit transformation


  //!!! Response function calculation
  cout << "***begin response function calculation***"<< endl;
  response1d *R=nullptr;
  // use different derivative class based on response flag
  if (parameter.responseflag=="condon") {
    R= new condon1d(parameter,trans1); std::cout << "Doing 1d condon calculation "<< std::endl;}
  else if (parameter.responseflag=="non-condon") {
    R= new non_condon1d(parameter,trans1); std::cout << " Doing 1d non-condon calculation "<< std::endl;}
  else if (parameter.responseflag=="cummulant") {
    R= new cummulant1d(parameter,trans1); std::cout << " Doing 1d cummulant calculation "<< std::endl;}
  else {
    std::cerr << " wrong responseflag!!" << std::endl;return ;}

  R->initial();
  R->calculation();
  R->printr("1dresponse");
  //!!! End of response_function calculation

  //!!! Fourier transformation
  cout << "zero_padding numbers:" << parameter.nzp << endl;
  R->zero_padding();
  //R->printr ("tmp1_zp");

  cout << "***begin Fast Fourier Transformation***"<< endl;
  dft1d dft1(parameter,R->R_r);
  dft1.spectra();


  //!!! End of Fourier transformation
  delete R;
  return ;

}
 
void spectra1d::readfreq()
{
  ifstream freqfile(parameter.freqfile,ios::in); //freq-file
  if (!freqfile) {
    cerr << "open freqfile fail\n" ;
    cerr << "file:"<<parameter.freqfile<<"doesn't exist\n" << endl;
    return ;
  }
  parameter.n=0;                  //total time step
  double freq_t[1]; //tmple freq(t)
  string oneline; 
  getline(freqfile,oneline); // remove fisrt line
  while(freqfile.peek() != EOF ){
    parameter.n++;
    freqfile >> freq_t[0]  ;
    trans1.add_freq(freq_t[0]);
    getline(freqfile,oneline); // remove other character
    for (int i=0;i<parameter.rtgap-1;i++) getline(freqfile,oneline); // skip oneline
  }
  freqfile.close();
}

void spectra1d::readtd()
{
  ifstream tdfile(parameter.tdfile,ios::in);
  string oneline;
  int nd=0;
  double td[1];
  if (!tdfile) {
    cerr << "open tdfile fail\n" ;
    cerr << "file:"<<parameter.tdfile <<"doesn't exist\n" << endl;
    return ;
  }
  getline(tdfile,oneline); // remove fisrt line
  if (parameter.responseflag=="non-condon")
  {
    while(tdfile.peek() != EOF ){
      nd++;
      tdfile >> td[0]  ;
      trans1.add_td(td[0]);
      getline(tdfile,oneline); // remove other character
      for (int i=0;i<parameter.rtgap-1;i++) getline(tdfile,oneline); // skip oneline
    }
    if (nd<parameter.n) {
      cerr << "too less transition dipole\n" <<endl ;
      return ;
    }
  }
  else
  {
  //read only one transisition dipole
    tdfile >> td[0]  ;
    trans1.add_td(td[0]);
  }
  tdfile.close();
}
