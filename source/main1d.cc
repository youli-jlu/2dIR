#include <iostream>
#include <iomanip>
#include <cstring>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <algorithm>
#include "class.h"  // transition and parameter determination
#include "1dIR.h"

using namespace std;

const double pi=3.141592653589793;

int main()
{
  transition trans1;
  para parameter;


  //!!!reading part
  ifstream parafile("parameter-1d.dat",ios::in);     //para-file

  if (!parafile) {
    cerr << "open parafile fail\n";
    return 1;
  }

  parafile >> parameter.responseflag;
  parafile >> parameter.spectraflag;
  parafile >> parameter.freqfile;       //frequency filename
  parafile >> parameter.spectrafile;       //frequency filename
  parafile >> parameter.tdfile;       //transitionn dipole filename
  parafile >> parameter.dimension;       //frequency dimension
  parafile >> parameter.dt;       //time interval between adjacent step
  parafile >> parameter.rtgap;      // read gap
  parafile >> parameter.tgap;     //time interval between adjacent ensemble 
  parafile >> parameter.tmax;      // relaxtion time of transition2
  parafile >> parameter.nzp;      // zero-podding size
  parafile >> parameter.omega_min;//frequency range in 2d-IR simulation
  parafile >> parameter.omega_max;
  parafile >> parameter.trans1T1;  // relaxtion time of transition
  parafile >> parameter.trans2T1;  // relaxtion time of transition





  para* parameterin=&parameter;
  spectra1d spectra(parameterin);
  spectra.calculation();


  parameterin=nullptr;
  return 0;



}


