#ifndef SPECTRA1d
#define SPECTRA1d
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
#include "response1d.h"
#include "dft.h"

//For 1d-IR simulation. Written by You Li

using namespace std;

class spectra1d
{
  public:
    spectra1d(para* parameterin):parameter(*parameterin) {};
    ~spectra1d() {};
    para parameter;
    transition trans1;
    void calculation();
    void readfreq();
    void readtd();
    void readtp();
};


#endif
