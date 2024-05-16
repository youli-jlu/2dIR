#ifndef SPECTRA2d
#define SPECTRA2d
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
#include "response2d.h"
#include "dft.h"

//For 2d-IR simulation. Written by You Li
//referenc:J. Phys. Chem. B,113,13118 ((11),(12),(13),(14),(15),(16))
//I also make a note for the details.

using namespace std;

class spectra2d
{
  public:
    spectra2d(para* parameterin):parameter(*parameterin) {};
    ~spectra2d() {};
    para parameter;
    transition trans1;
    transition trans2;
    void calculation();
    void readfreq();
    void readtd(string,transition&);
};


#endif
