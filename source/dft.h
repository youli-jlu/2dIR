#ifndef DFT
#define DFT
#include <iostream>
#include "class.h"  // transition and parameter determination
#include <complex.h>
#include <cmath>
#include <fftw3.h>
//using namespace std;


class dft
{
  public:
    dft(para);
    virtual ~dft();
    para parameter;
    int tmax;
    int dimension;
    int ntotal;
};

class dft1d:public dft
{
  public:
    dft1d(para,complex<double>*);
    ~dft1d() {};
    complex<double>* r=nullptr;
    void spectra();
    double* time2freq(complex<double>*);
};

class dft2d:public dft
{
  public:
    dft2d(para,complex<double>*,complex<double>*);
    ~dft2d() {};
    complex<double>* r=nullptr;
    complex<double>* nr=nullptr;
    void spectra();
    int index(int,int);
    double* time2freq(complex<double>*);
};


#endif
