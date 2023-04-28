#ifndef RESPONSE1d
#define RESPONSE1d

#include "class.h"  // transition and parameter determination
#include <iostream>
#include <complex>
#include <cmath>

//using namespace std;
class response1d
{
  public:
    //response1d(para parameterin,transition transin1,transition transin2):
    //  parameter(parameterin),trans1(transin1),trans2(transin2) ;
    response1d(para ,transition ) ;
    virtual ~response1d() ;
    para parameter;
    transition trans1;
    int tmax;
    int ntotal;
    int dimension;
    double gamma1;
    complex<double> *R_r=nullptr;

    virtual void initial() ;
    void printr(string);
    complex<double>  &r(int) ;
    void zero_padding1d(complex<double>*&,int);
    void zero_padding();
    void initial_value();
    void trans_mean();
    virtual void calculation()=0 ;
};

class condon1d:public response1d
{
  public:
    condon1d(para ,transition ) ; 
    ~condon1d() { };
    double gamma(int);
    void calculation();
};

class non_condon1d:public response1d
{
  public:
    non_condon1d(para ,transition ) ; 
    ~non_condon1d() { };
    void calculation();
    double gamma(int);
    double td10(int,int);
};
class cummulant1d:public response1d
{
  public:
    cummulant1d(para ,transition ) ; 
    ~cummulant1d() { };
    void calculation();
    double gamma(int);
    double* correlation_function(string,transition,transition);
};



#endif
