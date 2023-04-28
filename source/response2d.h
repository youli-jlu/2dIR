#ifndef RESPONSE2d
#define RESPONSE2d

#include "class.h"  // transition and parameter determination
#include <iostream>
#include <complex>
#include <cmath>

//using namespace std;
class response2d
{
  public:
    //response2d(para parameterin,transition transin1,transition transin2):
    //  parameter(parameterin),trans1(transin1),trans2(transin2) ;
    response2d(para ,transition ,transition ) ;
    virtual ~response2d() ;
    para parameter;
    transition trans1,trans2;
    int tmax;
    int ntotal;
    int dimension;
    double gamma1;
    double gamma2;
    complex<double> *R_r=nullptr;
    complex<double> *R_nr=nullptr;

    virtual void initial() ;
    void printr(string);
    void printnr(string);
    complex<double>  &r(int,int) ;
    complex<double> &nr(int,int) ;
    void zero_padding1d(complex<double>*&,int);
    void trans_mean();
    void zero_padding() ;
    void initial_value() ;
    virtual void calculation()=0 ;
};

class condon2d:public response2d
{
  public:
    condon2d(para ,transition ,transition ) ; 
    ~condon2d() { };
    double gamma_TA(int,int);
    double gamma_SE(int,int);
    void calculation();
};

class non_condon2d:public response2d
{
  public:
    non_condon2d(para ,transition ,transition ) ; 
    ~non_condon2d() { };
    void calculation();
    double gamma_TA(int,int);
    double gamma_SE(int,int);
    double td10(int,int,int);
    double td21(int,int,int);
};
class cummulant2d:public response2d
{
  public:
    cummulant2d(para ,transition ,transition ) ; 
    ~cummulant2d() { };
    void calculation();
    double gamma_TA(int,int);
    double gamma_SE(int,int);
    double* correlation_function(string,transition,transition);
};



#endif
