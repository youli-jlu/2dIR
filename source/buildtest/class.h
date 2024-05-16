#ifndef IR_PARA
#define IR_PARA

#include <iostream>
#include <vector>
#include <numeric>
#include <fstream>
#include <string>
#include <cstring>
#include <filesystem>

using namespace std;

class transition
{
  public:
    vector< double > freq;
    vector<double> tdx,tdy,tdz;
    double T1;
    double freq_mean;
    double real_mean_omega; //real mean frequency 


    //*** method
    void add_freq(double freq_t);
    //void add_td(double td_t);
    double fmean_cal(void);
    void fmean_removal();
    double freq_max(void);
    double freq_min(void);
    void tunit(double ut) ;
    void printfreq(string filename);
    double ens_mean_freq(int np,int tmax,int tgap);

};

class para
{
  public:
    //spectra category and response function
    int dimension; // freq dimension
    string responseflag;
    string spectraflag;
    string freqfile;
    string tdfile1;
    string tdfile2;
    string spectrafile;
    string responsefile;

    //basic parameter
    int n;    //total time step
    double dt;//time interval between adjacent step
    double rtgap; //read gap to save time
    double omega_min,omega_max; //frequency range in 2d-IR simulation
    double momega1,momega2; //frequency range in 2d-IR simulation
    double dw; //frequency range in 2d-IR simulation
    double domega; //momega2-momega1, with unit trans
    double t2;  //t2 time in 2d-IR simulation
    double tgap;//time interval between adjacent ensemble
    double nzp; // number of zero padding
    double trans1T1;
    double trans2T1;
    // ensemble parameter
    int tmax;
    int np=0;
    //constant parameter for unit transformation (if needed)
    const double h2cm=219474.63;
    const double ps2cm=1.0/0.18826966;

    void tunit(double ut); 
    void get_np(para parameter);

};

class checkname
{
  public:
    checkname(string namein):name(namein) {};
    ~checkname() {};
    string name;
    string ifreplace();
    bool is_file_exist(string);
};
        
      
    

#endif
