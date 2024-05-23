#include "class.h"


//*** method
void transition::add_freq(double freq_t)
{
  freq.push_back(freq_t);
}


//void transition::add_td(double td_t)
//{
//  td.push_back(td_t);
//}

double transition::fmean_cal(void)
{
  freq_mean = accumulate( freq.begin(),freq.end(),0.0) /freq.size();
  return freq_mean;
}

void transition::fmean_removal() {for (int i=0;i<(int) freq.size();i++) freq[i]=freq[i]-freq_mean;}

double transition::freq_max(void)
{
  double max=freq[0];
  for (int i=0;i<(int) freq.size();i++)
  {
    if(max<freq[i]) max=freq[i];
  }
  return max;
}

double transition::freq_min(void)
{
  double min=freq[0];
  for (int i=0;i<(int) freq.size();i++)
  {
    if(min>freq[i]) min=freq[i];
  }
  return min;
}
void transition::tunit(double ut) 
{
  T1*=ut;
}
void transition::printfreq(string filename)
{

  checkname check(filename);
  string fname=check.ifreplace();
  ofstream output(fname,ios::out);
  for (int i=0;i<(int) freq.size();i++) output << freq[i] << endl;
  output.close();
}
double transition::ens_mean_freq(int np,int tmax,int tgap)
{
  double emfreq=0;
  int nt=0;
  for (int t=0;t<tmax;t++)
  {
    for (int n=0;n<np;n++)
    {
      emfreq+=freq[t+n*tgap];
      nt++;
    }
  }
  return emfreq/nt;
}

void para::tunit(double ut) 
{
  dt*=ut*1e-3;
  t2*=ut;
}
void para::get_np(para parameter)
{
  np=int((parameter.n-dimension*tmax-parameter.t2)/parameter.tgap);
  cout << "total time:" << n*dt << "ps" << endl;
  cout << "max response time:" << tmax*dt << "ps" << endl;
  cout << "time interval in time average:" << tgap*dt*1e3 << "fs" << endl;
  if (np<=0) np=1;
}


//c++17 has error when compiling cython whl
//string checkname::ifreplace()
//{
//  int i=0;
//  string name2;
//  name2=name;
//  while(filesystem::exists(name2))
//  {
//    i++;
//    name2=name+"_"+to_string(i);
//  }
//  cout << "saved in "+name2 << endl;
//  return name2;
//}
string checkname::ifreplace()
{
  int i=0;
  string name2;
  name2=name;
  while(is_file_exist(name2))
  {
    i++;
    name2=name+"_" + to_string(i);
  }
  cout << "saved in "+  name2 << endl;
  return name2;
}
bool checkname::is_file_exist(string fileName)
{
  std::ifstream infile(fileName);
  return infile.good();
}





