#include <iostream>
#include <fstream>
#include <vector>

#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/math/Histogram.h>

using namespace DGtal;
using namespace std;


/**
 * Compute the EMD between to set of values
 * computes first two histograms with same max/min
 *
 */
int main(int argc, char**argv)
{
  //Loading values  argv[1]==data1  argv[2] == data2
  ifstream values1 (argv[1], std::ifstream::in);
  ifstream values2 (argv[2], std::ifstream::in);
  
  double val;
  std::vector<double> stats1;
  std::vector<double> stats2;
  trace.info() << "Loading data..."<<std::endl;
  
 
  while (values1.good())
  {
    values1>> val;
    stats1.push_back(val);
  }
  
  while(values2.good())
  {
    values2>> val;
    stats2.push_back(val);
  }
  values1.close();
  values2.close();
  trace.info()<< "Nb values, stat1="<<stats1.size()<<" stat2="<<stats2.size()<<std::endl;
  
  double emd = 0.0;
  double cdf1 = 0.0;
  double cdf2 = 0.0;
  double x, xprev=0.0;
  double dx; //distance between two events
  
  double m1 = 1.0;//(double)stats1.size();
  double m2 = 1.0;//;(double)stats2.size();
  
  std::vector<double>::const_iterator it1 = stats1.begin();
  std::vector<double>::const_iterator it2 = stats2.begin();
  
  
  if (*it1 == std::min(*it1,*it2))
  {
    xprev = *it1;
    cdf1+= m1;
    ++it1;
  }
  else
  {
    xprev = *it2;
    cdf2+= m2;
    ++it2;
  }
  
  while ((it1 != stats1.end()) && (it2 !=stats2.end()))
  {
    //We get the next event
    if (std::min(*it1, *it2) == *it1)
    {
      //Data from (1)
      x = *it1;
      dx = x - xprev;
      emd  +=  std::abs( cdf1/(double)stats1.size() - cdf2/(double)stats2.size() )*dx;
      cdf1+= m1;
      xprev = x;
      ++it1;
    }
    else
    {
      //Data from (2)
      x = *it2;
      dx = x - xprev;
      emd  +=  std::abs( cdf1/(double)stats1.size() - cdf2/(double)stats2.size() )*dx;
      cdf2+= m2;
      xprev = x;
      ++it2;
    }
    //std::cout << cdf1<<" "<<cdf2<<" "<<emd<< std::endl;
  }
  
  //Extra work for stat1
  while (it1 != stats1.end())
  {
    x = *it1;
    dx = x - xprev;
    emd  +=  std::abs( cdf1/(double)stats1.size() - cdf2/(double)stats2.size()  )*dx ;
    cdf1+= m1;
    xprev = x;
    ++it1;
    //std::cout << cdf1<<" "<<cdf2<<" "<<emd<< std::endl;
  }
  
  
  //Extra work for stat2
  while (it2 != stats2.end())
  {
    x = *it2;
    dx = x - xprev;
    emd  +=  std::abs( cdf1/(double)stats1.size() - cdf2/(double)stats2.size()  )*dx ;
    cdf2+= m2;
    xprev = x;
    ++it2;
   // std::cout << cdf1<<" "<<cdf2<<" "<<emd<< std::endl;
    
  }
  
  std::cout << "EMD = "<<emd<<std::endl;
  
  
  
  return 0;
}

