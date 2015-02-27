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
  int cpt= 0;
  std::vector<double> stats1(true);
  std::vector<double> stats2(true);
  trace.info() << "Loading data..."<<std::endl;
 
  while (values1.good())
  {
    values1>> val;
    stats1.push_back(val);
    ++cpt;
  }
  
  while(values2.good())
  {
      values2>> val;
      stats2.push_back(val);
      ++cpt;
  }
  values1.close();
  values2.close();
  trace.info()<< "Nb values= "<<cpt<<std::endl;
  
  std::vector<double>::const_iterator it1 = stats1.begin();
  std::vector<double>::const_iterator it2 = stats2.begin();
  double emd = 0.0;
  double cdf1 = 0.0;
  double cdf2 = 0.0;
  double integralcdf1 = 0.0, integralcdf2=0.0;
  double x, xprev=0.0;
  double dx; //distance between two events
  
  {
    //We get the next event
    if (std::min(*it1, *it2) == *it1)
    {
      //Data from (1)
      x = *it1;
      dx = x - xprev;
      integralcdf1 += cdf1*dx;
      integralcdf2 += cdf2*dx;
      cdf1++;
      xprev = x;
      ++it1;
    }
    else
    {
      //Data from (2)
      x = *it2;
      dx = x - xprev;
      integralcdf1 += cdf1*dx;
      integralcdf2 += cdf2*dx;
      cdf2++;
      xprev = x;
      ++it2
    }
    emd  +=  std::abs( cdf1 - cdf2 );
   
    
    
  }
  
  
  
  return 0;
}

