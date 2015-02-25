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
  //Loading values  argv[1]==data1  argv[2] == data2 argv[3] == nbBins
  
  ifstream values1 (argv[1], std::ifstream::in);
  ifstream values2 (argv[2], std::ifstream::in);

  double val;
  int cpt= 0;
  int bins = atoi(argv[3]);
  Statistic<double> stats;
  Statistic<double> stats1(true);
  Statistic<double> stats2(true);
  trace.info() << "Loading data..."<<std::endl;
  while (values1.good() || values2.good())
  {
    if (values1.good())
    {
      values1>> val;
      stats.addValue(val);
      stats1.addValue(val);
      ++cpt;
    }
    if (values2.good())
    {
      values2>> val;
      stats.addValue(val);
      stats2.addValue(val);
      ++cpt;
    }
  }
  values1.close();
  values2.close();
  stats.terminate();
  stats1.terminate();
  stats2.terminate();
  trace.info()<< "Nb values= "<<cpt<<std::endl;
  trace.info() << "Global Statistic= "<<stats<<std::endl;
  trace.info() << "Statistic 1= "<<stats1<<std::endl;
  trace.info() << "Statistic 2= "<<stats2<<std::endl;
  
  Histogram<double> hist1;
  hist1.init( bins, stats1);
  hist1.addValues(stats1.begin(), stats1.end());
  hist1.terminate();
  trace.info()<< "Histogram 1= "<<hist1<<std::endl;
  
  Histogram<double> hist2;
  hist2.init( bins, stats2);
  hist2.addValues(stats2.begin(), stats2.end());
  hist2.terminate();
  trace.info()<< "Histogram 2= "<<hist2<<std::endl;
  
  //EMD computation
  double EMD=0;
  for(unsigned int i = 0 ; i< hist1.size(); ++i)
    EMD += std::abs(hist1.cdf(i) - hist2.cdf(i));
  
  std::cout <<  "EMD= "<< EMD<<std::endl;
  return 0;
}

