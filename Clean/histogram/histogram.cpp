#include <iostream>
#include <fstream>
#include <vector>

#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/math/Histogram.h>

using namespace DGtal;
using namespace std;


/**
 * Compute the histogram on a given number if bins of a set of data
 *
 */
int main(int argc, char**argv)
{
  //Loading values  argv[1]==values  argv[2] == nbBins
  ifstream values (argv[1], std::ifstream::in);
  double val;
  int cpt= 0;
  int bins = atoi(argv[2]);
  Histogram<double> hist;
  Statistic<double> stats(true);
  while (values.good())
  {
    values>> val;
    stats.addValue(val);
    ++cpt;
  }
  values.close();
  stats.terminate();
  trace.info()<< "Nb values= "<<cpt<<std::endl;
  trace.info() << "Statistic= "<<stats<<std::endl;
  
 
  hist.init( bins, stats);
  hist.addValues(stats.begin(), stats.end());
  hist.terminate();
  trace.info()<< "Histogram= "<<hist<<std::endl;
  std::string outputName = "histogram.dat";
  std::ofstream output(outputName,std::ios::out) ;
  
  for ( unsigned int i = 0; i < hist.size(); ++i )
    output << i << " " << hist.pdf( i ) << std::endl;

  output.close();
  return 0;
}

