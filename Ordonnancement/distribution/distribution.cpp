#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <utility>
#include <algorithm>

/*************************************\
|  Compute the repartition of radius  |
|     corresponding to the input      |
| (thickness normalized distribution) |
\*************************************/



// argv[1] = name of the file which contain the distribution data ( .txt )
// argv[2] = precision p ( intervals : [0,p] [p,2p] ... [1-1/p,1] )

std::vector<double> extractDistribution(char * filename) ;
std::vector<int> computeRepartition(std::vector<double>& thick, int precision );
int writeInFile(std::vector<int>& repartition, std::string dataFilename, int precision, int size, double max );
int writeScript(std::string dataFilename, int precision,double max);
std::pair<double,double> expectationAndVariance ( std::vector<int>& repartition, int precision, double max );

int main ( int argc , char * argv[] ) {
  if ( argc != 3 ) {
      std::cerr << "Not the good number of arguments" << std::endl ;
      return -1 ;
  }
  int precision = atoi(argv[2]) ;
  // Extract data
  std::vector<double> thick = extractDistribution(argv[1]) ; 
  // sort the values
  std::sort(thick.begin(),thick.end());
  // Compute the repartition
  std::vector<int> repartition = computeRepartition(thick,precision) ;
  std::string filename = argv[1] ;
  std::string dataFilename = filename.substr(0,filename.size()-4)+"-rep.txt" ;
  // Write the repartition in a file
  int size = thick.size() ;
  double max = thick[size-1] ;
  int writeF = writeInFile(repartition,dataFilename,precision,size,max); 
  if (writeF != 0) {
      return -1 ;
  }
  // Create a script to plot it
  int writeS = writeScript(dataFilename,precision,max) ;
  if (writeS != 0 ) {
      return -1 ;
  }
  // Compute Expectation and Variance
  std::pair<double,double> values = expectationAndVariance(repartition,precision,max) ;
  // Print
  std::cout << "Expectation : " << values.first << ", variance : " << values.second << std::endl ;
  return 0 ;
}

/********************\
| AUXILIAR FUNCTIONS |
\********************/
 

double stringToDouble ( std::string word ) {
    return atof(word.c_str()) ;
}

std::vector<double> extractDistribution(char * filename) {
  std::cout << "Read data ... " ;
  std::vector<double> thick ;
  std::ifstream input(filename) ;
  if (!input) {
      std::cerr << "Error when opening input file" << std::endl ;
      return thick ; // empty 
  }
  std::string word ;
  int isData = 0 ;
  while(input >> word) {
      if (isData) {
          thick.push_back(stringToDouble(word)) ;
      }
      isData = 1 - isData ;
  } 
  std::cout<< "Ok" << std::endl ;
  return thick ;
}    



std::vector<int> computeRepartition(std::vector<double>& thick, int precision ) {
  std::cout << "Compute repartition ... " ;
  std::vector<int> repartition(precision,0) ;
  int i = 0 ;
  int size = thick.size() ;
  double max = thick[size-1] ;
  for( i = 0 ; i < size; i++ ) {
       int indice = (int) ( ceil( thick[i] / max * ((double) precision))) - 1 ;
       repartition[std::max(indice,0)]++ ;
  } 
  std::cout << "Ok" << std::endl ; 
  return repartition ;
}

int writeInFile(std::vector<int>& repartition, std::string dataFilename, int precision, int size, double max) {
  std::cout << "Write result in file ..." ;
  std::ofstream file(dataFilename.c_str(), std::ios::out);
  if (!file){
      std::cerr << "open file error" << std::endl ;
      return -1 ; 
  } 
  int i ;
  for ( i = 0 ; i < precision ; i++ ) {
        double mediumValue = (((double) i)+0.5)/((double) precision)* max    ;
        file << mediumValue << " " << repartition[i]/((double) size) * 100 << std::endl ;
  }
  std::cout << "Ok" << std::endl ;
  return 0 ;
}

int writeScript(std::string dataFilename, int precision, double max) {
    std::cout << "Write gnuplot script ... " ;
    std::string scriptName = dataFilename.substr(0,dataFilename.size()-8)+"-script.p" ;
    std::ofstream scriptFile(scriptName.c_str(),std::ios::out) ;
    if ( !scriptFile ) {
        std::cerr << "Open script file error" << std::endl ;
        return -1 ;
    }  
    scriptFile << "reset " << std::endl ;
    scriptFile << std::endl ;
    scriptFile << "set term pdfcairo" << std::endl ;
    scriptFile << "set output \"repartition.pdf\"" << std::endl ; // PDF file name
    scriptFile << "set title \" Thickness normalized repartition \" " << std::endl ;
    scriptFile << "set xlabel \" Radius  \" " << std::endl ;
    scriptFile << "set ylabel \" Frequency (%) \" " << std::endl ;
    scriptFile << std::endl ;
    scriptFile << "set style data line " << std::endl ;
    scriptFile << "set pointsize 2 " << std::endl ; 
    scriptFile << "set key top left " << std::endl ;
    scriptFile << "set boxwidth " << max /((double) precision) << std::endl ; // histogramm
    scriptFile << "set style fill solid border 3" << std::endl ;
    scriptFile << "plot \"" << dataFilename << "\" with boxes title \"repartition\"" ;
    std::cout << "Ok" << std::endl ;
    return 0 ;
}

std::pair<double,double> expectationAndVariance ( std::vector<int>& repartition, int precision, double max ) {
    double expectation = 0 ;
    double variance = 0 ;    // we use Var(X) = E[X^2] - E[X]^2 
    double sample = 0 ;
    int i ;
    for ( i = 0 ; i < precision ; i++ ) {
        expectation += repartition[i] * (((double) i)+0.5)/((double) precision)  ;
        variance += (((double) i)+0.5)/((double) precision) * (((double) i)+0.5)/((double) precision) * repartition[i] ;
        sample += repartition[i] ;
    } 
    expectation /= sample ;
    variance /= sample ;
    variance -= expectation * expectation ;
    expectation *= max ;
    variance *= max * max ;
    std::pair<double,double> values(expectation,variance) ;
    return values ;
}


