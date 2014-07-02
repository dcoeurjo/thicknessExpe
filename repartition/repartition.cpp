#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <utility>

/***********************************\
| Compute the repartition of radius |
|    corresponding to the input     |
|    ( thickness distribution )     |
\***********************************/



// argv[1] = name of the file which contain the distribution data ( .txt )

std::vector<double> extractDistribution(char * filename) ;
std::vector<std::pair<double,int> > computeRepartition(std::vector<double>& thick );
int writeInFile(std::vector<std::pair<double,int> >& repartition, std::string dataFilename);
int writeScript(std::string dataFilename);
std::pair<double,double> expectationAndVariance ( std::vector<std::pair<double,int> >& repartition );

int main ( int argc , char * argv[] ) {
  // argv[1] = name of the file which contain the distribution data ( .txt )
  if ( argc != 2 ) {
      std::cerr << "Not the good number of arguments" << std::endl ;
      return -1 ;
  }
  // Extract data
  std::vector<double> thick = extractDistribution(argv[1]) ;
  // Compute the repartition
  std::vector<std::pair<double,int> > repartition = computeRepartition(thick) ;
  std::string filename = argv[1] ;
  std::string dataFilename = filename.substr(0,filename.size()-4)+"-rep.txt" ;
  // Write the repartition in a file
  int writeF = writeInFile(repartition,dataFilename); 
  if (writeF != 0) {
      return -1 ;
  }
  // Create a script to plot it
  int writeS = writeScript(dataFilename) ;
  if (writeS != 0 ) {
      return -1 ;
  }
  // Compute and plot Expectation and Variance
  std::pair<double,double> values = expectationAndVariance(repartition) ;
  std::cout << "Expectation : " << values.first << ", variance : " << values.second << std::endl ;
  return 0 ;
}

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



std::vector<std::pair<double,int> > computeRepartition(std::vector<double>& thick ) {
  std::cout << "Compute repartition ... " ;
  std::vector<std::pair<double,int> > repartition ;
  int i = 0 ;
  int size = thick.size() ;
  double readValue = thick[0] ;
  int frequence = 0 ;
  for( i = 0 ; i < size; i++ ) {
      if (thick[i] == readValue) {
          frequence++ ;
      } 
      else {
          std::pair<double,int> paire(readValue,frequence) ;
          repartition.push_back(paire) ;  
          frequence = 1 ;
          readValue = thick[i] ;
      }
  } 
  // Finish : write last value
  std::pair<double,int> paire(readValue,frequence) ;
  repartition.push_back(paire) ;
  
  std::cout << "Ok" << std::endl ;
  return repartition ;
}

int writeInFile(std::vector<std::pair<double,int> >& repartition, std::string dataFilename) {
  std::cout << "Write result in file ..." ;
  std::ofstream file(dataFilename.c_str(), std::ios::out);
  if (!file){
      std::cerr << "open file error" << std::endl ;
      return -1 ; 
  } 
  unsigned int i ;
  for ( i = 0 ; i < repartition.size() ; i++ ) {
      file << repartition[i].first << " " << repartition[i].second << std::endl ;
  }
  std::cout << "Ok" << std::endl ;
  return 0 ;
}

int writeScript(std::string dataFilename) {
    std::cout << "Write gnuplot script ... " ;
    std::string scriptName = dataFilename.substr(0,dataFilename.size()-4)+"-script.p" ;
    std::ofstream scriptFile(scriptName.c_str(),std::ios::out) ;
    if ( !scriptFile ) {
        std::cerr << "Open script file error" << std::endl ;
        return -1 ;
    }  
    scriptFile << "reset " << std::endl ;
    scriptFile << std::endl ;
    scriptFile << "set term pdfcairo" << std::endl ;
    scriptFile << "set output \"repartition.pdf\"" << std::endl ; // PDF file name
    scriptFile << "set title \" Thickness repartition \" " << std::endl ;
    scriptFile << "set xlabel \" Radius  \" " << std::endl ;
    scriptFile << "set ylabel \" Frequency \" " << std::endl ;
    scriptFile << std::endl ;
    scriptFile << "set style data line " << std::endl ;
    scriptFile << "set xrange [0:1.01]" << std::endl ;
    scriptFile << "set pointsize 2 " << std::endl ; 
    scriptFile << "set key bottom left " << std::endl ;
    scriptFile << "set boxwidth 0.01" << std::endl ; // histogramm
    scriptFile << "set style fill solid border 3" << std::endl ;
    scriptFile << "plot \"" << dataFilename << "\" with boxes title \"repartition\"" ;
    std::cout << "Ok" << std::endl ;
    return 0 ;
}

std::pair<double,double> expectationAndVariance ( std::vector<std::pair<double,int> >& repartition ) {
    double expectation = 0 ;
    double variance = 0 ;    // we use Var(X) = E[X^2] - E[X]^2 
    double sample = 0 ;
    unsigned int i ;
    for ( i = 0 ; i < repartition.size() ; i++ ) {
        expectation += repartition[i].first * repartition[i].second ;
        variance += repartition[i].first * repartition[i].first * repartition[i].second ;
        sample += repartition[i].second ;
    } 
    expectation /= sample ;
    variance /= sample ;
    variance -= expectation * expectation ;
    std::pair<double,double> values(expectation,variance) ;
    return values ;
}


