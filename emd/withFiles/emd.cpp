#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <utility>
#include <algorithm>

/*
  1) read sdf values
  2) compute EMD 
  3) return value
*/

// Input : the 2 sdf files and the precision of the histogramm

double stringToDouble ( std::string word );
std::vector<double> readSdfValues( char * filename );
std::vector<double> computeRepartition(std::vector<double>& thick, int precision, double max );
double computeEMD(std::vector<double> d1, std::vector<double> d2) ;

int main ( int argc , char * argv[] ) {
  // argv[1] = first sdf values
  // argv[2] = second sdf values
  // argv[3] = precision pour repartition rayons
  int precision = atoi(argv[3]) ;
  std::vector<double> sdf1 = readSdfValues(argv[1]) ;
  std::vector<double> sdf2 = readSdfValues(argv[2]) ;
  double max = std::max(sdf1[sdf1.size()-1],sdf2[sdf2.size()-1]) ;
  std::vector<double> d1 = computeRepartition(sdf1,precision,max) ;
  std::vector<double> d2 = computeRepartition(sdf2,precision,max) ;
  double emd = computeEMD(d1,d2) ;
  std::cout << "EMD = " << emd << std::endl ; 
  return 0 ;
}

double stringToDouble ( std::string word ) {
    return atof(word.c_str()) ;
}

std::vector<double> readSdfValues( char * filename ) {
  std::vector<double> values ;
    std::cout << "Read data ... " ;
  std::ifstream input(filename) ;
  if (!input) {
      std::cerr << "Error when opening input file" << std::endl ;
      return values ; // empty 
  }
  std::string word ;
  int isData = 0 ;
  while(input >> word) {
      if (isData) {
          values.push_back(stringToDouble(word)) ;
      }
      isData = 1 - isData ;
  } 
  std::cout<< "Ok" << std::endl ;
  return values ;
}

std::vector<double> computeRepartition(std::vector<double>& thick, int precision, double max ) {
  std::cout << "Compute repartition ... " ;
  std::vector<double> repartition(precision,0.0) ;
  int i = 0 ;
  int size = thick.size() ;
  for( i = 0 ; i < size; i++ ) {
       int indice = (int) ( ceil( thick[i] / max * ((double) precision))) - 1 ;
       repartition[std::max(indice,0)]++ ;
  } 
  std::cout << "Ok" << std:: endl ; 
  // Normalize
  for ( i = 0 ; i < precision ; i++ ) {
      repartition[i] = repartition[i] / size * 100.0 ;
  }
  return repartition ;
}


double computeEMD(std::vector<double> d1, std::vector<double> d2) {
  std::vector<double> EMDi ;
  int size = d1.size() ;
  if ( size != ((int) d2.size()) ) {
      std::cerr << "Error : the datas don't have the same length " << std::endl ;
      return -1.0 ;
  }
  int i ;
  EMDi.push_back(0.0) ;
  for ( i = 0 ; i < size ; i++ ) {
      EMDi.push_back(d1[i]-d2[i]+EMDi[i]) ;
  }
  double emd = 0.0 ;
  for ( i = 0 ; i < size + 1 ; i++ ) {
      // std::cout << std::abs(EMDi[i]) << " " ;
      emd += std::abs(EMDi[i]) ;
  }
  return emd ;
}

