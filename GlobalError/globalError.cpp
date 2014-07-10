#include <iostream>
#include <fstream>
#include <cstdlib> 
#include <sstream>

double stringToDouble ( const std::string &s ) {
    std::stringstream convert(s) ;
    double number ;
    convert >> number ;
    return number ;
}

int main ( int argc, char * argv[] ) {
  std::ifstream fichier(argv[1],std::ios::out) ;
  if ( !fichier ) {
      std::cerr << "Open file error" << std::endl ;
      return -1 ;
  } 
  std::string word;
  double globalError = 0 ;
  double size = 0 ;
  unsigned int realValue = 0 ;
  while (fichier >> word) {
      if (realValue == 0) { // skip
          realValue++ ;   
      }
      else if (realValue == 1) { // add
           globalError += stringToDouble(word) ;
           size++ ;
           realValue-- ;
      }
      else {
          std::cerr << "Error : not a good value of realValue " << std::endl ;
      }
  }
  globalError /= size ;
  std::cout << " Global error : " << globalError << std::endl ;
  return 0 ;
}
