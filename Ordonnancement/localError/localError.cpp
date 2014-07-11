#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>

/*-------------------------------------------*\
|      Plot the distribution of the local     |
|  error on meshes with different magnitudes  |
|                  of noise.                  | 
\*-------------------------------------------*/

// PROTOTYPES //

double stringToDouble ( const std::string &s );
std::vector<double> readSdfValues( char * filename ) ;
float error( float a, float b );
float squaredError ( float a, float b) ;


int main ( int argc, char * argv[] ) {
    // argv[1] = thickness values ( .TXT file ) 
    // argv[2] = thickness values ( noised ) ( .TXT file ) 
    if ( argc != 3 ) {
        std::cerr << "Not the good number of arguments" << std::endl ;
        return -1 ;
    }
    std::vector<double> sdfValues  = readSdfValues(argv[1]) ;
    std::vector<double> sdfNoisedValues  = readSdfValues(argv[2]) ; 
    int size = sdfValues.size() ; 
    std::vector<float> localError(size,0.0) ;
    // Important note : The facet numbering is the same for the different OFF files
    std::string noiseName = argv[2] ;
    std::string errorResults = noiseName.substr(0,noiseName.size()-4)+"-local-error.txt";
    std::ofstream errorFile(errorResults.c_str(), std::ios::out);
    if (!errorFile) { 
        std::cerr << " Error in creating file " <<  std::endl;
        return -1 ;
    } 
    // Calculations
    int j ;
    for ( j = 0 ; j < size ; j++ ) {
         localError[j] = error(sdfValues[j],sdfNoisedValues[j]);
    } ;
    // Sort the values 
    std::sort(localError.begin(),localError.end());
    // Put in file
    for ( j = 0 ; j < size ; j++) {
        errorFile << (j+1)*100.0/( (float) size) << " " << localError[j] << std::endl ;
    };
    return 0;
}

/********************\
| Auxilary functions |
\********************/

double stringToDouble ( const std::string &s ) {
    std::stringstream convert(s) ;
    double number ;
    convert >> number ;
    return number ;
}

std::vector<double> readSdfValues( char * filename ) {
  std::vector<double> values ;
  std::ifstream input(filename) ;
  if (!input) {
      std::cerr << "Error when opening input file" << std::endl ;
      return values ; // empty 
  }
  std::string word ;
  int isData = 0 ; // Read only the second column
  while(input >> word) {
      if (isData) {
          values.push_back(stringToDouble(word)) ;
      }
      isData = 1 - isData ;
  } 
  return values ;
}
  

float error( float a, float b ) {
  if (a == 0) { // Pathological case, never happens...
      return b ;
  } 
  else {
      return std::abs(a-b)/a ;
  }
}

float squaredError ( float a, float b) {
  return (a-b)*(a-b) ;
}

std::string convertDouble(double number) {
   std::ostringstream ss;
   ss << number;
   return ss.str();
}



