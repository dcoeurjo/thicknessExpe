#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

/*------------------------------------------*\
| Create the data used to plot the thickness |
|   diagram computed by the volumic method.  |
\*------------------------------------------*/


int main ( int argc, char * argv[] ) {
  // argv[1] = name of the file ( type : name-resolution )
  // argv[2] = number of faces ( in sdf computation ) TODO : Unused
  // argv[3] = maximal thickness in argv[1].off
  // argv[4] = number of voxels in argv[1].off
  // These 3 arguments are used in the change of scale
  // to plot all the distributions
  double max = atof(argv[3]) ;
  double factor = 100 / atof(argv[4]) ;
  string input = argv[1] ;
  string fileInput = input+"-data.txt" ; 
  string fileOutput = input+"-newdata.txt" ;
  ifstream file1( fileInput.c_str(), ios::in); 
  ofstream file2( fileOutput.c_str(), ios::out);
  if ( file1 ) {
      if ( file2 ) {
          string line ; 
          double i = 1 ;
          while ( getline(file1,line) ) {
              double value = atof(line.c_str())/ max ;
              double x = i * factor ;
              file2 << x << " " << value << endl ;
              i++ ;
          }
	  file1.close() ;
          file2.close() ;
          return 0 ;
      }
      else
      cerr << "write file error" << endl ;
      return -1 ;
  }
  else
  cerr << "open file error" << endl ;
  return -1 ;
}
