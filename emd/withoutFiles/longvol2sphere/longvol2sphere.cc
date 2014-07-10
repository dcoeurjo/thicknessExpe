#include <stdio.h>
#include <stdlib.h>
#include <longvol.h>
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <fstream>

int main(int argc,char **argv)
{
  //Parsing command line
  if (argc != 2) {
    std::cerr << "Convert a longvol to Geomview SPHERE file 'x y z val' \nInput: a LONGVOL file\nOutput: the LIST file\n Usage : \n\t\t longvol2vol  <sourcefile.longvol>\n";
    exit(1);
  }

  std::string sourceFile = argv[1] ;
  std::string destFile = sourceFile.substr(0,sourceFile.size()-8)+"-spheres.txt" ;

  Longvol lv(argv[1]);
  if (!lv.isOK()) {
      fprintf( stderr, "Couldn't load \"%s\" file !\n", sourceFile.c_str() );
      return 1;
  };
  
  std::ofstream sphereFile(destFile.c_str(),std::ios::out) ;
  if ( !sphereFile ) {
      fprintf(stderr,"Open spheres file error\n");
      exit(-1);
  };  
  
  for(int k= lv.minZ() ; k < lv.maxZ() ; k++) {
      for(int j= lv.minY(); j< lv.maxY() ; j++) {
          for(int i= lv.minX(); i<lv.maxX() ; i++) {
	      if (lv(i,j,k)!=0) {
                  sphereFile << i << " " << j << " " << k << " " << sqrt(lv(i,j,k)) << std::endl ;
              }
          }
      }
  }

  return 0;
}

