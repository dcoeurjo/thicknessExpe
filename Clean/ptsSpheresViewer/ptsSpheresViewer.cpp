#include <iostream>
#include <fstream>
#include <vector>

#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/io/viewers/Viewer3D.h>

using namespace DGtal;

int main(int argc, char**argv)
{
  //Loading point set  argv[1]==PTS  argv[2]==SPHERES
  ifstream mypts (argv[1], std::ifstream::in);
  ifstream myspheres (argv[2], std::ifstream::in);

  QApplication app(argc,argv);
  Viewer3D<> viewer;
  double x,y,z,rad;
  viewer.show();
  
  int cpt= 0;
  while (myspheres.good())
  {
    myspheres>> x;
    myspheres >>y;
    myspheres>> z;
    myspheres >> rad;
    viewer.addBall(Z3i::RealPoint(x,y,z), 0.008,5);
    ++cpt;
  }
  myspheres.close();
  trace.info()<< "Nb balls= "<<cpt<<std::endl;
  
  
  viewer << CustomColors3D(Color::Red, Color::Red);
  
  while (mypts.good())
  {
    mypts>> x;
    mypts >>y;
    mypts>> z;
    viewer.addBall(Z3i::RealPoint(x,y,z), 0.008,5);
  }
  mypts.close();


  
  
  return app.exec();
}

