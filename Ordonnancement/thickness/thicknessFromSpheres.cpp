#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Bbox_3.h>
#include <iostream>
#include <fstream>
#include <algorithm>


typedef CGAL::Exact_predicates_inexact_constructions_kernel            Kernel;
typedef Kernel::Point_3                                                Point;
typedef Kernel::Sphere_3                                               Sphere;
typedef std::vector<Sphere>                                            Spheres;
typedef Spheres::iterator                                              Iterator;   
typedef CGAL::Box_intersection_d::Box_with_handle_d<double,3,Iterator> Box ;

Spheres points ;
std::vector<double> maxRadius ;

// PROTOTYPES //

int whatOption(char * opt) ;
std::string intToString(int number);
double stringToDouble ( const std::string &s );
std::vector<Sphere> readSpheres ( std::string& fileName );
bool belongsTo ( const Sphere& s1, const Sphere& s2 );
int nonZeroValues( std::vector<double>& values ) ;	
int writeThicknessInFile ( std::vector<double>& thick , std::string& fileName, int beginning, int resolution );
//double rayMaxBall ( const Point P, std::vector<Sphere>& spheres );


// A class for the "CGAL::box_intersection_d" call

class Calcul {
public:

  // init points and maxRadius

  Calcul(int resolution,int option,int thirdArgument) {
    switch(option) {
    case 0 : { // all voxels
      double i, j, k ;
      for ( k = 0 ; k < resolution ; k++ ) { 
        for ( j = 0 ; j < resolution ; j++ ) {
          for ( i = 0 ; i < resolution ; i++) {
            points.push_back(Sphere(Point(i,j,k),0)) ;
            maxRadius.push_back(0.0) ;
          }
        }
      }
      break ;
    }
    case 1 : { // by discretisation
      // thirdArgument = step of the discretisation
      srand(time(0)) ;
      int i,j,k ;
      for ( i = 0; i < thirdArgument; i++ ) {
        for ( j = 0; j < thirdArgument; j++ ) {
          for ( k = 0; k < thirdArgument; k++ ) {
            double x = ((double) i)/((double) thirdArgument) * resolution ;
            double y = ((double) j)/((double) thirdArgument) * resolution ;
            double z = ((double) k)/((double) thirdArgument) * resolution ;
            Point P(x,y,z) ;
            points.push_back(Sphere(Point(x,y,z),0)) ;
            maxRadius.push_back(0.0) ;
          }
        }
      }
      break ;
    }     
    case 2 : { // by Monte-Carlo
      // thirdArgument = number of points
      srand(time(0)) ;
      int nb = 0 ;
      while ( nb < thirdArgument ) {  
        double x = ((double) rand()) /((double) RAND_MAX ) * resolution ;
        double y = ((double) rand()) /((double) RAND_MAX ) * resolution ;
        double z = ((double) rand()) /((double) RAND_MAX ) * resolution ;
        Point P(x,y,z) ;
        points.push_back(Sphere(Point(x,y,z),0)) ;
        maxRadius.push_back(0.0) ;
        nb++ ;
      }
      break ;
    }
    default : {
      std::cerr << "Not a good option ! " << std::endl ;
    }
    }
  }

  // update the radius values  

  void operator()(const Box& a, const Box& b) {
    // b is a point
    if ( belongsTo( *(a.handle()), *(b.handle())) 
         && ((a.handle())->squared_radius() > maxRadius[b.handle()-points.begin()])) {
      maxRadius[b.handle()-points.begin()] = (a.handle())->squared_radius() ;
    }
  }

};


int main( int argc, char * argv[] ) {
  /* 
     argv[1] = .dat file containing the balls
     argv[2] = resolution for voxelisation
     argv[3] = option for the choosen points ( only voxels, by discretisation or randomly ( Monte-Carlo )
     argv[4] = number of points (Monte-Carlo) or the step ( Discretisation )
  */
  std::string filename = argv[1] ;
  int resolution = atoi(argv[2]) ;
  // Extract option and third argument if it exists
  int option = whatOption(argv[3]) ;
  int thirdArgument = 0 ;
  if ( argc == 5 ) {
    thirdArgument = atoi(argv[4]) ;
  }
  if ( argc != 5 && argc != 4 ) {
    std::cerr << " Not the good number of arguments " << std::endl ;
    return -1 ;
  }
  
  //Loading the ball set
  Spheres spheres = readSpheres(filename);
  // Init points and radius
  Calcul calc(resolution,option,thirdArgument) ; 
  std::vector<Box> boxes, pointBoxes ;
  for ( Iterator i = spheres.begin(); i != spheres.end(); ++i ) {
    boxes.push_back( Box(i->bbox(),i)) ;
  }
  for ( Iterator j = points.begin(); j != points.end(); ++j) {
    pointBoxes.push_back( Box(j->bbox(),j)) ;
  }
  CGAL::box_intersection_d( boxes.begin(), boxes.end(), pointBoxes.begin(), pointBoxes.end(), calc);
  sort(maxRadius.begin(),maxRadius.end()) ;
  int beginning = nonZeroValues(maxRadius) ;
  std::string dataFileName = filename.substr(0,filename.size()-4)+"-thick.txt" ;
  int w = writeThicknessInFile(maxRadius,dataFileName,beginning,resolution) ;
  if ( w != 0 ) {
    std::cout << " Open thickness file error" << std::endl ;
    return -1 ;
  }
  return 0; 
}

/************************\ 
| Fonctions auxilliaires |
\************************/

int whatOption(char * opt) {
  // opt = -v or -c or -mc 
  int option ;
  switch ( opt[1] ) {
  case 'v' :
    option = 0 ;
    break ;
  case 'c' : 
    option = 1 ;
    break ;
  case 'm' :
    option = 2 ;
    break ;
  default :
    option = -1 ;
  }
  return option ;
}
 

std::string intToString(int number) {
  std::ostringstream ss;
  ss << number;
  return ss.str();
}

double stringToDouble ( const std::string &s ) {
  std::stringstream convert(s) ;
  double number ;
  convert >> number ;
  return number ;
}

// Read Spheres from a .TXT file

Spheres readSpheres ( std::string& filename ) {
  Spheres spheres ;
  std::ifstream input(filename.c_str()) ;
  if (!input) {
    std::cerr << " Open sphere file error " << std::endl ;
  }
  std::string number ;
  while(input >> number) {
    double x = stringToDouble(number) ;
    if ( input >> number ) {
      double y = stringToDouble(number) ;
      if ( input >> number ) {        
        double z = stringToDouble(number) ;
        if ( input >> number ) {
          double r = stringToDouble(number) ;
          Point P(x,y,z) ;
          Sphere S(P,r*r) ;
          spheres.push_back(S) ;
        }
        else std::cerr << " Incorrect sphere file : cannot read R " << std::endl ;
      }
      else std::cerr << " Incorrect sphere file : cannot read z " << std::endl ;            
    }
    else std::cerr << " Incorrect sphere file : cannot read y " << std::endl ;
  }
  return spheres ;
}

// Test if the point represented by s2 belongs to s1.

bool belongsTo ( const Sphere& s1, const Sphere& s2 ) {
  Point p = s1.center() ;
  Point q = s2.center() ;
  double x = p.x() - q.x() ;
  double y = p.y() - q.y() ;
  double z = p.z() - q.z() ;
  return ( x*x+y*y+z*z <= s1.squared_radius() ) ;
}


int nonZeroValues( std::vector<double>& values ) {
  int i = 0 ;
  int size = values.size() ;
  while( i < size && values[i] == 0) {
    i++ ;
  }	
  return i ;
}
   


int writeThicknessInFile ( std::vector<double>& thick , std::string& fileName, int beginning, int resolution ) {
  std::ofstream thickFile(fileName.c_str(),std::ios::out) ;
  if ( !thickFile ) {
    return -1;
  }
  int i ;
  int size = thick.size()-beginning ;
  std::cout << "Number of points : " << size << std::endl ;
  double max = sqrt(thick[thick.size()-1]) ;
  double min = sqrt(thick[beginning]) ;
  std::cout << " Thickness min value : " << min << ", thickness max value : " << max << std::endl ;
  for ( i = beginning ; i < thick.size() ; i++ ) {
    double x = ((double) (i+1-beginning)) / ((double) size) * 100 ;
    thickFile << x << " " << sqrt(thick[i]) / resolution  << std::endl ;
  }
  return 0 ;
}

