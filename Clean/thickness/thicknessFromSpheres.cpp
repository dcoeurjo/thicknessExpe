#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Bbox_3.h>
#include <iostream>
#include <fstream>
#include <algorithm>

#include <CGAL/Min_sphere_of_spheres_d.h>

#include <CGAL/Cartesian_d.h>


#include <CGAL/bounding_box.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel            Kernel;


typedef Kernel::Point_3                                                Point;
typedef Kernel::Sphere_3                                               Sphere;


typedef std::vector<Sphere>                                            Spheres;
typedef Spheres::iterator                                              Iterator;   
typedef CGAL::Box_intersection_d::Box_with_handle_d<double,3,Iterator> Box ;

typedef Kernel::Iso_cuboid_3   BBox;


//Set of probing points
Spheres points ;

//Points index -> maxRadius (0 if outside)
std::vector<double> maxRadius ;
std::vector<double> maxRadiusClean ;

// PROTOTYPES //

int whatOption(char * opt) ;
std::string intToString(int number);
double stringToDouble ( const std::string &s );
std::vector<Sphere> readSpheres ( std::string& fileName );
bool belongsTo ( const Sphere& s1, const Sphere& s2 );
int nonZeroValues( std::vector<double>& values ) ;	
int writeThicknessInFile ( std::vector<double>& thick , std::string& fileName, int beginning, BBox bbox );
//double rayMaxBall ( const Point P, std::vector<Sphere>& spheres );




BBox findBoundingBox (  Spheres & S )
{
  std::vector<Point> pts;
  for(Iterator it = S.begin(), itend=S.end(); it!= itend; ++it)
  {
    Point center = it->center();
    double rad = sqrt(it->squared_radius());
    
    //std::cout<< "rad= "<< rad<<std::endl;
    
    pts.push_back(Point(center.x() , center.y(), center.z() + rad));
    pts.push_back(Point(center.x() , center.y(), center.z() - rad));
    pts.push_back(Point(center.x() +rad, center.y(), center.z() ));
    pts.push_back(Point(center.x() -rad, center.y(), center.z() ));
    pts.push_back(Point(center.x() , center.y()+rad, center.z() ));
    pts.push_back(Point(center.x() , center.y()-rad, center.z() ));
  }
  
  return CGAL::bounding_box(pts.begin(), pts.end());
  
}

double widthFromBBox(BBox &bbox)
{
  double values[] = {(bbox.xmax() - bbox.xmin()), (bbox.ymax() - bbox.ymin()),(bbox.zmax() - bbox.zmin())};
  return *std::max_element(values, values+3) / 2.0;
}


// A class for the "CGAL::box_intersection_d" call

class Calcul {
public:

  // init points and maxRadius

  Calcul(BBox bbox,int option,int thirdArgument) {
    switch(option) {
    case 0 :  { // regular
      // thirdArgument = width of the discretisation
      srand(time(0)) ;
      double i,j,k ;
      for ( i = bbox.xmin(); i < bbox.xmax();  )
      {
        
        for ( j = bbox.ymin(); j < bbox.ymax(); )
        {
          
          for ( k =  bbox.zmin(); k < bbox.zmax();  )
          {
            Point P(i,j,k) ;
            points.push_back(Sphere(P,0)) ;
            maxRadius.push_back(0.0) ;
            
            i += (bbox.xmax() - bbox.xmin()) / (double) thirdArgument;
          }
          j += (bbox.ymax() - bbox.ymin()) / (double) thirdArgument;
        }
        k += (bbox.zmax() - bbox.zmin()) / (double) thirdArgument;
        } 
      break ;
    }     
    case 1 : { // by Monte-Carlo
      // thirdArgument = number of points
      srand(time(0)) ;
      int nb = 0 ;
      while ( nb < thirdArgument )
      {
        double x = bbox.xmin() + ((double) rand()) /((double) RAND_MAX ) * (bbox.xmax() - bbox.xmin()) ;
        double y = bbox.ymin() + ((double) rand()) /((double) RAND_MAX ) * (bbox.ymax() - bbox.ymin()) ;
        double z = bbox.zmin() + ((double) rand()) /((double) RAND_MAX ) * (bbox.zmax() - bbox.zmin()) ;
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

  void operator()(const Box& a, const Box& b)
  {
    // b is a point
    if ( belongsTo( *(a.handle()), *(b.handle()))
         && ((a.handle())->squared_radius() > maxRadius[b.handle()-points.begin()]))
    {
      maxRadius[b.handle()-points.begin()] = (a.handle())->squared_radius() ;
     }

  }
};


/*
 * Read a set of spheres ("x y z radius") and estimate the thickness distribution using sampling (regular or uniform)
 *
 *  The output is a file (<input>-thick.txt) with
 *
 *   rank  thickness
 *
 * with rank [0:100] and thickness in [0:1]  (thickness normalized by the half-width of the bbox)
 *
 */

int main( int argc, char * argv[] ) {
  /* 
     argv[1] = .dat file containing the balls
     argv[2] = option for the choosen points ( only voxels, by discretisation or randomly ( Monte-Carlo )
     argv[3] = number of points (Monte-Carlo) or the step ( Discretisation )
  */
  std::string filename = argv[1] ;
  // Extract option and third argument if it exists
  int option = whatOption(argv[2]) ;
  int thirdArgument = 0 ;
  if ( argc == 4 ) {
    thirdArgument = atoi(argv[3]) ;
  }
  if ( argc != 4 && argc != 3 ) {
    std::cerr << " Not the good number of arguments " << std::endl ;
    return -1 ;
  }
  
  
  
  //Loading the ball set
  Spheres spheres = readSpheres(filename);
  
  std::cout << "Number of spheres= "<< spheres.size()  <<std::endl;
  std::cout << "Parameter (resolution or number of points)= "<< thirdArgument  <<std::endl;
  std::cout << "Samples= "<< 3*thirdArgument  <<std::endl;
  
  //Get bboxwith
  BBox bbox = findBoundingBox( spheres );
  std::cout << "BBox= "<<bbox  <<std::endl;
  
  // Init points and radius
  Calcul calc(bbox,option,3*thirdArgument) ;

  std::vector<Box> boxes, pointBoxes ;
  for ( Iterator i = spheres.begin(); i != spheres.end(); ++i ) {
    boxes.push_back( Box(i->bbox(),i)) ;
  }
  
  for ( Iterator j = points.begin(); j != points.end(); ++j) {
    pointBoxes.push_back( Box(j->bbox(),j)) ;
  }
  
  CGAL::box_intersection_d( boxes.begin(), boxes.end(), pointBoxes.begin(), pointBoxes.end(), calc);
  
  //we keep the same number of points with thick!=0 in MC
  int cpt=0;
  while ((cpt < maxRadius.size()) && (maxRadiusClean.size() < thirdArgument))
  {
    if (maxRadius[cpt] != 0)
    {
      maxRadiusClean.push_back(maxRadius[cpt]);
    }
    cpt++;
  }
  
  if (maxRadiusClean.size() < thirdArgument)
    std::cout<< "WARNING !!!!!! not enough point in the shape"<<std::endl;
  
  sort(maxRadiusClean.begin(),maxRadiusClean.end()) ;
  
  int firstNonZeroValue = nonZeroValues(maxRadiusClean) ;
  std::string dataFileName = filename.substr(0,filename.size()-4)+"-thick.txt" ;
  int w = writeThicknessInFile(maxRadiusClean,dataFileName,firstNonZeroValue,bbox) ;
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
  // opt = -r or -mc
  int option ;
  switch ( opt[1] ) {
  case 'r' :
    option = 0 ;
    break ;
  case 'm' :
    option = 1 ;
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
// "x y z radius" 
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
   


int writeThicknessInFile ( std::vector<double>& thick , std::string& fileName, int firstNonZeroValue, BBox bbox )
{
  std::ofstream thickFile(fileName.c_str(),std::ios::out) ;
  if ( !thickFile ) {
    return -1;
  }
  int i ;
  int size = thick.size()-firstNonZeroValue ;
  std::cout << "Number of points : " << size << std::endl ;
  double max = sqrt(thick[thick.size()-1]) ;
  double min = sqrt(thick[firstNonZeroValue]) ;
  std::cout << "Thickness min value : " << min << ", thickness max value : " << max << std::endl ;
  
  for ( i = firstNonZeroValue ; i < thick.size() ; i++ )
  {
    thickFile << sqrt(thick[i])   << std::endl ;
  }
  return 0 ;
}

