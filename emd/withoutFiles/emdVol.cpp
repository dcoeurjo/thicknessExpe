#include <iostream>
#include <fstream>
#include <list>
#include <vector>
//#include <math.h>
#include <cmath>
#include <cstdlib> 
#include <algorithm>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Sphere_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/mesh_segmentation.h>
#include <CGAL/property_map.h>

#include "spheresTools.h"

// Other prototypes

std::vector<double> computeRepartition(std::vector<double>& thick, int precision, double max );
double computeEMD(std::vector<double> d1, std::vector<double> d2);
void createScript( std::string& offName, std::string& fileResults );
void writeInFile(std::vector<double>& repartition, std::string dataFilename, int precision, int size, double max);


typedef CGAL::Exact_predicates_inexact_constructions_kernel    Kernel;
typedef CGAL::Polyhedron_3<Kernel>                             Polyhedron;
typedef Kernel::Vector_3                                       Vector;
typedef Kernel::Point_3                                        Point;
typedef Polyhedron::Vertex                                     Vertex;
typedef Kernel::Sphere_3                                       Sphere;
typedef std::vector<Sphere>                                    Spheres;
typedef Spheres::iterator                                      Iterator; 
typedef CGAL::Box_intersection_d::Box_with_handle_d<double,3,Iterator>  Box ;

typedef Polyhedron::Facet_iterator                             Facet_iterator;
typedef Polyhedron::Vertex_iterator                            Vertex_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator           HF_circulator;
typedef Polyhedron::Halfedge_around_vertex_const_circulator    HV_circulator;
 


/****************************************************\
|  Volumetric method robustness against noise (EMD)  |
\****************************************************/

int main ( int argc, char * argv[] ) {
    // argv[1] = filename ( type .off )
    // argv[2] = magnitudes of noise
    // argv[3] = type of noise (1,2,3)
    // argv[4] = resolution
     
    std::string offName = argv[1] ; 
    // get the magnitudes
    std::vector<double> magnitudes = getMagnitudes(argv[2]) ;
    // create EMD file and gnuplot script
    std::string fileResults =  offName.substr(0,offName.size()-4)+"-emd.txt" ;
    std::ofstream emdFile(fileResults.c_str(), std::ios::out);
    if (!emdFile) { 
        std::cout << " Error in creating EMD file " << fileResults << std::endl;
        return -1 ;
    }
    createScript(offName,fileResults) ;

    int resolution = atoi(argv[4]) ;
    unsigned int i ;
    for ( i = 0 ; i < magnitudes.size() ; i++ ) {
        // build the off noised mesh
        std::string noisedOffName = offName.substr(0,offName.size()-4)+"-s="+doubleToString(magnitudes[i])+"%.off" ;
        int s = createNoiseOff(argv[1], noisedOffName.c_str(),magnitudes[i],atoi(argv[3])) ;
        if ( s != 0 ) {
            std::cerr << "Error when creating the noised mesh " << std::endl ;
            return -1 ;
        }
        // compute thickness 
        Calcul calc(resolution) ; 
        std::vector<double> thickWithoutNoise = offToThick(offName,resolution,calc) ;
        initMaxRadius(resolution) ; // maxRadius = 0 for all points
        std::vector<double> thickWithNoise = offToThick(noisedOffName.c_str(),resolution,calc) ;
        freeArrays() ;
        std::string thickFilename = noisedOffName.substr(0,noisedOffName.size()-4)+"-thick.txt";
        std::string repFilename = offName.substr(0,noisedOffName.size()-4)+"-rep.txt";
        std::string repNoisedFilename = noisedOffName.substr(0,noisedOffName.size()-4)+"-rep.txt";
        std::ofstream thickFile(thickFilename.c_str(),std::ios::out);
        std::ofstream repFile(repFilename.c_str(),std::ios::out);
        std::ofstream repNoisedFile(repNoisedFilename.c_str(),std::ios::out);
        if( !thickFile || !repFile || !repNoisedFile ){ 
            std::cout << " Error in creating file " << fileResults << std::endl;
            return -1 ;
        } 

        int j = 0 ;
        int size = resolution * resolution * resolution ;
        std::vector<double> realValues ;
        std::vector<double> realValuesNoised ;
        for(j=0;j<size;j++) {
             // We only consider the points which belong to the two meshes
             if ((thickWithoutNoise[j] != 0) && (thickWithNoise[j] != 0)) {
                 realValues.push_back(thickWithoutNoise[j]) ;
                 realValuesNoised.push_back(thickWithNoise[j]) ;
             }
        }; 
        std::cout << "Sort ... " ;
        // Sort the values 
        std::sort(realValues.begin(),realValues.end());
        std::sort(realValuesNoised.begin(),realValuesNoised.end());
        // Normalize
        int taille = realValues.size() ;  
        for(j=0;j<taille;j++) {
            realValues[j] /= resolution ;
            realValuesNoised[j] /= resolution  ;

        };      
        // Compute EMD 
        std::cout << "Ok" << std::endl << "Calculations... " ;

        double max = std::max(realValues[taille-1],realValuesNoised[taille-1]) ;
        std::vector<double> d1 = computeRepartition(realValues,100,max) ;
        std::vector<double> d2 = computeRepartition(realValuesNoised,100,max) ;
        writeInFile(d1,repFilename,100,taille,max);
        writeInFile(d2,repNoisedFilename,100,taille,max);
        double emd = computeEMD(d1,d2) ;
        // Put in file
        std::cout << "Ok" << std::endl << "Write in files... " ;
        for(j=0;j<taille;j++) {
           thickFile << (j+1)*100.0/( (float) taille) << " " << realValuesNoised[j] << std::endl ;
        };
        emdFile << magnitudes[i] << " " << emd << std::endl ; 
        std::cout << "Ok" << std::endl;
    };
    return 0 ;
}



/********************\
| AUXILIAR FUNCTIONS |
\********************/


/************\
| Converters |
\************/

std::string doubleToString(double number) {
   std::ostringstream ss;
   ss << number;
   return ss.str();
}


std::string intToString(int number) {
   std::ostringstream ss;
   ss << number;
   return ss.str();
}



int stringToInt ( const std::string &s ) {
    std::stringstream convert(s) ;
    int number ;
    convert >> number ;
    return number ;
}

int stringToDouble ( const std::string &s ) {
    std::stringstream convert(s) ;
    double number ;
    convert >> number ;
    return number ;
}



/**************\
| Sphere tools |
\**************/


std::vector<Sphere> readSpheres ( const char * filename ) {
    std::vector<Sphere> spheres ;
    std::ifstream input(filename) ;
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



/***********************\
| Tools to create noise |
\***********************/

std::vector<double> findBoxBorder ( Polyhedron P ) {
  Vertex_iterator v = P.vertices_begin() ;
  std::vector<double> borders(6) ;
    borders[0] = v->point().x() ;
    borders[1] = v->point().x() ;
    borders[2] = v->point().y() ;
    borders[3] = v->point().y() ;
    borders[4] = v->point().z() ;
    borders[5] = v->point().z() ;
  do {
      double x = v->point().x() ;
      double y = v->point().y() ;
      double z = v->point().z() ; 
      if ( x < borders[0] ) {
          borders[0] = x ;
      }
      else if ( x > borders[1] ) {
          borders[1] = x ;
      }
      if ( y < borders[2] ) {
          borders[2] = y ;
      }
      else if ( y > borders[3] ) {
          borders[3] = y ;
      }
      if ( z < borders[4] ) {
          borders[4] = z ;
      }
      else if ( z > borders[5] ) {
          borders[5] = z ;
      }
  } while ( v++ != P.vertices_end() ) ;
   return borders ;
}


double longestDiagonalLength ( Polyhedron P ) {
  std::vector<double> borders = findBoxBorder(P) ;
  double deltaX = borders[1]-borders[0] ;
  double deltaY = borders[3]-borders[2] ;
  double deltaZ = borders[5]-borders[4] ;
  return sqrt(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ) ;
}

int createNoiseOff ( char * fileName, const char * outputName, double magnitude, int option) {    
    // Create and read Polyhedron
    Polyhedron mesh ;
    std::ifstream input(fileName) ;
    if ( !input || !(input >> mesh) || mesh.empty() ) {
        std::cerr << "Not a valid .off file" << std::endl; 
    };
    // Initialize generators
    std::default_random_engine generator;
    double s = magnitude/100*longestDiagonalLength(mesh);
    std::normal_distribution<double> distribution(0.0,s); 
    // apply noise forall vertices
    switch ( option ) {
        case 1 : { // "cubic noise"
            Vertex_iterator v = mesh.vertices_begin() ;
            do {
                double x = distribution(generator);
                double y = distribution(generator);
                double z = distribution(generator);
                Vector vec(x,y,z);
                v->point() = v->point() + vec ;
                } while ( v++ != mesh.vertices_end()) ;
        break ;
        }
        case 2 : { // "spheric noise"
            Vertex_iterator v = mesh.vertices_begin() ;
            do {
                bool inCircle = false ;
                double x, y, z ;
                while(!inCircle) {
                     x = distribution(generator);
                     y = distribution(generator);
                     z = distribution(generator);
                     if ( x*x + y*y + z*z <= s*s) {
                         inCircle = true ;
                     }
                }
                Vector vec(x,y,z);
                v->point() = v->point() + vec ;
                } while ( v++ != mesh.vertices_end()) ;
        break ;
        }
        case 3 : { // "normal noise"
            Vertex_iterator v = mesh.vertices_begin() ;
            do { 
                HV_circulator h = v->vertex_begin() ;
                Vector vec(0.0,0.0,0.0) ;
                do {
                    vec = vec + CGAL::cross_product(  
                                h->next()->vertex()->point() - h->vertex()->point(),
                                h->next()->next()->vertex()->point() - h->next()->vertex()->point()) ;
                    ++h ;
                    CGAL_assertion( h != v.vertex_begin()); // even degree guaranteed 
                    ++h ;
                   } while ( h != v->vertex_begin());
               double norme = sqrt(vec.x()*vec.x() + vec.y()*vec.y() + vec.z()*vec.z()) ;
               double factor = distribution(generator) ; 
               Vector vectorNoise( vec.x()*factor/norme, vec.y()*factor/norme, vec.z()*factor/norme) ;
               v->point() = v->point() + vectorNoise ;
               v++;
               } while ( v != mesh.vertices_end()) ;
        break ;
        }
        default : {
          std::cerr << " Erreur de saisie " << std::endl ;
        }
    }; 
    // Write polyhedron in .off format.
    std::ofstream offOutput(outputName,std::ios::out) ;
    CGAL::set_ascii_mode( offOutput ) ;
    offOutput << "OFF" << std::endl ;
    offOutput << mesh.size_of_vertices() << ' ' << mesh.size_of_facets() << ' ' << mesh.size_of_halfedges() << std::endl ; 
    std::copy( mesh.points_begin(), mesh.points_end(), std::ostream_iterator<Point>( offOutput, "\n"));
    for ( Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i) {
         HF_circulator j = i->facet_begin();
	 // facets in polyhedral surfaces are at least triangles.
	 CGAL_assertion( CGAL::circulator_size(j) >= 3);
	 offOutput << CGAL::circulator_size(j) << ' ' << std::distance(mesh.vertices_begin(), 
         i->facet_begin()->vertex()) ;
         while ( ++j != i->facet_begin()) {
	     offOutput << ' ' << std::distance(mesh.vertices_begin(), j->vertex());
         } 
	 offOutput << std::endl;
    }
    return 0;
}

/************************\
| Thickness calculations |
\************************/

std::vector<double> getMagnitudes ( char * magnitudeFile ) {
  std::vector<double> magnitudes ; 
  std::ifstream file(magnitudeFile,std::ios::in);
  if (!file) {
      std::cerr << "Error: open resolution file" << std::endl ;
  }
  else {
      std::string line ;
      while (getline(file,line)) {
          double magnitude = atof(line.c_str()) ;
          magnitudes.push_back(magnitude) ;
      }
  }
  file.close() ;
  sort(magnitudes.begin(),magnitudes.end()) ;
  return magnitudes ;
}

std::vector<double> offToThick(std::string filename, int resolution, Calcul& calc ) {
    std::string offName = filename ;
    std::string systemCall = "bash off2spheres.sh "+offName+" "+intToString(resolution) ;
    int s = system(systemCall.c_str()) ;
    if ( s != 0 ) {
        std::cerr << "Bash script error" << std::endl ;
    }
    std::string spheresFile = filename.substr(0,filename.size()-4)+"-"+intToString(resolution)+"-spheres.txt" ;
    std::vector<double> thick = spheresToThick(spheresFile.c_str(),resolution,0,calc) ;
    return thick ;
}


std::vector<double> spheresToThick( const char * filename , int resolution , int option , Calcul& calc ) {
    std::vector<Sphere> spheres = readSpheres(filename) ;
    std::vector<Box> boxes, pointBoxes ;
    std::cout << "Building boxes ... " ;
    for ( Iterator i = spheres.begin(); i != spheres.end(); ++i ) {
        boxes.push_back( Box(i->bbox(),i)) ;
    }
    for ( Iterator j = points.begin(); j != points.end(); ++j) {
      pointBoxes.push_back( Box(j->bbox(),j)) ;
    }
    std::cout << "Ok" << std::endl; 
    std::cout << "Compute intersections ... " ;
    CGAL::box_intersection_d( boxes.begin(), boxes.end(), pointBoxes.begin(), pointBoxes.end(), calc);
    std::cout << "Ok" << std::endl ;
    return maxRadius ;
}

/******************\
| EMD Calculations |
\******************/


std::vector<double> computeRepartition(std::vector<double>& thick, int precision, double max ) {
  std::cout << "Compute repartition ... " ;
  std::vector<double> repartition(precision,0.0) ;
  int i = 0 ;
  int size = thick.size() ;
  for( i = 0 ; i < size; i++ ) {
       int indice = (int) ( ceil( thick[i] / max * ((double) precision))) - 1 ;
       repartition[std::max(indice,0)]++ ;
  } 
  // normalize
  for( i = 0 ; i < precision ; i++ ) {
       repartition[i] /= size ;
       repartition[i] *= 100 ; 
  } 
  std::cout << "Ok" << std::endl ; 
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
      std::cout << std::abs(EMDi[i]) << " " ;
      emd += std::abs(EMDi[i]) ;
  }
  std::cout << std::endl ;
  return emd ;
}


/******************\
| Arrays functions |
\******************/


void initMaxRadius(int resolution) {
    int i ;
    for ( i = 0 ; i < resolution * resolution * resolution ; i++) {
        maxRadius[i] = 0 ;
    }
}

void freeArrays( void ) {
    maxRadius.resize(0) ; 
    points.resize(0) ;
}

/***************\
| Create script |
\***************/

void createScript( std::string& offName, std::string& fileResults ) {

    std::string scriptName = offName.substr(0,offName.size()-4)+"-emd.p" ;
    std::ofstream scriptFile(scriptName.c_str(),std::ios::out) ;
    if ( !scriptFile ) {
        std::cerr << "Open script file error" << std::endl ;
    }  
    scriptFile << "reset " << std::endl ;
    scriptFile << std::endl ;
    scriptFile << "set term pdfcairo" << std::endl ;
    scriptFile << "set output \"" << offName << "-emd.pdf\"" << std::endl ;
    scriptFile << "set title \" Earth Mover Distance \" " << std::endl ;
    scriptFile << "set xlabel \" Noise (%) \" " << std::endl ;
    scriptFile << "set ylabel \" EMD \" " << std::endl ;
    scriptFile << std::endl ;
    scriptFile << "set style data line " << std::endl ;
    scriptFile << "set pointsize 2 " << std::endl ; 
    scriptFile << "set key bottom right " << std::endl ;
    scriptFile << std::endl ;
    scriptFile << "plot \""+fileResults+"\" using 1:2 title \"EMD\"" ;
    std::cout << "Script ok !" << std::endl ;
}

void writeInFile(std::vector<double>& repartition, std::string dataFilename, int precision, int size, double max) {
  std::cout << "Write result in file ..." ;
  std::ofstream file(dataFilename.c_str(), std::ios::out);
  if (!file){
      std::cerr << "open file error" << std::endl ;
  } 
  int i ;
  for ( i = 0 ; i < precision ; i++ ) {
        double mediumValue = (((double) i)+0.5)/((double) precision)* max    ;
        file << mediumValue << " " << repartition[i] << std::endl ;
  }
  std::cout << "Ok" << std::endl ;
}





