#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/mesh_segmentation.h>
#include <CGAL/property_map.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel    Kernel;
typedef CGAL::Polyhedron_3<Kernel>                             Polyhedron;
typedef Kernel::Vector_3                                       Vector;
typedef Kernel::Point_3                                        Point;
typedef Polyhedron::Vertex                                     Vertex;
typedef Polyhedron::Facet_iterator                             Facet_iterator;
typedef Polyhedron::Vertex_iterator                            Vertex_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator           HF_circulator;
typedef Polyhedron::Halfedge_around_vertex_const_circulator    HV_circulator;

/*-------------------------------------------*\
|      Plot the distribution of the local     |
|  error on meshes with different magnitudes  |
|                  of noise.                  | 
\*-------------------------------------------*/

// PROTOTYPES //

std::vector<double> findMagnitudes ( char * magnitudeFile ) ;
float error( float a, float b );
float squaredError(float a, float b);
std::string convertDouble(double number);
std::vector<double> findBoxBorder ( Polyhedron P );
double longestDiagonalLength ( Polyhedron P );
float computeGlobalError( std::vector<float> localError );
int createNoiseOff ( char * fileName, const char * outputName, double magnitude, int option);
double findBoxDimension ( Polyhedron P );
std::vector<float> sdf_values ( const char * fileName );   	
void createGlobalErrorScript(std::string fileName, std::vector<double> globalError, std::vector<double> noiseMagnitudes );


int main ( int argc, char * argv[] ) {
    // argv[1] = filename
    // argv[2] = noise magnitudes file 
    // argv[3] = option ( for noise )
    if ( argc != 4 ) {
        std::cerr << "Not the good number of arguments" << std::endl ;
        return -1 ;
    }
    std::string fileName = argv[1] ;
    std::string sdfName = fileName.substr(0,fileName.size()-4)+"-sdf.txt";
    std::vector<float> sdfValues = sdf_values(argv[1]) ;  
    int size = sdfValues.size() ;
    // Ligne a enlever plus tard
    std::sort(sdfValues.begin(),sdfValues.end());
    std::ofstream file(sdfName.c_str(), std::ios::out);
    if( !file ){ 
        std::cerr << " Error in creating file " <<  std::endl;
        return -1 ;
    }
    int k = 0 ;
    // Put in file
    for(k=0;k<size;k++) {
        file << (k+1)*100.0/( (float) size) << " " << sdfValues[k] << std::endl ;
        };
    std::vector<double> noiseMagnitudes = findMagnitudes(argv[2]);
    unsigned int nbValues = noiseMagnitudes.size() ;
    //std::vector<double> globalError(nbValues) ;
    // Init script file header

    std::string scriptName = fileName.substr(0,fileName.size()-4)+"-script.p" ;
    std::ofstream scriptFile(scriptName.c_str(),std::ios::out) ;
    if ( !scriptFile ) {
        std::cerr << "Open script file error" << std::endl ;
        return -1 ;
    }  
    scriptFile << "reset " << std::endl ;
    scriptFile << std::endl ;
    scriptFile << "set term pdfcairo" << std::endl ;
    scriptFile << "set output \"" << argv[1] << "-local-error.pdf\"" << std::endl ;
    scriptFile << "set title \" Distribution of the local error \" " << std::endl ;
    scriptFile << "set xlabel \" Facets (%) \" " << std::endl ;
    scriptFile << "set ylabel \" Local error (%) \" " << std::endl ;
    scriptFile << std::endl ;
    scriptFile << "set logscale y 10" << std::endl ;
    scriptFile << "set style data line " << std::endl ;
    scriptFile << "set pointsize 2 " << std::endl ; 
    scriptFile << "set key bottom right " << std::endl ;
    scriptFile << std::endl ;
    scriptFile << "plot " ;
    
    unsigned int i ;
    for ( i = 0 ; i < nbValues ; i++ ) {
        // Compute the local error with the different noise magnitudes 
        std::string noiseName = fileName.substr(0,fileName.size()-4)+"s="+convertDouble(noiseMagnitudes[i])+"%.off" ;
        std::cout << "Create noised mesh... " ;	
        int create = createNoiseOff(argv[1],noiseName.c_str(),noiseMagnitudes[i],atoi(argv[3]));
        if ( create != 0) {
            std::cerr << " Error when create noised mesh " << std::endl ;
            return -1 ;  
        }
        std::cout << "Ok" << std::endl << "SDF VALUES ..." << std::endl ;
        std::vector<float> sdfNoisedValues = sdf_values(noiseName.c_str()) ; 
        if ( size == 0) {
            return -1 ;
        }
        
        std::vector<float> localError(size,0.0) ;
        // Important note : The facet numbering is the same for the different OFF files
        std::string errorResults = noiseName.substr(0,noiseName.size()-4)+"-local-error.txt";
        
        std::string sdfResults = noiseName.substr(0,noiseName.size()-4)+"-sdf.txt";
       // std::ofstream errorFile(errorResults.c_str(), std::ios::out);
        std::ofstream sdfFile(sdfResults.c_str(), std::ios::out);
        if( /*!errorFile ||*/ !sdfFile ){ 
            std::cerr << " Error in creating file " <<  std::endl;
            return -1 ;
        }
        // Calculations
        
        std::cout << "Compute local error... " ;
        
        int j = 0 ;
        
        for(j=0;j<size;j++) {
             localError[j] = error(sdfValues[j],sdfNoisedValues[j]);
        };
        // Sort the values 
        std::sort(localError.begin(),localError.end());
        
        std::sort(sdfNoisedValues.begin(),sdfNoisedValues.end());
        std::cout << "Ok" << std::endl << "Write in file ... " ;
        // Put in file
        for(j=0;j<size;j++) {
            errorFile << (j+1)*100.0/( (float) size) << " " << localError[j] << std::endl ;
            sdfFile << (j+1)*100.0/( (float) size) << " " << sdfNoisedValues[j] << std::endl ;
        };
        std::cout << "Ok" << std::endl ;
        scriptFile << "\""+errorResults+"\" using 1:2 title \"s="+convertDouble(noiseMagnitudes[i])+"%\", \\" << std::endl ;
        globalError[i] = computeGlobalError(localError) ; 
    } ;
    createGlobalErrorScript(fileName, globalError, noiseMagnitudes) ;
    return 0;
}

/********************\
| Auxilary functions |
\********************/

std::vector<double> findMagnitudes ( char * magnitudeFile ) {
  std::vector<double> magnitudes ;
  std::ifstream file(magnitudeFile,std::ios::in);
  if (!file) {
      std::cerr << "Error: open magnitude file" << std::endl ;
  }
  else {
      std::string word ;
      while (file >> word) {
          double magn = atof(word.c_str()) ;
          magnitudes.push_back(magn) ;
      }
  }
  file.close() ;
  return magnitudes ;
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


// return the global error if the function "error" is used ,
// or the mean squared error if the function "squaredError" is used.

float computeGlobalError( std::vector<float> localError ) {
  int size = localError.size() ;
  float moyenne = 0.0 ;
  for ( int i = 0; i < size; i++ ) {
      moyenne += localError[i] ;
  } 
  moyenne /= size ;
  return moyenne ;
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
        case 1 : {
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
        case 2 : {
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
        case 3 : {
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

double findBoxDimension ( Polyhedron P ) {
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
   double values[] = {borders[1]-borders[0],borders[3]-borders[2],borders[5]-borders[4] } ;
   return *std::max_element(values,values+3) ; 
}


std::vector<float> sdf_values ( const char * fileName ) {
  // create and read Polyhedron
  Polyhedron mesh;
  std::ifstream input(fileName);
  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr << "Not a valid .off file." << std::endl;
  }
  // To normalize
  std::cout << "Find box dimensions... " ;
  double factor = findBoxDimension(mesh) ;
  std::cout << "Ok" << std::endl ;
  // create a property-map
  typedef std::map<Polyhedron::Facet_const_handle, double> Facet_double_map;
  Facet_double_map internal_map;
  boost::associative_property_map<Facet_double_map> sdf_property_map(internal_map);
  // compute SDF values
  std::cout << "Compute SDF values ... " ;
  std::pair<double, double> min_max_sdf = CGAL::sdf_values(mesh, sdf_property_map);
  std::cout << "Ok" << std::endl ;
  int size = mesh.size_of_facets() ;
  // put SDF values in an array
  std::vector<float> values(size) ;
  int j = 0 ;
  std::cout << " Normalize ... " ;
  for(Polyhedron::Facet_const_iterator facet_it = mesh.facets_begin();
      facet_it != mesh.facets_end(); ++facet_it) {
      values[j] = ((min_max_sdf.second - min_max_sdf.first) * sdf_property_map[facet_it] + min_max_sdf.first) / (2*factor) ;
      j++;
  }
  std::cout << "Ok" << std::endl ;
  return values ;
}


void createGlobalErrorScript(std::string fileName, std::vector<double> globalError, std::vector<double> noiseMagnitudes ) {
    // Create file with results
    std::string resultsName = fileName.substr(0,fileName.size()-4)+"-global-error.txt" ; 
    std::ofstream resultsFile(resultsName.c_str(),std::ios::out) ;
    if ( !resultsFile ) {
        std::cerr << "Open global errors file error" << std::endl ;
    }
    unsigned int i ;
    for ( i = 0 ; i < globalError.size() ; i++ ) {
        resultsFile << noiseMagnitudes[i] << ' ' << globalError[i] << std::endl ;
    } 
    // Create Script    
    std::string scriptName = fileName.substr(0,fileName.size()-4)+"-global-error.p" ;
    std::ofstream scriptFile(scriptName.c_str(),std::ios::out) ;
    if ( !scriptFile ) {
        std::cerr << "Open script file error" << std::endl ;
    }     
    scriptFile << "reset " << std::endl ;
    scriptFile << std::endl ;
    scriptFile << "set term pdfcairo" << std::endl ;
    scriptFile << "set output \"" << fileName.substr(0,fileName.size()-4) << "-global-error.pdf\"" << std::endl ;
    scriptFile << "set title \" Global error for increasing noise magnitudes \" " << std::endl ;
    scriptFile << "set xlabel \" Noise magnitude (%) \" " << std::endl ;
    scriptFile << "set ylabel \" Global error (%) \" " << std::endl ;
    scriptFile << std::endl ;
    scriptFile << "set logscale x 10" << std::endl ;
    scriptFile << "set style data line " << std::endl ;
    scriptFile << "set pointsize 2 " << std::endl ; 
    scriptFile << "set key left top " << std::endl ;
    scriptFile << std::endl ;
    scriptFile << "plot \"" << resultsName << "\" using 1:2 " ;
}

