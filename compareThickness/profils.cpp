#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/mesh_segmentation.h>

#include <CGAL/property_map.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::Vertex_iterator Vertex_iterator;
 
std::vector<int> findResolutions ( char * resolutionFile ) ;
void afficheAide( void );
double findBoxDimension ( Polyhedron P );

int main ( int argc, char * argv[] ) {
    if ( argc == 1 ) {
        afficheAide() ;
        return 0 ;
    }
    if ( argc != 3 ) {
        std::cerr << " Not the good number of arguments " << std::endl ;
        return -1 ;
    }
    std::string fileName= argv[1] ;

    // Create gnuplot scripr
    std::ofstream scriptFile("script.p",std::ios::out) ;
    if ( !scriptFile ) {
        std::cerr << "Open script file error" << std::endl ;
        return -1 ;
    }  
    scriptFile << "reset " << std::endl ;
    scriptFile << std::endl ;
    // If you want/don't want to save the curves, uncomment/comment these two lines :
    scriptFile << "set term pdfcairo" << std::endl ;
    scriptFile << "set output \"" << argv[1] << "-compare.pdf\"" << std::endl ;
    scriptFile << "set title \" Thickness distribution \" " << std::endl ;
    scriptFile << "set xlabel \" Cumulative distribution  \" " << std::endl ;
    scriptFile << "set ylabel \" Radius \" " << std::endl ;
    scriptFile << std::endl ;
    scriptFile << "set style data line " << std::endl ;
    scriptFile << "set pointsize 2 " << std::endl ; 
    scriptFile << "set key bottom right " << std::endl ;
    scriptFile << std::endl ;
    scriptFile << "plot " ;

    // SDF computation //

    // create and read Polyhedron

  Polyhedron mesh;
  std::ifstream input(argv[1]);
  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr << "Not a valid .off file." << std::endl;
    return -1;
  }  
  // create a property-map

  typedef std::map<Polyhedron::Facet_const_handle, double> Facet_double_map;
  Facet_double_map internal_map;
  boost::associative_property_map<Facet_double_map> sdf_property_map(internal_map);
  // compute SDF values

  std::pair<double, double> min_max_sdf = CGAL::sdf_values(mesh, sdf_property_map);
  std::cout << " SDF min value : " << min_max_sdf.first << ", SDF max value : " << min_max_sdf.second << std::endl ; 

  int size = mesh.size_of_facets() ;
  // put SDF values in an array

  std::vector<float> values(size) ;
  int j = 0 ;
  double factor = findBoxDimension(mesh) ;
  for(Polyhedron::Facet_const_iterator facet_it = mesh.facets_begin();
      facet_it != mesh.facets_end(); ++facet_it) {
      // Normalize ( real values normalized by the size of the bounding box )
      values[j] = ((min_max_sdf.second - min_max_sdf.first) * sdf_property_map[facet_it] + min_max_sdf.first) / (2*factor) ;
      j++;
  }
  // sort the SDF values
  sort(values.begin(),values.end());
  // Write in files
  std::string fileResults = fileName.substr(0,fileName.size()-4)+"-sdf.txt";
  std::ofstream fichier(fileResults.c_str(), std::ios::out);
  if( fichier ){ 
      for(j=0;j<size;j++) {
          fichier << ((double) (j+1)) / ((double) size) * 100 << " " << values[j] << std::endl ;
      }
      fichier.close() ;
      scriptFile << "\""+fileResults+"\" using 1:2 title \"SDF\", \\" << std::endl ;
  }
  else {
  std::cout << "open file error" << std::endl ;
  return -1 ; 
  } ;

    // Volumic calculations //

    std::vector<int> resolutions = findResolutions(argv[2]) ;
    unsigned int taille = resolutions[0] ; // see findResolutions
    if ( taille == 0) {
        std::cerr << "Empty array" << std::endl ;
        return -1 ;
    }

    // Calculations
    unsigned int i ;
    for (i=0; i<taille; i++) {
        std::ostringstream oss;  
        oss << resolutions[i+1] ;
        std::string resol = oss.str();
        std::string systemCall = "bash profil_vol.sh "+fileName+" "+resol;
        int s = system(systemCall.c_str()) ;
        if ( s != 0 ) {
            std::cerr << "profil_vol call failed" << std::endl ;
            return -1 ;
        } ;
        // Fill the script file
        scriptFile << "\""+fileName.substr(0,fileName.size()-4)+"-"+resol+"-newdata.txt\" using 1:2 title \"resol-"+resol+"\", \\" << std::endl ;  
        std::cout << "ok" << std::endl ;        
    }
    scriptFile.close() ;
    return 0;
}

/********************\
| AUXILIAR FUNCTIONS |
\********************/

std::vector<int> findResolutions ( char * resolutionFile ) {
  std::vector<int> resolutions(1,0) ; // first = size 
  std::ifstream file(resolutionFile,std::ios::in);
  if (!file) {
      std::cerr << "Error: open resolution file" << std::endl ;
  }
  else {
      std::string line ;
      while (getline(file,line)) {
          int resol = atoi(line.c_str()) ;
          resolutions.push_back(resol) ;
          resolutions[0]++ ;
      }
  }
  file.close() ;
  return resolutions ;
}

void afficheAide( void ) {
    std::cout << " Compare the thickness distributions ( by the Shape Diameter Function and with a volumetric method, with different resolutions contained in a text file ) of an object in .off format." << std::endl ;
    std::cout << "Usage: ./profils <object file> <resolution file>" << std::endl ;
    std::cout << "Output: a gnuplot script and a pdf figure with the distributions." << std::endl ;
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

  // Note //
  
  // It is possible to compute the raw SDF values and post-process them using
  // the following lines: ( instead of line 71 )
  // const std::size_t number_of_rays = 25;  // cast 25 rays per facet
  // const double cone_angle = 2.0 / 3.0 * CGAL_PI; // set cone opening-angle
  // CGAL::sdf_values(mesh, sdf_property_map, cone_angle, number_of_rays, false);
  // std::pair<double, double> min_max_sdf =
  //  CGAL::sdf_values_postprocessing(mesh, sdf_property_map);


  
 


