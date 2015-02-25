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
#include <CGAL/bounding_box.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::Vertex_iterator Vertex_iterator;
typedef Kernel::Point_3 Point;
typedef std::map<Polyhedron::Facet_const_handle, double> Facet_double_map;

// Prototypes

std::vector<int> findResolutions ( char * resolutionFile ) ;
void afficheAide( void );
double findBoxDimension ( Polyhedron P );

int main ( int argc, char * argv[]) {
  if ( argc == 1 ) {
      afficheAide() ;
      return 0 ;
  }
  if ( argc != 2 ) {
      std::cerr << " Not the good number of arguments " << std::endl ;
      return -1 ;
  }
  std::string fileName= argv[1] ;

  // create and read Polyhedron

  Polyhedron mesh;
  std::ifstream input(argv[1]);
  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr << "Not a valid .off file." << std::endl;
    return -1;
  } 
 
  std::cout << "SDF Off loaded"<<std::endl;
  
  // create a property-map

  Facet_double_map internal_map;
  boost::associative_property_map<Facet_double_map> sdf_property_map(internal_map);

  // compute SDF values
  std::cout << "Computing sdf..."<<std::endl;

  const std::size_t number_of_rays = 25;  // cast 25 rays per facet
  const double cone_angle = 2.0 / 3.0 * CGAL_PI; // set cone opening-angle
  CGAL::sdf_values(mesh, sdf_property_map, cone_angle, number_of_rays, false);
  std::cout<<"Post processing..."<<std::endl;
  std::pair<double, double> min_max_sdf = CGAL::sdf_values_postprocessing(mesh, sdf_property_map);

  std::cout << "Exporting sdf..."<<std::endl;

  
  //SDF values are now postprocessed and thus in [0,1]
  //We use the minmax to reset the normalization
  // and scale from  the bounding box
  int size = mesh.size_of_facets() ;
  std::vector<double> values(size) ;
  int j = 0 ;
  double factor = findBoxDimension(mesh)  ;
  
  std::cout << "Scale factor (longest BBox edge)= "<<factor << std::endl;
  for(Polyhedron::Facet_const_iterator facet_it = mesh.facets_begin();
      facet_it != mesh.facets_end(); ++facet_it) {
      // Normalize ( real values normalized by the size of the bounding box + half-length)
      values[j] = ((min_max_sdf.second - min_max_sdf.first) * sdf_property_map[facet_it] + min_max_sdf.first) / (2.0*factor) ;
      j++;
  }
  
  //Sorting the sdf values
  sort(values.begin(), values.end());
  
  std::cout << "Number of samples= "<< size<< std::endl;
  std::cout << "Max sdf value= "<< values[size -1]<< std::endl;
  std::cout << "Min sdf value= "<< values[0]<< std::endl;
  std::cout << "Factor= "<< factor << std::endl;

  
  // Write in files
  std::string base_filename = fileName.substr(fileName.find_last_of("/\\") + 1);
  std::string::size_type const p(base_filename.find_last_of('.'));
  std::string file_without_extension = base_filename.substr(0, p);
  std::string fileResults = file_without_extension + "-sdf.txt" ;
  std::cout << "filename= "<<fileResults<<std::endl;
  
  std::ofstream fichier(fileResults.c_str(), std::ios::out);
  if (fichier) { 
      for(j=0;j<size;j++) {
          fichier << values[j] << std::endl ;
      }
      fichier.close() ;
  }
  else {
      std::cout << "open file error" << std::endl ;
      return -1 ; 
  }
  
  return 0 ;
}

/********************\
| AUXILIAR FUNCTIONS |
\********************/

void afficheAide( void ) {
    std::cout << " Compare the thickness distributions ( by the Shape Diameter Function and with a volumetric method, with different resolutions contained in a text file ) of an object in .off format." << std::endl ;
    std::cout << "Usage: ./sdf <off file> " << std::endl ;
    std::cout << "Output: all sorted sdf values (smoothed and normaliazed using longest edge of BBox)" << std::endl ;
}

double findBoxDimension ( Polyhedron P )
{
  std::vector<double> res;
  std::vector<Point> pts;
  for(Vertex_iterator it = P.vertices_begin(), itend=P.vertices_end(); it!= itend; ++it)
  {
    Point p = it->point();
    pts.push_back(p);
  }
  
  Kernel::Iso_cuboid_3 bbox=CGAL::bounding_box(pts.begin(), pts.end());
  
  double values[] = {(bbox.xmax() - bbox.xmin()) ,
                     (bbox.ymax() - bbox.ymin()),
    (bbox.zmax() - bbox.zmin())};
  
  return *std::max_element(values,values+3) ;
}


  
 
