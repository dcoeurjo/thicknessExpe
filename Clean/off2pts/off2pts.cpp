#include <iostream>
#include <fstream>
#include <random>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel    Kernel;
typedef CGAL::Polyhedron_3<Kernel>                             Polyhedron;
typedef Kernel::Vector_3                                       Vector;
typedef Kernel::Point_3                                        Point;
typedef Polyhedron::Vertex                                     Vertex;
typedef Polyhedron::Facet_iterator                             Facet_iterator;
typedef Polyhedron::Vertex_iterator                            Vertex_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator           HF_circulator;
typedef Polyhedron::Halfedge_around_vertex_const_circulator    HV_circulator;

std::string convertNumber(int number);
std::string convertNumber(double number);
std::pair<int,std::string> whatOption ( std::string& fileName, char * argv, double magnitude );
int createNoiseOff ( char * fileName, const char * outputName, int option, double s );
std::vector<double> findBoxBorder ( Polyhedron P );
double longestDiagonalLength ( Polyhedron P );


int main ( int argc, char * argv[] ) {
  // argv[1] = filename ( .off )
  
  std::string fileName = argv[1] ;

  std::cout <<  "Loading OFF..."<<std::endl;
  Polyhedron mesh ;
  std::ifstream input(fileName) ;
  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr << "Not a valid .off file" << std::endl;
  };
  std::cout << "Exporting..."<<std::endl;
  std::string outputName = fileName.substr(0,fileName.size()-4)+".pts" ;
  std::ofstream output(outputName,std::ios::out) ;
  
  for(Vertex_iterator v = mesh.vertices_begin(), vend=mesh.vertices_end();
      v!=vend; ++v)
    {
      output << v->point().x()<<" "<<v->point().y()<< " "<<v->point().z()<<std::endl;
    }
  
  return 0 ;
}

