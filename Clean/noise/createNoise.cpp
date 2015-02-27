#include <iostream>
#include <fstream>
#include <random>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/bounding_box.h>

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




int main ( int argc, char * argv[] ) {
  // argv[1] = filename ( .off )
  // argv[2] = type of noise
  // argv[3] = magnitude of the noise
  
  std::string fileName = argv[1] ;
  double magnitude = atof(argv[3]) ;
  // Parse option
  std::pair<int,std::string> option = whatOption(fileName,argv[2],magnitude) ;
  if ( option.first == 0 ) {
    std::cerr << "Error : incorrect option ! " << std::endl ;
    std::cout << "-c : cubic" << std::endl << "-s : spheric" << std::endl << "-n : normal" << std::endl ;
    return -1 ;
  }
  
  // Create the noised OFF files
  int res = createNoiseOff ( argv[1], (option.second).c_str(), option.first, magnitude ) ;
  if ( res != 0 ) {
    std::cerr << " Error when create noise " << std::endl ;
    return -1 ;
  }
  
  std::cout<<option.second << std::endl;
  return 0 ;
}

/*-----------------------*\
 | Fonctions auxilliaires  |
 \*-----------------------*/

std::string convertNumber(int number) {
  std::stringstream ss;
  ss << number;
  return ss.str();
}

std::string convertNumber(double number) {
  std::stringstream ss;
  ss <<  std::setprecision(11)<<number;
  return ss.str();
}

/*--------------------------------*\
 | Choice of the noise : option     |
 |  - 1 : gaussian on a cube        |
 |  - 2 : gaussian on a sphere      |
 |  - 3 : gaussian along the normal |
 \*--------------------------------*/

std::pair<int,std::string> whatOption ( std::string& fileName , char * argv, double magnitude ) {
  int option ;
  
  std::string base_filename = fileName.substr(fileName.find_last_of("/\\") + 1);
  std::string::size_type const p(base_filename.find_last_of('.'));
  std::string file_without_extension = base_filename.substr(0, p);
  std::string noisedName = file_without_extension ;
  
  switch ( argv[1] ) {
    case 'c' :
      option = 1 ;
      noisedName = noisedName+"-cubic.off" ;
      break ;
    case 's' :
      option = 2 ;
      noisedName = noisedName+"-s="+convertNumber(magnitude)+"-spheric.off" ;
      break ;
    case 'n' :
      option = 3 ;
      noisedName = noisedName+"-s="+convertNumber(magnitude)+"-normal.off" ;
      break ;
    default :
      option = 0 ;
      noisedName = "echec" ;
  }
  return std::pair<int,std::string>(option,noisedName) ;
}

int createNoiseOff ( char * fileName, const char * outputName, int option, double magnitude ) {
  // Create and read Polyhedron
  std::cout << "Loading Input mesh..."<<std::endl;
  
  Polyhedron mesh ;
  std::ifstream input(fileName) ;
  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr << "Not a valid .off file" << std::endl;
  };
  std::cout << "Noisification..."<<std::endl;
  
  // Init generator
  double s = magnitude/100*findBoxDimension(mesh);
  std::cout<< "longest bbox edge= "<<findBoxDimension(mesh)<<std::endl;
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0.0,s);
  // apply different noises
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
    case 2 : {      // in a sphere
      Vertex_iterator v = mesh.vertices_begin() ;
      do {
        bool inSphere = false ;
        double x, y, z ;
        while(!inSphere) {
          x = distribution(generator);
          y = distribution(generator);
          z = distribution(generator);
          if ( x*x + y*y + z*z <= s*s) {
            inSphere = true ;
          }
        }
        Vector vec(x,y,z);
        v->point() = v->point() + vec ;
      } while ( v++ != mesh.vertices_end()) ;
      break ;
    }
    case 3 : {     // along the normal
      Vertex_iterator v = mesh.vertices_begin() ;
      do {
        HV_circulator h = v->vertex_begin() ;
        Vector vec(0.0,0.0,0.0) ;
        do {
          vec = vec + CGAL::cross_product(
                                          h->next()->vertex()->point() - h->vertex()->point(),
                                          h->next()->next()->vertex()->point() - h->next()->vertex()->point()) ;
          ++h ;
          CGAL_assertion( h != v->vertex_begin()); // even degree guaranteed
          ++h ;
        } while ( h != v->vertex_begin());
        double norme = sqrt(vec.x()*vec.x() + vec.y()*vec.y() + vec.z()*vec.z()) ;
        double factor = distribution(generator) ;
        Vector vectorNoise( vec.x()*factor/norme, vec.y()*factor/norme, vec.z()*factor/norme) ;
        v->point() = v->point() + vectorNoise ;
        v++ ;
      } while ( v != mesh.vertices_end()) ;
      break ;
    }
    default :
      std::cerr << " Erreur de saisie " << std::endl ;
  };
  // Write polyhedron in .off format.
  std::cout << "Exporting..."<<std::endl;
  std::ofstream offOutput(outputName,std::ios::out) ;
  CGAL::set_ascii_mode( offOutput ) ;
  offOutput << mesh;
  
  /*
  offOutput << "OFF" << std::endl ;
  offOutput << mesh.size_of_vertices() << ' ' << mesh.size_of_facets() << ' ' << mesh.size_of_halfedges() << std::endl ;
  std::cout << "    vertices.."<<std::endl;
  std::copy( mesh.points_begin(), mesh.points_end(), std::ostream_iterator<Point>( offOutput, "\n"));
  std::cout << "    facets.."<<std::endl;

  
  
  
  
  for ( Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
  {
    HF_circulator j = i->facet_begin();
    // Facets in polyhedral surfaces are at least triangles.
    CGAL_assertion( CGAL::circulator_size(j) >= 3);
    offOutput << CGAL::circulator_size(j) << ' ' << std::distance(mesh.vertices_begin(),
                                                                  i->facet_begin()->vertex()) ;
    while ( ++j != i->facet_begin())
    {
	     offOutput << ' ' << std::distance(mesh.vertices_begin(), j->vertex());
    }
    offOutput << std::endl;
  }*/
  std::cout << "Done"<<std::endl;

  return 0;
}

