#ifndef SPHERETOOLS_H_INCLUDED
#define SPHERETOOLS_H_INCLUDED

#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <math.h>
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

Spheres points ;
std::vector<double> maxRadius ;

class Calcul {
public:

    // init points and maxRadius

  Calcul(int resolution) {
      double i, j, k ;
      for ( k = 0 ; k < resolution ; k++ ) { 
          for ( j = 0 ; j < resolution ; j++ ) {
              for ( i = 0 ; i < resolution ; i++) {
                  points.push_back(Sphere(Point(i,j,k),0)) ;
                  maxRadius.push_back(0.0) ;
              }
          }
      }
  }


  bool belongsTo ( const Sphere& s1, const Sphere& s2 ) {
      Point p = s1.center() ;
      Point q = s2.center() ;
      double x = p.x() - q.x() ;
      double y = p.y() - q.y() ;
      double z = p.z() - q.z() ;
      return ( x*x+y*y+z*z <= s1.squared_radius() ) ;
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


/************\
| Converters |
\************/

std::string intToString(int number);
std::string doubleToString(double number);
int stringToInt ( const std::string &s );
int stringToDouble ( const std::string &s );

/**************\
| Sphere tools |
\**************/

std::vector<Sphere> readSpheres ( const char * filename );

/***********************\
| Tools to create noise |
\***********************/

std::vector<double> findBoxBorder ( Polyhedron P );
double longestDiagonalLength ( Polyhedron P );
int createNoiseOff ( char * fileName, const char * outputName, double magnitude, int option);

/****************************************\
| thickness and local error calculations |
\****************************************/

std::vector<double> getMagnitudes ( char * resolutionFile );
float error( float a, float b );
std::vector<double> offToThick(std::string filename, int resolution, Calcul& calc );
std::vector<double> spheresToThick( const char * filename , int resolution , int option, Calcul& calc );

/******************\
| Arrays functions |
\******************/

void initMaxRadius(int resolution);
void freeArrays( void ) ;


#endif 
