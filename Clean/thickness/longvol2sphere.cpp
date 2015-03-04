/**
 * @file longvol2sphere.cpp
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2012/05/01
 *
 *
 * This file is part of the DGtal library.
 */

#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/io/readers/LongvolReader.h>
#include <DGtal/images/Image.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/base/BasicFunctors.h>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace DGtal;
using namespace Z3i;

namespace po = boost::program_options;



/**
 * Missing parameter error message.
 *
 * @param param
 */
void missingParam ( std::string param )
{
  trace.error() <<" Parameter: "<<param<<" is required..";
  trace.info() <<std::endl;
  exit ( 1 );
}

int main(int argc, char**argv)
{
  
  // parse command line ----------------------------------------------
  po::options_description general_opt ( "Allowed options are: " );
  general_opt.add_options()
  ( "help,h", "display this message." )
  ( "input,i", po::value<std::string>(), "Input longvol filename." )
  ;

  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
  }
  po::notify ( vm );
  if (!parseOK || vm.count ( "help" ) ||argc<=1 )
  {
    trace.info() << "Convert a longvol (long int) to a set of spheres files."<<std::endl
    << std::endl << "Basic usage: "<<std::endl
    << "\tlongvol2vol --input <LongvolFileName> "<<std::endl
    << general_opt << "\n";
    return 0;
  }
  
  //Parse options
  if ( ! ( vm.count ( "input" ) ) ) missingParam ( "--input" );
  std::string filename = vm["input"].as<std::string>();
 
  
  //Main program
  typedef ImageContainerBySTLVector<Z3i::Domain, DGtal::uint64_t>  MyImageC;
  MyImageC  imageC = LongvolReader< MyImageC >::importLongvol ( filename );
  
  for(MyImageC::Domain::ConstIterator it = imageC.domain().begin(),
        itend = imageC.domain().end(); it != itend;
      it++)
    if ( imageC(*it) != 0)
      std::cout << (*it)[0]<< " "<< (*it)[1]<< " "<< (*it)[2]<< " "<< sqrt(imageC(*it))<<std::endl;
  
  return 0;
}
