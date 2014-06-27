/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

/**
 * @file voAddBorder.cpp
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
#include <DGtal/io/readers/VolReader.h>
#include <DGtal/io/writers/VolWriter.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/Image.h>
#include <DGtal/images/ImageContainerBySTLVector.h>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
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
    ( "input,i", po::value<std::string>(), "Input vol file." )
    ( "output,o", po::value<string>(),"Output filename." );
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< endl;
  }
  po::notify ( vm );
  if ( !parseOK || vm.count ( "help" ) ||argc<=1 )
    {
      trace.info() << "Fill the shape by filling-out the exterior (+complement)"<<std::endl
                   << std::endl << "Basic usage: "<<std::endl
                   << "\tvolFillExterior --input <volFileName> --o <volOutputFileName> "<<std::endl
                   << general_opt << "\n";
      return 0;
    }

  //Parse options
  if ( ! ( vm.count ( "input" ) ) ) missingParam ( "--input" );
  std::string filename = vm["input"].as<std::string>();
  if ( ! ( vm.count ( "output" ) ) ) missingParam ( "--output" );
  std::string outputFileName = vm["output"].as<std::string>();

  typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char>  MyImageC;

  Z3i::Point px(1,0,0);
  Z3i::Point py(0,1,0);
  Z3i::Point pz(0,0,1);

  trace.beginBlock("Loading..");
  MyImageC  imageL = VolReader< MyImageC >::importVol ( filename );
  trace.endBlock();

  trace.beginBlock("Dilating..");
  MyImageC  imageC(imageL.domain());

  for(MyImageC::Range::Iterator it = imageL.range().begin(),itC = imageC.range().begin(),
        itend = imageL.range().end(); it != itend; ++it, ++itC)
    *itC = *it;

  for(MyImageC::Domain::ConstIterator it = imageC.domain().begin(),
        itend = imageC.domain().end(); it != itend; ++it)
    {
      unsigned char val = imageC(*it);

      if (val ==0)
        {

          if ((imageC.domain().isInside(*it - px)))
            val = val + imageL(*it - px);
          if ((imageC.domain().isInside(*it + px)))
            val = val + imageL(*it + px);

          if ((imageC.domain().isInside(*it - py)))
            val = val + imageL(*it - py);
          if ((imageC.domain().isInside(*it + py)))
            val = val + imageL(*it + py);

          if ((imageC.domain().isInside(*it - pz)))
            val = val + imageL(*it - pz);
          if ((imageC.domain().isInside(*it + pz)))
            val = val + imageL(*it + pz);

          if (val >1)
            imageC.setValue(*it, 1);
        }
    }
  trace.endBlock();




  trace.beginBlock("Saving");
  VolWriter< MyImageC>::exportVol("epaiss.vol", imageC);
  trace.endBlock();



  std::stack<Z3i::Point> pstack;
  pstack.push(*(imageC.domain().begin()));
  trace.beginBlock("filling");
  while (!pstack.empty())
    {
      Z3i::Point p= pstack.top();
      pstack.pop();
      imageC.setValue(p, 113);

      if ((imageC.domain().isInside(p - px)) &&
          (imageC( p - px) == 0))
        pstack.push( p -px);

      if ((imageC.domain().isInside(p - py)) &&
          (imageC( p - py) == 0))
        pstack.push( p -py);

      if ((imageC.domain().isInside(p - pz)) &&
          (imageC( p - pz) == 0))
        pstack.push( p -pz);

      if ((imageC.domain().isInside(p + px)) &&
          (imageC( p + px) == 0))
        pstack.push( p +px);

      if ((imageC.domain().isInside(p + py)) &&
          (imageC( p + py) == 0))
        pstack.push( p +py);

      if ((imageC.domain().isInside(p + pz)) &&
          (imageC( p + pz) == 0))
        pstack.push( p +pz);

    }
  trace.endBlock();

  trace.beginBlock("Complement");
   //Fastcopy
  for(MyImageC::Range::Iterator it = imageC.range().begin(),
        itend = imageC.range().end(); it != itend; ++it)
    if (*it == 113 )
      *it = 0;
    else
      *it = 128;
  trace.endBlock();

  trace.beginBlock("Saving");
  bool res =  VolWriter< MyImageC>::exportVol(outputFileName, imageC);
  trace.endBlock();

  if (res) return 0; else return 1;
}
