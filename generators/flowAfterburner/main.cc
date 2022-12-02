//
// Inspired by code from ATLAS.  Thanks!
//

#include "flowAfterburner.h"

#include <HepMC/GenEvent.h>
#include <HepMC/IO_BaseClass.h>
#include <HepMC/IO_GenEvent.h>

#include <CLHEP/Random/MTwistEngine.h>

#include <boost/version.hpp>  // to get BOOST_VERSION
#if (__GNUC__ == 4 && __GNUC_MINOR__ == 8 && (BOOST_VERSION == 106000 || BOOST_VERSION == 106700 || BOOST_VERSION == 107000))
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma message "ignoring bogus gcc warning in boost header ptree.hpp"
#include <boost/property_tree/ptree.hpp>
#pragma GCC diagnostic warning "-Wshadow"
#pragma GCC diagnostic warning "-Wunused-parameter"
#else
#include <boost/property_tree/ptree.hpp>
#endif

#include <boost/operators.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <algorithm>
#include <sstream>
#include <string>

namespace CLHEP
{
  class HepRandomEngine;
}

CLHEP::HepRandomEngine *engine;

int main()
{
  using namespace boost::property_tree;
  iptree pt;
  // These values (coded here or in the xml file are used only in the
  // flowAfterburner executable. If you want to adjust them in the flow
  // module of our simulations you need to change
  // generators/phhepmc/HepMCFlowAfterBurner.cc
  // for the default and/or the values set in the macro

  std::ifstream config_file("flowAfterburner.xml");

  if (config_file)
  {
    // Read XML configuration file.
    read_xml(config_file, pt);
  }
  long randomSeed = pt.get("FLOWAFTERBURNER.RANDOM.SEED", 11793);
  engine = new CLHEP::MTwistEngine(randomSeed);

  std::string input = pt.get("FLOWAFTERBURNER.INPUT", "sHijing.dat");
  std::string output = pt.get("FLOWAFTERBURNER.OUTPUT", "flowAfterburner.dat");

  float mineta = pt.get("FLOWAFTERBURNER.CUTS.MINETA", -5.0);
  float maxeta = pt.get("FLOWAFTERBURNER.CUTS.MAXETA", 5.0);

  float minpt = pt.get("FLOWAFTERBURNER.CUTS.MINPT", 0.0);
  float maxpt = pt.get("FLOWAFTERBURNER.CUTS.MAXPT", 100.0);

  std::string algorithmName = pt.get("FLOWAFTERBURNER.ALGORITHM", "MINBIAS");

  // Open input file.
  HepMC::IO_GenEvent ascii_in(input.c_str(), std::ios::in);
  HepMC::IO_GenEvent ascii_out(output.c_str(), std::ios::out);
  HepMC::GenEvent *evt;

  while (ascii_in >> evt)
  {
    flowAfterburner(evt, engine, algorithmName, mineta, maxeta, minpt, maxpt);

    ascii_out << evt;
    delete evt;
  }
}
