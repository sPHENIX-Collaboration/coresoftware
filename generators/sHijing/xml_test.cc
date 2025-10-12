// Allows the user to generate hijing events and store the results in
// a HepMC file.
//
// Inspired by code from ATLAS.  Thanks!
//
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <iostream>
#include <string>

#define f2cFortran
#define gFortran

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#include "cfortran.h"
#pragma GCC diagnostic pop

//using namespace boost;

float atl_ran(int * /*unused*/)
{
  return 0.0;
}
// This prevents cppcheck to flag the next line as error
// cppcheck-suppress *
FCALLSCFUN1(FLOAT, atl_ran, ATL_RAN, atl_ran, PINT)


int main() // NOLINT(bugprone-exception-escape)
{
  boost::property_tree::iptree pt;
  boost::property_tree::iptree null;

  std::ifstream config_file("xml_test.xml");
  if (config_file)
  {
    read_xml(config_file, pt);
  }

  boost::property_tree::iptree &it = pt.get_child("HIJING.FASTJET", null);
  BOOST_FOREACH (boost::property_tree::iptree::value_type &v, it)
  {
    if (boost::to_upper_copy(v.first) != "ALGORITHM")
    {
      continue;
    }
    std::string name = v.second.get("NAME", "ANTIKT");
    float R = v.second.get("R", 0.4);
    std::cout << name << " " << R << std::endl;
  }

  return 0;
}
