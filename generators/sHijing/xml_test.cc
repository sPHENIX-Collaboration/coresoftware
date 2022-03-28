// Allows the user to generate hijing events and store the results in
// a HepMC file.
//
// Inspired by code from ATLAS.  Thanks!
//
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <iostream>
#include <string>


#define f2cFortran
#define gFortran
#include "cfortran.h"

using namespace boost;
using namespace std;

float
atl_ran(int *)
{
  return 0.0;
}
// This prevents cppcheck to flag the next line as error
// cppcheck-suppress *
FCALLSCFUN1 (FLOAT, atl_ran, ATL_RAN, atl_ran, PINT)

using namespace boost::property_tree;

int
main ()
{
  iptree pt, null;

  std::ifstream config_file("xml_test.xml");
  if (config_file)
    {
      read_xml (config_file, pt);
    }

  iptree &it = pt.get_child("HIJING.FASTJET", null);
  BOOST_FOREACH(iptree::value_type &v, it)
    {
      if (to_upper_copy(v.first) != "ALGORITHM") continue;
      string name = v.second.get("NAME", "ANTIKT");
      float R = v.second.get("R",0.4);
      cout << name << " " << R << endl;
    }

  return 0;
}
