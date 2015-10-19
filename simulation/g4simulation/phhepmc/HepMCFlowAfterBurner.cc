//
// Inspired by code from ATLAS.  Thanks!
//
#include "HepMCFlowAfterBurner.h"
#include "PHHepMCGenEvent.h"

#include <flowAfterburner/flowAfterburner.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <iostream>
#include <cstdlib>
#include <string>
#include <memory>

#include <HepMC/GenEvent.h>

#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/MTwistEngine.h>
#include <CLHEP/Random/RandFlat.h>
#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Geometry/Point3D.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>


using namespace std;

CLHEP::HepRandomEngine * engine = NULL;

