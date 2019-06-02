// $$Id: PHG4RICHDetector.cc,v 1.1 2013/10/01 00:33:00 jinhuang Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.1 $$
 * \date $$Date: 2013/10/01 00:33:00 $$
 */

#include "PHG4RICHDetector.h"
#include "PHG4RICHSteppingAction.h"

#include <g4main/PHG4Detector.h>     // for PHG4Detector

#include <Geant4/G4Region.hh>        // for G4Region

#include <boost/foreach.hpp>

#include <map>
#include <sstream>

class PHCompositeNode;

using namespace std;
using namespace ePHENIXRICH;

PHG4RICHDetector::PHG4RICHDetector(PHG4RICHSubsystem *subsys, PHCompositeNode *Node, const RICH_Geometry &g)
  : PHG4Detector(Node)
  , ePHENIXRICHConstruction(subsys, g)
  , _region(nullptr)
{
}

PHG4RICHDetector::PHG4RICHDetector(PHG4RICHSubsystem *subsys, PHCompositeNode *Node)
  : PHG4Detector(Node)
  , ePHENIXRICHConstruction(subsys)
  , _region(nullptr)
{
}

void PHG4RICHDetector::Construct(G4LogicalVolume *logicWorld)
{
  _region = new G4Region("FCALREGION");
  _region->SetRegionalSteppingAction(new PHG4RICHSteppingAction(this));

  ePHENIXRICHConstruction::Construct_RICH(logicWorld);

  BOOST_FOREACH (map_log_vol_t::value_type &vol_pair, map_log_vol)
  {
    _region->AddRootLogicalVolume(vol_pair.second);
  }
}
