#include "PHG4SectorDetector.h"

#include <g4main/PHG4RegionInformation.h>
#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <Geant4/G4Material.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>

#include <Geant4/G4VisAttributes.hh>
#include <Geant4/G4Colour.hh>

#include <sstream>
#include <cassert>

using namespace std;
using namespace PHG4Sector;

//_______________________________________________________________
//note this inactive thickness is ~1.5% of a radiation length
PHG4SectorDetector::PHG4SectorDetector(PHCompositeNode *Node,
    const std::string &dnam) :
  PHG4Detector(Node, dnam), PHG4SectorConstructor(dnam), _region(NULL)
{
}

//_______________________________________________________________
//_______________________________________________________________
bool
PHG4SectorDetector::IsInSectorActive(G4VPhysicalVolume * volume)
{
  for (map_phy_vol_t::const_iterator it = map_active_phy_vol.begin();
      it != map_active_phy_vol.end(); ++it)
    {
      if (volume == (*it).second)
        {
          return true;
        }
    }

  return false;
}

//_______________________________________________________________
bool
PHG4SectorDetector::IsInSectorInactive(G4VPhysicalVolume * volume)
{
  return false;
}

//_______________________________________________________________
void
PHG4SectorDetector::Construct(G4LogicalVolume* logicWorld)
{

//  // vis staff has not been set
//  assert(DetectorVisAtt == NULL);
//
//  DetectorVisAtt = new G4VisAttributes();
//  PHG4Utils::SetColour(DetectorVisAtt, material);
//  DetectorVisAtt->SetVisibility(true);
//  DetectorVisAtt->SetForceSolid(true);
////  DetectorVisAtt->SetForce

  Construct_Sectors(logicWorld);

//  BOOST_FOREACH( map_log_vol_t::value_type &vol_pair, map_log_vol )
//    {
//      _region->AddRootLogicalVolume(vol_pair.second);
//    }


  for (map_log_vol_t::iterator it = map_log_vol.begin(); it != map_log_vol.end();
      ++it)
    {
      if ((*it).first != G4String(name_base + "_Log"))
        {
          // sub layers

          DetectorVisAtt = new G4VisAttributes();
          PHG4Utils::SetColour(DetectorVisAtt, (*it).second->GetMaterial()->GetName());
          DetectorVisAtt->SetVisibility(true);
          DetectorVisAtt->SetForceSolid(true);
          (*it).second->SetVisAttributes(DetectorVisAtt);
        }
    }
}
