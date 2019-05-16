#include "PHG4SectorDetector.h"
#include "PHG4SectorDisplayAction.h"
#include "PHG4SectorSubsystem.h"

#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <Geant4/G4Cons.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>

#include <Geant4/G4Colour.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cassert>
#include <sstream>

using namespace std;
using namespace PHG4Sector;

//_______________________________________________________________
//note this inactive thickness is ~1.5% of a radiation length
PHG4SectorDetector::PHG4SectorDetector(PHG4SectorSubsystem *subsys, PHCompositeNode *Node, const std::string &dnam)
  : PHG4Detector(Node, dnam)
  , PHG4SectorConstructor(dnam, subsys)
  , m_DisplayAction(dynamic_cast<PHG4SectorDisplayAction *>(subsys->GetDisplayAction()))
{
}

//_______________________________________________________________
//_______________________________________________________________
bool PHG4SectorDetector::IsInSectorActive(G4VPhysicalVolume *volume)
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
bool PHG4SectorDetector::IsInSectorInactive(G4VPhysicalVolume *volume)
{
  return false;
}

//_______________________________________________________________
void PHG4SectorDetector::Construct(G4LogicalVolume *logicWorld)
{
  Construct_Sectors(logicWorld);

  for (map_log_vol_t::iterator it = map_log_vol.begin(); it != map_log_vol.end();
       ++it)
  {
    if ((*it).first != G4String(name_base + "_Log"))
    {
      // sub layers
      m_DisplayAction->AddVolume((*it).second, "SectorDetector");
    }
  }
}
