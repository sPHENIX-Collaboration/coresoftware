#include "PHG4SectorDetector.h"
#include "PHG4SectorDisplayAction.h"

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>      // for PHG4Subsystem

#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4String.hh>  // for G4String

#include <map>  // for _Rb_tree_iterator, _Rb_tree_co...
#include <sstream>
#include <utility>  // for pair

class G4VPhysicalVolume;
class PHCompositeNode;

//_______________________________________________________________
//note this inactive thickness is ~1.5% of a radiation length
PHG4SectorDetector::PHG4SectorDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, const std::string &dnam)
  : PHG4Detector(subsys, Node, dnam)
  , PHG4SectorConstructor(dnam, subsys)
  , m_DisplayAction(dynamic_cast<PHG4SectorDisplayAction *>(subsys->GetDisplayAction()))
{
}

//_______________________________________________________________
//_______________________________________________________________
bool PHG4SectorDetector::IsInSectorActive(G4VPhysicalVolume *physvol)
{
  for (map_phy_vol_t::const_iterator it = map_active_phy_vol.begin();
       it != map_active_phy_vol.end(); ++it)
  {
    if (physvol == (*it).second)
    {
      return true;
    }
  }

  return false;
}

//_______________________________________________________________
void PHG4SectorDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  Construct_Sectors(logicWorld);

  for (auto &it : map_log_vol)
  {
    if (it.first != G4String(name_base + "_Log"))
    {
      // sub layers
      m_DisplayAction->AddVolume(it.second, "SectorDetector");
    }
  }
}
