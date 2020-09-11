//____________________________________________________________________________..
//
// This is a working template for the G4 Construct() method which needs to be implemented
// We wedge a method between the G4 Construct() to enable volume hierarchies on the macro
// so here it is called ConstructMe() but there is no functional difference
// Currently this installs a simple G4Box solid, creates a logical volume from it
// and places it. Put your own detector in place (just make sure all active volumes
// get inserted into the m_PhysicalVolumesSet)
// 
// Rather than using hardcoded values you should consider using the parameter class
// Parameter names and defaults are set in PHG4TpcEndCapSubsystem::SetDefaultParameters()
// Only parameters defined there can be used (also to override in the macro)
// to avoids typos.
// IMPORTANT: parameters have no inherent units, there is a convention (cm/deg)
// but in any case you need to multiply them here with the correct CLHEP/G4 unit 
// 
// The place where you put your own detector is marked with
// //begin implement your own here://
// //end implement your own here://
// Do not forget to include the G4 includes for your volumes
//____________________________________________________________________________..

#include "PHG4TpcEndCapDetector.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Color.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cmath>
#include <iostream>

class G4VSolid;
class PHCompositeNode;

using namespace std;

//____________________________________________________________________________..
PHG4TpcEndCapDetector::PHG4TpcEndCapDetector(PHG4Subsystem *subsys,
                                         PHCompositeNode *Node,
                                         PHParameters *parameters,
                                         const std::string &dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_Params(parameters)
{
}

//_______________________________________________________________
int PHG4TpcEndCapDetector::IsInDetector(G4VPhysicalVolume *volume) const
{
  set<G4VPhysicalVolume *>::const_iterator iter = m_PhysicalVolumesSet.find(volume);
  if (iter != m_PhysicalVolumesSet.end())
  {
    return 1;
  }
  return 0;
}

//_______________________________________________________________
void PHG4TpcEndCapDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
 //begin implement your own here://
 // Do not forget to multiply the parameters with their respective CLHEP/G4 unit !
  double xdim = m_Params->get_double_param("size_x") * cm;
  double ydim = m_Params->get_double_param("size_y") * cm;
  double zdim = m_Params->get_double_param("size_z") * cm;
  G4VSolid *solidbox = new G4Box("PHG4TpcEndCapSolid", xdim / 2., ydim / 2., zdim / 2.);
  G4LogicalVolume *logical = new G4LogicalVolume(solidbox, G4Material::GetMaterial(m_Params->get_string_param("material")), "PHG4TpcEndCapLogical");

  G4VisAttributes *vis = new G4VisAttributes(G4Color(G4Colour::Grey()));  // grey is good to see the tracks in the display
  vis->SetForceSolid(true);
  logical->SetVisAttributes(vis);
  G4RotationMatrix *rotm = new G4RotationMatrix();
  rotm->rotateX(m_Params->get_double_param("rot_x") * deg);
  rotm->rotateY(m_Params->get_double_param("rot_y") * deg);
  rotm->rotateZ(m_Params->get_double_param("rot_z") * deg);

  G4VPhysicalVolume *phy = new G4PVPlacement(
      rotm,
      G4ThreeVector(m_Params->get_double_param("place_x") * cm,
                    m_Params->get_double_param("place_y") * cm,
                    m_Params->get_double_param("place_z") * cm),
      logical, "PHG4TpcEndCap", logicWorld, 0, false, OverlapCheck());
  // add it to the list of placed volumes so the IsInDetector method
  // picks them up
  m_PhysicalVolumesSet.insert(phy);
 //end implement your own here://
  return;
}

//_______________________________________________________________
void PHG4TpcEndCapDetector::Print(const std::string &what) const
{
  cout << "PHG4TpcEndCap Detector:" << endl;
  if (what == "ALL" || what == "VOLUME")
  {
    cout << "Version 0.1" << endl;
    cout << "Parameters:" << endl;
    m_Params->Print();
  }
  return;
}
