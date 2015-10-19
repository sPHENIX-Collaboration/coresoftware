#include "PHG4CylinderDetector.h"
#include "PHG4CylinderRegionSteppingAction.h"

#include <g4main/PHG4PhenixDetector.h>
#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cmath>
#include <sstream>

using namespace std;

//_______________________________________________________________
//note this inactive thickness is ~1.5% of a radiation length
PHG4CylinderDetector::PHG4CylinderDetector( PHCompositeNode *Node, const std::string &dnam, const int lyr ): 
  PHG4Detector(Node,dnam),
  radius(100*cm),
  length(100*cm),
  xpos(0),
  ypos(0),
  zpos(0),
  active(0),
  layer(lyr),
  blackhole(0)
{
}

//_______________________________________________________________
//_______________________________________________________________
bool PHG4CylinderDetector::IsInCylinderActive(const G4VPhysicalVolume * volume) const
{
  if (active && volume == cylinder_physi)
    {
      return true;
    }
  return false;
}

//_______________________________________________________________
bool PHG4CylinderDetector::IsInCylinder(const G4VPhysicalVolume * volume) const
{
  if (volume == cylinder_physi)
    {
      return true;
    }
  return false;
}


//_______________________________________________________________
void PHG4CylinderDetector::Construct( G4LogicalVolume* logicWorld )
{
  TrackerMaterial = G4Material::GetMaterial(material);

//   _region = new G4Region(GetName().c_str());
//   _region->SetRegionalSteppingAction(new PHG4CylinderRegionSteppingAction(this));


  if ( ! TrackerMaterial)
    {
      std::cout << "Error: Can not set material" << std::endl;
      exit(-1);
    }

  G4VisAttributes* siliconVis= new G4VisAttributes();
  if (IsBlackHole())
    {
      PHG4Utils::SetColour(siliconVis, "BlackHole");
      siliconVis->SetVisibility(false);
      siliconVis->SetForceSolid(false);
    }
  else
    {
      PHG4Utils::SetColour(siliconVis, material);
      siliconVis->SetVisibility(true);
      siliconVis->SetForceSolid(true);
    }

  // determine length of cylinder using PHENIX's rapidity coverage if flag is true
  cylinder_solid = new G4Tubs(G4String(GetName().c_str()),
			      radius,
			      radius+TrackerThickness, 
			      length/2. ,0,twopi);

  cylinder_logic = new G4LogicalVolume(cylinder_solid, 
				       TrackerMaterial, 
				       G4String(GetName().c_str()),
				       0,0,0);
  cylinder_logic->SetVisAttributes(siliconVis);
  //  _region->AddRootLogicalVolume(cylinder_logic);
  cylinder_physi = new G4PVPlacement(0, G4ThreeVector(xpos,ypos,zpos), 
				     cylinder_logic, 
				     G4String(GetName().c_str()), 
				     logicWorld, 0, false, overlapcheck);

}
