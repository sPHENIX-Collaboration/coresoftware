#include "PHG4BlockDetector.h"
#include "PHG4BlockRegionSteppingAction.h"

#include <g4main/PHG4RegionInformation.h>
#include <g4main/PHG4Utils.h>


#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <Geant4/G4Material.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>

#include <Geant4/G4VisAttributes.hh>
#include <Geant4/G4Colour.hh>

#include <sstream>

using namespace std;

//_______________________________________________________________
//note this inactive thickness is ~1.5% of a radiation length
PHG4BlockDetector::PHG4BlockDetector( PHCompositeNode *Node, const std::string &dnam,const int lyr  ):
  PHG4Detector(Node, dnam),
  center_in_x(0*cm),
  center_in_y(0*cm),
  center_in_z(-200*cm),
  _region(NULL),
  active(0),
  layer(lyr),
  blackhole(0)
{
  //set the default radii
  for (int i = 0; i < 3; i++)
    {
      dimension[i] = 100 * cm;
    }

}

//_______________________________________________________________
bool PHG4BlockDetector::IsInBlock(G4VPhysicalVolume * volume) const
{
  if (active && volume == block_physi)
  {
    return true;
  }
  return false;
}

bool PHG4BlockDetector::IsInBlockActive(G4VPhysicalVolume * volume) const
{
  if (volume == block_physi)
  {
    return true;
  }
  return false;
}

//_______________________________________________________________
void PHG4BlockDetector::Construct( G4LogicalVolume* logicWorld )
{

  TrackerMaterial = G4Material::GetMaterial(material.c_str());


  if ( ! TrackerMaterial )
    {
      std::cout << "Error: Can not set material" << std::endl;
      exit(-1);
    }

  block_solid = new G4Box(G4String(GetName().c_str()),
                          dimension[0]/2., dimension[1]/2., dimension[2]/2.);

  block_logic = new G4LogicalVolume(block_solid,
                                    TrackerMaterial,
                                    G4String(GetName().c_str()),
                                    0, 0, 0);
  G4VisAttributes* matVis = new G4VisAttributes();
  if (IsBlackHole())
    {
      PHG4Utils::SetColour(matVis, "BlackHole");
      matVis->SetVisibility(false);
      matVis->SetForceSolid(false);
    }
  else
    {
      PHG4Utils::SetColour(matVis, material);
      matVis->SetVisibility(true);
      matVis->SetForceSolid(true);
    }
  block_logic->SetVisAttributes(matVis);

  G4RotationMatrix *rotm  = new G4RotationMatrix();
  rotm->rotateZ(z_rot);
  block_physi = new G4PVPlacement(rotm, G4ThreeVector(center_in_x, center_in_y, center_in_z),
                                  block_logic,
                                  G4String(GetName().c_str()),
                                  logicWorld, 0, false, overlapcheck);

}
