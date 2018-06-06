#include "PHG4ConeDetector.h"
#include "PHG4ConeRegionSteppingAction.h"

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

using namespace std;

//_______________________________________________________________
//note this inactive thickness is ~1.5% of a radiation length
PHG4ConeDetector::PHG4ConeDetector( PHCompositeNode *Node, const std::string &dnam,const int lyr  ):
  PHG4Detector(Node, dnam),
  TrackerMaterial(nullptr),
  InactiveMaterial(nullptr),
  block_solid(nullptr),
  block_logic(nullptr),
  block_physi(nullptr),
  place_in_x(0*cm),
  place_in_y(0*cm),
  place_in_z(300*cm),
  rMin1(5*cm),
  rMax1(100*cm),
  rMin2(5*cm),
  rMax2(200*cm),
  dZ(100*cm),
  sPhi(0),
  dPhi(2*M_PI),
  z_rot(0),
  _region(nullptr),
  active(0),
  layer(lyr)
{}

//_______________________________________________________________
//_______________________________________________________________
bool PHG4ConeDetector::IsInConeActive(G4VPhysicalVolume * volume)
{
  if (volume == block_physi)
    {
      return true;
    }
  return false;
}

//_______________________________________________________________
bool PHG4ConeDetector::IsInConeInactive(G4VPhysicalVolume * volume)
{
  return false;
}

//_______________________________________________________________
void PHG4ConeDetector::Construct( G4LogicalVolume* logicWorld )
{

  TrackerMaterial = G4Material::GetMaterial(material.c_str());


  if ( ! TrackerMaterial )
    {
      std::cout << "Error: Can not set material" << std::endl;
      exit(-1);
    }

  block_solid = new G4Cons(G4String(GetName().c_str()),
			   rMin1, rMax1, rMin2, rMax2, dZ, sPhi, dPhi);

  block_logic = new G4LogicalVolume(block_solid,
                                    TrackerMaterial,
                                    G4String(GetName().c_str()),
                                    0, 0, 0);
  G4VisAttributes* matVis = new G4VisAttributes();
  PHG4Utils::SetColour(matVis,material);
  matVis->SetVisibility(true);
  matVis->SetForceSolid(true);
  block_logic->SetVisAttributes(matVis);

  G4RotationMatrix *rotm  = new G4RotationMatrix();
  rotm->rotateZ(z_rot);
  block_physi = new G4PVPlacement(rotm, G4ThreeVector(place_in_x, place_in_y, place_in_z),
                                  block_logic,
                                  G4String(GetName().c_str()),
                                  logicWorld, 0, false, OverlapCheck());

}
