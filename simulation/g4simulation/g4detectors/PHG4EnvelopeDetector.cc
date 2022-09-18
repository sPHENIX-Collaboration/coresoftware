#include "PHG4EnvelopeDetector.h"

#include <g4main/PHG4Detector.h>  // for PHG4Detector

#include <Geant4/G4Colour.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4SystemOfUnits.hh>   // for mm, m
#include <Geant4/G4ThreeVector.hh>     // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>     // for G4Transform3D
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Types.hh>            // for G4double
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume
#include <Geant4/G4VisAttributes.hh>

#include <cmath>  // for M_PI
#include <iostream>

class G4Material;
class G4VSolid;
class PHCompositeNode;

using namespace std;

//___________________________________________________________________________________
PHG4EnvelopeDetector::PHG4EnvelopeDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, const std::string& dnam)
  : PHG4Detector(subsys, Node, dnam)
  , _placeInX(0.0 * mm)
  , _placeInY(0.0 * mm)
  , _placeInZ(2.0 * m)
  , _innerRadius(0.0 * m)
  , _outerRadius(1.0 * m)
  , _dZ(1.0 * mm)
  , _dZ_cyl(2.0 * m)
  , _sPhi(0)
  , _dPhi(2 * M_PI)
  , _materialCrystal("G4_PbWO4")
  , _active(1)
  , _layer(0)
  , _superdetector("NONE")
{
}

//_______________________________________________________________________________________
bool PHG4EnvelopeDetector::IsInEnvelope(G4VPhysicalVolume* volume) const
{
  if (volume->GetName().find("arbage") != string::npos)
  {
    return true;
  }
  else
  {
    return false;
  }
}

void PHG4EnvelopeDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  //***************
  //Backward Endcap
  //***************

  G4double placeInZ = _placeInZ;
  G4double placeInY = _placeInY;
  G4double placeInX = _placeInX;

  G4double rMin1 = _innerRadius;
  G4double rMax1 = _outerRadius;
  G4double rMin2 = rMin1;
  G4double rMax2 = rMax1;
  G4double dZ = _dZ;
  G4double sPhi = _sPhi;
  G4double dPhi = _dPhi;

  G4Material* material_crystal = GetDetectorMaterial("G4_PbWO4");

  G4VSolid* GarbageCollector_solid = new G4Cons("GarbageCollector_solid",
                                                rMin1, rMax1,
                                                rMin2, rMax2,
                                                dZ / 2.,
                                                sPhi, dPhi);

  G4LogicalVolume* GarbageCollector_logical = new G4LogicalVolume(GarbageCollector_solid, material_crystal, G4String("GarbageCollector"), nullptr, nullptr, nullptr);

  G4VisAttributes* ecalVisAtt = new G4VisAttributes();
  ecalVisAtt->SetVisibility(true);
  ecalVisAtt->SetForceSolid(false);
  ecalVisAtt->SetColour(G4Colour::Magenta());
  GarbageCollector_logical->SetVisAttributes(ecalVisAtt);

  G4RotationMatrix ecal_rotm;
  ecal_rotm.rotateX(0);
  ecal_rotm.rotateY(0);
  ecal_rotm.rotateZ(0);

  new G4PVPlacement(G4Transform3D(ecal_rotm, G4ThreeVector(placeInX, placeInY, placeInZ)),
                    GarbageCollector_logical,
                    "GarbageCollector",
                    logicWorld,
                    false,
                    false,
                    OverlapCheck());

  //**************
  //Forward Endcap
  //**************

  new G4PVPlacement(G4Transform3D(ecal_rotm, G4ThreeVector(placeInX, placeInY, -1 * placeInZ)),
                    GarbageCollector_logical,
                    "GarbageCollector_front",
                    logicWorld,
                    false,
                    false,
                    OverlapCheck());

  //*****************************
  //Cylinder Surrounding Detector
  //*****************************

  G4double cyl_place = 0 * mm;
  G4double r_min = _outerRadius + (1) * mm;
  G4double r_max = r_min + 1 * mm;
  G4double dZ_cyl = _dZ_cyl;

  G4VSolid* GarbageCollector_cyl_solid = new G4Tubs("GarbageCollector_cyl_solid",
                                                    r_min,
                                                    r_max,
                                                    dZ_cyl,
                                                    sPhi,
                                                    dPhi);

  G4LogicalVolume* GarbageCollector_cyl_logical = new G4LogicalVolume(GarbageCollector_cyl_solid,
                                                                      material_crystal,
                                                                      G4String("GarbageCollector_cyl"),
                                                                      nullptr,
                                                                      nullptr,
                                                                      nullptr);

  GarbageCollector_cyl_logical->SetVisAttributes(ecalVisAtt);

  new G4PVPlacement(G4Transform3D(ecal_rotm, G4ThreeVector(0, 0, cyl_place)),
                    GarbageCollector_cyl_logical,
                    "GarbageCollector_cyl",
                    logicWorld,
                    false,
                    false,
                    OverlapCheck());
}
