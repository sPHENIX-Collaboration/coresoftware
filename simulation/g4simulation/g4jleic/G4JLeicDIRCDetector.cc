#include "G4JLeicDIRCDetector.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>  // for PHG4Detector

#include <Geant4/G4Box.hh>
#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4VisAttributes.hh>

#include <cmath>
#include <iostream>  // for operator<<, endl, bas...
#include <utility>   // for pair

class PHCompositeNode;

using namespace std;

G4JLeicDIRCDetector::G4JLeicDIRCDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *params, const std::string &dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_Params(params)
{
}

//_______________________________________________________________
//_______________________________________________________________
int G4JLeicDIRCDetector::IsInDIRC(G4VPhysicalVolume *volume) const
{
  set<G4VPhysicalVolume *>::const_iterator iter = m_PhysicalVolumesSet.find(volume);

  if (iter != m_PhysicalVolumesSet.end())
  {
    return 1;
  }

  return 0;
}

void G4JLeicDIRCDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  cout << "constructing DIRC" << endl;
  double cb_DIRC_bars_DZ = 340 * cm;
  double cb_DIRC_bars_DY = 42. * cm;
  double cb_DIRC_bars_DX = 1.7 * cm;
  double dR = 83.65 * cm;
  double myL = 2 * M_PI * dR;
  int NUM = myL / cb_DIRC_bars_DY;
  double cb_DIRC_bars_deltaphi = 2 * M_PI / NUM;
  string solidname = "cb_DIRC_bars_Solid";
  G4VSolid *solid = new G4Box(solidname, cb_DIRC_bars_DX / 2., cb_DIRC_bars_DY / 2., cb_DIRC_bars_DZ / 2.);
  string logicname = "cb_DIRC_bars_Logic";
  G4LogicalVolume *logical = new G4LogicalVolume(solid, G4Material::GetMaterial("Quartz"), logicname);
  G4VisAttributes *vis = new G4VisAttributes(G4Color(0., 1., 0., 1.0));
  vis->SetForceSolid(true);
  logical->SetVisAttributes(vis);
  for (int ia = 0; ia < NUM; ia++)
  {
    double phi = (ia * (cb_DIRC_bars_deltaphi));
    double x = -dR * cos(phi);
    double y = -dR * sin(phi);
    G4RotationMatrix rot;
    rot.rotateZ(cb_DIRC_bars_deltaphi * ia);
    string physname = "cb_DIRC_bars_Phys_" + ia;
    G4VPhysicalVolume *phy = new G4PVPlacement(G4Transform3D(rot, G4ThreeVector(x, y, -400)),
                                               logical, physname,
                                               logicWorld, ia, false, OverlapCheck());
    m_PhysicalVolumesSet.insert(phy);
  }
  return;
}

void G4JLeicDIRCDetector::Print(const std::string &what) const
{
  cout << "JLeic DIRC Detector:" << endl;
  if (what == "ALL" || what == "VOLUME")
  {
    cout << "Version 0.1" << endl;
  }
  return;
}
