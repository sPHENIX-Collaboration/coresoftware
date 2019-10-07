#include "G4JLeicVTXDetector.h"

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <g4main/PHG4Detector.h>  // for PHG4Detector

#include <Geant4/G4Box.hh>
#include <Geant4/G4Color.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>              // for G4Transform3D
#include <Geant4/G4VisAttributes.hh>

#include <cmath>
#include <iostream>  // for operator<<, endl, bas...
#include <utility>   // for pair

class G4VSolid;
class PHCompositeNode;

using namespace std;

G4JLeicVTXDetector::G4JLeicVTXDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParametersContainer *params, const std::string &dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_ParamsContainer(params)
{
  const PHParameters *par = m_ParamsContainer->GetParameters(-1);
  m_IsActiveFlag = par->get_int_param("active");
  m_IsAbsorberActiveFlag = par->get_int_param("absorberactive");
  m_Layers = par->get_int_param("layers");
}

//_______________________________________________________________
//_______________________________________________________________
int G4JLeicVTXDetector::IsInVTX(G4VPhysicalVolume *volume) const
{
  map<G4VPhysicalVolume *, int>::const_iterator iter = m_PhysicalVolumesMap.find(volume);

  if (iter != m_PhysicalVolumesMap.end())
  {
    return iter->second;
  }

  return 0;
}

void G4JLeicVTXDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  cout << "constructing VTX" << endl;
  for (int ilayer = 0; ilayer < m_Layers; ilayer++)
  {
    // get parameters for this layer
    const PHParameters *par = m_ParamsContainer->GetParameters(ilayer);
    double cb_VTX_ladder_DZ = par->get_double_param("Dz") * cm;
    double cb_VTX_ladder_DY = par->get_double_param("Dy") * cm;
    double cb_VTX_ladder_Thickness = par->get_double_param("Dx") * cm;
    double dR = par->get_double_param("Rin") * cm;
    double myL = 2 * M_PI * dR;
    int NUM = myL / cb_VTX_ladder_DY;

    for (int i = 0; i < 2; i++)
    {
      double LN = cb_VTX_ladder_DY * NUM;
      double LN1 = cb_VTX_ladder_DY * (NUM + 1 + i);
      if (LN / LN1 > 0.8) NUM = NUM + 1;
    }

    double cb_VTX_ladder_deltaphi = 2 * M_PI / NUM;
    string solidname = "cb_VTX_ladder_Solid_" + to_string(ilayer);
    G4VSolid *solid = new G4Box(solidname, cb_VTX_ladder_Thickness / 2., cb_VTX_ladder_DY / 2., cb_VTX_ladder_DZ / 2.);
    string logical_name = "cb_VTX_ladder_Logic_" + to_string(ilayer);
    G4LogicalVolume *logical = new G4LogicalVolume(solid, G4Material::GetMaterial("G4_Si"), logical_name);
    G4VisAttributes *attr_cb_VTX_ladder = nullptr;
    switch (ilayer)
    {
    case 0:
    case 1:
      attr_cb_VTX_ladder = new G4VisAttributes(G4Color(0.0, 0.2, 0.8, 2.0));
      break;
    case 2:
      attr_cb_VTX_ladder = new G4VisAttributes(G4Color(0.0, 0.2, 0.8, 0.7));
      break;
    default:
      attr_cb_VTX_ladder = new G4VisAttributes(G4Color(0.0 + 0.1 * double(ilayer - 3), 1., 1. - 0.1 * double(ilayer - 3), 1.0));
      break;
    }
    attr_cb_VTX_ladder->SetForceSolid(true);
    logical->SetVisAttributes(attr_cb_VTX_ladder);
    for (int ia = 0; ia < NUM; ia++)
    {
      double phi = (ia * (cb_VTX_ladder_deltaphi));
      double x = -dR * cos(phi);
      double y = -dR * sin(phi);
      G4RotationMatrix rot;
      rot.rotateZ(cb_VTX_ladder_deltaphi * ia);
      rot.rotateZ(-7. * deg);
      ostringstream physname;
      physname << "cb_VTX_ladder_Phys_" << ilayer << "_" << ia;
      G4VPhysicalVolume *phy = new G4PVPlacement(G4Transform3D(rot, G4ThreeVector(x, y, -400)),
                                                 logical, physname.str(),
                                                 logicWorld, 0, false, OverlapCheck());
      // layer starts at zero but needs to be positive to mark active volume
      m_PhysicalVolumesMap.insert(make_pair(phy, ilayer + 1));
    }
  }

  return;
}

void G4JLeicVTXDetector::Print(const std::string &what) const
{
  cout << "JLeic VTX Detector:" << endl;
  if (what == "ALL" || what == "VOLUME")
  {
    cout << "Version 0.1" << endl;
  }
  return;
}
