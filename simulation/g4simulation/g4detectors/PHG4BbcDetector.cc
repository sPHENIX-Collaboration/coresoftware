#include "PHG4BbcDetector.h"

#include <phparameter/PHParameters.h>
#include <phool/phool.h>

#include <g4main/PHG4Detector.h>                // for PHG4Detector

#include <Geant4/G4Polyhedra.hh>
#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>           // for G4RotationMatrix
#include <Geant4/G4String.hh>                   // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>              // for G4ThreeVector
#include <Geant4/G4VisAttributes.hh>

#include <cmath>
#include <iostream>                             // for operator<<, endl, bas...
#include <utility>                              // for pair

class PHCompositeNode;

using namespace std;

PHG4BbcDetector::PHG4BbcDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *params, const std::string &dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_Params(params)
{
  //IsActive = par->get_int_param("active");
  //IsAbsorberActive = par->get_int_param("absorberactive");
}

//_______________________________________________________________
//_______________________________________________________________
int PHG4BbcDetector::IsInBbc(G4VPhysicalVolume *volume) const
{
  //cout << PHWHERE << " volume " << volume->GetName() << endl;

  set<G4VPhysicalVolume *>::const_iterator iter = m_PhysicalVolumesSet.find(volume);
  if (iter != m_PhysicalVolumesSet.end())
  {
    return 1;
  }

  return 0;
}

void PHG4BbcDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  cout << PHWHERE << " Constructing BBC" << endl;

  G4Material *Quartz = G4Material::GetMaterial("G4_SILICON_DIOXIDE");

  const G4double z_bbcq[] = {-1.5, 1.5};
  const G4double rInner_bbcq[] = {0.,0.};
  const G4double rOuter_bbcq[] = {1.27*cm,1.27*cm};
  G4Polyhedra *bbcq = new G4Polyhedra("bbcq",0.,2.0*M_PI,6,2,z_bbcq,rInner_bbcq,rOuter_bbcq); // bbc quartz radiator

  G4LogicalVolume *bbcq_log_vol = new G4LogicalVolume(bbcq, Quartz, G4String("Bbc_quartz"));
  G4VisAttributes *bbcqVisAtt = new G4VisAttributes();
  bbcqVisAtt->SetVisibility(true);
  bbcqVisAtt->SetForceSolid(true);
  bbcqVisAtt->SetColour(G4Colour::Red());
  bbcq_log_vol->SetVisAttributes(bbcqVisAtt);

  //m_bbcz = m_Params->get_double_param("z") * cm;
  m_bbcz = 250.0*cm;

  // Locations of the 64 PMT tubes in an arm.
  // Should prob move this to a geometry object...
  G4float xpos[] = {
    2.45951, -2.45951, 4.91902, 0, -4.91902, 7.37854, 2.45951, -2.45951, 
    -7.37854, 9.83805, 4.91902, 0, -4.91902, -9.83805, 7.37854, 2.45951, 
    -2.45951, -7.37854, 9.83805, 4.91902, -4.91902, -9.83805, 12.2976, 7.37854, 
    -7.37854, -12.2976, 9.83805, -9.83805, 12.2976, 7.37854, -7.37854, -12.2976, 
    12.2976, 7.37854, -7.37854, -12.2976, 9.83805, -9.83805, 12.2976, 7.37854, 
    -7.37854, -12.2976, 9.83805, 4.91902, -4.91902, -9.83805, 7.37854, 2.45951, 
    -2.45951, -7.37854, 9.83805, 4.91902, 0, -4.91902, -9.83805, 7.37854, 
    2.45951, -2.45951, -7.37854, 4.91902, 0, -4.91902, 2.45951, -2.45951, 
  };

  G4float ypos[] = {
    12.78, 12.78, 11.36, 11.36, 11.36, 9.94, 9.94, 9.94, 
    9.94, 8.52, 8.52, 8.52, 8.52, 8.52, 7.1, 7.1, 
    7.1, 7.1, 5.68, 5.68, 5.68, 5.68, 4.26, 4.26, 
    4.26, 4.26, 2.84, 2.84, 1.42, 1.42, 1.42, 1.42, 
    -1.42, -1.42, -1.42, -1.42, -2.84, -2.84, -4.26, -4.26, 
    -4.26, -4.26, -5.68, -5.68, -5.68, -5.68, -7.1, -7.1, 
    -7.1, -7.1, -8.52, -8.52, -8.52, -8.52, -8.52, -9.94, 
    -9.94, -9.94, -9.94, -11.36, -11.36, -11.36, -12.78, -12.78, 
  };

  const int NPMT = 64;  // No. PMTs per arm
  for ( int iarm = 0; iarm<2; iarm ++ )
  {
    G4float side = 1.0;
    if ( iarm==0 ) side = -1.0;

    for (int itube = 0; itube < NPMT; itube++)
    {
      G4VPhysicalVolume *vol = new G4PVPlacement(0, G4ThreeVector(side*xpos[itube]*cm, ypos[itube]*cm, side*m_bbcz),
          bbcq_log_vol, "BBC", logicWorld, false, iarm*NPMT + itube, OverlapCheck());

      m_PhysicalVolumesSet.insert(vol);
    }
  }

  return;
}

void PHG4BbcDetector::Print(const std::string &what) const
{
  cout << "Bbc Detector:" << endl;
  if (what == "ALL" || what == "VOLUME")
  {
    cout << "Version 0.1" << endl;
  }
  return;
}
