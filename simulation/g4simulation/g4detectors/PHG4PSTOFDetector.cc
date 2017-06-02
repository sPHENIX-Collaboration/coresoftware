#include "PHG4PSTOFDetector.h"
#include "PHG4Parameters.h"

#include <g4main/PHG4Utils.h>


#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4IntersectionSolid.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4ExtrudedSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4TwoVector.hh>
#include <Geant4/G4Trap.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4UserLimits.hh>

#include <Geant4/G4VisAttributes.hh>
#include <Geant4/G4Colour.hh>

#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Object.h>
#include <CGAL/Circular_kernel_intersections.h>
#include <CGAL/Boolean_set_operations_2.h>

#include <boost/math/special_functions/sign.hpp>

#include <cmath>
#include <sstream>

using namespace std;

PHG4PSTOFDetector::PHG4PSTOFDetector( PHCompositeNode *Node, PHG4Parameters *parames, const std::string &dnam):
  PHG4Detector(Node, dnam),
  active(1)
{
}

PHG4PSTOFDetector::~PHG4PSTOFDetector()
{
}

//_______________________________________________________________
//_______________________________________________________________
//int PHG4PSTOFDetector::IsInPSTOF(G4VPhysicalVolume * volume) const
int PHG4PSTOFDetector::IsInPSTOF(G4LogicalVolume * volume) const
{
  // G4AssemblyVolumes naming convention:
  if (active)
  {
    if ( volume == active_volume )
    {
      return 1;
    }
  }

  return 0;
}

void PHG4PSTOFDetector::Construct( G4LogicalVolume* logicWorld )
{
  cout << "In PHG4PSTOFDetector::Construct" << endl;

  G4Material* Glass = G4Material::GetMaterial("G4_GLASS_PLATE");
  G4Box *pstof_box = new G4Box("pstof_box",0.8*cm,6*cm,5*cm);

  G4LogicalVolume* pstof_log_vol =  new G4LogicalVolume(pstof_box, Glass, G4String("PSTOF_box"), 0, 0, 0);
  G4VisAttributes* pstofVisAtt = new G4VisAttributes();
  pstofVisAtt->SetVisibility(true);
  pstofVisAtt->SetForceSolid(true);
  pstofVisAtt->SetColour(G4Colour::Blue());
  pstof_log_vol->SetVisAttributes(pstofVisAtt);

  active_volume = pstof_log_vol;  // save active volume for stepping action

  const int NMOD = 21;   // number of modules in one row
  const int NROWS = 56;  // number of rows (in azimuth)
  //const int NROWS = 4;  // number of rows (in azimuth)
  // each module is 8 ch, 1x10cm strips, with strips aligned in phi direction

  // Locations of modules
  //Double_t tof_radius = 85.0*cm; // cm
  Double_t z_mod[2][NMOD] = { { -109.3, -96.66, -84.42, -72.55, -61.07, -49.97, -39.25, -28.72, -18.76, -9.191, 
    0, 9.191, 18.76, 28.72, 39.25, 49.97, 61.07, 72.55, 84.42, 96.66, 109.3 },
           { -107.2, -94.66, -82.52, -70.75, -59.37, -48.47, -37.85, -27.72, -18.76, -9.191, 
             0, 9.191, 18.76, 27.72, 37.85, 48.47, 59.37, 70.75, 82.52, 94.66, 107.2 } };
  Double_t r_mod[2][NMOD] = { { 85.6, 85.6, 85.6, 85.6, 86, 86.5, 86.5, 86.5, 85.5, 83.6, 
    87.5, 83.6, 85.5, 86.5, 86.5, 86.5, 86, 85.6, 85.6, 85.6, 85.6 },
           { 85.3, 85.2, 84.9, 84.8, 85.1, 85, 85, 84.8, 83.8, 81.9, 
             85.8, 81.9, 83.8, 84.8, 85, 85, 85.1, 84.8, 84.9, 85.2, 85.3 } };

  for (int irow=0; irow<NROWS; irow++)
  {
    int rowtype = irow%2; // odd or even row
    Double_t phi = irow*(2.0*M_PI/NROWS);

    for (int imod=0; imod<NMOD; imod++)
    {
      //int itof = abs(imod-NMOD/2);
      //Double_t z = itof*(9.0+2.1*itof/(NMOD/2))*cm;  // center of mrpc module
      Double_t z = z_mod[rowtype][imod]*cm;
      Double_t r = r_mod[rowtype][imod]*cm;

      // amount to rotate
      //Double_t theta = atan2(z+z_offset[rowtype][itof],tof_radius+y_offset[rowtype][itof]);
      Double_t theta = atan2(z,r);

      /*
      if (imod<NMOD/2)
      {
        z = -z;
        theta = -theta;
      }
      */

      G4RotationMatrix *rotm = new G4RotationMatrix();
      rotm->rotateZ( -phi );
      rotm->rotateY( theta );

      Double_t x = r*cos(phi);
      Double_t y = r*sin(phi);

      /*
      cout << irow << "\t" << imod << "\t" << itof << "\t"
        << x/cm << "\t" << y/cm << "\t" << z/cm << "\t" << r/cm << "\t" << theta/deg << endl;
      */

      int modnum = NMOD*irow + imod;
      new G4PVPlacement(rotm,G4ThreeVector(x,y,z), pstof_log_vol, "PSTOF", logicWorld, false, modnum, overlapcheck);
    }
  }

  return;
}

/*
int PHG4PSTOFDetector::DisplayVolume(G4VSolid *volume,  G4LogicalVolume* logvol,G4RotationMatrix *rotm )
{
  static int i = 0;
  G4LogicalVolume* checksolid = new G4LogicalVolume(volume,G4Material::GetMaterial("G4_POLYSTYRENE"),"DISPLAYLOGICAL", 0, 0, 0);
  G4VisAttributes* visattchk = new G4VisAttributes();
  visattchk->SetVisibility(true);
  visattchk->SetForceSolid(false);
  switch(i)
    {
    case 0:
      visattchk->SetColour(G4Colour::Red());
      i++;
      break;
    case 1:
      visattchk->SetColour(G4Colour::Magenta());
      i++;
      break;
    case 2:
      visattchk->SetColour(G4Colour::Yellow());
      i++;
      break;
    case 3:
      visattchk->SetColour(G4Colour::Blue());
      i++;
      break;
    case 4:
      visattchk->SetColour(G4Colour::Cyan());
      i++;
      break;
    default:
      visattchk->SetColour(G4Colour::Green());
      i=0;
      break;
    }

  checksolid->SetVisAttributes(visattchk);
  new G4PVPlacement(rotm,G4ThreeVector(0,0,0),checksolid,"DISPLAYVOL",logvol, 0, false, overlapcheck);
  return 0;
}
*/

void PHG4PSTOFDetector::Print(const string &what) const
{
  cout << "PSTOF Detector:" << endl;
  if (what == "ALL" || what == "VOLUME")
    {
      cout << "Version 0.1" << endl;
    }
  return;
}
