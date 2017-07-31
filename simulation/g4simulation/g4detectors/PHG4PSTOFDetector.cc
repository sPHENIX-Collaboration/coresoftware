#include "PHG4PSTOFDetector.h"
#include "PHG4Parameters.h"
#include "PHG4ParametersContainer.h"

#include <g4main/PHG4Utils.h>


#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <Geant4/G4Material.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4TwoVector.hh>
#include <Geant4/G4UserLimits.hh>

#include <Geant4/G4VisAttributes.hh>
#include <Geant4/G4Colour.hh>

#include <cmath>
#include <sstream>

using namespace std;

PHG4PSTOFDetector::PHG4PSTOFDetector( PHCompositeNode *Node, PHG4ParametersContainer *params, const std::string &dnam):
  PHG4Detector(Node, dnam)
  ,paramscontainer(params)
{
  const PHG4Parameters *par = paramscontainer->GetParameters(-1);
IsActive =  par->get_int_param("active");
IsAbsorberActive = par->get_int_param("absorberactive");
}

PHG4PSTOFDetector::~PHG4PSTOFDetector()
{
}

//_______________________________________________________________
//_______________________________________________________________
int PHG4PSTOFDetector::IsInPSTOF(G4VPhysicalVolume * volume) const
{
  // G4AssemblyVolumes naming convention:
  if ( active_phys_vols.find(volume) !=  active_phys_vols.end())
    {
      return 1;
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
  // each module is 8 ch, 1x10cm strips, with strips aligned in phi direction

  // Locations of modules (geometry version -1)
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
    double phi = irow*(2.0*M_PI/NROWS);

    for (int imod=0; imod<NMOD; imod++)
    {
      const PHG4Parameters *par = paramscontainer->GetParameters(imod);
      double z = NAN;
      double r = NAN;
      if (rowtype == 0)
      {
	z = par->get_double_param("z_mod_0")*cm;
        r = par->get_double_param("r_mod_0")*cm;
      }
      else
      {
	z = par->get_double_param("z_mod_1")*cm;
        r = par->get_double_param("r_mod_1")*cm;
      }

      //int itof = abs(imod-NMOD/2);
      //Double_t z = itof*(9.0+2.1*itof/(NMOD/2))*cm;  // center of mrpc module
      // Double_t z = z_mod[rowtype][imod]*cm;
      // Double_t r = r_mod[rowtype][imod]*cm;

      // amount to rotate
      //Double_t theta = atan2(z+z_offset[rowtype][itof],tof_radius+y_offset[rowtype][itof]);
      double theta = atan2(z,r);

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
      G4VPhysicalVolume *vol = new G4PVPlacement(rotm,G4ThreeVector(x,y,z), pstof_log_vol, "PSTOF", logicWorld, false, modnum, overlapcheck);
      if (IsActive)
      {
        active_phys_vols.insert(vol);
      }
    }
  }

  return;
}


void PHG4PSTOFDetector::Print(const std::string &what) const
{
  cout << "PSTOF Detector:" << endl;
  if (what == "ALL" || what == "VOLUME")
    {
      cout << "Version 0.1" << endl;
    }
  return;
}
