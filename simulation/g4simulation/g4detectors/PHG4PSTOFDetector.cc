#include "PHG4PSTOFDetector.h"

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SystemOfUnits.hh>

#include <Geant4/G4Colour.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cmath>

using namespace std;

PHG4PSTOFDetector::PHG4PSTOFDetector(PHCompositeNode *Node, PHParametersContainer *params, const std::string &dnam)
  : PHG4Detector(Node, dnam)
  , paramscontainer(params)
{
  const PHParameters *par = paramscontainer->GetParameters(-1);
  IsActive = par->get_int_param("active");
  IsAbsorberActive = par->get_int_param("absorberactive");
  nmod = par->get_int_param("modules");
  nrows = par->get_int_param("rows");
}

//_______________________________________________________________
//_______________________________________________________________
int PHG4PSTOFDetector::IsInPSTOF(G4VPhysicalVolume *volume) const
{
  // G4AssemblyVolumes naming convention:
  std::map<G4VPhysicalVolume *, int>::const_iterator iter = active_phys_vols.find(volume);

  if (iter != active_phys_vols.end())
  {
    return iter->second;
  }

  return 0;
}

void PHG4PSTOFDetector::Construct(G4LogicalVolume *logicWorld)
{
  G4Material *Glass = G4Material::GetMaterial("G4_GLASS_PLATE");
  G4Box *pstof_box = new G4Box("pstof_box", 0.8 * cm, 6 * cm, 5 * cm);

  G4LogicalVolume *pstof_log_vol = new G4LogicalVolume(pstof_box, Glass, G4String("PSTOF_box"), 0, 0, 0);
  G4VisAttributes *pstofVisAtt = new G4VisAttributes();
  pstofVisAtt->SetVisibility(true);
  pstofVisAtt->SetForceSolid(true);
  pstofVisAtt->SetColour(G4Colour::Blue());
  pstof_log_vol->SetVisAttributes(pstofVisAtt);

  for (int irow = 0; irow < nrows; irow++)
  {
    int rowtype = irow % 2;  // odd or even row
    double phi = irow * (2.0 * M_PI / nrows);

    for (int imod = 0; imod < nmod; imod++)
    {
      const PHParameters *par = paramscontainer->GetParameters(imod);
      double z = NAN;
      double r = NAN;
      if (rowtype == 0)
      {
        z = par->get_double_param("z_mod_0") * cm;
        r = par->get_double_param("r_mod_0") * cm;
      }
      else
      {
        z = par->get_double_param("z_mod_1") * cm;
        r = par->get_double_param("r_mod_1") * cm;
      }

      // amount to rotate
      //double theta = atan2(z+z_offset[rowtype][itof],tof_radius+y_offset[rowtype][itof]);
      double theta = atan2(z, r);

      G4RotationMatrix *rotm = new G4RotationMatrix();
      rotm->rotateZ(-phi);
      rotm->rotateY(theta);

      double x = r * cos(phi);
      double y = r * sin(phi);

      int modnum = nmod * irow + imod;
      G4VPhysicalVolume *vol = new G4PVPlacement(rotm, G4ThreeVector(x, y, z), pstof_log_vol, "PSTOF", logicWorld, false, modnum, OverlapCheck());
      if (IsActive)
      {
	active_phys_vols[vol] = modnum;
//        active_phys_vols.insert(vol);
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
