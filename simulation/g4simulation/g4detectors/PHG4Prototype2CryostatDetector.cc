#include "PHG4Prototype2CryostatDetector.h"
#include "PHG4Parameters.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeomv3.h"

#include <g4main/PHG4Utils.h>


#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4IntersectionSolid.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4ExtrudedSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4TwoVector.hh>
#include <Geant4/G4Trap.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4UserLimits.hh>

#include <Geant4/G4VisAttributes.hh>
#include <Geant4/G4Colour.hh>

#include <cmath>
#include <sstream>

using namespace std;

PHG4Prototype2CryostatDetector::PHG4Prototype2CryostatDetector( PHCompositeNode *Node, PHG4Parameters *parameters, const std::string &dnam  ):
  PHG4Detector(Node, dnam),
  params(parameters),
  cryostatassembly(NULL),
  alu_z(609.6*mm),
  n_alu_plates(3),
  active(params->get_int_param("active")),
  absorberactive(params->get_int_param("absorberactive")),
  layer(0)
{
  for (int i=0; i<3; i++)
    {
      volume_alu[i] = NAN;
    }
  double y_up = alu_z/2.; // it is a square in z/y
  double y_down = -alu_z/2.;
  alu_plate_corner_upper_left[0] = G4TwoVector(1419.6*mm,y_up);
  alu_plate_corner_upper_right[0] = G4TwoVector(1419.6*mm+9.5*mm,y_up);
  alu_plate_corner_lower_right[0] = G4TwoVector(1419.6*mm+9.5*mm,y_down);
  alu_plate_corner_lower_left[0] = G4TwoVector(1419.6*mm,y_down);

  alu_plate_corner_upper_left[1] = G4TwoVector(1507.2*mm,y_up);
  alu_plate_corner_upper_right[1] = G4TwoVector(1507.2*mm+88.9*mm,y_up);
  alu_plate_corner_lower_right[1] = G4TwoVector(1507.2*mm+88.9*mm,y_down);
  alu_plate_corner_lower_left[1] = G4TwoVector(1507.2*mm,y_down);

  alu_plate_corner_upper_left[2] = G4TwoVector(1739.3*mm,y_up);
  alu_plate_corner_upper_right[2] = G4TwoVector(1739.3*mm+25.4*mm,y_up);
  alu_plate_corner_lower_right[2] = G4TwoVector(1739.3*mm+25.4*mm,y_down);
  alu_plate_corner_lower_left[2] = G4TwoVector(1739.3*mm,y_down);

}

//_______________________________________________________________
//_______________________________________________________________
int
PHG4Prototype2CryostatDetector::IsInPrototype2Cryostat(G4VPhysicalVolume * volume) const
{
  // G4AssemblyVolumes naming convention:
  //     av_WWW_impr_XXX_YYY_ZZZ

  // where:

  //     WWW - assembly volume instance number
  //     XXX - assembly volume imprint number
  //     YYY - the name of the placed logical volume
  //     ZZZ - the logical volume index inside the assembly volume
  // e.g. av_1_impr_82_HcalInnerScinti_11_pv_11
  // 82 the number of the scintillator mother volume
  // HcalInnerScinti_11: name of scintillator slat
  // 11: number of scintillator slat logical volume
  if (absorberactive)
    {
      if (volume->GetName().find("CryostatAluPlate") != string::npos)
	{
	  return -1;
	}
    }
  if (active)
    {
      if (volume->GetName().find("InnerScinti") != string::npos)
	{
	  return 1;
	}
    }
  return 0;
}

G4LogicalVolume*
PHG4Prototype2CryostatDetector::ConstructAluPlate(G4LogicalVolume* hcalenvelope, const int n)
{
  if (n < 0 || n >= n_alu_plates)
    {
      cout << "invalid alu plate number " << n << endl;
      exit(1);
    }
  ostringstream name;
  G4VSolid* alu_plate;
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(alu_plate_corner_upper_left[n]);
  vertexes.push_back(alu_plate_corner_upper_right[n]);
  vertexes.push_back(alu_plate_corner_lower_right[n]);
  vertexes.push_back(alu_plate_corner_lower_left[n]);
  G4TwoVector zero(0, 0);
  name.str("");
  name << "CryostatAluPlateSolid_" << n;
  alu_plate =  new G4ExtrudedSolid(name.str(),
				     vertexes,
				     alu_z  / 2.0,
				     zero, 1.0,
				     zero, 1.0);

  volume_alu[n] = alu_plate->GetCubicVolume()*n_alu_plates;
  name.str("");
  name << "CryostatAluPlate_" << n;
  G4LogicalVolume *cryostataluplate = new G4LogicalVolume(alu_plate,G4Material::GetMaterial("G4_AL"),"CryostatAluPlate", 0, 0, 0);
  G4VisAttributes* visattchk = new G4VisAttributes();
  visattchk->SetVisibility(true);
  visattchk->SetForceSolid(false);
  visattchk->SetColour(G4Colour::Grey());
  cryostataluplate->SetVisAttributes(visattchk);
  return cryostataluplate;
}


// Construct the envelope and the call the
// actual inner hcal construction
void
PHG4Prototype2CryostatDetector::Construct( G4LogicalVolume* logicWorld )
{
  G4ThreeVector g4vec(0,0,0);
  G4RotationMatrix *Rot = new G4RotationMatrix();
  Rot->rotateX(params->get_double_param("rot_x")*deg);
  Rot->rotateY(params->get_double_param("rot_y")*deg);
  Rot->rotateZ(params->get_double_param("rot_z")*deg);
  cryostatassembly = new G4AssemblyVolume();
  ConstructCryostat(logicWorld);
  cryostatassembly->MakeImprint(logicWorld,g4vec,Rot,0,overlapcheck);
  return;
}

int
PHG4Prototype2CryostatDetector::ConstructCryostat(G4LogicalVolume* hcalenvelope)
{
  ostringstream name;
  for (int i = 0; i < n_alu_plates; i++)
    //      for (int i = 0; i < 2; i++)
    {
      G4LogicalVolume* alu_plate = ConstructAluPlate(hcalenvelope,i); // bottom alu plate
      name.str("");
      name << "CryostatAlu_" << i;
      G4RotationMatrix *Rot = new G4RotationMatrix();
      G4ThreeVector g4vec(0,0,0);
      cryostatassembly->AddPlacedVolume(alu_plate,g4vec,Rot);
    }
  return 0;
}


int
PHG4Prototype2CryostatDetector::DisplayVolume(G4VSolid *volume,  G4LogicalVolume* logvol, G4RotationMatrix *rotm )
{
  static int i = 0;
  G4LogicalVolume* checksolid = new G4LogicalVolume(volume, G4Material::GetMaterial("G4_POLYSTYRENE"), "DISPLAYLOGICAL", 0, 0, 0);
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
      i = 0;
      break;
    }

  checksolid->SetVisAttributes(visattchk);
  new G4PVPlacement(rotm, G4ThreeVector(0, 0, 0), checksolid, "DISPLAYVOL", logvol, 0, false, overlapcheck);
  //  new G4PVPlacement(rotm, G4ThreeVector(0, -460.3, 0), checksolid, "DISPLAYVOL", logvol, 0, false, overlapcheck);
  return 0;
}


void
PHG4Prototype2CryostatDetector::Print(const string &what) const
{
  cout << "Cryostat:" << endl;
  if (what == "ALL" || what == "VOLUME")
    {
      for (int i=0; i<3; i++)
	{
	  cout << "Volume Alu " << i << ": " << volume_alu[i]/cm/cm/cm << " cm^3" << endl;
	}
    }
  return;
}
