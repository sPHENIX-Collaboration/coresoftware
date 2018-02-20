#include "PHG4Prototype2OuterHcalDetector.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4Colour.hh>
#include <Geant4/G4ExtrudedSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4TwoVector.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cmath>
#include <sstream>

using namespace std;

double scinti_box_smaller = 0.02*mm;

PHG4Prototype2OuterHcalDetector::PHG4Prototype2OuterHcalDetector( PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam  ):
  PHG4Detector(Node, dnam),
  params(parameters),
  outerhcalsteelplate(NULL),
  outerhcalassembly(NULL),
  steel_plate_corner_upper_left(1777.6*mm,-433.5*mm),
  steel_plate_corner_upper_right(2600.4*mm,-417.4*mm), 
  steel_plate_corner_lower_right(2601.2*mm,-459.8*mm),
  steel_plate_corner_lower_left(1770.9*mm,-459.8*mm),

  scinti_u1_front_size(166.2*mm),
  scinti_u1_corner_upper_left(0*mm,0*mm),
  scinti_u1_corner_upper_right(828.9*mm,0*mm),
  scinti_u1_corner_lower_right(828.9*mm,-240.54*mm),
  scinti_u1_corner_lower_left(0*mm,-scinti_u1_front_size),

  scinti_u2_corner_upper_left(0*mm,0*mm),
  scinti_u2_corner_upper_right(828.9*mm,-74.3*mm),
  scinti_u2_corner_lower_right(828.9*mm,-320.44*mm),
  scinti_u2_corner_lower_left(0*mm,-171.0*mm),

  scinti_t9_distance_to_corner(0.86*mm),
  scinti_t9_front_size(241.5*mm),
  scinti_t9_corner_upper_left(0*mm,0*mm),
  scinti_t9_corner_upper_right(697.4*mm,-552.2*mm),
  scinti_t9_corner_lower_right(697.4*mm,-697.4*mm/tan(47.94/180.*M_PI)-scinti_t9_front_size),
  scinti_t9_corner_lower_left(0*mm,-scinti_t9_front_size),

  scinti_t10_front_size(241.4*mm),
  scinti_t10_corner_upper_left(0*mm,0*mm),
  scinti_t10_corner_upper_right(697.4*mm,-629.3*mm),
  scinti_t10_corner_lower_right(697.4*mm,-697.4*mm/tan(44.2/180.*M_PI)-scinti_t10_front_size),
  scinti_t10_corner_lower_left(0*mm,-scinti_t10_front_size),

  scinti_t11_front_size(241.4*mm),
  scinti_t11_corner_upper_left(0*mm,0*mm),
  scinti_t11_corner_upper_right(697.4*mm,-717.1*mm),
  scinti_t11_corner_lower_right(697.4*mm,-697.4*mm/tan(42.47/180.*M_PI)-scinti_t11_front_size),
  scinti_t11_corner_lower_left(0*mm,-scinti_t11_front_size),

  scinti_t12_front_size(312.7*mm),
  scinti_t12_corner_upper_left(0*mm,0*mm),
  scinti_t12_corner_upper_right(697.4*mm,-761.8*mm),
  scinti_t12_corner_lower_right(392.9*mm,-827.7),
  scinti_t12_corner_lower_left(0*mm,-scinti_t12_front_size),

  scinti_x(828.9),
  scinti_x_hi_eta(697.4*mm+121.09*mm),
  steel_z(1600.*mm),
  size_z(steel_z),
  scinti_tile_z(steel_z),
  scinti_tile_thickness(7*mm),
scinti_box_smaller(0.02*mm), // blargh - off by 20 microns bc scinti tilt angle, need to revisit at some point
  gap_between_tiles(1*mm),
  scinti_gap(8.5*mm),
  tilt_angle(12*deg),
  deltaphi(2*M_PI/320.),
  volume_steel(NAN),
  volume_scintillator(NAN),
  n_scinti_plates(20),
  n_steel_plates(n_scinti_plates+1),
  active(params->get_int_param("active")),
  absorberactive(params->get_int_param("absorberactive")),
  layer(0)
{}

//_______________________________________________________________
//_______________________________________________________________
int
PHG4Prototype2OuterHcalDetector::IsInPrototype2OuterHcal(G4VPhysicalVolume * volume) const
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
      if (volume->GetName().find("OuterHcalSteelPlate") != string::npos)
	{
	  return -1;
	}
    }
  if (active)
    {
      if (volume->GetName().find("OuterScinti") != string::npos)
	{
	  return 1;
	}
    }
  return 0;
}

G4LogicalVolume*
PHG4Prototype2OuterHcalDetector::ConstructSteelPlate(G4LogicalVolume* hcalenvelope)
{
  if (!outerhcalsteelplate)
    {
      G4VSolid* steel_plate;
      std::vector<G4TwoVector> vertexes;
      vertexes.push_back(steel_plate_corner_upper_left);
      vertexes.push_back(steel_plate_corner_upper_right);
      vertexes.push_back(steel_plate_corner_lower_right);
      vertexes.push_back(steel_plate_corner_lower_left);
      G4TwoVector zero(0, 0);
      steel_plate =  new G4ExtrudedSolid("OuterHcalSteelPlateSolid",
					 vertexes,
					 size_z  / 2.0,
					 zero, 1.0,
					 zero, 1.0);

      volume_steel = steel_plate->GetCubicVolume()*n_steel_plates;
      outerhcalsteelplate = new G4LogicalVolume(steel_plate,G4Material::GetMaterial("Steel_A36"),"OuterHcalSteelPlate", 0, 0, 0);
      G4VisAttributes* visattchk = new G4VisAttributes();
      visattchk->SetVisibility(true);
      visattchk->SetForceSolid(false);
      visattchk->SetColour(G4Colour::Blue());
      outerhcalsteelplate->SetVisAttributes(visattchk);
    }
  return outerhcalsteelplate;
}

G4LogicalVolume*
PHG4Prototype2OuterHcalDetector::ConstructScintillatorBox(G4LogicalVolume* hcalenvelope)
{ 
  G4VSolid* scintiboxsolid = new G4Box("OuterHcalScintiMother",scinti_x/2.,(scinti_gap-scinti_box_smaller)/2.,scinti_tile_z/2.);
  //  DisplayVolume(scintiboxsolid,hcalenvelope);
  G4LogicalVolume* scintiboxlogical = new G4LogicalVolume(scintiboxsolid,G4Material::GetMaterial("G4_AIR"),G4String("OuterHcalScintiMother"), 0, 0, 0);
  G4VisAttributes* hcalVisAtt = new G4VisAttributes();
  hcalVisAtt->SetVisibility(true);
  hcalVisAtt->SetForceSolid(false);
  hcalVisAtt->SetColour(G4Colour::Red());
  G4LogicalVolume *scintiu1_logic = ConstructScintiTileU1(hcalenvelope);
  scintiu1_logic->SetVisAttributes(hcalVisAtt);

  hcalVisAtt = new G4VisAttributes();
  hcalVisAtt->SetVisibility(true);
  hcalVisAtt->SetForceSolid(false);
  hcalVisAtt->SetColour(G4Colour::Cyan());
  G4LogicalVolume *scintiu2_logic = ConstructScintiTileU2(hcalenvelope);
  scintiu2_logic->SetVisAttributes(hcalVisAtt);
  G4RotationMatrix *Rot;  
  Rot = new G4RotationMatrix();  
  Rot->rotateX(-90*deg);
  new G4PVPlacement(Rot,G4ThreeVector(-scinti_x/2.,0,-scinti_u1_front_size-gap_between_tiles/2.-gap_between_tiles),scintiu2_logic,"OuterScinti_0", scintiboxlogical, false, 0, overlapcheck);

  Rot = new G4RotationMatrix();  
  Rot->rotateX(-90*deg);
  new G4PVPlacement(Rot,G4ThreeVector(-scinti_x/2.,0,-gap_between_tiles/2.),scintiu1_logic,"OuterScinti_1", scintiboxlogical, false, 0, overlapcheck);

  Rot = new G4RotationMatrix();  
  Rot->rotateX(90*deg);
  new G4PVPlacement(Rot,G4ThreeVector(-scinti_x/2.,0,gap_between_tiles/2.),scintiu1_logic,"OuterScinti_2", scintiboxlogical, false, 0, overlapcheck);

  Rot = new G4RotationMatrix();  
  Rot->rotateX(90*deg);
  new G4PVPlacement(Rot,G4ThreeVector(-scinti_x/2.,0,scinti_u1_front_size+gap_between_tiles/2.+gap_between_tiles),scintiu2_logic,"OuterScinti_3", scintiboxlogical, false, 0, overlapcheck);


  return scintiboxlogical;
}

G4LogicalVolume*
PHG4Prototype2OuterHcalDetector::ConstructScintiTileU1(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(scinti_u1_corner_upper_left);
  vertexes.push_back(scinti_u1_corner_upper_right);
  vertexes.push_back(scinti_u1_corner_lower_right);
  vertexes.push_back(scinti_u1_corner_lower_left);
  G4TwoVector zero(0, 0);
  G4VSolid *scintiu1 =  new G4ExtrudedSolid("OuterHcalScintiU1",
					    vertexes,
					    scinti_tile_thickness  / 2.0,
					    zero, 1.0,
					    zero, 1.0);

  G4LogicalVolume *scintiu1_logic = new G4LogicalVolume(scintiu1,G4Material::GetMaterial("G4_POLYSTYRENE"),"OuterHcalScintiU1", NULL, NULL, NULL);
  //   DisplayVolume(scintiu1,hcalenvelope);
  return scintiu1_logic;
}

G4LogicalVolume*
PHG4Prototype2OuterHcalDetector::ConstructScintiTileU2(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(scinti_u2_corner_upper_left);
  vertexes.push_back(scinti_u2_corner_upper_right);
  vertexes.push_back(scinti_u2_corner_lower_right);
  vertexes.push_back(scinti_u2_corner_lower_left);
  G4TwoVector zero(0, 0);
  G4VSolid *scintiu2 =  new G4ExtrudedSolid("OuterHcalScintiU2",
					    vertexes,
					    scinti_tile_thickness  / 2.0,
					    zero, 1.0,
					    zero, 1.0);

  G4LogicalVolume *scintiu2_logic = new G4LogicalVolume(scintiu2,G4Material::GetMaterial("G4_POLYSTYRENE"),"OuterHcalScintiU2", NULL, NULL, NULL);
  //   DisplayVolume(scintiu2,hcalenvelope);
  return scintiu2_logic;
}

G4LogicalVolume*
PHG4Prototype2OuterHcalDetector::ConstructScintillatorBoxHiEta(G4LogicalVolume* hcalenvelope)
{ 
  G4VSolid* scintiboxsolid = new G4Box("OuterHcalScintiMother",scinti_x/2.,(scinti_gap-scinti_box_smaller)/2.,scinti_tile_z/2.);
  //  DisplayVolume(scintiboxsolid,hcalenvelope);
  G4LogicalVolume* scintiboxlogical = new G4LogicalVolume(scintiboxsolid,G4Material::GetMaterial("G4_AIR"),G4String("OuterHcalScintiMother"), 0, 0, 0);

  G4VisAttributes* hcalVisAtt = new G4VisAttributes();
  hcalVisAtt->SetVisibility(true);
  hcalVisAtt->SetForceSolid(false);
  hcalVisAtt->SetColour(G4Colour::Magenta());
  G4LogicalVolume *scintit9_logic = ConstructScintiTile9(hcalenvelope);
  scintit9_logic->SetVisAttributes(hcalVisAtt);

  double distance_to_corner = -size_z/2.+scinti_t9_distance_to_corner;
  G4RotationMatrix *Rot;  
  Rot = new G4RotationMatrix();  
  Rot->rotateX(90*deg);
  new G4PVPlacement(Rot,G4ThreeVector(-scinti_x_hi_eta/2.,0,distance_to_corner),scintit9_logic,"OuterScinti_9", scintiboxlogical, false, 0, overlapcheck);

  hcalVisAtt = new G4VisAttributes();
  hcalVisAtt->SetVisibility(true);
  hcalVisAtt->SetForceSolid(false);
  hcalVisAtt->SetColour(G4Colour::Blue());
  G4LogicalVolume *scintit10_logic = ConstructScintiTile10(hcalenvelope);
  scintit10_logic->SetVisAttributes(hcalVisAtt);

  distance_to_corner += scinti_t9_front_size + gap_between_tiles;
  Rot = new G4RotationMatrix();  
  Rot->rotateX(90*deg);
  new G4PVPlacement(Rot,G4ThreeVector(-scinti_x_hi_eta/2.,0,distance_to_corner),scintit10_logic,"OuterScinti_10", scintiboxlogical, false, 0, overlapcheck);

  hcalVisAtt = new G4VisAttributes();
  hcalVisAtt->SetVisibility(true);
  hcalVisAtt->SetForceSolid(false);
  hcalVisAtt->SetColour(G4Colour::Yellow());
  G4LogicalVolume *scintit11_logic = ConstructScintiTile11(hcalenvelope);
  scintit11_logic->SetVisAttributes(hcalVisAtt);

  distance_to_corner += scinti_t10_front_size + gap_between_tiles;
  Rot = new G4RotationMatrix();  
  Rot->rotateX(90*deg);
  new G4PVPlacement(Rot,G4ThreeVector(-scinti_x_hi_eta/2.,0,distance_to_corner),scintit11_logic,"OuterScinti_11", scintiboxlogical, false, 0, overlapcheck);

  hcalVisAtt = new G4VisAttributes();
  hcalVisAtt->SetVisibility(true);
  hcalVisAtt->SetForceSolid(false);
  hcalVisAtt->SetColour(G4Colour::Cyan());
  G4LogicalVolume *scintit12_logic = ConstructScintiTile12(hcalenvelope);
  scintit12_logic->SetVisAttributes(hcalVisAtt);

  distance_to_corner += scinti_t11_front_size + gap_between_tiles;
  Rot = new G4RotationMatrix();  
  Rot->rotateX(90*deg);
  new G4PVPlacement(Rot,G4ThreeVector(-scinti_x_hi_eta/2.,0,distance_to_corner),scintit12_logic,"OuterScinti_12", scintiboxlogical, false, 0, overlapcheck);
  //DisplayVolume(scintiboxlogical,hcalenvelope);
  return scintiboxlogical;
}
G4LogicalVolume*
PHG4Prototype2OuterHcalDetector::ConstructScintiTile9(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(scinti_t9_corner_upper_left);
  vertexes.push_back(scinti_t9_corner_upper_right);
  vertexes.push_back(scinti_t9_corner_lower_right);
  vertexes.push_back(scinti_t9_corner_lower_left);
  G4TwoVector zero(0, 0);
  G4VSolid *scintit9 =  new G4ExtrudedSolid("OuterHcalScintiT9",
					    vertexes,
					    scinti_tile_thickness  / 2.0,
					    zero, 1.0,
					    zero, 1.0);

  G4LogicalVolume *scintit9_logic = new G4LogicalVolume(scintit9,G4Material::GetMaterial("G4_POLYSTYRENE"),"OuterHcalScintiT9", NULL, NULL, NULL);
  //     DisplayVolume(scintit9,hcalenvelope);
  return scintit9_logic;
}

G4LogicalVolume*
PHG4Prototype2OuterHcalDetector::ConstructScintiTile10(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(scinti_t10_corner_upper_left);
  vertexes.push_back(scinti_t10_corner_upper_right);
  vertexes.push_back(scinti_t10_corner_lower_right);
  vertexes.push_back(scinti_t10_corner_lower_left);
  G4TwoVector zero(0, 0);
  G4VSolid *scintit10 =  new G4ExtrudedSolid("OuterHcalScintiT10",
					    vertexes,
					    scinti_tile_thickness  / 2.0,
					    zero, 1.0,
					    zero, 1.0);

  G4LogicalVolume *scintit10_logic = new G4LogicalVolume(scintit10,G4Material::GetMaterial("G4_POLYSTYRENE"),"OuterHcalScintiT10", NULL, NULL, NULL);
  //     DisplayVolume(scintit10,hcalenvelope);
  return scintit10_logic;
}

G4LogicalVolume*
PHG4Prototype2OuterHcalDetector::ConstructScintiTile11(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(scinti_t11_corner_upper_left);
  vertexes.push_back(scinti_t11_corner_upper_right);
  vertexes.push_back(scinti_t11_corner_lower_right);
  vertexes.push_back(scinti_t11_corner_lower_left);
  G4TwoVector zero(0, 0);
  G4VSolid *scintit11 =  new G4ExtrudedSolid("OuterHcalScintiT11",
					    vertexes,
					    scinti_tile_thickness  / 2.0,
					    zero, 1.0,
					    zero, 1.0);

  G4LogicalVolume *scintit11_logic = new G4LogicalVolume(scintit11,G4Material::GetMaterial("G4_POLYSTYRENE"),"OuterHcalScintiT11", NULL, NULL, NULL);
  //     DisplayVolume(scintit11,hcalenvelope);
  return scintit11_logic;
}

G4LogicalVolume*
PHG4Prototype2OuterHcalDetector::ConstructScintiTile12(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(scinti_t12_corner_upper_left);
  vertexes.push_back(scinti_t12_corner_upper_right);
  vertexes.push_back(scinti_t12_corner_lower_right);
  vertexes.push_back(scinti_t12_corner_lower_left);
  G4TwoVector zero(0, 0);
  G4VSolid *scintit12 =  new G4ExtrudedSolid("OuterHcalScintiT12",
					    vertexes,
					    scinti_tile_thickness  / 2.0,
					    zero, 1.0,
					    zero, 1.0);

  G4LogicalVolume *scintit12_logic = new G4LogicalVolume(scintit12,G4Material::GetMaterial("G4_POLYSTYRENE"),"OuterHcalScintiT12", NULL, NULL, NULL);
  //     DisplayVolume(scintit12,hcalenvelope);
  return scintit12_logic;
}

// Construct the envelope and the call the
// actual outer hcal construction
void
PHG4Prototype2OuterHcalDetector::Construct( G4LogicalVolume* logicWorld )
{
  G4ThreeVector g4vec(params->get_double_param("place_x")*cm,
                      params->get_double_param("place_y")*cm,
		      params->get_double_param("place_z")*cm);
  G4RotationMatrix *Rot = new G4RotationMatrix();
  Rot->rotateX(params->get_double_param("rot_x")*deg);
  Rot->rotateY(params->get_double_param("rot_y")*deg);
  Rot->rotateZ(params->get_double_param("rot_z")*deg);
  //  ConstructScintillatorBoxHiEta(logicWorld);
  outerhcalassembly = new G4AssemblyVolume();
  //ConstructSteelPlate(hcal_envelope_log);
  // return;
  ConstructOuterHcal(logicWorld);
  outerhcalassembly->MakeImprint(logicWorld,g4vec,Rot,0,overlapcheck);
  return;
}

int
PHG4Prototype2OuterHcalDetector::ConstructOuterHcal(G4LogicalVolume* hcalenvelope)
{
  G4LogicalVolume* steel_plate = ConstructSteelPlate(hcalenvelope); // bottom steel plate
  G4LogicalVolume* scintibox = NULL;
  if (params->get_int_param("hi_eta"))
    {
      scintibox = ConstructScintillatorBoxHiEta(hcalenvelope);
    }
  else
    {
      scintibox = ConstructScintillatorBox(hcalenvelope);
    }
  double phi = 0.;
  double phislat = 0.;
  ostringstream name;
  // the coordinate of the center of the bottom of the bottom steel plate
  // to get the radius of the circle which is the center of the scintillator box
  double bottom_xmiddle_steel_tile = (steel_plate_corner_lower_right.x()-steel_plate_corner_lower_left.x())/2.+steel_plate_corner_lower_left.x();
  double bottom_ymiddle_steel_tile = steel_plate_corner_lower_right.y();
  double middlerad = sqrt(bottom_xmiddle_steel_tile*bottom_xmiddle_steel_tile + bottom_ymiddle_steel_tile * bottom_ymiddle_steel_tile);
  double philow = atan((bottom_ymiddle_steel_tile-scinti_gap/2.)/bottom_xmiddle_steel_tile);
  double scintiangle = GetScintiAngle();
  for (int i = 0; i < n_steel_plates; i++)
    //      for (int i = 0; i < 2; i++)
    {
      name.str("");
      name << "OuterHcalSteel_" << i;
      G4RotationMatrix *Rot = new G4RotationMatrix();
      Rot->rotateZ(phi*rad);
      G4ThreeVector g4vec(0,0,0);
      outerhcalassembly->AddPlacedVolume(steel_plate,g4vec,Rot);
      if (i > 0)
	{
	  double ypos = sin(phi+philow) * middlerad;
	  double xpos = cos(phi+philow) * middlerad;
	  // the center of the scintillator is not the center of the inner hcal
	  // but depends on the tilt angle. Therefore we need to shift
	  // the center from the mid point
	  ypos += sin((-tilt_angle)/rad - phi);
	  xpos -= cos((-tilt_angle)/rad - phi);
	  name.str("");
	  name << "OuterHcalScintiBox_" << i;
	  Rot = new G4RotationMatrix();
	  Rot->rotateZ(scintiangle+phislat);
	  G4ThreeVector g4vec(xpos, ypos, 0);

	  outerhcalassembly->AddPlacedVolume(scintibox,g4vec,Rot);
	  phislat += deltaphi;
	}
      phi += deltaphi;
    }
  return 0;
}

// calculate the angle of the bottom scintillator. It is the angle of the top edge
// of the steel plate
double
PHG4Prototype2OuterHcalDetector::GetScintiAngle()
{
  double xlen = steel_plate_corner_upper_right.x() - steel_plate_corner_upper_left.x();
  double ylen = steel_plate_corner_upper_right.y() - steel_plate_corner_upper_left.y();
  double angle =  atan(ylen/xlen);
  return angle;
}

int
PHG4Prototype2OuterHcalDetector::DisplayVolume(G4VSolid *volume,  G4LogicalVolume* logvol, G4RotationMatrix *rotm )
{
  G4LogicalVolume* checksolid = new G4LogicalVolume(volume, G4Material::GetMaterial("G4_POLYSTYRENE"), "DISPLAYLOGICAL", 0, 0, 0);
  DisplayVolume(checksolid,logvol,rotm);
  return 0;
}

int
PHG4Prototype2OuterHcalDetector::DisplayVolume(G4LogicalVolume *checksolid,  G4LogicalVolume* logvol, G4RotationMatrix *rotm )
{
  static int i = 0;
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
PHG4Prototype2OuterHcalDetector::Print(const string &what) const
{
  cout << "Outer Hcal Detector:" << endl;
  if (what == "ALL" || what == "VOLUME")
    {
      cout << "Volume Steel: " << volume_steel/cm/cm/cm << " cm^3" << endl;
      cout << "Volume Scintillator: " << volume_scintillator/cm/cm/cm << " cm^3" << endl;
    }
  return;
}

