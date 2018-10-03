#include "PHG4Prototype2InnerHcalDetector.h"

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
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4RotationMatrix.hh>

#include <boost/format.hpp>

#include <cmath>
#include <sstream>

using namespace std;

static const string scintimothername = "InnerHcalScintiMother";
static const string steelplatename = "InnerHcalSteelPlate";

PHG4Prototype2InnerHcalDetector::PHG4Prototype2InnerHcalDetector( PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam  ):
  PHG4Detector(Node, dnam),
  params(parameters),
  m_InnerHcalSteelPlate(nullptr),
  m_InnerHcalAssembly(nullptr),
  steel_plate_corner_upper_left(1157.5*mm,-151.44*mm),
  steel_plate_corner_upper_right(1308.5*mm,-286.96*mm), 
  steel_plate_corner_lower_right(1298.8*mm,-297.39*mm),
  steel_plate_corner_lower_left(1155.8*mm,-163.92*mm),

  scinti_u1_front_size(105.9*mm),
  scinti_u1_corner_upper_left(0*mm,0*mm),
  scinti_u1_corner_upper_right(198.1*mm,0*mm),
  scinti_u1_corner_lower_right(198.1*mm,-121.3*mm),
  scinti_u1_corner_lower_left(0*mm,-scinti_u1_front_size),

  scinti_u2_corner_upper_left(0*mm,0*mm),
  scinti_u2_corner_upper_right(198.1*mm,-15.4*mm),
  scinti_u2_corner_lower_right(198.1*mm,-141.5*mm),
  scinti_u2_corner_lower_left(0*mm,-110.59*mm),

  scinti_t9_distance_to_corner(26.44*mm),
  scinti_t9_front_size(140.3*mm),
  scinti_t9_corner_upper_left(0*mm,0*mm),
  scinti_t9_corner_upper_right(198.1*mm,-134.4*mm),
  scinti_t9_corner_lower_right(198.1*mm,-198.1*mm/tan(52.02/180.*M_PI)-scinti_t9_front_size),
  scinti_t9_corner_lower_left(0*mm,-scinti_t9_front_size),

  scinti_t10_front_size(149.2*mm),
  scinti_t10_corner_upper_left(0*mm,0*mm),
  scinti_t10_corner_upper_right(198.1*mm,-154.6*mm),
  scinti_t10_corner_lower_right(198.1*mm,-198.1*mm/tan(48.34/180.*M_PI)-scinti_t10_front_size),
  scinti_t10_corner_lower_left(0*mm,-scinti_t10_front_size),


  scinti_t11_front_size(144.3*mm),
  scinti_t11_corner_upper_left(0*mm,0*mm),
  scinti_t11_corner_upper_right(198.1*mm,-176.2*mm),
  scinti_t11_corner_lower_right(198.1*mm,-198.1*mm/tan(45.14/180.*M_PI)-scinti_t11_front_size),
  scinti_t11_corner_lower_left(0*mm,-scinti_t11_front_size),

  scinti_t12_front_size(186.6*mm),
  scinti_t12_corner_upper_left(0*mm,0*mm),
  scinti_t12_corner_upper_right(198.1*mm,-197.11*mm),
  scinti_t12_corner_lower_right(198.1*mm,-198.1*mm/tan(41.47/180.*M_PI)-scinti_t12_front_size),
  scinti_t12_corner_lower_left(0*mm,-scinti_t12_front_size),

  scinti_x(198.1),
  steel_z(901.7*mm),
  size_z(steel_z),
  scinti_tile_z(steel_z),
  scinti_tile_thickness(7*mm),
  scinti_box_smaller(0.02*mm), // blargh - off by 20 microns bc scinti tilt angle, need to revisit at some point
  gap_between_tiles(1*mm),
  scinti_gap(8.5*mm),
  deltaphi(2*M_PI/320.),
  volume_steel(NAN),
  volume_scintillator(NAN),
  n_scinti_plates(20),
  n_steel_plates(n_scinti_plates+1),
  m_ActiveFlag(params->get_int_param("active")),
  m_AbsorberActiveFlag(params->get_int_param("absorberactive")),
  layer(0)
{}

PHG4Prototype2InnerHcalDetector::~PHG4Prototype2InnerHcalDetector()
{
  delete m_InnerHcalAssembly;
}

//_______________________________________________________________
//_______________________________________________________________
int
PHG4Prototype2InnerHcalDetector::IsInPrototype2InnerHcal(G4VPhysicalVolume * volume) const
{
  G4LogicalVolume *logvol = volume->GetLogicalVolume();
  if (m_AbsorberActiveFlag && logvol == m_InnerHcalSteelPlate)
  {
    return -1;
  }
  if (m_ActiveFlag && m_ActiveVolumeSet.find(logvol) != m_ActiveVolumeSet.end())
  {
    return 1;
  }
  return 0;
}

G4LogicalVolume*
PHG4Prototype2InnerHcalDetector::ConstructSteelPlate(G4LogicalVolume* hcalenvelope)
{
  if (!m_InnerHcalSteelPlate)
    {
      G4VSolid* steel_plate;
      std::vector<G4TwoVector> vertexes;
      vertexes.push_back(steel_plate_corner_upper_left);
      vertexes.push_back(steel_plate_corner_upper_right);
      vertexes.push_back(steel_plate_corner_lower_right);
      vertexes.push_back(steel_plate_corner_lower_left);
      G4TwoVector zero(0, 0);
      steel_plate =  new G4ExtrudedSolid("InnerHcalSteelPlateSolid",
					 vertexes,
					 size_z  / 2.0,
					 zero, 1.0,
					 zero, 1.0);

      volume_steel = steel_plate->GetCubicVolume()*n_steel_plates;
      m_InnerHcalSteelPlate = new G4LogicalVolume(steel_plate,G4Material::GetMaterial("Steel_A36"),steelplatename, 0, 0, 0);
      G4VisAttributes visattchk;
      visattchk.SetVisibility(true);
      visattchk.SetForceSolid(false);
      visattchk.SetColour(G4Colour::Blue());
      m_InnerHcalSteelPlate->SetVisAttributes(visattchk);
    }
  return m_InnerHcalSteelPlate;
}

G4LogicalVolume*
PHG4Prototype2InnerHcalDetector::ConstructScintillatorBoxHiEta(G4LogicalVolume* hcalenvelope)
{
  int copynum = 9;
  G4VSolid* scintiboxsolid = new G4Box(scintimothername,scinti_x/2.,(scinti_gap-scinti_box_smaller)/2.,scinti_tile_z/2.);
  //  DisplayVolume(scintiboxsolid,hcalenvelope);

  G4LogicalVolume* scintiboxlogical = new G4LogicalVolume(scintiboxsolid,G4Material::GetMaterial("G4_AIR"),scintimothername, 0, 0, 0);
  G4VisAttributes hcalVisAtt;
  hcalVisAtt.SetVisibility(true);
  hcalVisAtt.SetForceSolid(false);
  hcalVisAtt.SetColour(G4Colour::Magenta());
  G4LogicalVolume *scintit9_logic = ConstructScintiTile9(hcalenvelope);
  scintit9_logic->SetVisAttributes(hcalVisAtt);

  double distance_to_corner = -size_z/2.+scinti_t9_distance_to_corner;
  G4RotationMatrix *Rot;  
  Rot = new G4RotationMatrix();  
  Rot->rotateX(90*deg);
  new G4PVPlacement(Rot,G4ThreeVector(-scinti_x/2.,0,distance_to_corner),scintit9_logic,(boost::format("InnerScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());

  hcalVisAtt.SetVisibility(true);
  hcalVisAtt.SetForceSolid(false);
  hcalVisAtt.SetColour(G4Colour::Blue());
  G4LogicalVolume *scintit10_logic = ConstructScintiTile10(hcalenvelope);
  scintit10_logic->SetVisAttributes(hcalVisAtt);

  distance_to_corner += scinti_t9_front_size + gap_between_tiles;
  Rot = new G4RotationMatrix();  
  Rot->rotateX(90*deg);
  copynum++;
  new G4PVPlacement(Rot,G4ThreeVector(-scinti_x/2.,0,distance_to_corner),scintit10_logic,(boost::format("InnerScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());

  hcalVisAtt.SetVisibility(true);
  hcalVisAtt.SetForceSolid(false);
  hcalVisAtt.SetColour(G4Colour::Yellow());
  G4LogicalVolume *scintit11_logic = ConstructScintiTile11(hcalenvelope);
  scintit11_logic->SetVisAttributes(hcalVisAtt);

  distance_to_corner += scinti_t10_front_size + gap_between_tiles;
  Rot = new G4RotationMatrix();  
  Rot->rotateX(90*deg);
  copynum++;
  new G4PVPlacement(Rot,G4ThreeVector(-scinti_x/2.,0,distance_to_corner),scintit11_logic,(boost::format("InnerScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());

  hcalVisAtt.SetVisibility(true);
  hcalVisAtt.SetForceSolid(false);
  hcalVisAtt.SetColour(G4Colour::Cyan());
  G4LogicalVolume *scintit12_logic = ConstructScintiTile12(hcalenvelope);
  scintit12_logic->SetVisAttributes(hcalVisAtt);

  distance_to_corner += scinti_t11_front_size + gap_between_tiles;
  Rot = new G4RotationMatrix();  
  Rot->rotateX(90*deg);
  copynum++;
  new G4PVPlacement(Rot,G4ThreeVector(-scinti_x/2.,0,distance_to_corner),scintit12_logic,(boost::format("InnerScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());

  //    DisplayVolume(scintiboxlogical,hcalenvelope);
  return scintiboxlogical;
}

G4LogicalVolume*
PHG4Prototype2InnerHcalDetector::ConstructScintillatorBox(G4LogicalVolume* hcalenvelope)
{ 
  int copynum = 0;
  G4VSolid* scintiboxsolid = new G4Box(scintimothername,scinti_x/2.,(scinti_gap-scinti_box_smaller)/2.,scinti_tile_z/2.);
  //    DisplayVolume(scintiboxsolid,hcalenvelope);
  G4LogicalVolume* scintiboxlogical = new G4LogicalVolume(scintiboxsolid,G4Material::GetMaterial("G4_AIR"),scintimothername, 0, 0, 0);
  G4VisAttributes hcalVisAtt;
  hcalVisAtt.SetVisibility(true);
  hcalVisAtt.SetForceSolid(false);
  hcalVisAtt.SetColour(G4Colour::Red());
  G4LogicalVolume *scintiu1_logic = ConstructScintiTileU1(hcalenvelope);
  scintiu1_logic->SetVisAttributes(hcalVisAtt);

  hcalVisAtt.SetVisibility(true);
  hcalVisAtt.SetForceSolid(false);
  hcalVisAtt.SetColour(G4Colour::Cyan());
  G4LogicalVolume *scintiu2_logic = ConstructScintiTileU2(hcalenvelope);
  scintiu2_logic->SetVisAttributes(hcalVisAtt);
  G4RotationMatrix *Rot;  
  Rot = new G4RotationMatrix();  
  Rot->rotateX(-90*deg);
  new G4PVPlacement(Rot,G4ThreeVector(-scinti_x/2.,0,-scinti_u1_front_size-gap_between_tiles/2.-gap_between_tiles),scintiu2_logic,(boost::format("InnerScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());

  Rot = new G4RotationMatrix();  
  Rot->rotateX(-90*deg);
  copynum++;
  new G4PVPlacement(Rot,G4ThreeVector(-scinti_x/2.,0,-gap_between_tiles/2.),scintiu1_logic,(boost::format("InnerScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());

  Rot = new G4RotationMatrix();  
  Rot->rotateX(90*deg);
  copynum++;
  new G4PVPlacement(Rot,G4ThreeVector(-scinti_x/2.,0,gap_between_tiles/2.),scintiu1_logic,(boost::format("InnerScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());

  Rot = new G4RotationMatrix();  
  Rot->rotateX(90*deg);
  copynum++;
  new G4PVPlacement(Rot,G4ThreeVector(-scinti_x/2.,0,scinti_u1_front_size+gap_between_tiles/2.+gap_between_tiles),scintiu2_logic,(boost::format("InnerScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());
  //  DisplayVolume(scintiboxlogical,hcalenvelope);
  return scintiboxlogical;
}

G4LogicalVolume*
PHG4Prototype2InnerHcalDetector::ConstructScintiTileU1(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(scinti_u1_corner_upper_left);
  vertexes.push_back(scinti_u1_corner_upper_right);
  vertexes.push_back(scinti_u1_corner_lower_right);
  vertexes.push_back(scinti_u1_corner_lower_left);
  G4TwoVector zero(0, 0);
  G4VSolid *scintiu1 =  new G4ExtrudedSolid("InnerHcalScintiU1",
					    vertexes,
					    scinti_tile_thickness  / 2.0,
					    zero, 1.0,
					    zero, 1.0);

  G4LogicalVolume *scintiu1_logic = new G4LogicalVolume(scintiu1,G4Material::GetMaterial("G4_POLYSTYRENE"),"InnerHcalScintiU1", nullptr, nullptr, nullptr);
  //     DisplayVolume(scintiu1,hcalenvelope);
  m_ActiveVolumeSet.insert(scintiu1_logic);
  return scintiu1_logic;
}

G4LogicalVolume*
PHG4Prototype2InnerHcalDetector::ConstructScintiTileU2(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(scinti_u2_corner_upper_left);
  vertexes.push_back(scinti_u2_corner_upper_right);
  vertexes.push_back(scinti_u2_corner_lower_right);
  vertexes.push_back(scinti_u2_corner_lower_left);
  G4TwoVector zero(0, 0);
  G4VSolid *scintiu2 =  new G4ExtrudedSolid("InnerHcalScintiU2",
					    vertexes,
					    scinti_tile_thickness  / 2.0,
					    zero, 1.0,
					    zero, 1.0);

  G4LogicalVolume *scintiu2_logic = new G4LogicalVolume(scintiu2,G4Material::GetMaterial("G4_POLYSTYRENE"),"InnerHcalScintiU2", nullptr, nullptr, nullptr);
  //   DisplayVolume(scintiu2,hcalenvelope);
  m_ActiveVolumeSet.insert(scintiu2_logic);
  return scintiu2_logic;
}

G4LogicalVolume*
PHG4Prototype2InnerHcalDetector::ConstructScintiTile9(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(scinti_t9_corner_upper_left);
  vertexes.push_back(scinti_t9_corner_upper_right);
  vertexes.push_back(scinti_t9_corner_lower_right);
  vertexes.push_back(scinti_t9_corner_lower_left);
  G4TwoVector zero(0, 0);
  G4VSolid *scintit9 =  new G4ExtrudedSolid("InnerHcalScintiT9",
					    vertexes,
					    scinti_tile_thickness  / 2.0,
					    zero, 1.0,
					    zero, 1.0);

  G4LogicalVolume *scintit9_logic = new G4LogicalVolume(scintit9,G4Material::GetMaterial("G4_POLYSTYRENE"),"InnerHcalScintiT9", nullptr, nullptr, nullptr);
  //     DisplayVolume(scintit9,hcalenvelope);
  m_ActiveVolumeSet.insert(scintit9_logic);
  return scintit9_logic;
}

G4LogicalVolume*
PHG4Prototype2InnerHcalDetector::ConstructScintiTile10(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(scinti_t10_corner_upper_left);
  vertexes.push_back(scinti_t10_corner_upper_right);
  vertexes.push_back(scinti_t10_corner_lower_right);
  vertexes.push_back(scinti_t10_corner_lower_left);
  G4TwoVector zero(0, 0);
  G4VSolid *scintit10 =  new G4ExtrudedSolid("InnerHcalScintiT10",
					    vertexes,
					    scinti_tile_thickness  / 2.0,
					    zero, 1.0,
					    zero, 1.0);

  G4LogicalVolume *scintit10_logic = new G4LogicalVolume(scintit10,G4Material::GetMaterial("G4_POLYSTYRENE"),"InnerHcalScintiT10", nullptr, nullptr, nullptr);
  //     DisplayVolume(scintit10,hcalenvelope);
  m_ActiveVolumeSet.insert(scintit10_logic);
  return scintit10_logic;
}

G4LogicalVolume*
PHG4Prototype2InnerHcalDetector::ConstructScintiTile11(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(scinti_t11_corner_upper_left);
  vertexes.push_back(scinti_t11_corner_upper_right);
  vertexes.push_back(scinti_t11_corner_lower_right);
  vertexes.push_back(scinti_t11_corner_lower_left);
  G4TwoVector zero(0, 0);
  G4VSolid *scintit11 =  new G4ExtrudedSolid("InnerHcalScintiT11",
					    vertexes,
					    scinti_tile_thickness  / 2.0,
					    zero, 1.0,
					    zero, 1.0);

  G4LogicalVolume *scintit11_logic = new G4LogicalVolume(scintit11,G4Material::GetMaterial("G4_POLYSTYRENE"),"InnerHcalScintiT11", nullptr, nullptr, nullptr);
  //     DisplayVolume(scintit11,hcalenvelope);
  m_ActiveVolumeSet.insert(scintit11_logic);
  return scintit11_logic;
}

G4LogicalVolume*
PHG4Prototype2InnerHcalDetector::ConstructScintiTile12(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(scinti_t12_corner_upper_left);
  vertexes.push_back(scinti_t12_corner_upper_right);
  vertexes.push_back(scinti_t12_corner_lower_right);
  vertexes.push_back(scinti_t12_corner_lower_left);
  G4TwoVector zero(0, 0);
  G4VSolid *scintit12 =  new G4ExtrudedSolid("InnerHcalScintiT12",
					    vertexes,
					    scinti_tile_thickness  / 2.0,
					    zero, 1.0,
					    zero, 1.0);

  G4LogicalVolume *scintit12_logic = new G4LogicalVolume(scintit12,G4Material::GetMaterial("G4_POLYSTYRENE"),"InnerHcalScintiT12", nullptr, nullptr, nullptr);
  //     DisplayVolume(scintit12,hcalenvelope);
  m_ActiveVolumeSet.insert(scintit12_logic);
  return scintit12_logic;
}

// Construct the envelope and the call the
// actual inner hcal construction
void
PHG4Prototype2InnerHcalDetector::Construct( G4LogicalVolume* logicWorld )
{
  G4ThreeVector g4vec(params->get_double_param("place_x")*cm,
                      params->get_double_param("place_y")*cm,
		      params->get_double_param("place_z")*cm);
  G4RotationMatrix Rot;
  Rot.rotateX(params->get_double_param("rot_x")*deg);
  Rot.rotateY(params->get_double_param("rot_y")*deg);
  Rot.rotateZ(params->get_double_param("rot_z")*deg);
  //  ConstructScintiTile9(logicWorld);
  //    ConstructScintillatorBoxHiEta(logicWorld);
  //ConstructScintillatorBox(logicWorld);
  //  return;
  m_InnerHcalAssembly = new G4AssemblyVolume();
  //ConstructSteelPlate(hcal_envelope_log);
  // return;
  ConstructInnerHcal(logicWorld);
  m_InnerHcalAssembly->MakeImprint(logicWorld,g4vec,&Rot,0,OverlapCheck());
// this is rather pathetic - there is no way to extract the name when a volume is added
// to the assembly. The only thing we can do is get an iterator over the placed volumes
// in the order in which they were placed. Since this code does not install the scintillators
// for the Al version, parsing the volume names to get the id does not work since it changes
// So now we loop over all volumes and store them in a map for fast lookup of the row
  int isteel = 0;
  int iscinti = 0;
  vector<G4VPhysicalVolume*>::iterator it = m_InnerHcalAssembly->GetVolumesIterator();
  for (unsigned int i=0; i<m_InnerHcalAssembly-> TotalImprintedVolumes();i++)
  {
    string volname = (*it)->GetName();
    if (volname.find(steelplatename) != string::npos)
    { 
      m_SteelPlateIdMap.insert(make_pair(volname,isteel));
      ++isteel;
    }
    else if (volname.find(scintimothername) != string::npos)
    {
      m_ScintillatorIdMap.insert(make_pair(volname,iscinti));
      ++iscinti;
    }
    ++it;
  }
// print out volume names and their assigned id
   // map<string,int>::const_iterator iter;
   // for (iter = m_SteelPlateIdMap.begin(); iter != m_SteelPlateIdMap.end(); ++iter)
   // {
   //   cout << iter->first << ", " << iter->second << endl;
   // }
   // for (iter = m_ScintillatorIdMap.begin(); iter != m_ScintillatorIdMap.end(); ++iter)
   // {
   //   cout << iter->first << ", " << iter->second << endl;
   // }
  return;
}

int
PHG4Prototype2InnerHcalDetector::ConstructInnerHcal(G4LogicalVolume* hcalenvelope)
{
  G4LogicalVolume* steel_plate = ConstructSteelPlate(hcalenvelope); // bottom steel plate
  G4LogicalVolume* scintibox = nullptr;
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
  double bottom_xmiddle_steel_tile = (steel_plate_corner_lower_right.x()+steel_plate_corner_lower_left.x())/2.;
  double bottom_ymiddle_steel_tile = (steel_plate_corner_lower_left.y()+steel_plate_corner_lower_right.y())/2.;
  double middlerad = sqrt(bottom_xmiddle_steel_tile*bottom_xmiddle_steel_tile + bottom_ymiddle_steel_tile * bottom_ymiddle_steel_tile);
  double philow = atan((bottom_ymiddle_steel_tile-scinti_gap/2.-0.87*mm)/bottom_xmiddle_steel_tile);
  double scintiangle = GetScintiAngle();
  for (int i = 0; i < n_steel_plates; i++)
    //      for (int i = 0; i < 2; i++)
    {
      name.str("");
      name << "InnerHcalSteel_" << i;
      G4RotationMatrix Rot;
      Rot.rotateZ(phi*rad);
      G4ThreeVector g4vec(0,0,0);
      m_InnerHcalAssembly->AddPlacedVolume(steel_plate,g4vec,&Rot);
      if (i > 0)
	{
	  double ypos = sin(phi+philow) * middlerad;
	  double xpos = cos(phi+philow) * middlerad;
	  name.str("");
	  name << "InnerHcalScintiBox_" << i;
	  G4RotationMatrix Rot1;
	  Rot1.rotateZ(scintiangle+phislat);
	  G4ThreeVector g4vecsc(xpos, ypos, 0);
	  m_InnerHcalAssembly->AddPlacedVolume(scintibox,g4vecsc,&Rot1);
	  phislat += deltaphi;
	}
      phi += deltaphi;
    }
  return 0;
}

// calculate the angle of the bottom scintillator. It is the angle of the top edge
// of the steel plate
double
PHG4Prototype2InnerHcalDetector::GetScintiAngle()
{
  double xlen = steel_plate_corner_upper_right.x() - steel_plate_corner_upper_left.x();
  double ylen = steel_plate_corner_upper_right.y() - steel_plate_corner_upper_left.y();
  double angle = atan(ylen/xlen);
  return angle;
}

int
PHG4Prototype2InnerHcalDetector::DisplayVolume(G4VSolid *volume,  G4LogicalVolume* logvol, G4RotationMatrix *rotm )
{
  G4LogicalVolume* checksolid = new G4LogicalVolume(volume, G4Material::GetMaterial("G4_POLYSTYRENE"), "DISPLAYLOGICAL", 0, 0, 0);
  DisplayVolume(checksolid, logvol, rotm);
  return 0;
}

int
PHG4Prototype2InnerHcalDetector::DisplayVolume(G4LogicalVolume *checksolid,  G4LogicalVolume* logvol, G4RotationMatrix *rotm )
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
  new G4PVPlacement(rotm, G4ThreeVector(0, 0, 0), checksolid, "DISPLAYVOL", logvol, 0, false, OverlapCheck());
  //  new G4PVPlacement(rotm, G4ThreeVector(0, -460.3, 0), checksolid, "DISPLAYVOL", logvol, 0, false, OverlapCheck());
  return 0;
}

void
PHG4Prototype2InnerHcalDetector::Print(const string &what) const
{
  cout << "Inner Hcal Detector:" << endl;
  if (what == "ALL" || what == "VOLUME")
    {
      cout << "Volume Steel: " << volume_steel/cm/cm/cm << " cm^3" << endl;
      cout << "Volume Scintillator: " << volume_scintillator/cm/cm/cm << " cm^3" << endl;
    }
  return;
}
int PHG4Prototype2InnerHcalDetector::get_scinti_row_id(const string &volname)
{
  int id=-9999;
  auto it = m_ScintillatorIdMap.find(volname);
  if (it != m_ScintillatorIdMap.end())
  {
    id = it->second;
  }
  else
  {
    cout << "unknown scintillator volume name: " << volname << endl;
  }

  return id;
}

int PHG4Prototype2InnerHcalDetector::get_steel_plate_id(const string &volname)
{
  int id=-9999;
  auto it = m_SteelPlateIdMap.find(volname);
  if (it != m_SteelPlateIdMap.end())
  {
    id = it->second;
  }
  else
  {
    cout << "unknown steel volume name: " << volname << endl;
  }
  return id;
}
