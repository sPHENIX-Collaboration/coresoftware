#include "PHG4Prototype3InnerHcalDetector.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Units.h>
#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <TSystem.h>

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

PHG4Prototype3InnerHcalDetector::PHG4Prototype3InnerHcalDetector(PHCompositeNode *Node, PHParameters *parameters, const std::string& dnam)
  : PHG4Detector(Node, dnam)
  , m_params(parameters)
  , m_InnerHcalSteelPlate(nullptr)
  , m_InnerHcalAssembly(nullptr)
  , m_scintibox(nullptr)
  , m_SteelPlateCornerUpperLeft(1154.49 * mm, -189.06 * mm)
  , m_SteelPlateCornerUpperRight(1297.94 * mm, -349.22 * mm)
  , m_SteelPlateCornerLowerRight(1288.18 * mm, -357.8 * mm)
  , m_SteelPlateCornerLowerLeft(1157.3 * mm, -205.56 * mm)
//
  , m_ScintiTile9DistanceToCorner(26.44 * mm)
  , m_ScintiTile9FrontSize(140.3 * mm)
  , m_ScintiTile9CornerUpperLeft(0 * mm, 0 * mm)
  , m_ScintiTile9CornerUpperRight(191. * mm, -134.4 * 191. / 198.1 * mm)
  , m_ScintiTile9CornerLowerRight(191. * mm, -191. * mm / tan(52.02 / 180. * M_PI) - m_ScintiTile9FrontSize)
  , m_ScintiTile9CornerLowerLeft(0 * mm, -m_ScintiTile9FrontSize)
//
  , m_ScintiTile10FrontSize(149.2 * mm)
  , m_ScintiTile10CornerUpperLeft(0 * mm, 0 * mm)
  , m_ScintiTile10CornerUpperRight(191. * mm, -154.6 * 191. / 198.1 * mm)
  , m_ScintiTile10CornerLowerRight(191. * mm, -191. * mm / tan(48.34 / 180. * M_PI) - m_ScintiTile10FrontSize)
  , m_ScintiTile10CornerLowerLeft(0 * mm, -m_ScintiTile10FrontSize)
//
  , m_ScintiTile11FrontSize(144.3 * mm)
  , m_ScintiTile11CornerUpperLeft(0 * mm, 0 * mm)
  , m_ScintiTile11CornerUpperRight(191. * mm, -176.2 * 191. / 198.1 * mm)
  , m_ScintiTile11CornerLowerRight(191. * mm, -191. * mm / tan(45.14 / 180. * M_PI) - m_ScintiTile11FrontSize)
  , m_ScintiTile11CornerLowerLeft(0 * mm, -m_ScintiTile11FrontSize)
//
  , m_ScintiTile12FrontSize(186.6 * mm)
  , m_ScintiTile12CornerUpperLeft(0 * mm, 0 * mm)
  , m_ScintiTile12CornerUpperRight(191. * mm, -197.11 * 191. / 198.1 * mm)
  , m_ScintiTile12CornerLowerRight(191. * mm, -191. * mm / tan(41.47 / 180. * M_PI) - m_ScintiTile12FrontSize)
  , m_ScintiTile12CornerLowerLeft(0 * mm, -m_ScintiTile12FrontSize)
//
  , m_ScintiX(198.1 * mm)
  , m_SteelZ(901.7 * mm)
  , m_ScintiTileZ(m_SteelZ)
  , m_ScintiTileThickness(7 * mm)
  , m_GapBetweenTiles(1 * mm)
  , m_ScintiGap(10 * mm)
  , m_DeltaPhi(2 * M_PI / 256.)
  , m_VolumeSteel(NAN)
  , m_VolumeScintillator(NAN)
  , m_ScintiCornerLowerLeft(45.573 * inch, -7.383 * inch)
  , m_ScintiCornerLowerRight(50.604 * inch, -12.972 * inch)
// leave this in in case we ever need those coordinates
// m_ScintiCornerUpperRight(50.809*inch,-12.787*inch),
// m_ScintiCornerUpperLeft(45.778*inch,-7.198*inch),
  , m_NumScintiPlates(16)
  , m_NumSteelPlates(m_NumScintiPlates + 1)
  , m_Active(m_params->get_int_param("active"))
  , m_AbsorberActive(m_params->get_int_param("absorberactive"))
  , m_Layer(0)
{
  m_DeltaPhi += 0.00125 * m_DeltaPhi;
}

//_______________________________________________________________
//_______________________________________________________________
int PHG4Prototype3InnerHcalDetector::IsInPrototype3InnerHcal(G4VPhysicalVolume *volume) const
{
  G4LogicalVolume *logvol = volume->GetLogicalVolume();
  if (m_AbsorberActive && logvol == m_InnerHcalSteelPlate)
  {
      return -1;
  }
  if (m_Active && m_ActiveVolumeSet.find(logvol) != m_ActiveVolumeSet.end())
  {
      return 1;
  }
  return 0;
}

G4LogicalVolume*
PHG4Prototype3InnerHcalDetector::ConstructSteelPlate(G4LogicalVolume *hcalenvelope)
{
  if (!m_InnerHcalSteelPlate)
  {
    G4VSolid *steel_plate;

    std::vector<G4TwoVector> vertexes;
    vertexes.push_back(m_SteelPlateCornerUpperLeft);
    vertexes.push_back(m_SteelPlateCornerUpperRight);
    vertexes.push_back(m_SteelPlateCornerLowerRight);
    vertexes.push_back(m_SteelPlateCornerLowerLeft);
    G4TwoVector zero(0, 0);
    steel_plate = new G4ExtrudedSolid("InnerHcalSteelPlateSolid",
                                      vertexes,
                                      m_SteelZ / 2.0,
                                      zero, 1.0,
                                      zero, 1.0);

    m_VolumeSteel = steel_plate->GetCubicVolume() * m_NumSteelPlates;
    m_InnerHcalSteelPlate = new G4LogicalVolume(steel_plate, G4Material::GetMaterial("Steel_A36"), "InnerHcalSteelPlate", 0, 0, 0);
    G4VisAttributes *visattchk = new G4VisAttributes();
    visattchk->SetVisibility(true);
    visattchk->SetForceSolid(false);
    visattchk->SetColour(G4Colour::Blue());
    m_InnerHcalSteelPlate->SetVisAttributes(visattchk);
  }
  return m_InnerHcalSteelPlate;
}

G4LogicalVolume*
PHG4Prototype3InnerHcalDetector::ConstructScintillatorBoxHiEta(G4LogicalVolume *hcalenvelope)
{
  G4VSolid *scintiboxsolid = new G4Box("InnerHcalScintiMother", m_ScintiX / 2., (m_ScintiGap) / 2., m_ScintiTileZ / 2.);
  //  DisplayVolume(scintiboxsolid,hcalenvelope);

  G4LogicalVolume *scintiboxlogical = new G4LogicalVolume(scintiboxsolid, G4Material::GetMaterial("G4_AIR"), G4String("InnerHcalScintiMother"), 0, 0, 0);
  G4VisAttributes *hcalVisAtt = new G4VisAttributes();
  hcalVisAtt->SetVisibility(true);
  hcalVisAtt->SetForceSolid(false);
  hcalVisAtt->SetColour(G4Colour::Magenta());
  G4LogicalVolume *scintit9_logic = ConstructScintiTile9(hcalenvelope);
  scintit9_logic->SetVisAttributes(hcalVisAtt);
  double distance_to_corner = -m_SteelZ / 2. + m_ScintiTile9DistanceToCorner;
  G4RotationMatrix *Rot;
  Rot = new G4RotationMatrix();
  Rot->rotateX(90 * deg);
  new G4PVPlacement(Rot, G4ThreeVector(-m_ScintiX / 2., 0, distance_to_corner), scintit9_logic, "InnerScinti_9", scintiboxlogical, false, 0, overlapcheck);

  hcalVisAtt = new G4VisAttributes();
  hcalVisAtt->SetVisibility(true);
  hcalVisAtt->SetForceSolid(false);
  hcalVisAtt->SetColour(G4Colour::Blue());
  G4LogicalVolume *scintit10_logic = ConstructScintiTile10(hcalenvelope);
  scintit10_logic->SetVisAttributes(hcalVisAtt);

  distance_to_corner += m_ScintiTile9FrontSize + m_GapBetweenTiles;
  Rot = new G4RotationMatrix();
  Rot->rotateX(90 * deg);
  new G4PVPlacement(Rot, G4ThreeVector(-m_ScintiX / 2., 0, distance_to_corner), scintit10_logic, "InnerScinti_10", scintiboxlogical, false, 0, overlapcheck);

  hcalVisAtt = new G4VisAttributes();
  hcalVisAtt->SetVisibility(true);
  hcalVisAtt->SetForceSolid(false);
  hcalVisAtt->SetColour(G4Colour::Yellow());
  G4LogicalVolume *scintit11_logic = ConstructScintiTile11(hcalenvelope);
  scintit11_logic->SetVisAttributes(hcalVisAtt);

  distance_to_corner += m_ScintiTile10FrontSize + m_GapBetweenTiles;
  Rot = new G4RotationMatrix();
  Rot->rotateX(90 * deg);
  new G4PVPlacement(Rot, G4ThreeVector(-m_ScintiX / 2., 0, distance_to_corner), scintit11_logic, "InnerScinti_11", scintiboxlogical, false, 0, overlapcheck);

  hcalVisAtt = new G4VisAttributes();
  hcalVisAtt->SetVisibility(true);
  hcalVisAtt->SetForceSolid(false);
  hcalVisAtt->SetColour(G4Colour::Cyan());
  G4LogicalVolume *scintit12_logic = ConstructScintiTile12(hcalenvelope);
  scintit12_logic->SetVisAttributes(hcalVisAtt);

  distance_to_corner += m_ScintiTile11FrontSize + m_GapBetweenTiles;
  Rot = new G4RotationMatrix();
  Rot->rotateX(90 * deg);
  new G4PVPlacement(Rot, G4ThreeVector(-m_ScintiX / 2., 0, distance_to_corner), scintit12_logic, "InnerScinti_12", scintiboxlogical, false, 0, overlapcheck);

  //    DisplayVolume(scintiboxlogical,hcalenvelope);
  return scintiboxlogical;
}

G4LogicalVolume*
PHG4Prototype3InnerHcalDetector::ConstructScintiTile9(G4LogicalVolume *hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(m_ScintiTile9CornerUpperLeft);
  vertexes.push_back(m_ScintiTile9CornerUpperRight);
  vertexes.push_back(m_ScintiTile9CornerLowerRight);
  vertexes.push_back(m_ScintiTile9CornerLowerLeft);
  G4TwoVector zero(0, 0);
  G4VSolid *scintit9 = new G4ExtrudedSolid("InnerHcalScintiT9",
                                           vertexes,
                                           m_ScintiTileThickness / 2.0,
                                           zero, 1.0,
                                           zero, 1.0);

  G4LogicalVolume *scintit9_logic = new G4LogicalVolume(scintit9, G4Material::GetMaterial("G4_POLYSTYRENE"), "InnerHcalScintiT9", nullptr, nullptr, nullptr);
  //     DisplayVolume(scintit9,hcalenvelope);
  m_ActiveVolumeSet.insert(scintit9_logic);
  return scintit9_logic;
}

G4LogicalVolume*
PHG4Prototype3InnerHcalDetector::ConstructScintiTile10(G4LogicalVolume *hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(m_ScintiTile10CornerUpperLeft);
  vertexes.push_back(m_ScintiTile10CornerUpperRight);
  vertexes.push_back(m_ScintiTile10CornerLowerRight);
  vertexes.push_back(m_ScintiTile10CornerLowerLeft);
  G4TwoVector zero(0, 0);
  G4VSolid *scintit10 = new G4ExtrudedSolid("InnerHcalScintiT10",
                                            vertexes,
                                            m_ScintiTileThickness / 2.0,
                                            zero, 1.0,
                                            zero, 1.0);

  G4LogicalVolume *scintit10_logic = new G4LogicalVolume(scintit10, G4Material::GetMaterial("G4_POLYSTYRENE"), "InnerHcalScintiT10", nullptr, nullptr, nullptr);
  //     DisplayVolume(scintit10,hcalenvelope);
  m_ActiveVolumeSet.insert(scintit10_logic);
  return scintit10_logic;
}

G4LogicalVolume*
PHG4Prototype3InnerHcalDetector::ConstructScintiTile11(G4LogicalVolume *hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(m_ScintiTile11CornerUpperLeft);
  vertexes.push_back(m_ScintiTile11CornerUpperRight);
  vertexes.push_back(m_ScintiTile11CornerLowerRight);
  vertexes.push_back(m_ScintiTile11CornerLowerLeft);
  G4TwoVector zero(0, 0);
  G4VSolid *scintit11 = new G4ExtrudedSolid("InnerHcalScintiT11",
                                            vertexes,
                                            m_ScintiTileThickness / 2.0,
                                            zero, 1.0,
                                            zero, 1.0);

  G4LogicalVolume *scintit11_logic = new G4LogicalVolume(scintit11, G4Material::GetMaterial("G4_POLYSTYRENE"), "InnerHcalScintiT11", nullptr, nullptr, nullptr);
  //     DisplayVolume(scintit11,hcalenvelope);
  m_ActiveVolumeSet.insert(scintit11_logic);
  return scintit11_logic;
}

G4LogicalVolume*
PHG4Prototype3InnerHcalDetector::ConstructScintiTile12(G4LogicalVolume *hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(m_ScintiTile12CornerUpperLeft);
  vertexes.push_back(m_ScintiTile12CornerUpperRight);
  vertexes.push_back(m_ScintiTile12CornerLowerRight);
  vertexes.push_back(m_ScintiTile12CornerLowerLeft);
  G4TwoVector zero(0, 0);
  G4VSolid *scintit12 = new G4ExtrudedSolid("InnerHcalScintiT12",
                                            vertexes,
                                            m_ScintiTileThickness / 2.0,
                                            zero, 1.0,
                                            zero, 1.0);

  G4LogicalVolume *scintit12_logic = new G4LogicalVolume(scintit12, G4Material::GetMaterial("G4_POLYSTYRENE"), "InnerHcalScintiT12", nullptr, nullptr, nullptr);
  //     DisplayVolume(scintit12,hcalenvelope);
  m_ActiveVolumeSet.insert(scintit12_logic);
  return scintit12_logic;
}

// Construct the envelope and the call the
// actual inner hcal construction
void PHG4Prototype3InnerHcalDetector::Construct(G4LogicalVolume *logicWorld)
{
  G4ThreeVector g4vec(m_params->get_double_param("place_x") * cm,
                      m_params->get_double_param("place_y") * cm,
                      m_params->get_double_param("place_z") * cm);
  G4RotationMatrix *Rot = new G4RotationMatrix();
  Rot->rotateX(m_params->get_double_param("rot_x") * deg);
  Rot->rotateY(m_params->get_double_param("rot_y") * deg);
  Rot->rotateZ(m_params->get_double_param("rot_z") * deg);
  m_InnerHcalAssembly = new G4AssemblyVolume();
  ConstructInnerHcal(logicWorld);
  m_InnerHcalAssembly->MakeImprint(logicWorld, g4vec, Rot, 0, overlapcheck);
  return;
}

int PHG4Prototype3InnerHcalDetector::ConstructInnerHcal(G4LogicalVolume *hcalenvelope)
{
  G4LogicalVolume *steel_plate = ConstructSteelPlate(hcalenvelope);  // bottom steel plate
  if (m_params->get_int_param("hi_eta"))
  {
    m_scintibox = ConstructScintillatorBoxHiEta(hcalenvelope);
  }
  else
  {
    cout << "midrapidity scintillator not implemented" << endl;
    gSystem->Exit(1);
  }
  double phi = 0.;
  double phislat = 0.;
  ostringstream name;
  // the coordinate of the center of the bottom of the bottom steel plate
  // to get the radius of the circle which is the center of the scintillator box

  double bottom_xmiddle_steel_tile = (m_SteelPlateCornerLowerRight.x() + m_SteelPlateCornerLowerLeft.x()) / 2.;
  double bottom_ymiddle_steel_tile = (m_SteelPlateCornerLowerLeft.y() + m_SteelPlateCornerLowerRight.y()) / 2.;
  // the math is not exact, need to move the middle radius by 14mm to
  // get the upper steel plate right
  double middlerad = sqrt(bottom_xmiddle_steel_tile * bottom_xmiddle_steel_tile + bottom_ymiddle_steel_tile * bottom_ymiddle_steel_tile) - 0.14 * cm;
// another fudge factor to get the upper steel location right
  double philow = atan((bottom_ymiddle_steel_tile - (m_ScintiGap * 25. / 32.)) / bottom_xmiddle_steel_tile);

  double scintiangle = GetScintiAngle();

  // need to shift in x to get the upper steel plate close to right
  double xstart = 0;
  double xoff = 0.015 * cm;
  for (int i = 0; i < m_NumSteelPlates; i++)
  //           for (int i = 0; i < 2; i++)
  {
    name.str("");
    name << "InnerHcalSteel_" << i;
    G4RotationMatrix *Rot = new G4RotationMatrix();
    Rot->rotateZ(phi * rad);
    G4ThreeVector g4vec(xstart, 0, 0);
    m_InnerHcalAssembly->AddPlacedVolume(steel_plate, g4vec, Rot);
    if (i > 0)
    {
      double ypos = sin(phi + philow) * middlerad;
      double xpos = cos(phi + philow) * middlerad;
      name.str("");
      name << "InnerHcalScintiBox_" << i;
      Rot = new G4RotationMatrix();
      Rot->rotateZ(scintiangle + phislat);
      G4ThreeVector g4vecsc(xpos + xstart, ypos, 0);
      m_InnerHcalAssembly->AddPlacedVolume(m_scintibox, g4vecsc, Rot);
      phislat += m_DeltaPhi;
    }
    phi += m_DeltaPhi;
    xstart += xoff;
  }
  return 0;
}

// calculate the angle of the bottom scintillator. It is the angle of the top edge
// of the steel plate
double
PHG4Prototype3InnerHcalDetector::GetScintiAngle()
{
  double xlen = m_ScintiCornerLowerRight.x() - m_ScintiCornerLowerLeft.x();
  double ylen = m_ScintiCornerLowerRight.y() - m_ScintiCornerLowerLeft.y();
  double angle = atan(ylen / xlen);
  return angle;
}

int PHG4Prototype3InnerHcalDetector::DisplayVolume(G4VSolid *volume, G4LogicalVolume *logvol, G4RotationMatrix *rotm)
{
  G4LogicalVolume *checksolid = new G4LogicalVolume(volume, G4Material::GetMaterial("G4_POLYSTYRENE"), "DISPLAYLOGICAL", 0, 0, 0);
  DisplayVolume(checksolid, logvol, rotm);
  return 0;
}

int PHG4Prototype3InnerHcalDetector::DisplayVolume(G4LogicalVolume *checksolid, G4LogicalVolume *logvol, G4RotationMatrix *rotm)
{
  static int i = 0;
  G4VisAttributes *visattchk = new G4VisAttributes();
  visattchk->SetVisibility(true);
  visattchk->SetForceSolid(false);
  switch (i)
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

void PHG4Prototype3InnerHcalDetector::Print(const string& what) const
{
  cout << "Inner Hcal Detector:" << endl;
  if (what == "ALL" || what == "VOLUME")
  {
    cout << "Volume Steel: " << m_VolumeSteel / cm3 << " cm^3" << endl;
    cout << "Volume Scintillator: " << m_VolumeScintillator / cm3 << " cm^3" << endl;
  }
  return;
}
