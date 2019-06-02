#include "PHG4Prototype3InnerHcalDetector.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>                   // for PHG4Detector
#include <g4main/PHG4Units.h>

#include <TSystem.h>

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4Colour.hh>
#include <Geant4/G4ExtrudedSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>              // for G4RotationMatrix
#include <Geant4/G4String.hh>                      // for G4String
#include <Geant4/G4SystemOfUnits.hh>               // for mm, deg, cm, cm3, rad
#include <Geant4/G4ThreeVector.hh>                 // for G4ThreeVector
#include <Geant4/G4TwoVector.hh>
#include <Geant4/G4VisAttributes.hh>
#include <Geant4/G4VPhysicalVolume.hh>             // for G4VPhysicalVolume
#include <Geant4/G4VSolid.hh>                      // for G4VSolid

#include <boost/format.hpp>

#include <cmath>
#include <iostream>                                // for operator<<, endl
#include <sstream>
#include <utility>                                 // for pair, make_pair
#include <vector>                                  // for vector, vector<>::...

class PHCompositeNode;

using namespace std;

static const string scintimothername = "InnerHcalScintiMother";
static const string steelplatename = "InnerHcalSteelPlate";

PHG4Prototype3InnerHcalDetector::PHG4Prototype3InnerHcalDetector(PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam)
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
  , m_VolumeScintillator(0)
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
// patch to get the upper steel plate location right within less than a mm
  m_DeltaPhi += 0.00125 * m_DeltaPhi;
}

PHG4Prototype3InnerHcalDetector::~PHG4Prototype3InnerHcalDetector()
{
  delete m_InnerHcalAssembly;
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

G4LogicalVolume *
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
    m_InnerHcalSteelPlate = new G4LogicalVolume(steel_plate, G4Material::GetMaterial(m_params->get_string_param("material")), steelplatename, 0, 0, 0);
    G4VisAttributes visattchk;
    visattchk.SetVisibility(true);
    visattchk.SetForceSolid(false);
    visattchk.SetColour(G4Colour::Blue());
    m_InnerHcalSteelPlate->SetVisAttributes(visattchk);
  }
  return m_InnerHcalSteelPlate;
}

G4LogicalVolume *
PHG4Prototype3InnerHcalDetector::ConstructScintillatorBoxHiEta(G4LogicalVolume *hcalenvelope)
{
  int copynum = 0;
  G4VSolid *scintiboxsolid = new G4Box(scintimothername, m_ScintiX / 2., (m_ScintiGap) / 2., m_ScintiTileZ / 2.);
  //DisplayVolume(scintiboxsolid,hcalenvelope);

  G4LogicalVolume *scintiboxlogical = new G4LogicalVolume(scintiboxsolid, G4Material::GetMaterial("G4_AIR"), G4String(scintimothername), 0, 0, 0);
  G4VisAttributes hcalVisAtt;
  hcalVisAtt.SetVisibility(true);
  hcalVisAtt.SetForceSolid(false);
  hcalVisAtt.SetColour(G4Colour::Magenta());
  G4LogicalVolume *scintit9_logic = ConstructScintiTile9(hcalenvelope);
  scintit9_logic->SetVisAttributes(hcalVisAtt);
  double distance_to_corner = -m_SteelZ / 2. + m_ScintiTile9DistanceToCorner;
  G4RotationMatrix *Rot;
  Rot = new G4RotationMatrix();
  Rot->rotateX(90 * deg);
  new G4PVPlacement(Rot, G4ThreeVector(-m_ScintiX / 2., 0, distance_to_corner), scintit9_logic, (boost::format("InnerScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());

  hcalVisAtt.SetVisibility(true);
  hcalVisAtt.SetForceSolid(false);
  hcalVisAtt.SetColour(G4Colour::Blue());
  G4LogicalVolume *scintit10_logic = ConstructScintiTile10(hcalenvelope);
  scintit10_logic->SetVisAttributes(hcalVisAtt);

  distance_to_corner += m_ScintiTile9FrontSize + m_GapBetweenTiles;
  Rot = new G4RotationMatrix();
  Rot->rotateX(90 * deg);
  copynum++;
  new G4PVPlacement(Rot, G4ThreeVector(-m_ScintiX / 2., 0, distance_to_corner), scintit10_logic, (boost::format("InnerScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());

  hcalVisAtt.SetVisibility(true);
  hcalVisAtt.SetForceSolid(false);
  hcalVisAtt.SetColour(G4Colour::Yellow());
  G4LogicalVolume *scintit11_logic = ConstructScintiTile11(hcalenvelope);
  scintit11_logic->SetVisAttributes(hcalVisAtt);

  distance_to_corner += m_ScintiTile10FrontSize + m_GapBetweenTiles;
  Rot = new G4RotationMatrix();
  Rot->rotateX(90 * deg);
  copynum++;
  new G4PVPlacement(Rot, G4ThreeVector(-m_ScintiX / 2., 0, distance_to_corner), scintit11_logic, (boost::format("InnerScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());

  hcalVisAtt.SetVisibility(true);
  hcalVisAtt.SetForceSolid(false);
  hcalVisAtt.SetColour(G4Colour::Cyan());
  G4LogicalVolume *scintit12_logic = ConstructScintiTile12(hcalenvelope);
  scintit12_logic->SetVisAttributes(hcalVisAtt);

  distance_to_corner += m_ScintiTile11FrontSize + m_GapBetweenTiles;
  Rot = new G4RotationMatrix();
  Rot->rotateX(90 * deg);
  copynum++;
  new G4PVPlacement(Rot, G4ThreeVector(-m_ScintiX / 2., 0, distance_to_corner), scintit12_logic, (boost::format("InnerScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());
  //    DisplayVolume(scintiboxlogical,hcalenvelope);
  return scintiboxlogical;
}

G4LogicalVolume *
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

  m_VolumeScintillator += scintit9->GetCubicVolume() * m_NumScintiPlates;
  G4LogicalVolume *scintit9_logic = new G4LogicalVolume(scintit9, G4Material::GetMaterial("G4_POLYSTYRENE"), "InnerHcalScintiT9", nullptr, nullptr, nullptr);
  //     DisplayVolume(scintit9,hcalenvelope);
  m_ActiveVolumeSet.insert(scintit9_logic);
  return scintit9_logic;
}

G4LogicalVolume *
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

  m_VolumeScintillator += scintit10->GetCubicVolume() * m_NumScintiPlates;
  G4LogicalVolume *scintit10_logic = new G4LogicalVolume(scintit10, G4Material::GetMaterial("G4_POLYSTYRENE"), "InnerHcalScintiT10", nullptr, nullptr, nullptr);
  //     DisplayVolume(scintit10,hcalenvelope);
  m_ActiveVolumeSet.insert(scintit10_logic);
  return scintit10_logic;
}

G4LogicalVolume *
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

  m_VolumeScintillator += scintit11->GetCubicVolume() * m_NumScintiPlates;
  G4LogicalVolume *scintit11_logic = new G4LogicalVolume(scintit11, G4Material::GetMaterial("G4_POLYSTYRENE"), "InnerHcalScintiT11", nullptr, nullptr, nullptr);
  //     DisplayVolume(scintit11,hcalenvelope);
  m_ActiveVolumeSet.insert(scintit11_logic);
  return scintit11_logic;
}

G4LogicalVolume *
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

  m_VolumeScintillator += scintit12->GetCubicVolume() * m_NumScintiPlates;
  G4LogicalVolume *scintit12_logic = new G4LogicalVolume(scintit12, G4Material::GetMaterial("G4_POLYSTYRENE"), "InnerHcalScintiT12", nullptr, nullptr, nullptr);
  //       DisplayVolume(scintit12,hcalenvelope);
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
  G4RotationMatrix Rot;
  Rot.rotateX(m_params->get_double_param("rot_x") * deg);
  Rot.rotateY(m_params->get_double_param("rot_y") * deg);
  Rot.rotateZ(m_params->get_double_param("rot_z") * deg);
  m_InnerHcalAssembly = new G4AssemblyVolume();
  ConstructInnerHcal(logicWorld);
  m_InnerHcalAssembly->MakeImprint(logicWorld, g4vec, &Rot, 0, OverlapCheck());
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

int PHG4Prototype3InnerHcalDetector::ConstructInnerHcal(G4LogicalVolume *hcalenvelope)
{
  G4LogicalVolume *steel_plate = ConstructSteelPlate(hcalenvelope);  // bottom steel plate
  if (m_params->get_int_param("scintillators"))
  {
    if (m_params->get_int_param("hi_eta"))
    {
      m_scintibox = ConstructScintillatorBoxHiEta(hcalenvelope);
    }
    else
    {
      cout << "midrapidity scintillator not implemented" << endl;
      gSystem->Exit(1);
    }
  }
  double phi = 0.;
  double phislat = 0.;
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
  {
    G4RotationMatrix Rot;
    Rot.rotateZ(phi * rad);
    G4ThreeVector g4vec(xstart, 0, 0);
    m_InnerHcalAssembly->AddPlacedVolume(steel_plate, g4vec, &Rot);
    if (m_scintibox && i > 0)
    {
      double ypos = sin(phi + philow) * middlerad;
      double xpos = cos(phi + philow) * middlerad;
      G4RotationMatrix Rot1;
      Rot1.rotateZ(scintiangle + phislat);
      G4ThreeVector g4vecsc(xpos + xstart, ypos, 0);
      m_InnerHcalAssembly->AddPlacedVolume(m_scintibox, g4vecsc, &Rot1);
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

void PHG4Prototype3InnerHcalDetector::Print(const string &what) const
{
  cout << "Inner Hcal Detector:" << endl;
  if (what == "ALL" || what == "VOLUME")
  {
    cout << "Volume Steel: " << m_VolumeSteel / cm3 << " cm^3" << endl;
    cout << "Volume Scintillator: " << m_VolumeScintillator / cm3 << " cm^3" << endl;
  }
  return;
}

int PHG4Prototype3InnerHcalDetector::get_scinti_row_id(const string &volname)
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

int PHG4Prototype3InnerHcalDetector::get_steel_plate_id(const string &volname)
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
