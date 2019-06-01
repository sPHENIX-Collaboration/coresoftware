#include "PHG4Prototype2OuterHcalDetector.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>                   // for PHG4Detector

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4Colour.hh>
#include <Geant4/G4ExtrudedSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>              // for G4RotationMatrix
#include <Geant4/G4String.hh>                      // for G4String
#include <Geant4/G4SystemOfUnits.hh>
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

static const string scintimothername = "OuterHcalScintiMother";
static const string steelplatename = "OuterHcalSteelPlate";

PHG4Prototype2OuterHcalDetector::PHG4Prototype2OuterHcalDetector(PHCompositeNode* Node, PHParameters* parameters, const std::string& dnam)
  : PHG4Detector(Node, dnam)
  , m_Params(parameters)
  , m_OuterHcalSteelPlate(nullptr)
  , m_OuterHcalAssembly(nullptr)
  , m_SteelPlateCornerUpperLeft(1777.6 * mm, -433.5 * mm)
  , m_SteelPlateCornerUpperRight(2600.4 * mm, -417.4 * mm)
  , m_SteelPlateCornerLowerRight(2601.2 * mm, -459.8 * mm)
  , m_SteelPlateCornerLowerLeft(1770.9 * mm, -459.8 * mm)
  ,

  m_ScintiUoneFrontSize(166.2 * mm)
  , m_ScintiUoneCornerUpperLeft(0 * mm, 0 * mm)
  , m_ScintiUoneCornerUpperRight(828.9 * mm, 0 * mm)
  , m_ScintiUoneCornerLowerRight(828.9 * mm, -240.54 * mm)
  , m_ScintiUoneCornerLowerLeft(0 * mm, -m_ScintiUoneFrontSize)
  ,

  m_ScintiU2CornerUpperLeft(0 * mm, 0 * mm)
  , m_ScintiU2CornerUpperRight(828.9 * mm, -74.3 * mm)
  , m_ScintiU2CornerLowerRight(828.9 * mm, -320.44 * mm)
  , m_ScintiU2CornerLowerLeft(0 * mm, -171.0 * mm)
  ,

  m_ScintiT9DistanceToCorner(0.86 * mm)
  , m_ScintiT9FrontSize(241.5 * mm)
  , m_ScintiT9CornerUpperLeft(0 * mm, 0 * mm)
  , m_ScintiT9CornerUpperRight(697.4 * mm, -552.2 * mm)
  , m_ScintiT9CornerLowerRight(697.4 * mm, -697.4 * mm / tan(47.94 / 180. * M_PI) - m_ScintiT9FrontSize)
  , m_ScintiT9CornerLowerLeft(0 * mm, -m_ScintiT9FrontSize)
  ,

  m_ScintiT10FrontSize(241.4 * mm)
  , m_ScintiT10CornerUpperLeft(0 * mm, 0 * mm)
  , m_ScintiT10CornerUpperRight(697.4 * mm, -629.3 * mm)
  , m_ScintiT10CornerLowerRight(697.4 * mm, -697.4 * mm / tan(44.2 / 180. * M_PI) - m_ScintiT10FrontSize)
  , m_ScintiT10CornerLowerLeft(0 * mm, -m_ScintiT10FrontSize)
  ,

  m_ScintiT11FrontSize(241.4 * mm)
  , m_ScintiT11CornerUpperLeft(0 * mm, 0 * mm)
  , m_ScintiT11CornerUpperRight(697.4 * mm, -717.1 * mm)
  , m_ScintiT11CornerLowerRight(697.4 * mm, -697.4 * mm / tan(42.47 / 180. * M_PI) - m_ScintiT11FrontSize)
  , m_ScintiT11CornerLowerLeft(0 * mm, -m_ScintiT11FrontSize)
  ,

  m_ScintiT12FrontSize(312.7 * mm)
  , m_ScintiT12CornerUpperLeft(0 * mm, 0 * mm)
  , m_ScintiT12CornerUpperRight(697.4 * mm, -761.8 * mm)
  , m_ScintiT12CornerLowerRight(392.9 * mm, -827.7)
  , m_ScintiT12CornerLowerLeft(0 * mm, -m_ScintiT12FrontSize)
  ,

  m_ScintiX(828.9)
  , m_ScintiXHiEta(697.4 * mm + 121.09 * mm)
  , m_SteelZ(1600. * mm)
  , m_SizeZ(m_SteelZ)
  , m_ScintiTileZ(m_SteelZ)
  , m_ScintiTileThickness(7 * mm)
  , m_ScintiBoxSmaller(0.02 * mm)
  ,  // blargh - off by 20 microns bc scinti tilt angle, need to revisit at some point
  m_GapBetweenTiles(1 * mm)
  , m_ScintiGap(8.5 * mm)
  , m_TiltAngle(12 * deg)
  , m_DeltaPhi(2 * M_PI / 320.)
  , m_VolumeSteel(NAN)
  , m_VolumeScintillator(NAN)
  , m_NScintiPlates(20)
  , m_NSteelPlates(m_NScintiPlates + 1)
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_AbsorberActiveFlag(m_Params->get_int_param("absorberactive"))
  , m_Layer(0)
{
}

PHG4Prototype2OuterHcalDetector::~PHG4Prototype2OuterHcalDetector()
{
  delete m_OuterHcalAssembly;
}

//_______________________________________________________________
//_______________________________________________________________
int PHG4Prototype2OuterHcalDetector::IsInPrototype2OuterHcal(G4VPhysicalVolume* volume) const
{
  G4LogicalVolume* logvol = volume->GetLogicalVolume();
  if (m_AbsorberActiveFlag && logvol == m_OuterHcalSteelPlate)
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
PHG4Prototype2OuterHcalDetector::ConstructSteelPlate(G4LogicalVolume* hcalenvelope)
{
  if (!m_OuterHcalSteelPlate)
  {
    G4VSolid* steel_plate;
    std::vector<G4TwoVector> vertexes;
    vertexes.push_back(m_SteelPlateCornerUpperLeft);
    vertexes.push_back(m_SteelPlateCornerUpperRight);
    vertexes.push_back(m_SteelPlateCornerLowerRight);
    vertexes.push_back(m_SteelPlateCornerLowerLeft);
    G4TwoVector zero(0, 0);
    steel_plate = new G4ExtrudedSolid("OuterHcalSteelPlateSolid",
                                      vertexes,
                                      m_SizeZ / 2.0,
                                      zero, 1.0,
                                      zero, 1.0);

    m_VolumeSteel = steel_plate->GetCubicVolume() * m_NSteelPlates;
    m_OuterHcalSteelPlate = new G4LogicalVolume(steel_plate, G4Material::GetMaterial("Steel_A36"), "OuterHcalSteelPlate", 0, 0, 0);
    G4VisAttributes visattchk;
    visattchk.SetVisibility(true);
    visattchk.SetForceSolid(false);
    visattchk.SetColour(G4Colour::Blue());
    m_OuterHcalSteelPlate->SetVisAttributes(visattchk);
  }
  return m_OuterHcalSteelPlate;
}

G4LogicalVolume*
PHG4Prototype2OuterHcalDetector::ConstructScintillatorBox(G4LogicalVolume* hcalenvelope)
{
  int copynum = 0;
  G4VSolid* scintiboxsolid = new G4Box(scintimothername, m_ScintiX / 2., (m_ScintiGap - m_ScintiBoxSmaller) / 2., m_ScintiTileZ / 2.);
  //  DisplayVolume(scintiboxsolid,hcalenvelope);
  G4LogicalVolume* scintiboxlogical = new G4LogicalVolume(scintiboxsolid, G4Material::GetMaterial("G4_AIR"), G4String(scintimothername), 0, 0, 0);
  G4VisAttributes* hcalVisAtt = new G4VisAttributes();
  hcalVisAtt->SetVisibility(true);
  hcalVisAtt->SetForceSolid(false);
  hcalVisAtt->SetColour(G4Colour::Red());
  G4LogicalVolume* scintiu1_logic = ConstructScintiTileU1(hcalenvelope);
  scintiu1_logic->SetVisAttributes(hcalVisAtt);

  hcalVisAtt = new G4VisAttributes();
  hcalVisAtt->SetVisibility(true);
  hcalVisAtt->SetForceSolid(false);
  hcalVisAtt->SetColour(G4Colour::Cyan());
  G4LogicalVolume* scintiu2_logic = ConstructScintiTileU2(hcalenvelope);
  scintiu2_logic->SetVisAttributes(hcalVisAtt);
  G4RotationMatrix* Rot;
  Rot = new G4RotationMatrix();
  Rot->rotateX(-90 * deg);
  new G4PVPlacement(Rot, G4ThreeVector(-m_ScintiX / 2., 0, -m_ScintiUoneFrontSize - m_GapBetweenTiles / 2. - m_GapBetweenTiles), scintiu2_logic, (boost::format("OuterScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());

  Rot = new G4RotationMatrix();
  Rot->rotateX(-90 * deg);
  copynum++;
  new G4PVPlacement(Rot, G4ThreeVector(-m_ScintiX / 2., 0, -m_GapBetweenTiles / 2.), scintiu1_logic, (boost::format("OuterScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());

  Rot = new G4RotationMatrix();
  Rot->rotateX(90 * deg);
  copynum++;
  new G4PVPlacement(Rot, G4ThreeVector(-m_ScintiX / 2., 0, m_GapBetweenTiles / 2.), scintiu1_logic, (boost::format("OuterScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());

  Rot = new G4RotationMatrix();
  Rot->rotateX(90 * deg);
  copynum++;
  new G4PVPlacement(Rot, G4ThreeVector(-m_ScintiX / 2., 0, m_ScintiUoneFrontSize + m_GapBetweenTiles / 2. + m_GapBetweenTiles), scintiu2_logic, (boost::format("OuterScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());

  return scintiboxlogical;
}

G4LogicalVolume*
PHG4Prototype2OuterHcalDetector::ConstructScintiTileU1(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(m_ScintiUoneCornerUpperLeft);
  vertexes.push_back(m_ScintiUoneCornerUpperRight);
  vertexes.push_back(m_ScintiUoneCornerLowerRight);
  vertexes.push_back(m_ScintiUoneCornerLowerLeft);
  G4TwoVector zero(0, 0);
  G4VSolid* scintiu1 = new G4ExtrudedSolid("OuterHcalScintiU1",
                                           vertexes,
                                           m_ScintiTileThickness / 2.0,
                                           zero, 1.0,
                                           zero, 1.0);

  G4LogicalVolume* scintiu1_logic = new G4LogicalVolume(scintiu1, G4Material::GetMaterial("G4_POLYSTYRENE"), "OuterHcalScintiU1", nullptr, nullptr, nullptr);
  //   DisplayVolume(scintiu1,hcalenvelope);
  m_ActiveVolumeSet.insert(scintiu1_logic);
  return scintiu1_logic;
}

G4LogicalVolume*
PHG4Prototype2OuterHcalDetector::ConstructScintiTileU2(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(m_ScintiU2CornerUpperLeft);
  vertexes.push_back(m_ScintiU2CornerUpperRight);
  vertexes.push_back(m_ScintiU2CornerLowerRight);
  vertexes.push_back(m_ScintiU2CornerLowerLeft);
  G4TwoVector zero(0, 0);
  G4VSolid* scintiu2 = new G4ExtrudedSolid("OuterHcalScintiU2",
                                           vertexes,
                                           m_ScintiTileThickness / 2.0,
                                           zero, 1.0,
                                           zero, 1.0);

  G4LogicalVolume* scintiu2_logic = new G4LogicalVolume(scintiu2, G4Material::GetMaterial("G4_POLYSTYRENE"), "OuterHcalScintiU2", nullptr, nullptr, nullptr);
  //   DisplayVolume(scintiu2,hcalenvelope);
  m_ActiveVolumeSet.insert(scintiu2_logic);
  return scintiu2_logic;
}

G4LogicalVolume*
PHG4Prototype2OuterHcalDetector::ConstructScintillatorBoxHiEta(G4LogicalVolume* hcalenvelope)
{
  int copynum = 0;
  G4VSolid* scintiboxsolid = new G4Box(scintimothername, m_ScintiX / 2., (m_ScintiGap - m_ScintiBoxSmaller) / 2., m_ScintiTileZ / 2.);
  //  DisplayVolume(scintiboxsolid,hcalenvelope);
  G4LogicalVolume* scintiboxlogical = new G4LogicalVolume(scintiboxsolid, G4Material::GetMaterial("G4_AIR"), G4String(scintimothername), 0, 0, 0);

  G4VisAttributes hcalVisAtt;
  hcalVisAtt.SetVisibility(true);
  hcalVisAtt.SetForceSolid(false);
  hcalVisAtt.SetColour(G4Colour::Magenta());
  G4LogicalVolume* scintit9_logic = ConstructScintiTile9(hcalenvelope);
  scintit9_logic->SetVisAttributes(hcalVisAtt);

  double distance_to_corner = -m_SizeZ / 2. + m_ScintiT9DistanceToCorner;
  G4RotationMatrix* Rot;
  Rot = new G4RotationMatrix();
  Rot->rotateX(90 * deg);
  new G4PVPlacement(Rot, G4ThreeVector(-m_ScintiXHiEta / 2., 0, distance_to_corner), scintit9_logic, (boost::format("OuterScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());

  hcalVisAtt.SetVisibility(true);
  hcalVisAtt.SetForceSolid(false);
  hcalVisAtt.SetColour(G4Colour::Blue());
  G4LogicalVolume* scintit10_logic = ConstructScintiTile10(hcalenvelope);
  scintit10_logic->SetVisAttributes(hcalVisAtt);

  distance_to_corner += m_ScintiT9FrontSize + m_GapBetweenTiles;
  Rot = new G4RotationMatrix();
  Rot->rotateX(90 * deg);
  copynum++;
  new G4PVPlacement(Rot, G4ThreeVector(-m_ScintiXHiEta / 2., 0, distance_to_corner), scintit10_logic, (boost::format("OuterScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());

  hcalVisAtt.SetVisibility(true);
  hcalVisAtt.SetForceSolid(false);
  hcalVisAtt.SetColour(G4Colour::Yellow());
  G4LogicalVolume* scintit11_logic = ConstructScintiTile11(hcalenvelope);
  scintit11_logic->SetVisAttributes(hcalVisAtt);

  distance_to_corner += m_ScintiT10FrontSize + m_GapBetweenTiles;
  Rot = new G4RotationMatrix();
  Rot->rotateX(90 * deg);
  copynum++;
  new G4PVPlacement(Rot, G4ThreeVector(-m_ScintiXHiEta / 2., 0, distance_to_corner), scintit11_logic, (boost::format("OuterScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());

  hcalVisAtt.SetVisibility(true);
  hcalVisAtt.SetForceSolid(false);
  hcalVisAtt.SetColour(G4Colour::Cyan());
  G4LogicalVolume* scintit12_logic = ConstructScintiTile12(hcalenvelope);
  scintit12_logic->SetVisAttributes(hcalVisAtt);

  distance_to_corner += m_ScintiT11FrontSize + m_GapBetweenTiles;
  Rot = new G4RotationMatrix();
  Rot->rotateX(90 * deg);
  copynum++;
  new G4PVPlacement(Rot, G4ThreeVector(-m_ScintiXHiEta / 2., 0, distance_to_corner), scintit12_logic, (boost::format("OuterScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());
  //DisplayVolume(scintiboxlogical,hcalenvelope);
  return scintiboxlogical;
}
G4LogicalVolume*
PHG4Prototype2OuterHcalDetector::ConstructScintiTile9(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(m_ScintiT9CornerUpperLeft);
  vertexes.push_back(m_ScintiT9CornerUpperRight);
  vertexes.push_back(m_ScintiT9CornerLowerRight);
  vertexes.push_back(m_ScintiT9CornerLowerLeft);
  G4TwoVector zero(0, 0);
  G4VSolid* scintit9 = new G4ExtrudedSolid("OuterHcalScintiT9",
                                           vertexes,
                                           m_ScintiTileThickness / 2.0,
                                           zero, 1.0,
                                           zero, 1.0);

  G4LogicalVolume* scintit9_logic = new G4LogicalVolume(scintit9, G4Material::GetMaterial("G4_POLYSTYRENE"), "OuterHcalScintiT9", nullptr, nullptr, nullptr);
  //     DisplayVolume(scintit9,hcalenvelope);
  m_ActiveVolumeSet.insert(scintit9_logic);
  return scintit9_logic;
}

G4LogicalVolume*
PHG4Prototype2OuterHcalDetector::ConstructScintiTile10(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(m_ScintiT10CornerUpperLeft);
  vertexes.push_back(m_ScintiT10CornerUpperRight);
  vertexes.push_back(m_ScintiT10CornerLowerRight);
  vertexes.push_back(m_ScintiT10CornerLowerLeft);
  G4TwoVector zero(0, 0);
  G4VSolid* scintit10 = new G4ExtrudedSolid("OuterHcalScintiT10",
                                            vertexes,
                                            m_ScintiTileThickness / 2.0,
                                            zero, 1.0,
                                            zero, 1.0);

  G4LogicalVolume* scintit10_logic = new G4LogicalVolume(scintit10, G4Material::GetMaterial("G4_POLYSTYRENE"), "OuterHcalScintiT10", nullptr, nullptr, nullptr);
  //     DisplayVolume(scintit10,hcalenvelope);
  m_ActiveVolumeSet.insert(scintit10_logic);
  return scintit10_logic;
}

G4LogicalVolume*
PHG4Prototype2OuterHcalDetector::ConstructScintiTile11(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(m_ScintiT11CornerUpperLeft);
  vertexes.push_back(m_ScintiT11CornerUpperRight);
  vertexes.push_back(m_ScintiT11CornerLowerRight);
  vertexes.push_back(m_ScintiT11CornerLowerLeft);
  G4TwoVector zero(0, 0);
  G4VSolid* scintit11 = new G4ExtrudedSolid("OuterHcalScintiT11",
                                            vertexes,
                                            m_ScintiTileThickness / 2.0,
                                            zero, 1.0,
                                            zero, 1.0);

  G4LogicalVolume* scintit11_logic = new G4LogicalVolume(scintit11, G4Material::GetMaterial("G4_POLYSTYRENE"), "OuterHcalScintiT11", nullptr, nullptr, nullptr);
  //     DisplayVolume(scintit11,hcalenvelope);
  m_ActiveVolumeSet.insert(scintit11_logic);
  return scintit11_logic;
}

G4LogicalVolume*
PHG4Prototype2OuterHcalDetector::ConstructScintiTile12(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(m_ScintiT12CornerUpperLeft);
  vertexes.push_back(m_ScintiT12CornerUpperRight);
  vertexes.push_back(m_ScintiT12CornerLowerRight);
  vertexes.push_back(m_ScintiT12CornerLowerLeft);
  G4TwoVector zero(0, 0);
  G4VSolid* scintit12 = new G4ExtrudedSolid("OuterHcalScintiT12",
                                            vertexes,
                                            m_ScintiTileThickness / 2.0,
                                            zero, 1.0,
                                            zero, 1.0);

  G4LogicalVolume* scintit12_logic = new G4LogicalVolume(scintit12, G4Material::GetMaterial("G4_POLYSTYRENE"), "OuterHcalScintiT12", nullptr, nullptr, nullptr);
  //     DisplayVolume(scintit12,hcalenvelope);
  m_ActiveVolumeSet.insert(scintit12_logic);
  return scintit12_logic;
}

// Construct the envelope and the call the
// actual outer hcal construction
void PHG4Prototype2OuterHcalDetector::Construct(G4LogicalVolume* logicWorld)
{
  G4ThreeVector g4vec(m_Params->get_double_param("place_x") * cm,
                      m_Params->get_double_param("place_y") * cm,
                      m_Params->get_double_param("place_z") * cm);
  G4RotationMatrix Rot;
  Rot.rotateX(m_Params->get_double_param("rot_x") * deg);
  Rot.rotateY(m_Params->get_double_param("rot_y") * deg);
  Rot.rotateZ(m_Params->get_double_param("rot_z") * deg);
  //  ConstructScintillatorBoxHiEta(logicWorld);
  m_OuterHcalAssembly = new G4AssemblyVolume();
  //ConstructSteelPlate(hcal_envelope_log);
  // return;
  ConstructOuterHcal(logicWorld);
  m_OuterHcalAssembly->MakeImprint(logicWorld, g4vec, &Rot, 0, OverlapCheck());
  // this is rather pathetic - there is no way to extract the name when a volume is added
  // to the assembly. The only thing we can do is get an iterator over the placed volumes
  // in the order in which they were placed. Since this code does not install the scintillators
  // for the Al version, parsing the volume names to get the id does not work since it changes
  // So now we loop over all volumes and store them in a map for fast lookup of the row
  int isteel = 0;
  int iscinti = 0;
  vector<G4VPhysicalVolume*>::iterator it = m_OuterHcalAssembly->GetVolumesIterator();
  for (unsigned int i = 0; i < m_OuterHcalAssembly->TotalImprintedVolumes(); i++)
  {
    string volname = (*it)->GetName();
    if (volname.find(steelplatename) != string::npos)
    {
      m_SteelPlateIdMap.insert(make_pair(volname, isteel));
      ++isteel;
    }
    else if (volname.find(scintimothername) != string::npos)
    {
      m_ScintillatorIdMap.insert(make_pair(volname, iscinti));
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

int PHG4Prototype2OuterHcalDetector::ConstructOuterHcal(G4LogicalVolume* hcalenvelope)
{
  G4LogicalVolume* steel_plate = ConstructSteelPlate(hcalenvelope);  // bottom steel plate
  G4LogicalVolume* scintibox = nullptr;
  if (m_Params->get_int_param("hi_eta"))
  {
    scintibox = ConstructScintillatorBoxHiEta(hcalenvelope);
  }
  else
  {
    scintibox = ConstructScintillatorBox(hcalenvelope);
  }
  double phi = 0.;
  double phislat = 0.;
  // the coordinate of the center of the bottom of the bottom steel plate
  // to get the radius of the circle which is the center of the scintillator box
  double bottom_xmiddle_steel_tile = (m_SteelPlateCornerLowerRight.x() - m_SteelPlateCornerLowerLeft.x()) / 2. + m_SteelPlateCornerLowerLeft.x();
  double bottom_ymiddle_steel_tile = m_SteelPlateCornerLowerRight.y();
  double middlerad = sqrt(bottom_xmiddle_steel_tile * bottom_xmiddle_steel_tile + bottom_ymiddle_steel_tile * bottom_ymiddle_steel_tile);
  double philow = atan((bottom_ymiddle_steel_tile - m_ScintiGap / 2.) / bottom_xmiddle_steel_tile);
  double scintiangle = GetScintiAngle();
  for (int i = 0; i < m_NSteelPlates; i++)
  {
    G4RotationMatrix Rot;
    Rot.rotateZ(phi * rad);
    G4ThreeVector g4vec(0, 0, 0);
    m_OuterHcalAssembly->AddPlacedVolume(steel_plate, g4vec, &Rot);
    if (i > 0)
    {
      double ypos = sin(phi + philow) * middlerad;
      double xpos = cos(phi + philow) * middlerad;
      // the center of the scintillator is not the center of the inner hcal
      // but depends on the tilt angle. Therefore we need to shift
      // the center from the mid point
      ypos += sin((-m_TiltAngle) / rad - phi);
      xpos -= cos((-m_TiltAngle) / rad - phi);
      G4RotationMatrix Rota;
      Rota.rotateZ(scintiangle + phislat);
      G4ThreeVector g4vec(xpos, ypos, 0);

      m_OuterHcalAssembly->AddPlacedVolume(scintibox, g4vec, &Rota);
      phislat += m_DeltaPhi;
    }
    phi += m_DeltaPhi;
  }
  return 0;
}

// calculate the angle of the bottom scintillator. It is the angle of the top edge
// of the steel plate
double
PHG4Prototype2OuterHcalDetector::GetScintiAngle()
{
  double xlen = m_SteelPlateCornerUpperRight.x() - m_SteelPlateCornerUpperLeft.x();
  double ylen = m_SteelPlateCornerUpperRight.y() - m_SteelPlateCornerUpperLeft.y();
  double angle = atan(ylen / xlen);
  return angle;
}

void PHG4Prototype2OuterHcalDetector::Print(const string& what) const
{
  cout << "Outer Hcal Detector:" << endl;
  if (what == "ALL" || what == "VOLUME")
  {
    cout << "Volume Steel: " << m_VolumeSteel / cm3 << " cm^3" << endl;
    cout << "Volume Scintillator: " << m_VolumeScintillator / cm3 << " cm^3" << endl;
  }
  return;
}

int PHG4Prototype2OuterHcalDetector::get_scinti_row_id(const string& volname)
{
  int id = -9999;
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

int PHG4Prototype2OuterHcalDetector::get_steel_plate_id(const string& volname)
{
  int id = -9999;
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
