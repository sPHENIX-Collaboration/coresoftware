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
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4TwoVector.hh>
#include <Geant4/G4VisAttributes.hh>

#include <boost/format.hpp>

#include <cmath>
#include <sstream>

using namespace std;

static const string scintimothername = "InnerHcalScintiMother";
static const string steelplatename = "InnerHcalSteelPlate";

PHG4Prototype2InnerHcalDetector::PHG4Prototype2InnerHcalDetector(PHCompositeNode* Node, PHParameters* parameters, const std::string& dnam)
  : PHG4Detector(Node, dnam)
  , m_Params(parameters)
  , m_InnerHcalSteelPlate(nullptr)
  , m_InnerHcalAssembly(nullptr)
  , m_SteelPlateCornerUpperLeft(1157.5 * mm, -151.44 * mm)
  , m_SteelPlateCornerUpperRight(1308.5 * mm, -286.96 * mm)
  , m_SteelPlateCornerLowerRight(1298.8 * mm, -297.39 * mm)
  , m_SteelPlateCornerLowerLeft(1155.8 * mm, -163.92 * mm)
  ,

  m_ScintiUoneFrontSize(105.9 * mm)
  , m_ScintiUoneCornerUpperLeft(0 * mm, 0 * mm)
  , m_ScintiUoneCornerUpperRight(198.1 * mm, 0 * mm)
  , m_ScintiUoneCornerLowerRight(198.1 * mm, -121.3 * mm)
  , m_ScintiUoneCornerLowerLeft(0 * mm, -m_ScintiUoneFrontSize)
  ,

  m_ScintiU2CornerUpperLeft(0 * mm, 0 * mm)
  , m_ScintiU2CornerUpperRight(198.1 * mm, -15.4 * mm)
  , m_ScintiU2CornerLowerRight(198.1 * mm, -141.5 * mm)
  , m_ScintiU2CornerLowerLeft(0 * mm, -110.59 * mm)
  ,

  m_ScintiT9DistanceToCorner(26.44 * mm)
  , m_ScintiT9FrontSize(140.3 * mm)
  , m_ScintiT9CornerUpperLeft(0 * mm, 0 * mm)
  , m_ScintiT9CornerUpperRight(198.1 * mm, -134.4 * mm)
  , m_ScintiT9CornerLowerRight(198.1 * mm, -198.1 * mm / tan(52.02 / 180. * M_PI) - m_ScintiT9FrontSize)
  , m_ScintiT9CornerLowerLeft(0 * mm, -m_ScintiT9FrontSize)
  ,

  m_ScintiT10FrontSize(149.2 * mm)
  , m_ScintiT10CornerUpperLeft(0 * mm, 0 * mm)
  , m_ScintiT10CornerUpperRight(198.1 * mm, -154.6 * mm)
  , m_ScintiT10CornerLowerRight(198.1 * mm, -198.1 * mm / tan(48.34 / 180. * M_PI) - m_ScintiT10FrontSize)
  , m_ScintiT10CornerLowerLeft(0 * mm, -m_ScintiT10FrontSize)
  ,

  m_ScintiT11FrontSize(144.3 * mm)
  , m_ScintiT11CornerUpperLeft(0 * mm, 0 * mm)
  , m_ScintiT11CornerUpperRight(198.1 * mm, -176.2 * mm)
  , m_ScintiT11CornerLowerRight(198.1 * mm, -198.1 * mm / tan(45.14 / 180. * M_PI) - m_ScintiT11FrontSize)
  , m_ScintiT11CornerLowerLeft(0 * mm, -m_ScintiT11FrontSize)
  ,

  m_ScintiT12FrontSize(186.6 * mm)
  , m_ScintiT12CornerUpperLeft(0 * mm, 0 * mm)
  , m_ScintiT12CornerUpperRight(198.1 * mm, -197.11 * mm)
  , m_ScintiT12CornerLowerRight(198.1 * mm, -198.1 * mm / tan(41.47 / 180. * M_PI) - m_ScintiT12FrontSize)
  , m_ScintiT12CornerLowerLeft(0 * mm, -m_ScintiT12FrontSize)
  ,

  m_ScintiX(198.1)
  , m_SteelZ(901.7 * mm)
  , m_SizeZ(m_SteelZ)
  , m_ScintiTileZ(m_SteelZ)
  , m_ScintiTileThickness(7 * mm)
  , m_ScintiBoxSmaller(0.02 * mm)
  ,  // blargh - off by 20 microns bc scinti tilt angle, need to revisit at some point
  m_GapBetweenTiles(1 * mm)
  , m_ScintiGap(8.5 * mm)
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

PHG4Prototype2InnerHcalDetector::~PHG4Prototype2InnerHcalDetector()
{
  delete m_InnerHcalAssembly;
}

//_______________________________________________________________
//_______________________________________________________________
int PHG4Prototype2InnerHcalDetector::IsInPrototype2InnerHcal(G4VPhysicalVolume* volume) const
{
  G4LogicalVolume* logvol = volume->GetLogicalVolume();
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
    vertexes.push_back(m_SteelPlateCornerUpperLeft);
    vertexes.push_back(m_SteelPlateCornerUpperRight);
    vertexes.push_back(m_SteelPlateCornerLowerRight);
    vertexes.push_back(m_SteelPlateCornerLowerLeft);
    G4TwoVector zero(0, 0);
    steel_plate = new G4ExtrudedSolid("InnerHcalSteelPlateSolid",
                                      vertexes,
                                      m_SizeZ / 2.0,
                                      zero, 1.0,
                                      zero, 1.0);

    m_VolumeSteel = steel_plate->GetCubicVolume() * m_NSteelPlates;
    m_InnerHcalSteelPlate = new G4LogicalVolume(steel_plate, G4Material::GetMaterial("Steel_A36"), steelplatename, 0, 0, 0);
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
  int copynum = 0;
  G4VSolid* scintiboxsolid = new G4Box(scintimothername, m_ScintiX / 2., (m_ScintiGap - m_ScintiBoxSmaller) / 2., m_ScintiTileZ / 2.);
  //  DisplayVolume(scintiboxsolid,hcalenvelope);

  G4LogicalVolume* scintiboxlogical = new G4LogicalVolume(scintiboxsolid, G4Material::GetMaterial("G4_AIR"), scintimothername, 0, 0, 0);
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
  new G4PVPlacement(Rot, G4ThreeVector(-m_ScintiX / 2., 0, distance_to_corner), scintit9_logic, (boost::format("InnerScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());

  hcalVisAtt.SetVisibility(true);
  hcalVisAtt.SetForceSolid(false);
  hcalVisAtt.SetColour(G4Colour::Blue());
  G4LogicalVolume* scintit10_logic = ConstructScintiTile10(hcalenvelope);
  scintit10_logic->SetVisAttributes(hcalVisAtt);

  distance_to_corner += m_ScintiT9FrontSize + m_GapBetweenTiles;
  Rot = new G4RotationMatrix();
  Rot->rotateX(90 * deg);
  copynum++;
  new G4PVPlacement(Rot, G4ThreeVector(-m_ScintiX / 2., 0, distance_to_corner), scintit10_logic, (boost::format("InnerScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());

  hcalVisAtt.SetVisibility(true);
  hcalVisAtt.SetForceSolid(false);
  hcalVisAtt.SetColour(G4Colour::Yellow());
  G4LogicalVolume* scintit11_logic = ConstructScintiTile11(hcalenvelope);
  scintit11_logic->SetVisAttributes(hcalVisAtt);

  distance_to_corner += m_ScintiT10FrontSize + m_GapBetweenTiles;
  Rot = new G4RotationMatrix();
  Rot->rotateX(90 * deg);
  copynum++;
  new G4PVPlacement(Rot, G4ThreeVector(-m_ScintiX / 2., 0, distance_to_corner), scintit11_logic, (boost::format("InnerScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());

  hcalVisAtt.SetVisibility(true);
  hcalVisAtt.SetForceSolid(false);
  hcalVisAtt.SetColour(G4Colour::Cyan());
  G4LogicalVolume* scintit12_logic = ConstructScintiTile12(hcalenvelope);
  scintit12_logic->SetVisAttributes(hcalVisAtt);

  distance_to_corner += m_ScintiT11FrontSize + m_GapBetweenTiles;
  Rot = new G4RotationMatrix();
  Rot->rotateX(90 * deg);
  copynum++;
  new G4PVPlacement(Rot, G4ThreeVector(-m_ScintiX / 2., 0, distance_to_corner), scintit12_logic, (boost::format("InnerScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());

  //    DisplayVolume(scintiboxlogical,hcalenvelope);
  return scintiboxlogical;
}

G4LogicalVolume*
PHG4Prototype2InnerHcalDetector::ConstructScintillatorBox(G4LogicalVolume* hcalenvelope)
{
  int copynum = 0;
  G4VSolid* scintiboxsolid = new G4Box(scintimothername, m_ScintiX / 2., (m_ScintiGap - m_ScintiBoxSmaller) / 2., m_ScintiTileZ / 2.);
  //    DisplayVolume(scintiboxsolid,hcalenvelope);
  G4LogicalVolume* scintiboxlogical = new G4LogicalVolume(scintiboxsolid, G4Material::GetMaterial("G4_AIR"), scintimothername, 0, 0, 0);
  G4VisAttributes hcalVisAtt;
  hcalVisAtt.SetVisibility(true);
  hcalVisAtt.SetForceSolid(false);
  hcalVisAtt.SetColour(G4Colour::Red());
  G4LogicalVolume* scintiu1_logic = ConstructScintiTileU1(hcalenvelope);
  scintiu1_logic->SetVisAttributes(hcalVisAtt);

  hcalVisAtt.SetVisibility(true);
  hcalVisAtt.SetForceSolid(false);
  hcalVisAtt.SetColour(G4Colour::Cyan());
  G4LogicalVolume* scintiu2_logic = ConstructScintiTileU2(hcalenvelope);
  scintiu2_logic->SetVisAttributes(hcalVisAtt);
  G4RotationMatrix* Rot;
  Rot = new G4RotationMatrix();
  Rot->rotateX(-90 * deg);
  new G4PVPlacement(Rot, G4ThreeVector(-m_ScintiX / 2., 0, -m_ScintiUoneFrontSize - m_GapBetweenTiles / 2. - m_GapBetweenTiles), scintiu2_logic, (boost::format("InnerScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());

  Rot = new G4RotationMatrix();
  Rot->rotateX(-90 * deg);
  copynum++;
  new G4PVPlacement(Rot, G4ThreeVector(-m_ScintiX / 2., 0, -m_GapBetweenTiles / 2.), scintiu1_logic, (boost::format("InnerScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());

  Rot = new G4RotationMatrix();
  Rot->rotateX(90 * deg);
  copynum++;
  new G4PVPlacement(Rot, G4ThreeVector(-m_ScintiX / 2., 0, m_GapBetweenTiles / 2.), scintiu1_logic, (boost::format("InnerScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());

  Rot = new G4RotationMatrix();
  Rot->rotateX(90 * deg);
  copynum++;
  new G4PVPlacement(Rot, G4ThreeVector(-m_ScintiX / 2., 0, m_ScintiUoneFrontSize + m_GapBetweenTiles / 2. + m_GapBetweenTiles), scintiu2_logic, (boost::format("InnerScinti_%d") % copynum).str(), scintiboxlogical, false, copynum, OverlapCheck());
  //  DisplayVolume(scintiboxlogical,hcalenvelope);
  return scintiboxlogical;
}

G4LogicalVolume*
PHG4Prototype2InnerHcalDetector::ConstructScintiTileU1(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(m_ScintiUoneCornerUpperLeft);
  vertexes.push_back(m_ScintiUoneCornerUpperRight);
  vertexes.push_back(m_ScintiUoneCornerLowerRight);
  vertexes.push_back(m_ScintiUoneCornerLowerLeft);
  G4TwoVector zero(0, 0);
  G4VSolid* scintiu1 = new G4ExtrudedSolid("InnerHcalScintiU1",
                                           vertexes,
                                           m_ScintiTileThickness / 2.0,
                                           zero, 1.0,
                                           zero, 1.0);

  G4LogicalVolume* scintiu1_logic = new G4LogicalVolume(scintiu1, G4Material::GetMaterial("G4_POLYSTYRENE"), "InnerHcalScintiU1", nullptr, nullptr, nullptr);
  //     DisplayVolume(scintiu1,hcalenvelope);
  m_ActiveVolumeSet.insert(scintiu1_logic);
  return scintiu1_logic;
}

G4LogicalVolume*
PHG4Prototype2InnerHcalDetector::ConstructScintiTileU2(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(m_ScintiU2CornerUpperLeft);
  vertexes.push_back(m_ScintiU2CornerUpperRight);
  vertexes.push_back(m_ScintiU2CornerLowerRight);
  vertexes.push_back(m_ScintiU2CornerLowerLeft);
  G4TwoVector zero(0, 0);
  G4VSolid* scintiu2 = new G4ExtrudedSolid("InnerHcalScintiU2",
                                           vertexes,
                                           m_ScintiTileThickness / 2.0,
                                           zero, 1.0,
                                           zero, 1.0);

  G4LogicalVolume* scintiu2_logic = new G4LogicalVolume(scintiu2, G4Material::GetMaterial("G4_POLYSTYRENE"), "InnerHcalScintiU2", nullptr, nullptr, nullptr);
  //   DisplayVolume(scintiu2,hcalenvelope);
  m_ActiveVolumeSet.insert(scintiu2_logic);
  return scintiu2_logic;
}

G4LogicalVolume*
PHG4Prototype2InnerHcalDetector::ConstructScintiTile9(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(m_ScintiT9CornerUpperLeft);
  vertexes.push_back(m_ScintiT9CornerUpperRight);
  vertexes.push_back(m_ScintiT9CornerLowerRight);
  vertexes.push_back(m_ScintiT9CornerLowerLeft);
  G4TwoVector zero(0, 0);
  G4VSolid* scintit9 = new G4ExtrudedSolid("InnerHcalScintiT9",
                                           vertexes,
                                           m_ScintiTileThickness / 2.0,
                                           zero, 1.0,
                                           zero, 1.0);

  G4LogicalVolume* scintit9_logic = new G4LogicalVolume(scintit9, G4Material::GetMaterial("G4_POLYSTYRENE"), "InnerHcalScintiT9", nullptr, nullptr, nullptr);
  //     DisplayVolume(scintit9,hcalenvelope);
  m_ActiveVolumeSet.insert(scintit9_logic);
  return scintit9_logic;
}

G4LogicalVolume*
PHG4Prototype2InnerHcalDetector::ConstructScintiTile10(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(m_ScintiT10CornerUpperLeft);
  vertexes.push_back(m_ScintiT10CornerUpperRight);
  vertexes.push_back(m_ScintiT10CornerLowerRight);
  vertexes.push_back(m_ScintiT10CornerLowerLeft);
  G4TwoVector zero(0, 0);
  G4VSolid* scintit10 = new G4ExtrudedSolid("InnerHcalScintiT10",
                                            vertexes,
                                            m_ScintiTileThickness / 2.0,
                                            zero, 1.0,
                                            zero, 1.0);

  G4LogicalVolume* scintit10_logic = new G4LogicalVolume(scintit10, G4Material::GetMaterial("G4_POLYSTYRENE"), "InnerHcalScintiT10", nullptr, nullptr, nullptr);
  //     DisplayVolume(scintit10,hcalenvelope);
  m_ActiveVolumeSet.insert(scintit10_logic);
  return scintit10_logic;
}

G4LogicalVolume*
PHG4Prototype2InnerHcalDetector::ConstructScintiTile11(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(m_ScintiT11CornerUpperLeft);
  vertexes.push_back(m_ScintiT11CornerUpperRight);
  vertexes.push_back(m_ScintiT11CornerLowerRight);
  vertexes.push_back(m_ScintiT11CornerLowerLeft);
  G4TwoVector zero(0, 0);
  G4VSolid* scintit11 = new G4ExtrudedSolid("InnerHcalScintiT11",
                                            vertexes,
                                            m_ScintiTileThickness / 2.0,
                                            zero, 1.0,
                                            zero, 1.0);

  G4LogicalVolume* scintit11_logic = new G4LogicalVolume(scintit11, G4Material::GetMaterial("G4_POLYSTYRENE"), "InnerHcalScintiT11", nullptr, nullptr, nullptr);
  //     DisplayVolume(scintit11,hcalenvelope);
  m_ActiveVolumeSet.insert(scintit11_logic);
  return scintit11_logic;
}

G4LogicalVolume*
PHG4Prototype2InnerHcalDetector::ConstructScintiTile12(G4LogicalVolume* hcalenvelope)
{
  std::vector<G4TwoVector> vertexes;
  vertexes.push_back(m_ScintiT12CornerUpperLeft);
  vertexes.push_back(m_ScintiT12CornerUpperRight);
  vertexes.push_back(m_ScintiT12CornerLowerRight);
  vertexes.push_back(m_ScintiT12CornerLowerLeft);
  G4TwoVector zero(0, 0);
  G4VSolid* scintit12 = new G4ExtrudedSolid("InnerHcalScintiT12",
                                            vertexes,
                                            m_ScintiTileThickness / 2.0,
                                            zero, 1.0,
                                            zero, 1.0);

  G4LogicalVolume* scintit12_logic = new G4LogicalVolume(scintit12, G4Material::GetMaterial("G4_POLYSTYRENE"), "InnerHcalScintiT12", nullptr, nullptr, nullptr);
  //     DisplayVolume(scintit12,hcalenvelope);
  m_ActiveVolumeSet.insert(scintit12_logic);
  return scintit12_logic;
}

// Construct the envelope and the call the
// actual inner hcal construction
void PHG4Prototype2InnerHcalDetector::Construct(G4LogicalVolume* logicWorld)
{
  G4ThreeVector g4vec(m_Params->get_double_param("place_x") * cm,
                      m_Params->get_double_param("place_y") * cm,
                      m_Params->get_double_param("place_z") * cm);
  G4RotationMatrix Rot;
  Rot.rotateX(m_Params->get_double_param("rot_x") * deg);
  Rot.rotateY(m_Params->get_double_param("rot_y") * deg);
  Rot.rotateZ(m_Params->get_double_param("rot_z") * deg);
  //  ConstructScintiTile9(logicWorld);
  //    ConstructScintillatorBoxHiEta(logicWorld);
  //ConstructScintillatorBox(logicWorld);
  //  return;
  m_InnerHcalAssembly = new G4AssemblyVolume();
  //ConstructSteelPlate(hcal_envelope_log);
  // return;
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
  for (unsigned int i = 0; i < m_InnerHcalAssembly->TotalImprintedVolumes(); i++)
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

int PHG4Prototype2InnerHcalDetector::ConstructInnerHcal(G4LogicalVolume* hcalenvelope)
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
  double bottom_xmiddle_steel_tile = (m_SteelPlateCornerLowerRight.x() + m_SteelPlateCornerLowerLeft.x()) / 2.;
  double bottom_ymiddle_steel_tile = (m_SteelPlateCornerLowerLeft.y() + m_SteelPlateCornerLowerRight.y()) / 2.;
  double middlerad = sqrt(bottom_xmiddle_steel_tile * bottom_xmiddle_steel_tile + bottom_ymiddle_steel_tile * bottom_ymiddle_steel_tile);
  double philow = atan((bottom_ymiddle_steel_tile - m_ScintiGap / 2. - 0.87 * mm) / bottom_xmiddle_steel_tile);
  double scintiangle = GetScintiAngle();
  for (int i = 0; i < m_NSteelPlates; i++)
  {
    G4RotationMatrix Rot;
    Rot.rotateZ(phi * rad);
    G4ThreeVector g4vec(0, 0, 0);
    m_InnerHcalAssembly->AddPlacedVolume(steel_plate, g4vec, &Rot);
    if (i > 0)
    {
      double ypos = sin(phi + philow) * middlerad;
      double xpos = cos(phi + philow) * middlerad;
      G4RotationMatrix Rot1;
      Rot1.rotateZ(scintiangle + phislat);
      G4ThreeVector g4vecsc(xpos, ypos, 0);
      m_InnerHcalAssembly->AddPlacedVolume(scintibox, g4vecsc, &Rot1);
      phislat += m_DeltaPhi;
    }
    phi += m_DeltaPhi;
  }
  return 0;
}

// calculate the angle of the bottom scintillator. It is the angle of the top edge
// of the steel plate
double
PHG4Prototype2InnerHcalDetector::GetScintiAngle()
{
  double xlen = m_SteelPlateCornerUpperRight.x() - m_SteelPlateCornerUpperLeft.x();
  double ylen = m_SteelPlateCornerUpperRight.y() - m_SteelPlateCornerUpperLeft.y();
  double angle = atan(ylen / xlen);
  return angle;
}

void PHG4Prototype2InnerHcalDetector::Print(const string& what) const
{
  cout << "Inner Hcal Detector:" << endl;
  if (what == "ALL" || what == "VOLUME")
  {
    cout << "Volume Steel: " << m_VolumeSteel / cm3 << " cm^3" << endl;
    cout << "Volume Scintillator: " << m_VolumeScintillator / cm3 << " cm^3" << endl;
  }
  return;
}
int PHG4Prototype2InnerHcalDetector::get_scinti_row_id(const string& volname)
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

int PHG4Prototype2InnerHcalDetector::get_steel_plate_id(const string& volname)
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
