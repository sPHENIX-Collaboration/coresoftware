#include "PHG4ForwardEcalDetector.h"

#include "PHG4ForwardEcalDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4gdml/PHG4GDMLConfig.hh>
#include <g4gdml/PHG4GDMLUtility.hh>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>  // for G4Transform3D
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Types.hh>            // for G4double, G4int
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume

#include <TSystem.h>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>  // for pair, make_pair

class G4VSolid;
class PHCompositeNode;

using namespace std;

//_______________________________________________________________________
PHG4ForwardEcalDetector::PHG4ForwardEcalDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters* parameters, const std::string& dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4ForwardEcalDisplayAction*>(subsys->GetDisplayAction()))
  , m_Params(parameters)
  , m_GdmlConfig(PHG4GDMLUtility::GetOrMakeConfigNode(Node))
  , m_XRot(0.0)
  , m_YRot(0.0)
  , m_ZRot(0.0)
  , m_PlaceX(0.0 * mm)
  , m_PlaceY(0.0 * mm)
  , m_PlaceZ(3150.0 * mm)
  , m_dZ(170 * mm)
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_AbsorberActiveFlag(m_Params->get_int_param("absorberactive"))
  , m_Layer(0)
  , m_SuperDetector("NONE")
  , m_TowerLogicNamePrefix("hEcalTower")
{
  m_Params->Print();
  for (int i = 0; i < 3; i++)
  {
    m_TowerDx[i] = 30 * mm;
    m_TowerDy[i] = 30 * mm;
    m_TowerDz[i] = 170.0 * mm;
  }
  for (int i = 3; i < 7; i++)
  {
    m_TowerDx[i] = NAN;
    m_TowerDy[i] = NAN;
    m_TowerDz[i] = NAN;
  }
  m_RMin[0] = 110 * mm;
  m_RMax[0] = 2250 * mm;
  m_RMin[1] = 120 * mm;
  m_RMax[1] = 2460 * mm;
  assert(m_GdmlConfig);
}

//_______________________________________________________________________
int PHG4ForwardEcalDetector::IsInForwardEcal(G4VPhysicalVolume* volume) const
{
  G4LogicalVolume* mylogvol = volume->GetLogicalVolume();
  if (m_ActiveFlag)
  {
    if (m_ScintiLogicalVolSet.find(mylogvol) != m_ScintiLogicalVolSet.end())
    {
      return 1;
    }
  }
  if (m_AbsorberActiveFlag)
  {
    if (m_AbsorberLogicalVolSet.find(mylogvol) != m_AbsorberLogicalVolSet.end())
    {
      return -1;
    }
  }
  return 0;
}

//_______________________________________________________________________
void PHG4ForwardEcalDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  //if ( Verbosity() > 0 )
  {
    cout << "PHG4ForwardEcalDetector: Begin Construction" << endl;
  }

  /* Read parameters for detector construction and mappign from file */
  ParseParametersFromTable();

  /* Create the cone envelope = 'world volume' for the crystal calorimeter */
  G4Material* Air = G4Material::GetMaterial("G4_AIR");

  G4VSolid* ecal_envelope_solid = new G4Cons("hEcal_envelope_solid",
                                             m_RMin[0], m_RMax[0],
                                             m_RMin[1], m_RMax[1],
                                             m_dZ / 2.,
                                             0, 2 * M_PI);

  G4LogicalVolume* ecal_envelope_log = new G4LogicalVolume(ecal_envelope_solid, Air, G4String("hEcal_envelope"), 0, 0, 0);

  /* Define visualization attributes for envelope cone */
  GetDisplayAction()->AddVolume(ecal_envelope_log, "Envelope");

  /* Define rotation attributes for envelope cone */
  G4RotationMatrix ecal_rotm;
  ecal_rotm.rotateX(m_XRot);
  ecal_rotm.rotateY(m_YRot);
  ecal_rotm.rotateZ(m_ZRot);

  /* Place envelope cone in simulation */
  string name_envelope = m_TowerLogicNamePrefix + "_envelope";

  new G4PVPlacement(G4Transform3D(ecal_rotm, G4ThreeVector(m_PlaceX, m_PlaceY, m_PlaceZ)),
                    ecal_envelope_log, name_envelope, logicWorld, 0, false, OverlapCheck());

  /* Construct single calorimeter towers */
  G4LogicalVolume* singletower[7] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  typedef std::map<std::string, towerposition>::iterator it_type;
  for (it_type iterator = m_TowerPositionMap.begin(); iterator != m_TowerPositionMap.end(); ++iterator)
  {
    for (int i = 0; i < 7; i++)
    {
      if (iterator->second.type == i)
      {
        singletower[i] = ConstructTower(i);
      }
    }
  }

  if (Verbosity() > 1)
  {
    cout << singletower << endl;
  }
  /* Place calorimeter towers within envelope */
  PlaceTower(ecal_envelope_log, singletower);

  return;
}

//_______________________________________________________________________
G4LogicalVolume*
PHG4ForwardEcalDetector::ConstructTower(int type)
{
  if (Verbosity() > 0)
  {
    cout << "PHG4ForwardEcalDetector: Build logical volume for single tower, type = " << type << endl;
  }
  assert(type >= 0 && type <= 6);
  // This method allows construction of Type 0,1 tower (PbGl or PbW04).
  // Call a separate routine to generate Type 2 towers (PbSc)
  // Call a separate routine to generate Type 3-6 towers (E864 Pb-Scifi)

  if (type == 2) return ConstructTowerType2();
  if ((type == 3) || (type == 4) || (type == 5) || (type == 6)) return ConstructTowerType3_4_5_6(type);

  /* create logical volume for single tower */
  G4Material* material_air = G4Material::GetMaterial("G4_AIR");

  G4Material* material_scintillator;
  double tower_dx = m_TowerDx[type];
  double tower_dy = m_TowerDy[type];
  double tower_dz = m_TowerDz[type];
  cout << "building type " << type << " towers" << endl;
  if (type == 0)
  {
    material_scintillator = G4Material::GetMaterial("G4_LEAD_OXIDE");
  }
  else if (type == 1)
  {
    material_scintillator = G4Material::GetMaterial("G4_PbWO4");
  }
  else
  {
    cout << "PHG4ForwardEcalDetector::ConstructTower invalid type = " << type << endl;
    material_scintillator = nullptr;
  }

  ostringstream single_tower_solid_name;
  single_tower_solid_name << m_TowerLogicNamePrefix << "_single_scintillator_type" << type;

  G4VSolid* single_tower_solid = new G4Box(G4String(single_tower_solid_name.str().c_str()),
                                           tower_dx / 2.0,
                                           tower_dy / 2.0,
                                           tower_dz / 2.0);

  ostringstream single_tower_logic_name;
  single_tower_logic_name << "single_tower_logic_type" << type;

  G4LogicalVolume* single_tower_logic = new G4LogicalVolume(single_tower_solid,
                                                            material_air,
                                                            single_tower_logic_name.str().c_str(),
                                                            0, 0, 0);

  ostringstream single_scintillator_name;
  single_scintillator_name << "single_scintillator_type" << type;

  G4VSolid* solid_scintillator = new G4Box(G4String(single_scintillator_name.str().c_str()),
                                           tower_dx / 2.0,
                                           tower_dy / 2.0,
                                           tower_dz / 2.0);

  ostringstream hEcal_scintillator_plate_logic_name;
  hEcal_scintillator_plate_logic_name << "hEcal_scintillator_plate_logic_type" << type;

  G4LogicalVolume* logic_scint = new G4LogicalVolume(solid_scintillator,
                                                     material_scintillator,
                                                     hEcal_scintillator_plate_logic_name.str(),
                                                     0, 0, 0);

  GetDisplayAction()->AddVolume(logic_scint, "Scintillator");

  /* place physical volumes for scintillator */

  string name_scintillator = m_TowerLogicNamePrefix + "_single_plate_scintillator";

  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0),
                    logic_scint,
                    name_scintillator,
                    single_tower_logic,
                    0, 0, OverlapCheck());

  GetDisplayAction()->AddVolume(single_tower_logic, "ScintillatorSingleTower");

  if (Verbosity() > 0)
  {
    cout << "PHG4ForwardEcalDetector: Building logical volume for single tower done, type = " << type << endl;
  }

  return single_tower_logic;
}

G4LogicalVolume*
PHG4ForwardEcalDetector::ConstructTowerType2()
{
  if (Verbosity() > 0)
  {
    cout << "PHG4ForwardEcalDetector: Build logical volume for single tower type 2..." << endl;
  }
  /* create logical volume for single tower */
  G4Material* material_air = G4Material::GetMaterial("G4_AIR");

  G4VSolid* single_tower_solid = new G4Box(G4String("single_tower_solid2"),
                                           m_TowerDx[2] / 2.0,
                                           m_TowerDy[2] / 2.0,
                                           m_TowerDz[2] / 2.0);

  G4LogicalVolume* single_tower_logic = new G4LogicalVolume(single_tower_solid,
                                                            material_air,
                                                            "single_tower_logic2",
                                                            0, 0, 0);

  /* create geometry volumes for scintillator and absorber plates to place inside single_tower */
  // PHENIX EMCal JGL 3/27/2016
  G4int nlayers = 66;
  G4double thickness_layer = m_TowerDz[2] / (float) nlayers;
  G4double thickness_absorber = thickness_layer * 0.23;      // 27% absorber by length
  G4double thickness_scintillator = thickness_layer * 0.73;  // 73% scintillator by length

  G4VSolid* solid_absorber = new G4Box(G4String("single_plate_absorber_solid2"),
                                       m_TowerDx[2] / 2.0,
                                       m_TowerDy[2] / 2.0,
                                       thickness_absorber / 2.0);

  G4VSolid* solid_scintillator = new G4Box(G4String("single_plate_scintillator2"),
                                           m_TowerDx[2] / 2.0,
                                           m_TowerDy[2] / 2.0,
                                           thickness_scintillator / 2.0);

  /* create logical volumes for scintillator and absorber plates to place inside single_tower */
  G4Material* material_scintillator = G4Material::GetMaterial("G4_POLYSTYRENE");
  G4Material* material_absorber = G4Material::GetMaterial("G4_Pb");

  G4LogicalVolume* logic_absorber = new G4LogicalVolume(solid_absorber,
                                                        material_absorber,
                                                        "single_plate_absorber_logic2",
                                                        0, 0, 0);

  G4LogicalVolume* logic_scint = new G4LogicalVolume(solid_scintillator,
                                                     material_scintillator,
                                                     "hEcal_scintillator_plate_logic2",
                                                     0, 0, 0);

  m_AbsorberLogicalVolSet.insert(logic_absorber);
  m_ScintiLogicalVolSet.insert(logic_scint);
  GetDisplayAction()->AddVolume(logic_absorber, "Absorber");
  GetDisplayAction()->AddVolume(logic_scint, "Scintillator");

  /* place physical volumes for absorber and scintillator plates */
  G4double xpos_i = 0;
  G4double ypos_i = 0;
  G4double zpos_i = (-1 * m_TowerDz[2] / 2.0) + thickness_absorber / 2.0;

  string name_absorber = m_TowerLogicNamePrefix + "_single_plate_absorber2";
  string name_scintillator = m_TowerLogicNamePrefix + "_single_plate_scintillator2";
  for (int i = 1; i <= nlayers; i++)
  {
    new G4PVPlacement(0, G4ThreeVector(xpos_i, ypos_i, zpos_i),
                      logic_absorber,
                      name_absorber,
                      single_tower_logic,
                      0, 0, OverlapCheck());

    zpos_i += (thickness_absorber / 2. + thickness_scintillator / 2.);

    new G4PVPlacement(0, G4ThreeVector(xpos_i, ypos_i, zpos_i),
                      logic_scint,
                      name_scintillator,
                      single_tower_logic,
                      0, 0, OverlapCheck());

    zpos_i += (thickness_absorber / 2. + thickness_scintillator / 2.);
  }

  GetDisplayAction()->AddVolume(single_tower_logic, "ScintillatorSingleTower");

  if (Verbosity() > 0)
  {
    cout << "PHG4ForwardEcalDetector: Building logical volume for single tower done." << endl;
  }

  return single_tower_logic;
}

G4LogicalVolume*
PHG4ForwardEcalDetector::ConstructTowerType3_4_5_6(int type)
{
  if (Verbosity() > 0)
  {
    cout << "PHG4ForwardEcalDetector: Build logical volume for single tower type ..." << type << endl;
  }

  double tower_dx, tower_dy, tower_dz;
  int num_fibers_x, num_fibers_y;
  tower_dx = m_TowerDx[type];
  tower_dy = m_TowerDy[type];
  tower_dz = m_TowerDz[type];
  switch (type)
  {
  case 3:
    num_fibers_x = 10;
    num_fibers_y = 10;
    break;
  case 4:
    num_fibers_x = 9;
    num_fibers_y = 10;
    break;
  case 5:
    num_fibers_x = 10;
    num_fibers_y = 9;
    break;
  case 6:
    num_fibers_x = 9;
    num_fibers_y = 9;
    break;
  default:
    cout << "PHG4ForwardEcalDetector: Invalid tower type in ConstructTowerType3_4_5_6, stopping..." << endl;
    return nullptr;
  }

  /* create logical volume for single tower */
  G4Material* material_air = G4Material::GetMaterial("G4_AIR");

  G4String solidName = "single_tower_solid";
  solidName += type;
  G4VSolid* single_tower_solid = new G4Box(solidName,
                                           tower_dx / 2.0,
                                           tower_dy / 2.0,
                                           tower_dz / 2.0);

  ostringstream name_single_tower_logic;
  name_single_tower_logic << "single_tower_logic" << type;

  G4LogicalVolume* single_tower_logic = new G4LogicalVolume(single_tower_solid,
                                                            material_air,
                                                            name_single_tower_logic.str(),
                                                            0, 0, 0);

  // Now the absorber and then the fibers:

  G4String absorberName = "single_absorber_solid";
  absorberName += type;
  G4VSolid* single_absorber_solid = new G4Box(absorberName,
                                              tower_dx / 2.0,
                                              tower_dy / 2.0,
                                              tower_dz / 2.0);

  G4String absorberLogicName = "single_absorber_logic";
  absorberLogicName += type;
  // E864 Pb-Scifi calorimeter
  // E864 Calorimeter is 99% Pb, 1% Antimony
  G4LogicalVolume* single_absorber_logic = new G4LogicalVolume(single_absorber_solid,
                                                               G4Material::GetMaterial("E864_Absorber"),

                                                               absorberLogicName,
                                                               0, 0, 0);

  /* create geometry volumes for scintillator and place inside single_tower */
  // 1.1mm fibers

  G4String fiberName = "single_fiber_scintillator_solid";
  fiberName += type;
  G4VSolid* single_scintillator_solid = new G4Tubs(fiberName,
                                                   0.0, 0.055 * cm, (tower_dz / 2.0), 0.0, CLHEP::twopi);

  /* create logical volumes for scintillator and absorber plates to place inside single_tower */
  G4Material* material_scintillator = G4Material::GetMaterial("G4_POLYSTYRENE");

  G4String fiberLogicName = "hEcal_scintillator_fiber_logic";
  fiberLogicName += type;
  G4LogicalVolume* single_scintillator_logic = new G4LogicalVolume(single_scintillator_solid,
                                                                   material_scintillator,
                                                                   fiberLogicName,
                                                                   0, 0, 0);
  m_AbsorberLogicalVolSet.insert(single_absorber_logic);
  m_ScintiLogicalVolSet.insert(single_scintillator_logic);
  GetDisplayAction()->AddVolume(single_absorber_logic, "Absorber");
  GetDisplayAction()->AddVolume(single_scintillator_logic, "Fiber");

  // place array of fibers inside absorber

  G4double fiber_unit_cell = 10.0 * cm / 47.0;
  G4double xpos_i = -(tower_dx / 2.0) + (fiber_unit_cell / 2.0);
  G4double ypos_i = -(tower_dy / 2.0) + (fiber_unit_cell / 2.0);
  G4double zpos_i = 0.0;

  ostringstream name_scintillator;
  name_scintillator << m_TowerLogicNamePrefix << "_single_fiber_scintillator" << type;

  for (int i = 0; i < num_fibers_x; i++)
  {
    for (int j = 0; j < num_fibers_y; j++)
    {
      new G4PVPlacement(0, G4ThreeVector(xpos_i + i * fiber_unit_cell, ypos_i + j * fiber_unit_cell, zpos_i),
                        single_scintillator_logic,
                        name_scintillator.str(),
                        single_absorber_logic,
                        0, 0, OverlapCheck());
    }
  }

  // Place the absorber inside the envelope

  ostringstream name_absorber;
  name_absorber << m_TowerLogicNamePrefix << "_single_absorber" << type;

  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0),
                    single_absorber_logic,
                    name_absorber.str().c_str(),
                    single_tower_logic,
                    0, 0, OverlapCheck());
  GetDisplayAction()->AddVolume(single_tower_logic, "ScintillatorSingleTower");

  if (Verbosity() > 0)
  {
    cout << "PHG4ForwardEcalDetector: Building logical volume for single tower done." << endl;
  }

  return single_tower_logic;
}

int PHG4ForwardEcalDetector::PlaceTower(G4LogicalVolume* ecalenvelope, G4LogicalVolume* singletowerIn[7])
{
  /* Loop over all tower positions in vector and place tower */
  typedef std::map<std::string, towerposition>::iterator it_type;

  for (it_type iterator = m_TowerPositionMap.begin(); iterator != m_TowerPositionMap.end(); ++iterator)
  {
    if (Verbosity() > 0)
    {
      cout << "PHG4ForwardEcalDetector: Place tower " << iterator->first
           << " at x = " << iterator->second.x << " , y = " << iterator->second.y << " , z = " << iterator->second.z << endl;
    }

    assert(iterator->second.type >= 0 && iterator->second.type <= 6);
    G4LogicalVolume* singletower = singletowerIn[iterator->second.type];

    G4PVPlacement* tower_placement =
        new G4PVPlacement(0, G4ThreeVector(iterator->second.x, iterator->second.y, iterator->second.z),
                          singletower,
                          iterator->first,
                          ecalenvelope,
                          0, 0, OverlapCheck());

    m_GdmlConfig->exclude_physical_vol(tower_placement);
  }

  return 0;
}

int PHG4ForwardEcalDetector::ParseParametersFromTable()
{
  /* Open the datafile, if it won't open return an error */
  ifstream istream_mapping;
  istream_mapping.open(m_Params->get_string_param("mapping_file"));
  if (!istream_mapping.is_open())
  {
    cout << "ERROR in PHG4ForwardEcalDetector: Failed to open mapping file " << m_Params->get_string_param("mapping_file") << endl;
    gSystem->Exit(1);
  }

  /* loop over lines in file */
  string line_mapping;
  while (getline(istream_mapping, line_mapping))
  {
    /* Skip lines starting with / including a '#' */
    if (line_mapping.find("#") != string::npos)
    {
      if (Verbosity() > 0)
      {
        cout << "PHG4ForwardEcalDetector: SKIPPING line in mapping file: " << line_mapping << endl;
      }
      continue;
    }

    istringstream iss(line_mapping);

    /* If line starts with keyword Tower, add to tower positions */
    if (line_mapping.find("Tower ") != string::npos)
    {
      unsigned idx_j, idx_k, idx_l;
      G4double pos_x, pos_y, pos_z;
      G4double size_x, size_y, size_z;
      G4double rot_x, rot_y, rot_z;
      int type;
      string dummys;

      /* read string- break if error */
      if (!(iss >> dummys >> type >> idx_j >> idx_k >> idx_l >> pos_x >> pos_y >> pos_z >> size_x >> size_y >> size_z >> rot_x >> rot_y >> rot_z))
      {
        cout << "ERROR in PHG4ForwardEcalDetector: Failed to read line in mapping file " << m_Params->get_string_param("mapping_file") << endl;
        gSystem->Exit(1);
      }

      /* Construct unique name for tower */
      /* Mapping file uses cm, this class uses mm for length */
      ostringstream towername;
      towername << m_TowerLogicNamePrefix << "_t_" << type << "_j_" << idx_j << "_k_" << idx_k;
      /* Add Geant4 units */
      pos_x = pos_x * cm;
      pos_y = pos_y * cm;
      pos_z = pos_z * cm;

      /* insert tower into tower map */
      towerposition tower_new;
      tower_new.x = pos_x;
      tower_new.y = pos_y;
      tower_new.z = pos_z;
      tower_new.type = type;
      m_TowerPositionMap.insert(make_pair(towername.str(), tower_new));
    }
    else
    {
      /* If this line is not a comment and not a tower, save parameter as string / value. */
      string parname;
      double parval;

      /* read string- break if error */
      if (!(iss >> parname >> parval))
      {
        cerr << "ERROR in PHG4ForwardEcalDetector: Failed to read line in mapping file " << m_Params->get_string_param("mapping_file") << endl;
        gSystem->Exit(1);
      }

      m_GlobalParameterMap.insert(make_pair(parname, parval));
    }
  }
  /* Update member variables for global parameters based on parsed parameter file */
  std::map<string, double>::iterator parit;
  ostringstream twr;
  for (int i = 0; i < 7; i++)
  {
    twr.str("");
    twr << "Gtower" << i << "_dx";
    parit = m_GlobalParameterMap.find(twr.str());
    m_TowerDx[i] = parit->second * cm;
    twr.str("");
    twr << "Gtower" << i << "_dy";
    parit = m_GlobalParameterMap.find(twr.str());
    m_TowerDy[i] = parit->second * cm;
    twr.str("");
    twr << "Gtower" << i << "_dz";
    parit = m_GlobalParameterMap.find(twr.str());
    m_TowerDz[i] = parit->second * cm;
  }
  ostringstream rad;
  for (int i = 0; i < 2; i++)
  {
    int index = i + 1;
    rad.str("");
    rad << "Gr" << index << "_inner";
    parit = m_GlobalParameterMap.find(rad.str());
    if (parit != m_GlobalParameterMap.end())
    {
      m_RMin[i] = parit->second * cm;
    }
    rad.str("");
    rad << "Gr" << index << "_outer";
    parit = m_GlobalParameterMap.find(rad.str());
    if (parit != m_GlobalParameterMap.end())
    {
      m_RMax[i] = parit->second * cm;
    }
  }
  parit = m_GlobalParameterMap.find("Gdz");
  if (parit != m_GlobalParameterMap.end())
  {
    m_dZ = parit->second * cm;
  }
  parit = m_GlobalParameterMap.find("Gx0");
  if (parit != m_GlobalParameterMap.end())
  {
    m_PlaceX = parit->second * cm;
  }
  parit = m_GlobalParameterMap.find("Gy0");
  if (parit != m_GlobalParameterMap.end())
  {
    m_PlaceY = parit->second * cm;
  }
  parit = m_GlobalParameterMap.find("Gz0");
  if (parit != m_GlobalParameterMap.end())
  {
    m_PlaceZ = parit->second * cm;
  }
  parit = m_GlobalParameterMap.find("Grot_x");
  if (parit != m_GlobalParameterMap.end())
  {
    m_XRot = parit->second;
  }
  parit = m_GlobalParameterMap.find("Grot_y");
  if (parit != m_GlobalParameterMap.end())
  {
    m_YRot = parit->second;
  }
  parit = m_GlobalParameterMap.find("Grot_z");
  if (parit != m_GlobalParameterMap.end())
  {
    m_ZRot = parit->second;
  }
  return 0;
}

void PHG4ForwardEcalDetector::SetTowerDimensions(double dx, double dy, double dz, int type)
{
  assert(type >= 0 && type <= 6);
  m_TowerDx[type] = dx;
  m_TowerDy[type] = dy;
  m_TowerDz[type] = dz;
}
