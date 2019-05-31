#include "PHG4CrystalCalorimeterDetector.h"
#include "PHG4CrystalCalorimeterDisplayAction.h"
#include "PHG4CrystalCalorimeterSubsystem.h"

#include <g4main/PHG4Detector.h>                  // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>             // for PHG4DisplayAction

#include <Geant4/G4Box.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4Element.hh>                    // for G4Element
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>             // for G4RotationMatrix
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>                // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>                // for G4Transform3D
#include <Geant4/G4Types.hh>                      // for G4double
#include <Geant4/G4VPhysicalVolume.hh>            // for G4VPhysicalVolume

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <utility>                                // for pair, make_pair

class PHCompositeNode;

using namespace std;

//_______________________________________________________________________
PHG4CrystalCalorimeterDetector::PHG4CrystalCalorimeterDetector(PHG4CrystalCalorimeterSubsystem* subsys, PHCompositeNode* Node, const std::string& dnam)
  : PHG4Detector(Node, dnam)
  , _place_in_x(0.0 * mm)
  , _place_in_y(0.0 * mm)
  , _place_in_z(-1080.0 * mm)
  , _rot_in_x(0.0)
  , _rot_in_y(M_PI)
  , _rot_in_z(0.0)
  , _rMin1(22 * mm)
  , _rMax1(656 * mm)
  , _rMin2(26 * mm)
  , _rMax2(775 * mm)
  , _dZ(180 * mm)
  , _sPhi(0)
  , _dPhi(2 * M_PI)
  , _crystal_dx(20 * mm)
  , _crystal_dy(20 * mm)
  , _crystal_dz(180.0 * mm)
  , _materialCrystal("G4_PbWO4")
  , _active(1)
  , _absorberactive(0)
  , _layer(0)
  , _blackhole(0)
  , _superdetector("NONE")
  , _mapping_tower_file("")
  , m_DisplayAction(dynamic_cast<PHG4CrystalCalorimeterDisplayAction*>(subsys->GetDisplayAction()))
  , _towerlogicnameprefix("CrystalCalorimeterTower")
{
}

//_______________________________________________________________________
int PHG4CrystalCalorimeterDetector::IsInCrystalCalorimeter(G4VPhysicalVolume* volume) const
{
  if (volume->GetName().find(_towerlogicnameprefix) != string::npos)
  {
    if (volume->GetName().find("crystal") != string::npos)
    {
      return 1;
    }
    /* only record energy in actual absorber- drop energy lost in air gaps inside envelope */
    else if (volume->GetName().find("absorber") != string::npos)
    {
      return -1;
    }
    else if (volume->GetName().find("envelope") != string::npos)
    {
      return 0;
    }
  }

  return 0;
}

//_______________________________________________________________________
void PHG4CrystalCalorimeterDetector::Construct(G4LogicalVolume* logicWorld)
{
  if (Verbosity() > 0)
  {
    cout << "PHG4CrystalCalorimeterDetector: Begin Construction" << endl;
  }

  if (_mapping_tower_file.empty())
  {
    cout << "ERROR in PHG4CrystalCalorimeterDetector: No tower mapping file specified. Abort detector construction." << endl;
    cout << "Please run SetTowerMappingFile( std::string filename ) first." << endl;
    exit(1);
  }

  /* Read parameters for detector construction and mappign from file */
  ParseParametersFromTable();

  /* Create the cone envelope = 'world volume' for the crystal calorimeter */
  G4Material* Air = G4Material::GetMaterial("G4_AIR");

  G4VSolid* eemc_envelope_solid = new G4Cons("eemc_envelope_solid",
                                             _rMin1, _rMax1,
                                             _rMin2, _rMax2,
                                             _dZ / 2.,
                                             _sPhi, _dPhi);

  G4LogicalVolume* eemc_envelope_log = new G4LogicalVolume(eemc_envelope_solid, Air, G4String("eemc_envelope"), 0, 0, 0);

  GetDisplayAction()->AddVolume(eemc_envelope_log, "Envelope");
  /* Define rotation attributes for envelope cone */
  G4RotationMatrix eemc_rotm;
  eemc_rotm.rotateX(_rot_in_x);
  eemc_rotm.rotateY(_rot_in_y);
  eemc_rotm.rotateZ(_rot_in_z);

  /* Place envelope cone in simulation */
  ostringstream name_envelope;
  name_envelope.str("");
  name_envelope << _towerlogicnameprefix << "_envelope" << endl;

  new G4PVPlacement(G4Transform3D(eemc_rotm, G4ThreeVector(_place_in_x, _place_in_y, _place_in_z)),
                    eemc_envelope_log, name_envelope.str().c_str(), logicWorld, 0, false, OverlapCheck());

  /* Construct single calorimeter tower */
  G4LogicalVolume* singletower = ConstructTower();

  /* Place calorimeter tower within envelope */
  PlaceTower(eemc_envelope_log, singletower);

  return;
}

//_______________________________________________________________________
G4LogicalVolume*
PHG4CrystalCalorimeterDetector::ConstructTower()
{
  if (Verbosity() > 0)
  {
    cout << "PHG4CrystalCalorimeterDetector: Build logical volume for single tower..." << endl;
  }

  G4double carbon_thickness = 0.009 * cm;
  G4double airgap_crystal_carbon = 0.012 * cm;

  /* dimesnions of full tower */
  G4double tower_dx = _crystal_dx + 2 * (carbon_thickness + airgap_crystal_carbon);
  G4double tower_dy = _crystal_dy + 2 * (carbon_thickness + airgap_crystal_carbon);
  G4double tower_dz = _crystal_dz;

  /* create logical volume for single tower */
  G4Material* material_air = G4Material::GetMaterial("G4_AIR");

  G4VSolid* single_tower_solid = new G4Box(G4String("single_tower_solid"),
                                           tower_dx / 2.0,
                                           tower_dy / 2.0,
                                           tower_dz / 2.0);

  G4LogicalVolume* single_tower_logic = new G4LogicalVolume(single_tower_solid,
                                                            material_air,
                                                            "single_tower_logic",
                                                            0, 0, 0);

  /* create geometry volume for crystal inside single_tower */
  G4VSolid* solid_crystal = new G4Box(G4String("single_crystal_solid"),
                                      _crystal_dx / 2.0,
                                      _crystal_dy / 2.0,
                                      _crystal_dz / 2.0);

  /* create geometry volume for frame (carbon fober shell) inside single_tower */
  G4VSolid* Carbon_hunk_solid = new G4Box(G4String("Carbon_hunk_solid"),
                                          tower_dx / 2.0,
                                          tower_dy / 2.0,
                                          ((tower_dz / 2.0) - 1 * mm));

  G4double lead_dx = tower_dx - 2.0 * carbon_thickness;
  G4double lead_dy = tower_dy - 2.0 * carbon_thickness;

  G4VSolid* lead_solid = new G4Box(G4String("lead_solid"),
                                   lead_dx / 2.0,
                                   lead_dy / 2.0,
                                   tower_dz / 2.0);

  G4SubtractionSolid* Carbon_shell_solid = new G4SubtractionSolid(G4String("Carbon_Shell_solid"),
                                                                  Carbon_hunk_solid,
                                                                  lead_solid,
                                                                  0,
                                                                  G4ThreeVector(0.00 * mm, 0.00 * mm, 0.00 * mm));

  /* create logical volumes for crystal inside single_tower */
  G4Material* material_crystal = G4Material::GetMaterial(_materialCrystal.c_str());

  G4LogicalVolume* logic_crystal = new G4LogicalVolume(solid_crystal,
                                                       material_crystal,
                                                       "single_crystal_logic",
                                                       0, 0, 0);

  GetDisplayAction()->AddVolume(logic_crystal, "Crystal");

  /* create logical volumes for structural frame */
  //Carbon Fiber
  G4double a = 12.01 * g / mole;
  G4Element* elC = new G4Element("Carbon", "C", 6., a);
  G4double density_carbon_fiber = 0.144 * g / cm3;
  G4Material* CarbonFiber = new G4Material("CarbonFiber", density_carbon_fiber, 1);
  CarbonFiber->AddElement(elC, 1);

  G4Material* material_shell = CarbonFiber;

  G4LogicalVolume* logic_shell = new G4LogicalVolume(Carbon_shell_solid,
                                                     material_shell,
                                                     "single_absorber_logic",
                                                     0, 0, 0);

  GetDisplayAction()->AddVolume(logic_shell, "CarbonShell");

  /* Place structural frame in logical tower volume */
  ostringstream name_shell;
  name_shell.str("");
  name_shell << _towerlogicnameprefix << "_single_absorber" << endl;

  new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                    logic_shell,
                    name_shell.str().c_str(),
                    single_tower_logic,
                    0, 0, OverlapCheck());

  /* Place crystal in logical tower volume */
  ostringstream name_crystal;
  name_crystal.str("");
  name_crystal << _towerlogicnameprefix << "_single_crystal" << endl;

  new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                    logic_crystal,
                    name_crystal.str().c_str(),
                    single_tower_logic,
                    0, 0, OverlapCheck());

  if (Verbosity() > 0)
  {
    cout << "PHG4CrystalCalorimeterDetector: Building logical volume for single tower done." << endl;
  }

  return single_tower_logic;
}

int PHG4CrystalCalorimeterDetector::PlaceTower(G4LogicalVolume* eemcenvelope, G4LogicalVolume* singletower)
{
  /* Loop over all tower positions in vector and place tower */
  typedef std::map<std::string, towerposition>::iterator it_type;

  for (it_type iterator = _map_tower.begin(); iterator != _map_tower.end(); ++iterator)
  {
    if (Verbosity() > 0)
    {
      cout << "PHG4CrystalCalorimeterDetector: Place tower " << iterator->first
           << " at x = " << iterator->second.x << " , y = " << iterator->second.y << " , z = " << iterator->second.z << endl;
    }

    new G4PVPlacement(0, G4ThreeVector(iterator->second.x, iterator->second.y, iterator->second.z),
                      singletower,
                      iterator->first.c_str(),
                      eemcenvelope,
                      0, 0, OverlapCheck());
  }

  return 0;
}

int PHG4CrystalCalorimeterDetector::ParseParametersFromTable()
{
  /* Open the datafile, if it won't open return an error */
  ifstream istream_mapping;
  if (!istream_mapping.is_open())
  {
    istream_mapping.open(_mapping_tower_file.c_str());
    if (!istream_mapping)
    {
      cerr << "ERROR in PHG4CrystalCalorimeterDetector: Failed to open mapping file " << _mapping_tower_file << endl;
      exit(1);
    }
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
        cout << "PHG4CrystalCalorimeterDetector: SKIPPING line in mapping file: " << line_mapping << endl;
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
      G4double dummy;
      string dummys;

      /* read string- break if error */
      if (!(iss >> dummys >> dummy >> idx_j >> idx_k >> idx_l >> pos_x >> pos_y >> pos_z >> size_x >> size_y >> size_z >> rot_x >> rot_y >> rot_z))
      {
        cerr << "ERROR in PHG4CrystalCalorimeterDetector: Failed to read line in mapping file " << _mapping_tower_file << endl;
        exit(1);
      }

      /* Construct unique name for tower */
      /* Mapping file uses cm, this class uses mm for length */
      ostringstream towername;
      towername.str("");
      towername << _towerlogicnameprefix << "_j_" << idx_j << "_k_" << idx_k;

      /* Add Geant4 units */
      pos_x = pos_x * cm;
      pos_y = pos_y * cm;
      pos_z = pos_z * cm;

      /* insert tower into tower map */
      towerposition tower_new;
      tower_new.x = pos_x;
      tower_new.y = pos_y;
      tower_new.z = pos_z;
      _map_tower.insert(make_pair(towername.str(), tower_new));
    }
    else
    {
      /* If this line is not a comment and not a tower, save parameter as string / value. */
      string parname;
      G4double parval;

      /* read string- break if error */
      if (!(iss >> parname >> parval))
      {
        cerr << "ERROR in PHG4CrystalCalorimeterDetector: Failed to read line in mapping file " << _mapping_tower_file << endl;
        exit(1);
      }

      _map_global_parameter.insert(make_pair(parname, parval));
    }
  }

  /* Update member variables for global parameters based on parsed parameter file */
  std::map<string, G4double>::iterator parit;

  parit = _map_global_parameter.find("Gcrystal_dx");
  if (parit != _map_global_parameter.end())
    _crystal_dx = parit->second * cm;

  parit = _map_global_parameter.find("Gcrystal_dy");
  if (parit != _map_global_parameter.end())
    _crystal_dy = parit->second * cm;

  parit = _map_global_parameter.find("Gcrystal_dz");
  if (parit != _map_global_parameter.end())
    _crystal_dz = parit->second * cm;

  parit = _map_global_parameter.find("Gr1_inner");
  if (parit != _map_global_parameter.end())
    _rMin1 = parit->second * cm;

  parit = _map_global_parameter.find("Gr1_outer");
  if (parit != _map_global_parameter.end())
    _rMax1 = parit->second * cm;

  parit = _map_global_parameter.find("Gr2_inner");
  if (parit != _map_global_parameter.end())
    _rMin2 = parit->second * cm;

  parit = _map_global_parameter.find("Gr2_outer");
  if (parit != _map_global_parameter.end())
    _rMax2 = parit->second * cm;

  parit = _map_global_parameter.find("Gdz");
  if (parit != _map_global_parameter.end())
    _dZ = parit->second * cm;

  parit = _map_global_parameter.find("Gx0");
  if (parit != _map_global_parameter.end())
    _place_in_x = parit->second * cm;

  parit = _map_global_parameter.find("Gy0");
  if (parit != _map_global_parameter.end())
    _place_in_y = parit->second * cm;

  parit = _map_global_parameter.find("Gz0");
  if (parit != _map_global_parameter.end())
    _place_in_z = parit->second * cm;

  parit = _map_global_parameter.find("Grot_x");
  if (parit != _map_global_parameter.end())
    _rot_in_x = parit->second;

  parit = _map_global_parameter.find("Grot_y");
  if (parit != _map_global_parameter.end())
    _rot_in_y = parit->second;

  parit = _map_global_parameter.find("Grot_z");
  if (parit != _map_global_parameter.end())
    _rot_in_z = parit->second;

  return 0;
}
