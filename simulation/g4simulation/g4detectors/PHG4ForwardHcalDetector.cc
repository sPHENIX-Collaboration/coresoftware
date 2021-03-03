#include "PHG4ForwardHcalDetector.h"
#include "PHG4ForwardHcalDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <phool/recoConsts.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>      // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>      // for G4Transform3D
#include <Geant4/G4Types.hh>            // for G4double, G4int
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume

#include <TSystem.h>

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
PHG4ForwardHcalDetector::PHG4ForwardHcalDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters* parameters, const std::string& dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4ForwardHcalDisplayAction*>(subsys->GetDisplayAction()))
  , m_Params(parameters)
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_AbsorberActiveFlag(m_Params->get_int_param("absorberactive"))
  , m_TowerLogicNamePrefix("hHcalTower")
  , m_SuperDetector("NONE")
{
}
//_______________________________________________________________________
int PHG4ForwardHcalDetector::IsInForwardHcal(G4VPhysicalVolume* volume) const
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
void PHG4ForwardHcalDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  if (Verbosity() > 0)
  {
    std::cout << "PHG4ForwardHcalDetector: Begin Construction" << std::endl;
  }

  if (m_Params->get_string_param("mapping_file").empty())
  {
    cout << "ERROR in PHG4ForwardHcalDetector: No mapping file specified. Abort detector construction." << endl;
    cout << "Please run set_string_param(\"mapping_file\", std::string filename ) first." << endl;
    gSystem->Exit(1);
  }

  /* Read parameters for detector construction and mappign from file */
  ParseParametersFromTable();

  /* Create the cone envelope = 'world volume' for the crystal calorimeter */
  recoConsts* rc = recoConsts::instance();
  G4Material* WorldMaterial = G4Material::GetMaterial(rc->get_StringFlag("WorldMaterial"));

  G4VSolid* hcal_envelope_solid = new G4Cons("hHcal_envelope_solid",
                                             m_Params->get_double_param("rMin1") * cm,
                                             m_Params->get_double_param("rMax1") * cm,
                                             m_Params->get_double_param("rMin2") * cm,
                                             m_Params->get_double_param("rMax2") * cm,
                                             m_Params->get_double_param("dz") * cm / 2.,
                                             0., 2. * M_PI);

  G4LogicalVolume* hcal_envelope_log = new G4LogicalVolume(hcal_envelope_solid, WorldMaterial, "hHcal_envelope", 0, 0, 0);

  m_DisplayAction->AddVolume(hcal_envelope_log, "FHcalEnvelope");

  /* Define rotation attributes for envelope cone */
  G4RotationMatrix hcal_rotm;
  hcal_rotm.rotateX(m_Params->get_double_param("rot_x") * deg);
  hcal_rotm.rotateY(m_Params->get_double_param("rot_y") * deg);
  hcal_rotm.rotateZ(m_Params->get_double_param("rot_z") * deg);

  /* Place envelope cone in simulation */
  string name_envelope = m_TowerLogicNamePrefix + "_envelope";

  new G4PVPlacement(G4Transform3D(hcal_rotm, G4ThreeVector(m_Params->get_double_param("place_x") * cm,
                                                           m_Params->get_double_param("place_y") * cm,
                                                           m_Params->get_double_param("place_z") * cm)),
                    hcal_envelope_log, name_envelope, logicWorld, 0, false, OverlapCheck());

  /* Construct single calorimeter tower */
  G4LogicalVolume* singletower = ConstructTower();

  /* Place calorimeter tower within envelope */
  PlaceTower(hcal_envelope_log, singletower);

  return;
}

//_______________________________________________________________________
G4LogicalVolume*
PHG4ForwardHcalDetector::ConstructTower()
{
  if (Verbosity() > 0)
  {
    std::cout << "PHG4ForwardHcalDetector: Build logical volume for single tower..." << std::endl;
  }

  /* create logical volume for single tower */
  recoConsts* rc = recoConsts::instance();
  G4Material* WorldMaterial = G4Material::GetMaterial(rc->get_StringFlag("WorldMaterial"));
  double TowerDx = m_Params->get_double_param("tower_dx") * cm;
  double TowerDy = m_Params->get_double_param("tower_dy") * cm;
  double TowerDz = m_Params->get_double_param("tower_dz") * cm;
  double WlsDw = m_Params->get_double_param("wls_dw") * cm;
  double SupportDw = m_Params->get_double_param("support_dw") * cm;
  G4VSolid* single_tower_solid = new G4Box("single_tower_solid",
                                           TowerDx / 2.0,
                                           TowerDy / 2.0,
                                           TowerDz / 2.0);

  G4LogicalVolume* single_tower_logic = new G4LogicalVolume(single_tower_solid,
                                                            WorldMaterial,
                                                            "single_tower_logic",
                                                            0, 0, 0);

  /* create geometry volumes for scintillator and absorber plates to place inside single_tower */
  // based on STAR forward upgrade design: https://drupal.star.bnl.gov/STAR/files/ForwardUpgrade.v20.pdf
  G4double thickness_absorber = m_Params->get_double_param("thickness_absorber") * cm;
  G4double thickness_scintillator = m_Params->get_double_param("thickness_scintillator") * cm;
  G4int nlayers = TowerDz / (thickness_absorber + thickness_scintillator);

  G4VSolid* solid_absorber = new G4Box("single_plate_absorber_solid",
                                       (TowerDx - WlsDw) / 2.0,
                                       (TowerDy - SupportDw) / 2.0,
                                       thickness_absorber / 2.0);

  G4VSolid* solid_scintillator = new G4Box("single_plate_scintillator",
                                           (TowerDx - WlsDw) / 2.0,
                                           (TowerDy - SupportDw) / 2.0,
                                           thickness_scintillator / 2.0);

  G4VSolid* solid_WLS_plate = new G4Box("single_plate_wls",
                                        (WlsDw) / 2.0,
                                        (TowerDy - SupportDw) / 2.0,
                                        TowerDz / 2.0);

  G4VSolid* solid_support_plate = new G4Box("single_plate_support",
                                            (TowerDx) / 2.0,
                                            (SupportDw) / 2.0,
                                            TowerDz / 2.0);

  /* create logical volumes for scintillator and absorber plates to place inside single_tower */
  G4Material* material_scintillator = G4Material::GetMaterial(m_Params->get_string_param("scintillator"));
  G4Material* material_absorber = G4Material::GetMaterial(m_Params->get_string_param("absorber"));
  G4Material* material_wls = G4Material::GetMaterial(m_Params->get_string_param("scintillator"));
  G4Material* material_support = G4Material::GetMaterial(m_Params->get_string_param("support"));

  G4LogicalVolume* logic_absorber = new G4LogicalVolume(solid_absorber,
                                                        material_absorber,
                                                        "single_plate_absorber_logic",
                                                        0, 0, 0);

  m_AbsorberLogicalVolSet.insert(logic_absorber);

  G4LogicalVolume* logic_scint = new G4LogicalVolume(solid_scintillator,
                                                     material_scintillator,
                                                     "hHcal_scintillator_plate_logic",
                                                     0, 0, 0);
  m_ScintiLogicalVolSet.insert(logic_scint);

  G4LogicalVolume* logic_wls = new G4LogicalVolume(solid_WLS_plate,
                                                   material_wls,
                                                   "hHcal_wls_plate_logic",
                                                   0, 0, 0);

  m_AbsorberLogicalVolSet.insert(logic_wls);
  G4LogicalVolume* logic_support = new G4LogicalVolume(solid_support_plate,
                                                       material_support,
                                                       "hHcal_support_plate_logic",
                                                       0, 0, 0);

  m_AbsorberLogicalVolSet.insert(logic_support);
  m_DisplayAction->AddVolume(logic_absorber, "Absorber");
  m_DisplayAction->AddVolume(logic_scint, "Scintillator");
  m_DisplayAction->AddVolume(logic_wls, "WLSplate");
  m_DisplayAction->AddVolume(logic_support, "SupportPlate");

  /* place physical volumes for absorber and scintillator plates */
  G4double xpos_i = -WlsDw / 2.0;
  G4double ypos_i = -SupportDw / 2.0;
  G4double zpos_i = (-1 * TowerDz / 2.0) + thickness_absorber / 2.0;

  string name_absorber = m_TowerLogicNamePrefix + "_single_plate_absorber";

  string name_scintillator = m_TowerLogicNamePrefix + "_single_plate_scintillator";

  string name_wls = m_TowerLogicNamePrefix + "_single_plate_wls";

  string name_support = m_TowerLogicNamePrefix + "_single_plate_support";

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
  new G4PVPlacement(0, G4ThreeVector(0, (TowerDy / 2) - SupportDw / 2, 0),
                    logic_support,
                    name_support,
                    single_tower_logic,
                    0, 0, OverlapCheck());

  new G4PVPlacement(0, G4ThreeVector((TowerDx / 2) - WlsDw / 2, -SupportDw / 2, 0),
                    logic_wls,
                    name_wls,
                    single_tower_logic,
                    0, 0, OverlapCheck());
  m_DisplayAction->AddVolume(single_tower_logic, "SingleScintillator");

  if (Verbosity() > 0)
  {
    std::cout << "PHG4ForwardHcalDetector: Building logical volume for single tower done." << std::endl;
  }

  return single_tower_logic;
}

int PHG4ForwardHcalDetector::PlaceTower(G4LogicalVolume* hcalenvelope, G4LogicalVolume* singletower)
{
  /* Loop over all tower positions in vector and place tower */
  for (std::map<std::string, towerposition>::iterator iterator = m_TowerPostionMap.begin(); iterator != m_TowerPostionMap.end(); ++iterator)
  {
    if (Verbosity() > 0)
    {
      std::cout << "PHG4ForwardHcalDetector: Place tower " << iterator->first
                << " idx_j = " << iterator->second.idx_j << ", idx_k = " << iterator->second.idx_k
                << " at x = " << iterator->second.x << " , y = " << iterator->second.y << " , z = " << iterator->second.z << std::endl;
    }

    int copyno = (iterator->second.idx_j << 16) + iterator->second.idx_k;
    new G4PVPlacement(0, G4ThreeVector(iterator->second.x, iterator->second.y, iterator->second.z),
                      singletower,
                      iterator->first,
                      hcalenvelope,
                      0, copyno, OverlapCheck());
  }

  return 0;
}

int PHG4ForwardHcalDetector::ParseParametersFromTable()
{
  /* Open the datafile, if it won't open return an error */
  ifstream istream_mapping;
  istream_mapping.open(m_Params->get_string_param("mapping_file"));
  if (!istream_mapping.is_open())
  {
    std::cout << "ERROR in PHG4ForwardHcalDetector: Failed to open mapping file " << m_Params->get_string_param("mapping_file") << std::endl;
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
        std::cout << "PHG4ForwardHcalDetector: SKIPPING line in mapping file: " << line_mapping << std::endl;
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
        cout << "ERROR in PHG4ForwardHcalDetector: Failed to read line in mapping file " << m_Params->get_string_param("mapping_file") << std::endl;
        gSystem->Exit(1);
      }

      /* Construct unique name for tower */
      /* Mapping file uses cm, this class uses mm for length */
      ostringstream towername;
      towername.str("");
      towername << m_TowerLogicNamePrefix << "_j_" << idx_j << "_k_" << idx_k;

      /* Add Geant4 units */
      pos_x = pos_x * cm;
      pos_y = pos_y * cm;
      pos_z = pos_z * cm;

      /* insert tower into tower map */
      towerposition tower_new;
      tower_new.x = pos_x;
      tower_new.y = pos_y;
      tower_new.z = pos_z;
      tower_new.idx_j = idx_j;
      tower_new.idx_k = idx_k;
      m_TowerPostionMap.insert(make_pair(towername.str(), tower_new));
    }
    else
    {
      /* If this line is not a comment and not a tower, save parameter as string / value. */
      string parname;
      G4double parval;

      /* read string- break if error */
      if (!(iss >> parname >> parval))
      {
        cout << "ERROR in PHG4ForwardHcalDetector: Failed to read line in mapping file " << m_Params->get_string_param("mapping_file") << std::endl;
        gSystem->Exit(1);
      }

      m_GlobalParameterMap.insert(make_pair(parname, parval));
    }
  }

  /* Update member variables for global parameters based on parsed parameter file */
  std::map<string, G4double>::iterator parit;

  parit = m_GlobalParameterMap.find("Gtower_dx");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("tower_dx", parit->second);  // in cm
  }

  parit = m_GlobalParameterMap.find("Gtower_dy");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("tower_dy", parit->second);  // in cm
  }

  parit = m_GlobalParameterMap.find("Gtower_dz");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("tower_dz", parit->second);  // in cm
  }

  parit = m_GlobalParameterMap.find("Gr1_inner");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("rMin1", parit->second);
  }

  parit = m_GlobalParameterMap.find("Gr1_outer");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("rMax1", parit->second);
  }

  parit = m_GlobalParameterMap.find("Gr2_inner");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("rMin2", parit->second);
  }

  parit = m_GlobalParameterMap.find("Gr2_outer");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("rMax2", parit->second);
  }

  parit = m_GlobalParameterMap.find("Gdz");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("dZ", parit->second);
  }

  parit = m_GlobalParameterMap.find("Gx0");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("place_x", parit->second);
  }

  parit = m_GlobalParameterMap.find("Gy0");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("place_y", parit->second);
  }

  parit = m_GlobalParameterMap.find("Gz0");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("place_z", parit->second);
  }

  parit = m_GlobalParameterMap.find("Grot_x");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("rot_x", parit->second * rad / deg);
  }

  parit = m_GlobalParameterMap.find("Grot_y");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("rot_y", parit->second * rad / deg);
  }

  parit = m_GlobalParameterMap.find("Grot_z");
  if (parit != m_GlobalParameterMap.end())
  {
    m_Params->set_double_param("rot_z", parit->second * rad / deg);
  }

  return 0;
}
