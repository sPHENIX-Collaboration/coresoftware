#include "PHG4TpcDetector.h"
#include "PHG4TpcDefs.h"
#include "PHG4TpcDisplayAction.h"

#include <g4detectors/PHG4CellDefs.h> 
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>
#include <g4detectors/PHG4TpcCylinderGeom.h>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <phparameter/PHParameters.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>       
#include <phool/recoConsts.h>

#include <TSystem.h>

#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4String.hh>  // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4UserLimits.hh>
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume

#include <algorithm>  // for max, copy
#include <cassert>
#include <cmath>
#include <cstdlib>   // for exit
#include <iostream>  // for basic_ostream::operator<<
#include <map>       // for map
#include <sstream>

class G4VSolid;
class PHCompositeNode;

//_______________________________________________________________
PHG4TpcDetector::PHG4TpcDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4TpcDisplayAction *>(subsys->GetDisplayAction()))
  , m_Params(parameters)
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_AbsorberActiveFlag(m_Params->get_int_param("absorberactive"))
  , m_InnerCageRadius(m_Params->get_double_param("gas_inner_radius") * cm - m_Params->get_double_param("cage_layer_9_thickness") * cm - m_Params->get_double_param("cage_layer_8_thickness") * cm - m_Params->get_double_param("cage_layer_7_thickness") * cm - m_Params->get_double_param("cage_layer_6_thickness") * cm - m_Params->get_double_param("cage_layer_5_thickness") * cm - m_Params->get_double_param("cage_layer_4_thickness") * cm - m_Params->get_double_param("cage_layer_3_thickness") * cm - m_Params->get_double_param("cage_layer_2_thickness") * cm - m_Params->get_double_param("cage_layer_1_thickness") * cm)
  , m_OuterCageRadius(m_Params->get_double_param("gas_outer_radius") * cm + m_Params->get_double_param("cage_layer_9_thickness") * cm + m_Params->get_double_param("cage_layer_8_thickness") * cm + m_Params->get_double_param("cage_layer_7_thickness") * cm + m_Params->get_double_param("cage_layer_6_thickness") * cm + m_Params->get_double_param("cage_layer_5_thickness") * cm + m_Params->get_double_param("cage_layer_4_thickness") * cm + m_Params->get_double_param("cage_layer_3_thickness") * cm + m_Params->get_double_param("cage_layer_2_thickness") * cm + m_Params->get_double_param("cage_layer_1_thickness") * cm)
{
}

//_______________________________________________________________
int PHG4TpcDetector::IsInTpc(G4VPhysicalVolume *volume) const
{
  if (m_ActiveFlag)
  {
    if (m_ActiveVolumeSet.find(volume) != m_ActiveVolumeSet.end())
    {
      return 1;
    }
  }
  if (m_AbsorberActiveFlag)
  {
    if (m_AbsorberVolumeSet.find(volume) != m_AbsorberVolumeSet.end())
    {
      return -1;
    }
  }
  return 0;
}

//_______________________________________________________________
void PHG4TpcDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  // create Tpc envelope
  // tpc consists of (from inside to gas volume, outside is reversed up to now)
  // 1st layer cu
  // 2nd layer FR4
  // 3rd layer HoneyComb
  // 4th layer cu
  // 5th layer FR4
  // 6th layer Kapton
  // 7th layer cu
  // 8th layer Kapton
  // 9th layer cu

  double steplimits = m_Params->get_double_param("steplimits") * cm;
  if (std::isfinite(steplimits))
  {
    m_G4UserLimits = new G4UserLimits(steplimits);
  }

  G4VSolid *tpc_envelope = new G4Tubs("tpc_envelope", m_InnerCageRadius, m_OuterCageRadius, m_Params->get_double_param("tpc_length") * cm / 2., 0., 2 * M_PI);

  recoConsts *rc = recoConsts::instance();
  G4LogicalVolume *tpc_envelope_logic = new G4LogicalVolume(tpc_envelope,
                                                            GetDetectorMaterial(rc->get_StringFlag("WorldMaterial")),
                                                            "tpc_envelope");
  m_DisplayAction->AddVolume(tpc_envelope_logic, "TpcEnvelope");

  ConstructTpcExternalSupports(logicWorld);

  ConstructTpcCageVolume(tpc_envelope_logic);
  ConstructTpcGasVolume(tpc_envelope_logic);

  new G4PVPlacement(0, G4ThreeVector(m_Params->get_double_param("place_x") * cm, m_Params->get_double_param("place_y") * cm, m_Params->get_double_param("place_z") * cm),
    tpc_envelope_logic, "tpc_envelope",
    logicWorld, 0, false, OverlapCheck());
                  
  // geometry node
  add_geometry_node();
  
                  
}

int PHG4TpcDetector::ConstructTpcGasVolume(G4LogicalVolume *tpc_envelope)
{
  static std::map<int, std::string> tpcgasvolname =
      {{PHG4TpcDefs::North, "tpc_gas_north"},
       {PHG4TpcDefs::South, "tpc_gas_south"}};

  // Window / central membrane
  double tpc_window_thickness = m_Params->get_double_param("window_thickness") * cm;
  double tpc_half_length = (m_Params->get_double_param("tpc_length") * cm - tpc_window_thickness) / 2.;

  //'window' (modernly called central membrane only) material is ENIG, not Copper:
  //thickness in this recipe are just a ratio.  we set the usual thickness below.
  std::vector<double> thickness;
  std::vector<std::string> material;
  material.push_back("G4_Ni");
  thickness.push_back(.240 * cm);
  material.push_back("G4_Au");
  thickness.push_back(.008 * cm);
  G4Material *temp = nullptr;
  temp = GetDetectorMaterial("ENIG", false);
  if (temp == nullptr)
  {
    CreateCompositeMaterial("ENIG", material, thickness);  //see new function below
  }

  G4VSolid *tpc_window = new G4Tubs("tpc_window", m_Params->get_double_param("gas_inner_radius") * cm, m_Params->get_double_param("gas_outer_radius") * cm, tpc_window_thickness / 2., 0., 2 * M_PI);
  //we build our CM surface:
  G4LogicalVolume *tpc_window_logic = new G4LogicalVolume(tpc_window,
                                                          GetDetectorMaterial("ENIG"),
                                                          "tpc_window");
  //previously:                                                            GetDetectorMaterial(m_Params->get_string_param("window_surface1_material")),

  //  G4VisAttributes *visatt = new G4VisAttributes();
  //  visatt->SetVisibility(true);
  //  visatt->SetForceSolid(true);
  //  visatt->SetColor(PHG4TPCColorDefs::tpc_cu_color);
  //  tpc_window_logic->SetVisAttributes(visatt);

  m_DisplayAction->AddVolume(tpc_window_logic, "TpcWindow");
  G4VPhysicalVolume *tpc_window_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                                                         tpc_window_logic, "tpc_window",
                                                         tpc_envelope, false, PHG4TpcDefs::Window, OverlapCheck());

  m_AbsorberVolumeSet.insert(tpc_window_phys);

  //now build the FR4 layer beneath that:
  // Window / central membrane core
  double tpc_window_surface1_thickness = m_Params->get_double_param("window_surface1_thickness") * cm;
  double tpc_window_surface2_thickness = m_Params->get_double_param("window_surface2_thickness") * cm;
  double tpc_window_surface2_core_thickness = tpc_window_thickness - 2 * tpc_window_surface1_thickness;
  double tpc_window_core_thickness = tpc_window_surface2_core_thickness - 2 * (tpc_window_surface2_thickness);

  G4VSolid *tpc_window_surface2_core =
      new G4Tubs("tpc_window_surface2_core", m_Params->get_double_param("gas_inner_radius") * cm, m_Params->get_double_param("gas_outer_radius") * cm,
                 tpc_window_surface2_core_thickness / 2., 0., 2 * M_PI);
  G4LogicalVolume *tpc_window_surface2_core_logic = new G4LogicalVolume(tpc_window_surface2_core,
                                                                        GetDetectorMaterial(m_Params->get_string_param("window_surface2_material")),
                                                                        "tpc_window_surface2_core");

  m_DisplayAction->AddVolume(tpc_window_surface2_core_logic, "TpcWindow");
  //  visatt = new G4VisAttributes();
  //  visatt->SetVisibility(true);
  //  visatt->SetForceSolid(true);
  //  visatt->SetColor(PHG4TPCColorDefs::tpc_pcb_color);
  //  tpc_window_surface2_core_logic->SetVisAttributes(visatt);
  G4VPhysicalVolume *tpc_window_surface2_core_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                                                                       tpc_window_surface2_core_logic, "tpc_window_surface2_core",
                                                                       tpc_window_logic, false, PHG4TpcDefs::WindowCore, OverlapCheck());

  m_AbsorberVolumeSet.insert(tpc_window_surface2_core_phys);

  //and now the honeycomb core:
  G4VSolid *tpc_window_core =
      new G4Tubs("tpc_window", m_Params->get_double_param("gas_inner_radius") * cm, m_Params->get_double_param("gas_outer_radius") * cm,
                 tpc_window_core_thickness / 2., 0., 2 * M_PI);
  G4LogicalVolume *tpc_window_core_logic = new G4LogicalVolume(tpc_window_core,
                                                               GetDetectorMaterial(m_Params->get_string_param("window_core_material")),
                                                               "tpc_window_core");

  m_DisplayAction->AddVolume(tpc_window_core_logic, "TpcHoneyComb");
  G4VPhysicalVolume *tpc_window_core_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                                                              tpc_window_core_logic, "tpc_window_core",
                                                              tpc_window_surface2_core_logic, false, PHG4TpcDefs::WindowCore, OverlapCheck());

  m_AbsorberVolumeSet.insert(tpc_window_core_phys);

  // Gas
  G4VSolid *tpc_gas = new G4Tubs("tpc_gas", m_Params->get_double_param("gas_inner_radius") * cm, m_Params->get_double_param("gas_outer_radius") * cm, tpc_half_length / 2., 0., 2 * M_PI);

  G4LogicalVolume *tpc_gas_logic = new G4LogicalVolume(tpc_gas,
                                                       GetDetectorMaterial(m_Params->get_string_param("tpc_gas")),
                                                       "tpc_gas");

  tpc_gas_logic->SetUserLimits(m_G4UserLimits);
  m_DisplayAction->AddVolume(tpc_gas_logic, "TpcGas");
  G4VPhysicalVolume *tpc_gas_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, (tpc_half_length + tpc_window_thickness) / 2.),
                                                      tpc_gas_logic, tpcgasvolname[PHG4TpcDefs::North],
                                                      tpc_envelope, false, PHG4TpcDefs::North, OverlapCheck());
  std::cout << "north copy no: " << tpc_gas_phys->GetCopyNo() << std::endl;

  m_ActiveVolumeSet.insert(tpc_gas_phys);
  tpc_gas_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, -(tpc_half_length + tpc_window_thickness) / 2.),
                                   tpc_gas_logic, tpcgasvolname[PHG4TpcDefs::South],
                                   tpc_envelope, false, PHG4TpcDefs::South, OverlapCheck());

  std::cout << "south copy no: " << tpc_gas_phys->GetCopyNo() << std::endl;
  m_ActiveVolumeSet.insert(tpc_gas_phys);

#if G4VERSION_NUMBER >= 1033
  const G4RegionStore *theRegionStore = G4RegionStore::GetInstance();
  G4Region *tpcregion = theRegionStore->GetRegion("REGION_TPCGAS");
  tpc_gas_logic->SetRegion(tpcregion);
  tpcregion->AddRootLogicalVolume(tpc_gas_logic);
#endif
  return 0;
}

int PHG4TpcDetector::ConstructTpcExternalSupports(G4LogicalVolume *logicWorld)
{
  //note that these elements are outside the tpc logical volume!

  // Two two-inch diam. 304 Stainless Steel solid 'hanger beams' at 32.05" from beam center
  // at +/- 41.39 degrees left and right of vertical
  //stainless steel: 0.695 iron, 0.190 chromium, 0.095 nickel, 0.020 manganese.
  //if we're being pedantic, that is.  But store-bought stainless is probably okay.
  // G4Material *StainlessSteel = new G4Material("StainlessSteel",   density = 8.02*g/cm3, 5);
  //StainlessSteel->AddMaterial(matman->FindOrBuildMaterial("G4_Si"), 0.01);
  //StainlessSteel->AddMaterial(matman->FindOrBuildMaterial("G4_Mn"), 0.02);
  //StainlessSteel->AddMaterial(matman->FindOrBuildMaterial("G4_Cr"), 0.19);
  //StainlessSteel->AddMaterial(matman->FindOrBuildMaterial("G4_Ni"), 0.10);
  //StainlessSteel->AddMaterial(matman->FindOrBuildMaterial("G4_Fe"), 0.68);
  G4Material *stainlessSteel = GetDetectorMaterial("G4_STAINLESS-STEEL");
  double inch = 2.54 * cm;
  // double hangerAngle = 41.39 * M_PI / 180.;
  // double hangerRadius = 32.05 * inch;
  double hangerX = 21.193 *inch;
  double hangerY = 24.246* inch;
  double hangerDiameter = 2. * inch;
  double extra_length = 0.941*inch;




  G4VSolid *hangerSupport = new G4Tubs("tpc_hanger_support", 0, hangerDiameter / 2., (6*inch)+extra_length, 0., 2 * M_PI);
  G4LogicalVolume *hangerSupportLogic = new G4LogicalVolume(hangerSupport,
							    stainlessSteel,
							    "tpc_hanger_support");
  
  double tpc_steel_location = (90 * inch) / 2 - 3 *inch - extra_length / 2.0;

  m_DisplayAction->AddVolume(hangerSupportLogic, "TpcHangerSupport");
 G4VPhysicalVolume *tpc_hanger_support_phys[4] = {nullptr, nullptr,nullptr,nullptr};
  tpc_hanger_support_phys[0] = new G4PVPlacement(0, G4ThreeVector(hangerX, hangerY, tpc_steel_location),
                                              hangerSupportLogic, "tpc_hanger_support_northleft",
                                              logicWorld, false, 0, OverlapCheck());
  tpc_hanger_support_phys[1] = new G4PVPlacement(0, G4ThreeVector(-hangerX, hangerY, tpc_steel_location),
                                              hangerSupportLogic, "tpc_hanger_support_northright",
                                              logicWorld, false, 1, OverlapCheck());
  tpc_hanger_support_phys[2] = new G4PVPlacement(0, G4ThreeVector(hangerX, hangerY, -tpc_steel_location),
                                              hangerSupportLogic, "tpc_hanger_support_southleft",
                                              logicWorld, false, 2, OverlapCheck());
  tpc_hanger_support_phys[3] = new G4PVPlacement(0, G4ThreeVector(-hangerX, hangerY, -tpc_steel_location),
                                              hangerSupportLogic, "tpc_hanger_support_southright",
                                              logicWorld, false, 3, OverlapCheck());

  m_AbsorberVolumeSet.insert(tpc_hanger_support_phys[0]);
  m_AbsorberVolumeSet.insert(tpc_hanger_support_phys[1]);
  m_AbsorberVolumeSet.insert(tpc_hanger_support_phys[2]);
  m_AbsorberVolumeSet.insert(tpc_hanger_support_phys[3]);


  G4Material *carbonFiber = GetDetectorMaterial("CFRP_INTT");
  double hangerBeamInnerDiameter = 2.0 * inch;
  double hangerBeamOuterDiameter = 2.521 * inch;

  G4VSolid *hangerBeam = new G4Tubs("tpc_hanger_beam", hangerBeamInnerDiameter / 2., hangerBeamOuterDiameter / 2., (90 * inch) / 2., 0., 2 * M_PI);
  G4LogicalVolume *hangerBeamLogic = new G4LogicalVolume(hangerBeam,
                                                         carbonFiber,
                                                         "tpc_hanger_beam");

  m_DisplayAction->AddVolume(hangerBeamLogic, "TpcHangerBeam");
 G4VPhysicalVolume *tpc_hanger_beam_phys[2] = {nullptr, nullptr};
 tpc_hanger_beam_phys[0] = new G4PVPlacement(0, G4ThreeVector(hangerX, hangerY, 0),
                                              hangerBeamLogic, "tpc_hanger_beam_left",
                                              logicWorld, false, 0, OverlapCheck());
  tpc_hanger_beam_phys[1] = new G4PVPlacement(0, G4ThreeVector(-hangerX, hangerY, 0),
                                              hangerBeamLogic, "tpc_hanger_beam_right",
                                              logicWorld, false, 1, OverlapCheck());

  m_AbsorberVolumeSet.insert(tpc_hanger_beam_phys[0]);
  m_AbsorberVolumeSet.insert(tpc_hanger_beam_phys[1]);
  
  






  //Twelve one-inch diam carbon fiber rods of thickness 1/16" at 31.04" from beam center
  //borrowed from the INTT specification of carbon fiber
  //note that this defines a clocking!
  double rodAngleStart = M_PI / 12.;
  double rodAngularSpacing = 2 * M_PI / 12.;
  double rodRadius = 31.5 * inch;
  double rodWallThickness = 1. / 8. * inch;
  double rodDiameter = 3. / 4. * inch;
  G4VSolid *tieRod = new G4Tubs("tpc_tie_rod", rodDiameter / 2. - rodWallThickness, rodDiameter / 2., (m_Params->get_double_param("tpc_length") * cm) / 2., 0., 2 * M_PI);
  G4LogicalVolume *tieRodLogic = new G4LogicalVolume(tieRod,
                                                     carbonFiber,
                                                     "tpc_tie_rod");

  G4VPhysicalVolume *tpc_tie_rod_phys[12] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                                             nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};

  std::ostringstream name;
  for (int i = 0; i < 12; i++)
  {
    double ang = rodAngleStart + rodAngularSpacing * i;
    name.str("");
    name << "tpc_tie_rod_" << i;
    tpc_tie_rod_phys[i] = new G4PVPlacement(0, G4ThreeVector(rodRadius * cos(ang), rodRadius * sin(ang), 0),
                                            tieRodLogic, name.str(),
                                            logicWorld, false, i, OverlapCheck());
    m_AbsorberVolumeSet.insert(tpc_tie_rod_phys[i]);
  }

  //  G4VisAttributes *visatt = new G4VisAttributes();
  //  visatt->SetVisibility(true);
  //  visatt->SetForceSolid(true);
  //  visatt->SetColor(PHG4TPCColorDefs::tpc_cu_color);
  //  tpc_window_logic->SetVisAttributes(visatt);

  return 0;
}

int PHG4TpcDetector::ConstructTpcCageVolume(G4LogicalVolume *tpc_envelope)
{
  // 8th layer cu
  // 9th layer FR4
  // 10th layer HoneyComb
  // 11th layer FR4
  // 12th layer Kapton
  // 13th layer FR4
  static const int nlayers = 9;
  static const double thickness[nlayers] = {m_Params->get_double_param("cage_layer_1_thickness") * cm,
                                            m_Params->get_double_param("cage_layer_2_thickness") * cm,
                                            m_Params->get_double_param("cage_layer_3_thickness") * cm,
                                            m_Params->get_double_param("cage_layer_4_thickness") * cm,
                                            m_Params->get_double_param("cage_layer_5_thickness") * cm,
                                            m_Params->get_double_param("cage_layer_6_thickness") * cm,
                                            m_Params->get_double_param("cage_layer_7_thickness") * cm,
                                            m_Params->get_double_param("cage_layer_8_thickness") * cm,
                                            m_Params->get_double_param("cage_layer_9_thickness") * cm};

  static const std::string material[nlayers] = {m_Params->get_string_param("cage_layer_1_material"),
                                                m_Params->get_string_param("cage_layer_2_material"),
                                                m_Params->get_string_param("cage_layer_3_material"),
                                                m_Params->get_string_param("cage_layer_4_material"),
                                                m_Params->get_string_param("cage_layer_5_material"),
                                                m_Params->get_string_param("cage_layer_6_material"),
                                                m_Params->get_string_param("cage_layer_7_material"),
                                                m_Params->get_string_param("cage_layer_8_material"),
                                                m_Params->get_string_param("cage_layer_9_material")};

  double tpc_cage_radius = m_InnerCageRadius;
  std::ostringstream name;
  for (int i = 0; i < nlayers; i++)
  {
    name.str("");
    int layerno = i + 1;
    name << "tpc_cage_layer_" << layerno;
    G4VSolid *tpc_cage_layer = new G4Tubs(name.str(), tpc_cage_radius, tpc_cage_radius + thickness[i], m_Params->get_double_param("tpc_length") * cm / 2., 0., 2 * M_PI);
    G4LogicalVolume *tpc_cage_layer_logic = new G4LogicalVolume(tpc_cage_layer,
                                                                GetDetectorMaterial(material[i]),
                                                                name.str());
    m_DisplayAction->AddTpcInnerLayer(tpc_cage_layer_logic);
    G4VPhysicalVolume *tpc_cage_layer_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                                                               tpc_cage_layer_logic, name.str(),
                                                               tpc_envelope, false, layerno, OverlapCheck());
    m_AbsorberVolumeSet.insert(tpc_cage_layer_phys);
    tpc_cage_radius += thickness[i];
  }
  // outer cage

  tpc_cage_radius = m_OuterCageRadius;
  for (int i = 0; i < nlayers; i++)
  {
    tpc_cage_radius -= thickness[i];
    name.str("");
    int layerno = 10 + 1 + i;  // so the accompanying inner layer is layer - 10
    name << "tpc_cage_layer_" << layerno;
    G4VSolid *tpc_cage_layer = new G4Tubs(name.str(), tpc_cage_radius, tpc_cage_radius + thickness[i], m_Params->get_double_param("tpc_length") * cm / 2., 0., 2 * M_PI);
    G4LogicalVolume *tpc_cage_layer_logic = new G4LogicalVolume(tpc_cage_layer,
                                                                GetDetectorMaterial(material[i]),
                                                                name.str());
    m_DisplayAction->AddTpcOuterLayer(tpc_cage_layer_logic);
    G4VPhysicalVolume *tpc_cage_layer_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                                                               tpc_cage_layer_logic, name.str(),
                                                               tpc_envelope, false, layerno, OverlapCheck());
    m_AbsorberVolumeSet.insert(tpc_cage_layer_phys);
  }

  return 0;
}

void PHG4TpcDetector ::CreateCompositeMaterial(
    std::string compositeName,
    std::vector<std::string> materialName,
    std::vector<double> thickness)
{
  //takes in a list of material names known to Geant already, and thicknesses, and creates a new material called compositeName.

  //check that desired material name doesn't already exist
  G4Material *tempmat = GetDetectorMaterial(compositeName, false);

  if (tempmat != nullptr)
  {
    std::cout << __PRETTY_FUNCTION__ << " Fatal Error: composite material " << compositeName << " already exists" << std::endl;
    assert(!tempmat);
  }

  //check that both arrays have the same depth
  assert(materialName.size() == thickness.size());

  //sum up the areal density and total thickness so we can divvy it out
  double totalArealDensity = 0, totalThickness = 0;
  for (std::vector<double>::size_type i = 0; i < thickness.size(); i++)
  {
    tempmat = GetDetectorMaterial(materialName[i]);
    if (tempmat == nullptr)
    {
      std::cout << __PRETTY_FUNCTION__ << " Fatal Error: component material " << materialName[i] << " does not exist." << std::endl;
      gSystem->Exit(1);
      exit(1);
    }
    totalArealDensity += tempmat->GetDensity() * thickness[i];
    totalThickness += thickness[i];
  }

  //register a new material with the average density of the whole:
  double compositeDensity = totalArealDensity / totalThickness;
  G4Material *composite = new G4Material(compositeName, compositeDensity, thickness.size());

  //now calculate the fraction due to each material, and register those
  for (std::vector<double>::size_type i = 0; i < thickness.size(); i++)
  {
    tempmat = GetDetectorMaterial(materialName[i]);  //don't need to check this, since we did in the previous loop.
    composite->AddMaterial(tempmat, thickness[i] * tempmat->GetDensity() / totalArealDensity);
  }

  //how to register our finished material?
  return;
}


//_______________________________________________________________
void PHG4TpcDetector::add_geometry_node()
{

  // create PHG4TpcCylinderGeomContainer and put on node tree
  const std::string geonode_name = "CYLINDERCELLGEOM_SVTX";  
  auto geonode = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode(), geonode_name);
  if (!geonode)
  {
    geonode = new PHG4TpcCylinderGeomContainer;
    PHNodeIterator iter(topNode());
    auto runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
    auto newNode = new PHIODataNode<PHObject>(geonode, geonode_name, "PHObject");
    runNode->addNode(newNode);
  }

  const std::array<int, 3> NTpcLayers =
  {{
    m_Params->get_int_param("ntpc_layers_inner"),
    m_Params->get_int_param("ntpc_layers_mid"),
    m_Params->get_int_param("ntpc_layers_outer")
  }};
    
  std::array<double, 3> MinLayer;
  MinLayer[0] = m_Params->get_int_param("tpc_minlayer_inner");
  MinLayer[1] = MinLayer[0] + NTpcLayers[0];
  MinLayer[2] = MinLayer[1] + NTpcLayers[1];
    
  const std::array<double, 3> MinRadius =
  {{
    m_Params->get_double_param("tpc_minradius_inner"),
    m_Params->get_double_param("tpc_minradius_mid"),
    m_Params->get_double_param("tpc_minradius_outer")
  }};

  const std::array<double, 3> MaxRadius = 
  {{
    m_Params->get_double_param("tpc_maxradius_inner"),
    m_Params->get_double_param("tpc_maxradius_mid"),
    m_Params->get_double_param("tpc_maxradius_outer"),
  }};
  
  const std::array<double, 5> Thickness =
  {{
    0.687,
    1.012,     
    1.088,     
    0.534,     
    0.595,     
  }}; 

  const double drift_velocity = m_Params->get_double_param("drift_velocity");
  const double tpc_adc_clock = m_Params->get_double_param("tpc_adc_clock");
  const double MaxZ = m_Params->get_double_param("maxdriftlength");
  const double TBinWidth = tpc_adc_clock;
  const double extended_readout_time = m_Params->get_double_param("extended_readout_time");
  const double MaxT = extended_readout_time + 2.* MaxZ / drift_velocity;  // allows for extended time readout
  const double MinT = 0;
  const int NTBins = (int) ((MaxT - MinT) / TBinWidth) + 1;

  std::cout << PHWHERE << "MaxT " << MaxT << " TBinWidth " << TBinWidth << " extended readout time " 
	    << extended_readout_time  << " NTBins = " << NTBins << " drift velocity " << drift_velocity << std::endl;

  const std::array<double, 3> SectorPhi =
  {{
    m_Params->get_double_param("tpc_sector_phi_inner"),
    m_Params->get_double_param("tpc_sector_phi_mid"),
    m_Params->get_double_param("tpc_sector_phi_outer")
  }};
  
  const std::array<int, 3> NPhiBins =
  {{
    m_Params->get_int_param("ntpc_phibins_inner"),
    m_Params->get_int_param("ntpc_phibins_mid"),
    m_Params->get_int_param("ntpc_phibins_outer")
  }};

  const std::array<double, 3> PhiBinWidth =
  {{
    SectorPhi[0] * 12 / (double) NPhiBins[0],
    SectorPhi[1] * 12 / (double) NPhiBins[1],
    SectorPhi[2] * 12 / (double) NPhiBins[2]
  }};
  
  // should move to a common file 
  static constexpr int NSides = 2;
  static constexpr int NSectors = 12;

  std::array<std::vector<double>, NSides > sector_R_bias;
  std::array<std::vector<double>, NSides > sector_Phi_bias;
  std::array<std::vector<double>, NSides > sector_min_Phi;
  std::array<std::vector<double>, NSides > sector_max_Phi;
  

  for (int iregion = 0; iregion < 3; ++iregion)
  {
    //int zside = 0;
    for (int zside = 0; zside < 2;zside++)
    {
      sector_R_bias[zside].clear();
      sector_Phi_bias[zside].clear();
      sector_min_Phi[zside].clear();
      sector_max_Phi[zside].clear();
      //int eff_layer = 0;
      for (int isector = 0; isector < NSectors; ++isector)//12 sectors
      {

        // no bias per default
        // TODO: confirm with what is in PHG4TpcPadPlane Readout
        sector_R_bias[zside].push_back(0);
        sector_Phi_bias[zside].push_back(0);
        
        double sec_gap = (2*M_PI - SectorPhi[iregion]*12)/12;
        double sec_max_phi = M_PI - SectorPhi[iregion]/2 - sec_gap - 2 * M_PI / 12 * isector;// * (isector+1) ;
        double sec_min_phi = sec_max_phi - SectorPhi[iregion];
        sector_min_Phi[zside].push_back(sec_min_phi);
        sector_max_Phi[zside].push_back(sec_max_phi);
      }// isector
    }
    
    double sum_r = 0;
    for (int layer = MinLayer[iregion]; layer < MinLayer[iregion] + NTpcLayers[iregion]; ++layer)
    {          
      double r_length = Thickness[iregion];
      if(iregion == 0 && layer>0){
        if(layer%2==0) r_length = Thickness[4];
        else r_length = Thickness[3];
      }
      sum_r += r_length;
    }    
    const double pad_space = (MaxRadius[iregion] - MinRadius[iregion] - sum_r)/(NTpcLayers[iregion]-1);
    double current_r = MinRadius[iregion];
    
    for (int layer = MinLayer[iregion]; layer < MinLayer[iregion] + NTpcLayers[iregion]; ++layer)
    {
      if (Verbosity())
      {
        std::cout << " layer " << layer << " MinLayer " << MinLayer[iregion] << " region " << iregion
          << " radius " << MinRadius[iregion] + ((double) (layer - MinLayer[iregion]) + 0.5) * Thickness[iregion]
          << " thickness " << Thickness[iregion]
          << " NTBins " << NTBins << " tmin " << MinT << " tstep " << TBinWidth
          << " phibins " << NPhiBins[iregion] << " phistep " << PhiBinWidth[iregion] << std::endl;
      }
      
      auto layerseggeo = new PHG4TpcCylinderGeom;
      layerseggeo->set_layer(layer);
            
      double r_length = Thickness[iregion];
      if(iregion == 0 && layer>0)
      {
        if(layer%2==0) r_length = Thickness[4];
        else r_length = Thickness[3];
      }
      layerseggeo->set_thickness(r_length);
      layerseggeo->set_radius(current_r+r_length/2);
      layerseggeo->set_binning(PHG4CellDefs::sizebinning);
      layerseggeo->set_zbins(NTBins);
      layerseggeo->set_zmin(MinT);
      layerseggeo->set_zstep(TBinWidth);
      layerseggeo->set_phibins(NPhiBins[iregion]);
      layerseggeo->set_phistep(PhiBinWidth[iregion]);
      layerseggeo->set_r_bias(sector_R_bias);
      layerseggeo->set_phi_bias(sector_Phi_bias);
      layerseggeo->set_sector_min_phi(sector_min_Phi);
      layerseggeo->set_sector_max_phi(sector_max_Phi);
      
      // Chris Pinkenburg: greater causes huge memory growth which causes problems
      // on our farm. If you need to increase this - TALK TO ME first
      if (NPhiBins[iregion] * NTBins > 5100000)
      {
        std::cout << "increase Tpc cellsize, number of cells "
          << NPhiBins[iregion] * NTBins << " for layer " << layer
          << " exceed 5.1M limit" << std::endl;
        gSystem->Exit(1);
      }
      geonode->AddLayerCellGeom(layerseggeo);
      
      current_r += r_length + pad_space;
    }
  }
 }
