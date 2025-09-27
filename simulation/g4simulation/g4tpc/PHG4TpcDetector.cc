#include "PHG4TpcDetector.h"
#include "PHG4TpcDefs.h"
#include "PHG4TpcDisplayAction.h"

#include <g4detectors/PHG4CellDefs.h>
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <cdbobjects/CDBTTree.h>

#include <ffamodules/CDBInterface.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>
#include <phool/sphenix_constants.h>

#include <TSystem.h>

#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4NistManager.hh>
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
PHG4TpcDetector::~PHG4TpcDetector()
{
  delete m_cdbttree;
  delete m_G4UserLimits;
}
//_______________________________________________________________
int PHG4TpcDetector::IsInTpc(G4VPhysicalVolume *volume) const
{
  if (m_ActiveFlag)
  {
    if (m_ActiveVolumeSet.contains(volume))
    {
      return 1;
    }
  }
  if (m_AbsorberActiveFlag)
  {
    if (m_AbsorberVolumeSet.contains(volume))
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
  if (std::isfinite(steplimits) && steplimits > 0)
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

  new G4PVPlacement(nullptr, G4ThreeVector(m_Params->get_double_param("place_x") * cm, m_Params->get_double_param("place_y") * cm, m_Params->get_double_param("place_z") * cm),
                    tpc_envelope_logic, "tpc_envelope",
                    logicWorld, false, false, OverlapCheck());

  // geometry node
  add_geometry_node();
}

void PHG4TpcDetector::CreateTpcGasMixture()
{
  G4double density;
  G4int ncomponents;
  G4int natoms;

  G4double tpcGasTemperature = (273.15 + m_Params->get_double_param("TPC_gas_temperature")) * kelvin;
  G4double tpcGasPressure = m_Params->get_double_param("TPC_gas_pressure") * atmosphere;

  G4Material *CF4 = new G4Material("CF4", density = sphenix_constants::CF4_density * mg / cm3, ncomponents = 2, kStateGas, tpcGasTemperature, tpcGasPressure);
  CF4->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), natoms = 1);
  CF4->AddElement(G4NistManager::Instance()->FindOrBuildElement("F"), natoms = 4);

  G4Material *N2 = new G4Material("N2", density = 1.165 * mg / cm3, ncomponents = 1, kStateGas, tpcGasTemperature, tpcGasPressure);
  N2->AddElement(G4NistManager::Instance()->FindOrBuildElement("N"), natoms = 2);

  // Create isobutane as only butane is in the standard G4Material list (they have different densities)
  // https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/Appendix/materialNames.html
  G4Material *isobutane = new G4Material("isobutane", density = 2.59 * mg / cm3, ncomponents = 2, kStateGas, tpcGasTemperature, tpcGasPressure);
  isobutane->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), natoms = 4);  // Could use AddElement(PHG4Detector::GetDetectorElement("C", true), 4); instead?
  isobutane->AddElement(G4NistManager::Instance()->FindOrBuildElement("H"), natoms = 10);

  double Ne_frac = m_Params->get_double_param("Ne_frac");
  double Ar_frac = m_Params->get_double_param("Ar_frac");
  double CF4_frac = m_Params->get_double_param("CF4_frac");
  double N2_frac = m_Params->get_double_param("N2_frac");
  double isobutane_frac = m_Params->get_double_param("isobutane_frac");

  const double Ne_den = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ne")->GetDensity();
  const double Ar_den = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ar")->GetDensity();
  const double CF4_den = CF4->GetDensity();
  const double N2_den = N2->GetDensity();
  const double isobutane_den = isobutane->GetDensity();

  const double sphenix_tpc_gas_den = (Ne_den * Ne_frac) + (Ar_den * Ar_frac) + (CF4_den * CF4_frac) + (N2_den * N2_frac) + (isobutane_den * isobutane_frac);

  G4Material *sPHENIX_tpc_gas = new G4Material("sPHENIX_TPC_Gas", sphenix_tpc_gas_den, ncomponents = 5, kStateGas);  //, tpcGasTemperature, tpcGasPressure);
  sPHENIX_tpc_gas->AddMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_Ne"), (Ne_den * Ne_frac) / sphenix_tpc_gas_den);
  sPHENIX_tpc_gas->AddMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_Ar"), (Ar_den * Ar_frac) / sphenix_tpc_gas_den);
  sPHENIX_tpc_gas->AddMaterial(CF4, (CF4_den * CF4_frac) / sphenix_tpc_gas_den);
  sPHENIX_tpc_gas->AddMaterial(N2, (N2_den * N2_frac) / sphenix_tpc_gas_den);
  sPHENIX_tpc_gas->AddMaterial(isobutane, (isobutane_den * isobutane_frac) / sphenix_tpc_gas_den);
}

int PHG4TpcDetector::ConstructTpcGasVolume(G4LogicalVolume *tpc_envelope)
{
  static std::map<int, std::string> tpcgasvolname =
      {{PHG4TpcDefs::North, "tpc_gas_north"},
       {PHG4TpcDefs::South, "tpc_gas_south"}};

  CreateTpcGasMixture();

  // Window / central membrane
  double tpc_window_thickness = m_Params->get_double_param("window_thickness") * cm;
  double tpc_half_length = (m_Params->get_double_param("tpc_length") * cm - tpc_window_thickness) / 2.;

  //'window' (modernly called central membrane only) material is ENIG, not Copper:
  // thickness in this recipe are just a ratio.  we set the usual thickness below.
  std::vector<double> thickness;
  std::vector<std::string> material;
  material.emplace_back("G4_Ni");
  thickness.push_back(.240 * cm);
  material.emplace_back("G4_Au");
  thickness.push_back(.008 * cm);
  G4Material *temp = nullptr;
  temp = GetDetectorMaterial("ENIG", false);
  if (temp == nullptr)
  {
    CreateCompositeMaterial("ENIG", material, thickness);  // see new function below
  }

  G4VSolid *tpc_window = new G4Tubs("tpc_window", m_Params->get_double_param("gas_inner_radius") * cm, m_Params->get_double_param("gas_outer_radius") * cm, tpc_window_thickness / 2., 0., 2 * M_PI);
  // we build our CM surface:
  G4LogicalVolume *tpc_window_logic = new G4LogicalVolume(tpc_window,
                                                          GetDetectorMaterial("ENIG"),
                                                          "tpc_window");
  // previously:                                                            GetDetectorMaterial(m_Params->get_string_param("window_surface1_material")),

  //  G4VisAttributes *visatt = new G4VisAttributes();
  //  visatt->SetVisibility(true);
  //  visatt->SetForceSolid(true);
  //  visatt->SetColor(PHG4TPCColorDefs::tpc_cu_color);
  //  tpc_window_logic->SetVisAttributes(visatt);

  m_DisplayAction->AddVolume(tpc_window_logic, "TpcWindow");
  G4VPhysicalVolume *tpc_window_phys = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0),
                                                         tpc_window_logic, "tpc_window",
                                                         tpc_envelope, false, PHG4TpcDefs::Window, OverlapCheck());

  m_AbsorberVolumeSet.insert(tpc_window_phys);

  // now build the FR4 layer beneath that:
  //  Window / central membrane core
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
  G4VPhysicalVolume *tpc_window_surface2_core_phys = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0),
                                                                       tpc_window_surface2_core_logic, "tpc_window_surface2_core",
                                                                       tpc_window_logic, false, PHG4TpcDefs::WindowCore, OverlapCheck());

  m_AbsorberVolumeSet.insert(tpc_window_surface2_core_phys);

  // and now the honeycomb core:
  G4VSolid *tpc_window_core =
      new G4Tubs("tpc_window", m_Params->get_double_param("gas_inner_radius") * cm, m_Params->get_double_param("gas_outer_radius") * cm,
                 tpc_window_core_thickness / 2., 0., 2 * M_PI);
  G4LogicalVolume *tpc_window_core_logic = new G4LogicalVolume(tpc_window_core,
                                                               GetDetectorMaterial(m_Params->get_string_param("window_core_material")),
                                                               "tpc_window_core");

  m_DisplayAction->AddVolume(tpc_window_core_logic, "TpcHoneyComb");
  G4VPhysicalVolume *tpc_window_core_phys = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0),
                                                              tpc_window_core_logic, "tpc_window_core",
                                                              tpc_window_surface2_core_logic, false, PHG4TpcDefs::WindowCore, OverlapCheck());

  m_AbsorberVolumeSet.insert(tpc_window_core_phys);

  // Gas
  G4VSolid *tpc_gas = new G4Tubs("tpc_gas", m_Params->get_double_param("gas_inner_radius") * cm, m_Params->get_double_param("gas_outer_radius") * cm, tpc_half_length / 2., 0., 2 * M_PI);

  G4LogicalVolume *tpc_gas_logic = new G4LogicalVolume(tpc_gas,
                                                       GetDetectorMaterial("sPHENIX_TPC_Gas"),
                                                       "tpc_gas");

  tpc_gas_logic->SetUserLimits(m_G4UserLimits);
  m_DisplayAction->AddVolume(tpc_gas_logic, "TpcGas");
  G4VPhysicalVolume *tpc_gas_phys = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, (tpc_half_length + tpc_window_thickness) / 2.),
                                                      tpc_gas_logic, tpcgasvolname[PHG4TpcDefs::North],
                                                      tpc_envelope, false, PHG4TpcDefs::North, OverlapCheck());
  //  std::cout << "north copy no: " << tpc_gas_phys->GetCopyNo() << std::endl;

  m_ActiveVolumeSet.insert(tpc_gas_phys);
  tpc_gas_phys = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -(tpc_half_length + tpc_window_thickness) / 2.),
                                   tpc_gas_logic, tpcgasvolname[PHG4TpcDefs::South],
                                   tpc_envelope, false, PHG4TpcDefs::South, OverlapCheck());

  //  std::cout << "south copy no: " << tpc_gas_phys->GetCopyNo() << std::endl;
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
  // note that these elements are outside the tpc logical volume!

  // Two two-inch diam. 304 Stainless Steel solid 'hanger beams' at 32.05" from beam center
  // at +/- 41.39 degrees left and right of vertical
  // stainless steel: 0.695 iron, 0.190 chromium, 0.095 nickel, 0.020 manganese.
  // if we're being pedantic, that is.  But store-bought stainless is probably okay.
  // G4Material *StainlessSteel = new G4Material("StainlessSteel",   density = 8.02*g/cm3, 5);
  // StainlessSteel->AddMaterial(matman->FindOrBuildMaterial("G4_Si"), 0.01);
  // StainlessSteel->AddMaterial(matman->FindOrBuildMaterial("G4_Mn"), 0.02);
  // StainlessSteel->AddMaterial(matman->FindOrBuildMaterial("G4_Cr"), 0.19);
  // StainlessSteel->AddMaterial(matman->FindOrBuildMaterial("G4_Ni"), 0.10);
  // StainlessSteel->AddMaterial(matman->FindOrBuildMaterial("G4_Fe"), 0.68);
  G4Material *stainlessSteel = GetDetectorMaterial("G4_STAINLESS-STEEL");
  double inch = 2.54 * cm;
  // double hangerAngle = 41.39 * M_PI / 180.;
  // double hangerRadius = 32.05 * inch;
  double hangerX = 21.193 * inch;
  double hangerY = 24.246 * inch;
  double hangerDiameter = 2. * inch;
  double extra_length = 0.941 * inch;

  G4VSolid *hangerSupport = new G4Tubs("tpc_hanger_support", 0, hangerDiameter / 2., (6 * inch) + extra_length, 0., 2 * M_PI);
  G4LogicalVolume *hangerSupportLogic = new G4LogicalVolume(hangerSupport,
                                                            stainlessSteel,
                                                            "tpc_hanger_support");

  double tpc_steel_location = (90 * inch) / 2 - 3 * inch - extra_length / 2.0;

  m_DisplayAction->AddVolume(hangerSupportLogic, "TpcHangerSupport");
  G4VPhysicalVolume *tpc_hanger_support_phys[4] = {nullptr, nullptr, nullptr, nullptr};
  tpc_hanger_support_phys[0] = new G4PVPlacement(nullptr, G4ThreeVector(hangerX, hangerY, tpc_steel_location),
                                                 hangerSupportLogic, "tpc_hanger_support_northleft",
                                                 logicWorld, false, 0, OverlapCheck());
  tpc_hanger_support_phys[1] = new G4PVPlacement(nullptr, G4ThreeVector(-hangerX, hangerY, tpc_steel_location),
                                                 hangerSupportLogic, "tpc_hanger_support_northright",
                                                 logicWorld, false, 1, OverlapCheck());
  tpc_hanger_support_phys[2] = new G4PVPlacement(nullptr, G4ThreeVector(hangerX, hangerY, -tpc_steel_location),
                                                 hangerSupportLogic, "tpc_hanger_support_southleft",
                                                 logicWorld, false, 2, OverlapCheck());
  tpc_hanger_support_phys[3] = new G4PVPlacement(nullptr, G4ThreeVector(-hangerX, hangerY, -tpc_steel_location),
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
  tpc_hanger_beam_phys[0] = new G4PVPlacement(nullptr, G4ThreeVector(hangerX, hangerY, 0),
                                              hangerBeamLogic, "tpc_hanger_beam_left",
                                              logicWorld, false, 0, OverlapCheck());
  tpc_hanger_beam_phys[1] = new G4PVPlacement(nullptr, G4ThreeVector(-hangerX, hangerY, 0),
                                              hangerBeamLogic, "tpc_hanger_beam_right",
                                              logicWorld, false, 1, OverlapCheck());

  m_AbsorberVolumeSet.insert(tpc_hanger_beam_phys[0]);
  m_AbsorberVolumeSet.insert(tpc_hanger_beam_phys[1]);

  // Twelve one-inch diam carbon fiber rods of thickness 1/16" at 31.04" from beam center
  // borrowed from the INTT specification of carbon fiber
  // note that this defines a clocking!
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
    tpc_tie_rod_phys[i] = new G4PVPlacement(nullptr, G4ThreeVector(rodRadius * cos(ang), rodRadius * sin(ang), 0),
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
    G4VPhysicalVolume *tpc_cage_layer_phys = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0),
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
    G4VPhysicalVolume *tpc_cage_layer_phys = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0),
                                                               tpc_cage_layer_logic, name.str(),
                                                               tpc_envelope, false, layerno, OverlapCheck());
    m_AbsorberVolumeSet.insert(tpc_cage_layer_phys);
  }

  return 0;
}

void PHG4TpcDetector::CreateCompositeMaterial(
    const std::string &compositeName,
    std::vector<std::string> materialName,
    const std::vector<double> &thickness)
{
  // takes in a list of material names known to Geant already, and thicknesses, and creates a new material called compositeName.

  // check that desired material name doesn't already exist
  G4Material *tempmat = GetDetectorMaterial(compositeName, false);

  if (tempmat != nullptr)
  {
    std::cout << __PRETTY_FUNCTION__ << " Fatal Error: composite material " << compositeName << " already exists" << std::endl;
    assert(!tempmat);
  }

  // check that both arrays have the same depth
  assert(materialName.size() == thickness.size());

  // sum up the areal density and total thickness so we can divvy it out
  double totalArealDensity = 0;
  double totalThickness = 0;
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

  // register a new material with the average density of the whole:
  double compositeDensity = totalArealDensity / totalThickness;
  G4Material *composite = new G4Material(compositeName, compositeDensity, thickness.size());

  // now calculate the fraction due to each material, and register those
  for (std::vector<double>::size_type i = 0; i < thickness.size(); i++)
  {
    tempmat = GetDetectorMaterial(materialName[i]);  // don't need to check this, since we did in the previous loop.
    composite->AddMaterial(tempmat, thickness[i] * tempmat->GetDensity() / totalArealDensity);
  }

  // how to register our finished material?
  return;
}

//_______________________________________________________________
void PHG4TpcDetector::add_geometry_node()
{
  // create PHG4TpcCylinderGeomContainer and put on node tree
  const std::string geonode_name = "CYLINDERCELLGEOM_SVTX";
  auto *geonode = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode(), geonode_name);
  if (!geonode)
  {
    geonode = new PHG4TpcCylinderGeomContainer;
    PHNodeIterator iter(topNode());
    auto *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
    auto *newNode = new PHIODataNode<PHObject>(geonode, geonode_name, "PHObject");
    runNode->addNode(newNode);
  }

  m_cdb = CDBInterface::instance();
  std::string calibdir = m_cdb->getUrl("TPC_FEE_CHANNEL_MAP");

  if (!calibdir.empty())
  {
    m_cdbttree = new CDBTTree(calibdir);
    m_cdbttree->LoadCalibrations();
  }
  else
  {
    std::cout << "PHG4TpcPadPlaneReadout::InitRun no TPC_FEE_CHANNEL_MAP found" << std::endl;
    exit(1);
  }

  const std::array<int, 3> NTpcLayers =
      {{m_Params->get_int_param("ntpc_layers_inner"),
        m_Params->get_int_param("ntpc_layers_mid"),
        m_Params->get_int_param("ntpc_layers_outer")}};

  std::array<double, 3> MinLayer{};
  MinLayer[0] = m_Params->get_int_param("tpc_minlayer_inner");
  MinLayer[1] = MinLayer[0] + NTpcLayers[0];
  MinLayer[2] = MinLayer[1] + NTpcLayers[1];

  const std::array<double, 5> Thickness =
      {{
          0.56598621677629212,
          1.0206889851687158,
          1.0970475085472556,
          0.5630547309825637,
          0.56891770257002054,
      }};

  const double drift_velocity = m_Params->get_double_param("drift_velocity");
  const double tpc_adc_clock = m_Params->get_double_param("tpc_adc_clock");
  const double MaxZ = m_Params->get_double_param("maxdriftlength");
  const double TBinWidth = tpc_adc_clock;
  const double extended_readout_time = m_Params->get_double_param("extended_readout_time");
  const double MaxT = extended_readout_time + 2. * MaxZ / drift_velocity;  // allows for extended time readout
  const double MinT = 0;
  const int NTBins = (int) ((MaxT - MinT) / TBinWidth) + 1;

  std::cout << PHWHERE << "MaxT " << MaxT << " TBinWidth " << TBinWidth << " extended readout time "
            << extended_readout_time << " NTBins = " << NTBins << " drift velocity " << drift_velocity << std::endl;

  const std::array<int, 3> NPhiBins =
      {{m_Params->get_int_param("ntpc_phibins_inner"),
        m_Params->get_int_param("ntpc_phibins_mid"),
        m_Params->get_int_param("ntpc_phibins_outer")}};

  // should move to a common file
  static constexpr int NSides = 2;
  static constexpr int NSectors = 12;
  static constexpr int NLayers = 16 * 3;

  std::array<std::vector<double>, NSides> sector_R_bias;
  std::array<std::vector<double>, NSides> sector_Phi_bias;
  std::array<std::vector<double>, NSides> sector_min_Phi;
  std::array<std::vector<double>, NSides> sector_max_Phi;
  std::array<std::vector<double>, NLayers> pad_phi;
  std::array<std::vector<double>, NLayers> pad_R;
  std::array<double, NLayers> layer_radius{};
  std::array<double, NLayers> phi_bin_width_cdb{};
  std::array<std::array<std::array<double, 3>, NSectors>, NSides> sec_max_phi{};
  std::array<std::array<std::array<double, 3>, NSectors>, NSides> sec_min_phi{};
  int Nfee = 26;
  int Nch = 256;

  for (int f = 0; f < Nfee; f++)
  {
    for (int ch = 0; ch < Nch; ch++)
    {
      unsigned int key = 256 * (f) + ch;
      std::string varname = "layer";
      int l = m_cdbttree->GetIntValue(key, varname);
      if (l > 6)
      {
        int v_layer = l - 7;
        std::string phiname = "phi";
        pad_phi[v_layer].push_back(m_cdbttree->GetDoubleValue(key, phiname));
        std::string rname = "R";
        pad_R[v_layer].push_back(m_cdbttree->GetDoubleValue(key, rname));
      }
    }
  }

  for (size_t layer = 0; layer < NLayers; layer++)
  {
    layer_radius[layer] = 0;
    for (int pad = 0; pad < (int) pad_R[(int) layer].size(); pad++)
    {
      layer_radius[(int) layer] += pad_R[(int) layer][pad];
    }

    layer_radius[(int) layer] = layer_radius[(int) layer] / pad_R[(int) layer].size();
    layer_radius[(int) layer] = layer_radius[(int) layer] / 10.;

    auto min_phi_iter = std::min_element(pad_phi[layer].begin(), pad_phi[layer].end());
    auto max_phi_iter = std::max_element(pad_phi[layer].begin(), pad_phi[layer].end());
    double min_phi = static_cast<double>(*min_phi_iter);
    double max_phi = static_cast<double>(*max_phi_iter);
    min_phi = min_phi - M_PI / 2.;
    max_phi = max_phi - M_PI / 2.;
    phi_bin_width_cdb[layer] = std::abs(max_phi - min_phi) / (NPhiBins[(int) (layer / 16)] / 12 - 1);  // NOLINT(bugprone-integer-division)
    double SectorPhi = std::abs(max_phi - min_phi) + phi_bin_width_cdb[layer];
    for (int zside = 0; zside < 2; zside++)
    {
      for (int isector = 0; isector < NSectors; isector++)  // 12 sectors
      {
        if (zside == 0)
        {
          sec_min_phi[zside][isector][(int) layer / 16] = M_PI - 2 * M_PI / 12 * (isector + 1) + (-(max_phi) -phi_bin_width_cdb[layer] / 2.);
          sec_max_phi[zside][isector][(int) layer / 16] = sec_min_phi[zside][isector][(int) layer / 16] + SectorPhi;
        }
        if (zside == 1)
        {
          sec_max_phi[zside][isector][(int) layer / 16] = M_PI - 2 * M_PI / 12 * (isector + 1) + (max_phi + phi_bin_width_cdb[layer] / 2.);
          sec_min_phi[zside][isector][(int) layer / 16] = sec_max_phi[zside][isector][(int) layer / 16] - SectorPhi;
        }
      }
    }
  }

  for (int iregion = 0; iregion < 3; ++iregion)
  {
    // int zside = 0;
    for (int zside = 0; zside < 2; zside++)
    {
      sector_R_bias[zside].clear();
      sector_Phi_bias[zside].clear();
      sector_min_Phi[zside].clear();
      sector_max_Phi[zside].clear();
      // int eff_layer = 0;
      for (int isector = 0; isector < NSectors; ++isector)  // 12 sectors
      {
        // no bias per default
        // TODO: confirm with what is in PHG4TpcPadPlane Readout
        sector_R_bias[zside].push_back(0);
        sector_Phi_bias[zside].push_back(0);

        sector_min_Phi[zside].push_back(sec_min_phi[zside][isector][iregion]);
        sector_max_Phi[zside].push_back(sec_max_phi[zside][isector][iregion]);
      }  // isector
    }

    for (int layer = MinLayer[iregion]; layer < MinLayer[iregion] + NTpcLayers[iregion]; ++layer)
    {
      if (Verbosity())
      {
        std::cout << " layer " << layer << " MinLayer " << MinLayer[iregion] << " region " << iregion
                  << " radius " << layer_radius[layer - 7]
                  << " thickness " << Thickness[iregion]
                  << " NTBins " << NTBins << " tmin " << MinT << " tstep " << TBinWidth
                  << " phibins " << NPhiBins[iregion] << " phistep " << phi_bin_width_cdb[layer] << std::endl;
      }

      auto *layerseggeo = new PHG4TpcCylinderGeom;
      layerseggeo->set_layer(layer);

      double r_length = Thickness[iregion];
      if (iregion == 0)
      {
        if (layer % 2 == 0)
        {
          r_length = Thickness[4];
        }
        else
        {
          r_length = Thickness[3];
        }
      }
      int v_layer = layer - 7;
      if (v_layer >= 0)
      {
        layerseggeo->set_thickness(r_length);
        layerseggeo->set_radius(layer_radius[v_layer]);
        layerseggeo->set_binning(PHG4CellDefs::sizebinning);
        layerseggeo->set_zbins(NTBins);
        layerseggeo->set_zmin(MinT);
        layerseggeo->set_zstep(TBinWidth);
        layerseggeo->set_phibins(NPhiBins[iregion]);
        layerseggeo->set_phistep(phi_bin_width_cdb[v_layer]);
        layerseggeo->set_r_bias(sector_R_bias);
        layerseggeo->set_phi_bias(sector_Phi_bias);
        layerseggeo->set_sector_min_phi(sector_min_Phi);
        layerseggeo->set_sector_max_phi(sector_max_Phi);
      }

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
    }
  }
}
