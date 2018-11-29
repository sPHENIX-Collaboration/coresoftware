#include "PHG4TPCDetector.h"
#include "PHG4TPCDefs.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4ColorDefs.h>
#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4Region.hh>
#include <Geant4/G4RegionStore.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4UserLimits.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cmath>
#include <sstream>

using namespace std;

//_______________________________________________________________
PHG4TPCDetector::PHG4TPCDetector(PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam)
  : PHG4Detector(Node, dnam)
  , params(parameters)
  , g4userlimits(nullptr)
  , active(params->get_int_param("active"))
  , absorberactive(params->get_int_param("absorberactive"))
  , inner_cage_radius(params->get_double_param("gas_inner_radius") * cm - params->get_double_param("cage_layer_9_thickness") * cm - params->get_double_param("cage_layer_8_thickness") * cm - params->get_double_param("cage_layer_7_thickness") * cm - params->get_double_param("cage_layer_6_thickness") * cm - params->get_double_param("cage_layer_5_thickness") * cm - params->get_double_param("cage_layer_4_thickness") * cm - params->get_double_param("cage_layer_3_thickness") * cm - params->get_double_param("cage_layer_2_thickness") * cm - params->get_double_param("cage_layer_1_thickness") * cm)
  , outer_cage_radius(params->get_double_param("gas_outer_radius") * cm + params->get_double_param("cage_layer_9_thickness") * cm + params->get_double_param("cage_layer_8_thickness") * cm + params->get_double_param("cage_layer_7_thickness") * cm + params->get_double_param("cage_layer_6_thickness") * cm + params->get_double_param("cage_layer_5_thickness") * cm + params->get_double_param("cage_layer_4_thickness") * cm + params->get_double_param("cage_layer_3_thickness") * cm + params->get_double_param("cage_layer_2_thickness") * cm + params->get_double_param("cage_layer_1_thickness") * cm)

{
}

//_______________________________________________________________
int PHG4TPCDetector::IsInTPC(G4VPhysicalVolume *volume) const
{
  if (active)
  {
    if (activevols.find(volume) != activevols.end())
    {
      return 1;
    }
  }
  if (absorberactive)
  {
    if (absorbervols.find(volume) != absorbervols.end())
    {
      return -1;
    }
  }
  return 0;
}

//_______________________________________________________________
void PHG4TPCDetector::Construct(G4LogicalVolume *logicWorld)
{
  // create TPC envelope
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

  double steplimits = params->get_double_param("steplimits") * cm;
  if (isfinite(steplimits))
  {
    g4userlimits = new G4UserLimits(steplimits);
  }

  G4VSolid *tpc_envelope = new G4Tubs("tpc_envelope", inner_cage_radius, outer_cage_radius, params->get_double_param("tpc_length") * cm / 2., 0., 2 * M_PI);

  G4LogicalVolume *tpc_envelope_logic = new G4LogicalVolume(tpc_envelope,
                                                            G4Material::GetMaterial("G4_AIR"),
                                                            "tpc_envelope");
  G4VisAttributes *visatt = new G4VisAttributes();
  visatt->SetVisibility(false);
  tpc_envelope_logic->SetVisAttributes(visatt);

  ConstructTPCCageVolume(tpc_envelope_logic);
  ConstructTPCGasVolume(tpc_envelope_logic);

  new G4PVPlacement(0, G4ThreeVector(params->get_double_param("place_x") * cm,
                                     params->get_double_param("place_y") * cm,
                                     params->get_double_param("place_z") * cm),
                    tpc_envelope_logic, "tpc_envelope",
                    logicWorld, 0, false, OverlapCheck());
}

int PHG4TPCDetector::ConstructTPCGasVolume(G4LogicalVolume *tpc_envelope)
{
  static map<int, string> tpcgasvolname = 
    {{PHG4TPCDefs::North, "tpc_gas_north"}, 
     {PHG4TPCDefs::South,"tpc_gas_south"}};
  double tpc_window_thickness = params->get_double_param("window_thickness")*cm;
  double tpc_half_length = (params->get_double_param("tpc_length")*cm-tpc_window_thickness)/2.;
  G4VSolid *tpc_window = new G4Tubs("tpc_window",params->get_double_param("gas_inner_radius") * cm, params->get_double_param("gas_outer_radius") * cm, tpc_window_thickness/2., 0., 2 * M_PI);

  G4VSolid *tpc_gas = new G4Tubs("tpc_gas", params->get_double_param("gas_inner_radius") * cm, params->get_double_param("gas_outer_radius") * cm,   tpc_half_length / 2., 0., 2 * M_PI);
  G4LogicalVolume *tpc_window_logic = new G4LogicalVolume(tpc_window, 
                                                          G4Material::GetMaterial(params->get_string_param("window_material")),
                                                          "tpc_window");

  G4VisAttributes *visatt = new G4VisAttributes();
  visatt->SetVisibility(true);
  visatt->SetForceSolid(true);
  visatt->SetColor(PHG4TPCColorDefs::tpc_kapton_color);
  tpc_window_logic->SetVisAttributes(visatt);
  G4VPhysicalVolume *tpc_window_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                                                      tpc_window_logic, "tpc_window",
							 tpc_envelope, false, PHG4TPCDefs::Window, OverlapCheck());

  absorbervols.insert(tpc_window_phys);
  G4LogicalVolume *tpc_gas_logic = new G4LogicalVolume(tpc_gas,
                                                       G4Material::GetMaterial(params->get_string_param("tpc_gas")),
                                                       "tpc_gas");

  tpc_gas_logic->SetUserLimits(g4userlimits);
  visatt = new G4VisAttributes();
  visatt->SetVisibility(true);
  visatt->SetForceSolid(true);
  visatt->SetColor(PHG4TPCColorDefs::tpc_gas_color);
  tpc_gas_logic->SetVisAttributes(visatt);
  G4VPhysicalVolume *tpc_gas_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, (tpc_half_length+tpc_window_thickness)/2.),
                                                      tpc_gas_logic, tpcgasvolname[PHG4TPCDefs::North],
                                                      tpc_envelope, false, PHG4TPCDefs::North, OverlapCheck());
  cout << "north copy no: " << tpc_gas_phys->GetCopyNo() << endl;

  activevols.insert(tpc_gas_phys);
  tpc_gas_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, -(tpc_half_length+tpc_window_thickness)/2.),
                                                      tpc_gas_logic, tpcgasvolname[PHG4TPCDefs::South],
                                                      tpc_envelope, false, PHG4TPCDefs::South, OverlapCheck());

  cout << "south copy no: " << tpc_gas_phys->GetCopyNo() << endl;
  activevols.insert(tpc_gas_phys);


#if G4VERSION_NUMBER >= 1033
  const G4RegionStore* theRegionStore = G4RegionStore::GetInstance();
  G4Region *tpcregion = theRegionStore->GetRegion("REGION_TPCGAS");
  tpc_gas_logic->SetRegion(tpcregion);
  tpcregion->AddRootLogicalVolume(tpc_gas_logic);
#endif
  return 0;
}

int PHG4TPCDetector::ConstructTPCCageVolume(G4LogicalVolume *tpc_envelope)
{
  // 8th layer cu
  // 9th layer FR4
  // 10th layer HoneyComb
  // 11th layer FR4
  // 12th layer Kapton
  // 13th layer FR4
  static const int nlayers = 9;
  static const double thickness[nlayers] = {params->get_double_param("cage_layer_1_thickness") * cm,
                                            params->get_double_param("cage_layer_2_thickness") * cm,
                                            params->get_double_param("cage_layer_3_thickness") * cm,
                                            params->get_double_param("cage_layer_4_thickness") * cm,
                                            params->get_double_param("cage_layer_5_thickness") * cm,
                                            params->get_double_param("cage_layer_6_thickness") * cm,
                                            params->get_double_param("cage_layer_7_thickness") * cm,
                                            params->get_double_param("cage_layer_8_thickness") * cm,
                                            params->get_double_param("cage_layer_9_thickness") * cm};

  static const string material[nlayers] = {params->get_string_param("cage_layer_1_material"),
                                           params->get_string_param("cage_layer_2_material"),
                                           params->get_string_param("cage_layer_3_material"),
                                           params->get_string_param("cage_layer_4_material"),
                                           params->get_string_param("cage_layer_5_material"),
                                           params->get_string_param("cage_layer_6_material"),
                                           params->get_string_param("cage_layer_7_material"),
                                           params->get_string_param("cage_layer_8_material"),
                                           params->get_string_param("cage_layer_9_material")};

  static const G4Colour color[nlayers] = {PHG4TPCColorDefs::tpc_cu_color,
                                          PHG4TPCColorDefs::tpc_pcb_color,
                                          PHG4TPCColorDefs::tpc_honeycomb_color,
                                          PHG4TPCColorDefs::tpc_cu_color,
                                          PHG4TPCColorDefs::tpc_pcb_color,
                                          PHG4TPCColorDefs::tpc_kapton_color,
                                          PHG4TPCColorDefs::tpc_cu_color,
                                          PHG4TPCColorDefs::tpc_kapton_color,
                                          PHG4TPCColorDefs::tpc_cu_color};
  double tpc_cage_radius = inner_cage_radius;
  ostringstream name;
  for (int i = 0; i < nlayers; i++)
  {
    name.str("");
    int layerno = i + 1;
    name << "tpc_cage_layer_" << layerno;
    G4VSolid *tpc_cage_layer = new G4Tubs(name.str(), tpc_cage_radius, tpc_cage_radius + thickness[i], params->get_double_param("tpc_length") * cm / 2., 0., 2 * M_PI);
    G4LogicalVolume *tpc_cage_layer_logic = new G4LogicalVolume(tpc_cage_layer,
                                                                G4Material::GetMaterial(material[i]),
                                                                name.str());
    G4VisAttributes *visatt = new G4VisAttributes();
    visatt->SetVisibility(true);
    visatt->SetForceSolid(true);
    visatt->SetColor(color[i]);
    tpc_cage_layer_logic->SetVisAttributes(visatt);
    G4VPhysicalVolume *tpc_cage_layer_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                                                               tpc_cage_layer_logic, name.str(),
                                                               tpc_envelope, false, layerno, OverlapCheck());
    absorbervols.insert(tpc_cage_layer_phys);
    tpc_cage_radius += thickness[i];
  }
  // outer cage
  tpc_cage_radius = outer_cage_radius;
  for (int i = 0; i < nlayers; i++)
  {
    tpc_cage_radius -= thickness[i];
    name.str("");
    int layerno = 10 + 1 + i;  // so the accompanying inner layer is layer - 10
    name << "tpc_cage_layer_" << layerno;
    G4VSolid *tpc_cage_layer = new G4Tubs(name.str(), tpc_cage_radius, tpc_cage_radius + thickness[i], params->get_double_param("tpc_length") * cm / 2., 0., 2 * M_PI);
    G4LogicalVolume *tpc_cage_layer_logic = new G4LogicalVolume(tpc_cage_layer,
                                                                G4Material::GetMaterial(material[i]),
                                                                name.str());
    G4VisAttributes *visatt = new G4VisAttributes();
    visatt->SetVisibility(true);
    visatt->SetForceSolid(true);
    visatt->SetColor(color[i]);
    tpc_cage_layer_logic->SetVisAttributes(visatt);
    G4VPhysicalVolume *tpc_cage_layer_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                                                               tpc_cage_layer_logic, name.str(),
                                                               tpc_envelope, false, layerno, OverlapCheck());
    absorbervols.insert(tpc_cage_layer_phys);
  }

  return 0;
}

int PHG4TPCDetector::DisplayVolume(G4VSolid *volume, G4LogicalVolume *logvol, G4RotationMatrix *rotm)
{
  static int i = 0;
  G4LogicalVolume *checksolid = new G4LogicalVolume(volume, G4Material::GetMaterial("G4_POLYSTYRENE"), "DISPLAYLOGICAL", 0, 0, 0);
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
  new G4PVPlacement(rotm, G4ThreeVector(0, 0, 0), checksolid, "DISPLAYVOL", logvol, 0, false, OverlapCheck());
  return 0;
}
