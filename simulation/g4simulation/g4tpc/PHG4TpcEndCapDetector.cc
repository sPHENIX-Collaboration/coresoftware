#include "PHG4TpcEndCapDetector.h"

#include "PHG4TpcEndCapDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <TSystem.h>

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4ExtrudedSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4String.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4Transform3D.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4TwoVector.hh>
#include <Geant4/G4Types.hh>  // for G4double
#include <Geant4/G4VPhysicalVolume.hh>

#include <CLHEP/Vector/RotationZ.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <boost/format.hpp>
#pragma GCC diagnostic pop

#include <algorithm>  // for max, copy
#include <cassert>
#include <cmath>
#include <cstdlib>  // for exit
#include <iostream>

class G4VSolid;
class PHCompositeNode;

//____________________________________________________________________________..
PHG4TpcEndCapDetector::PHG4TpcEndCapDetector(PHG4Subsystem *subsys,
                                             PHCompositeNode *Node,
                                             PHParameters *parameters,
                                             const std::string &dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_Params(parameters)
  , m_DisplayAction(dynamic_cast<PHG4TpcEndCapDisplayAction *>(subsys->GetDisplayAction()))

{
  assert(subsys->GetDisplayAction());
  assert(m_DisplayAction);
  Verbosity(m_Params->get_int_param("construction_verbosity"));
}

PHG4TpcEndCapDetector::~PHG4TpcEndCapDetector()
{
  if (m_EndCapAssembly)
  {
    if (Verbosity())
    {
      std::cout << __PRETTY_FUNCTION__ << " delete m_EndCapAssembly" << std::endl;
    }

    delete m_EndCapAssembly;
  }
}

//_______________________________________________________________
int PHG4TpcEndCapDetector::IsInDetector(G4VPhysicalVolume *volume) const
{
  G4LogicalVolume *logvol = volume->GetLogicalVolume();
  std::set<G4LogicalVolume *>::const_iterator iter = m_LogicalVolumesSet.find(logvol);
  if (iter != m_LogicalVolumesSet.end())
  {
    return 1;
  }
  return 0;
}

//_______________________________________________________________
void PHG4TpcEndCapDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  assert(m_DisplayAction);

  assert(m_EndCapAssembly == nullptr);
  m_EndCapAssembly = ConstructEndCapAssembly();
  assert(m_EndCapAssembly);

  G4TranslateZ3D g4vec_front_z(m_Params->get_double_param("envelop_front_surface_z") * cm);

  G4RotateY3D rotm_otherside(180 * deg);

  G4ThreeVector g4vec_center(m_Params->get_double_param("place_x") * cm,
                             m_Params->get_double_param("place_y") * cm,
                             m_Params->get_double_param("place_z") * cm);
  G4RotationMatrix rotm_center;
  rotm_center.rotateX(m_Params->get_double_param("rot_x") * deg);
  rotm_center.rotateY(m_Params->get_double_param("rot_y") * deg);
  rotm_center.rotateZ(m_Params->get_double_param("rot_z") * deg);
  G4Transform3D transform_center(rotm_center, g4vec_center);

  int i = 0;
  //  G4Transform3D transform_side1 = g4vec_front_z * transform_center;
  G4Transform3D transform_side1 = transform_center * g4vec_front_z;
  m_EndCapAssembly->MakeImprint(logicWorld, transform_side1, i++, OverlapCheck());
  // the other side
  G4Transform3D transform_side2 = transform_center * rotm_otherside * g4vec_front_z;
  m_EndCapAssembly->MakeImprint(logicWorld, transform_side2, i++, OverlapCheck());

  return;
}

G4AssemblyVolume *PHG4TpcEndCapDetector::ConstructEndCapAssembly()
{
  G4AssemblyVolume *assemblyvol = new G4AssemblyVolume();
  G4double starting_z(0);

  // Internal HBD structure
  // From doi:10.1016/j.nima.2011.04.015
  // Component Material X0 (cm) Thickness (cm) Area (%) Rad. Length (%)
  //  Mesh SS 1.67 0.003 11.5 0.021  <- not used for GEMs trackers
  //  AddLayer("Mesh", "Steel",
  //          0.003 * cm, false, 11.5);

  //  //  GEM frames FR4 17.1 0.15x4 6.5 0.228 <- not used for GEMs trackers
  //  AddLayer("Frame0", "G10",
  //          0.15 * cm, false, 6.5);

  std::vector<double> thickness;
  std::vector<std::string> material;
  material.push_back("G4_Cu");
  thickness.push_back(0.0005 * 2. * cm);
  material.push_back("G4_KAPTON");
  thickness.push_back(0.005 * cm);
  material.push_back("sPHENIX_TPC_Gas");  // proper gas name, but should be pulled from params to match TpcSubsystem?
  thickness.push_back(0.2 * cm);
  G4Material *temp = GetDetectorMaterial("GEMeffective", false);
  if (temp == nullptr)
  {
    CreateCompositeMaterial("GEMeffective", material, thickness);  //see new function below
  }
  double totalThickness = 0;
  for (std::vector<double>::size_type i = 0; i < thickness.size(); i++)
  {
    totalThickness += thickness[i];
  }

  const int n_GEM_layers = m_Params->get_int_param("n_GEM_layers");

  //instead of building this layer-by-layer, we build a single block corresponding to all the gems that were previously handled in this fashion:
  totalThickness *= n_GEM_layers;
  AddLayer(assemblyvol, starting_z, G4String("GEMAllParts"), "GEMeffective", totalThickness, 64);  //note this slightly undercounts the gas because the gas fill should be 100%, and slightly mispositions the inner edge of the material because the way it is made <100% in AddLayer is by making it thinner than nominally requested but centering it in the region it would have occupied.

  // 16 layer readout plane by TTM
  // https://indico.bnl.gov/event/8307/contributions/36744/attachments/27646/42337/R3-Review.pptx
  const int n_PCB_layers(16);
  // 35 um / layer Cu
  AddLayer(assemblyvol, starting_z, G4String("PCBCu"), "G4_Cu", 0.0035 * cm * n_PCB_layers, 80);
  // 7 mil / layer board
  AddLayer(assemblyvol, starting_z, "PCBBase", "FR4", 0.00254 * cm * 7 * n_PCB_layers, 100);

  ConstructWagonWheel(assemblyvol, starting_z);
  ConstructElectronics(assemblyvol, starting_z);

  return assemblyvol;
}

void PHG4TpcEndCapDetector ::CreateCompositeMaterial(
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

void PHG4TpcEndCapDetector ::AddLayer(  //
    G4AssemblyVolume *assemblyvol,
    G4double &z_start,
    const std::string &_name,  //! name base for this layer
    std::string _material,     //! material name in G4
    G4double _depth,           //! depth in G4 units
    double _percentage_filled  //! percentage filled//
)
{
  z_start += _depth / 2.;
  G4ThreeVector g4vec(0, 0, z_start);
  z_start += _depth / 2.;

  std::string name_base = boost::str(boost::format("%1%_Layer_%2%") % GetName() % _name);

  G4VSolid *solid_layer = new G4Tubs(
      name_base,
      m_Params->get_double_param("envelop_r_min") * cm,
      m_Params->get_double_param("envelop_r_max") * cm,
      _depth * _percentage_filled / 100. / 2.,
      0, CLHEP::twopi);

  auto material = GetDetectorMaterial(_material);
  if (material == nullptr)
  {
    std::cout << __PRETTY_FUNCTION__ << " Fatal Error: missing material " << _material << std::endl;
    assert(material);
  }

  G4LogicalVolume *logical_layer = new G4LogicalVolume(solid_layer, material, name_base);
  m_LogicalVolumesSet.insert(logical_layer);

  assemblyvol->AddPlacedVolume(logical_layer, g4vec, nullptr);

  assert(m_DisplayAction);
  m_DisplayAction->AddVolume(logical_layer, _material);

  return;
}

void PHG4TpcEndCapDetector::ConstructWagonWheel(G4AssemblyVolume *assmeblyvol,
                                                G4double &z_start)  // careful z_start is modified and being used later
{
  const int n_sectors = m_Params->get_int_param("n_sectors");
  assert(n_sectors >= 1);
  const int n_radial_modules = m_Params->get_int_param("n_radial_modules");
  assert(n_radial_modules >= 1);

  const std::string material_name(m_Params->get_string_param("wagon_wheel_material"));
  auto material = GetDetectorMaterial(material_name);
  if (material == nullptr)
  {
    std::cout << __PRETTY_FUNCTION__ << " Fatal Error: missing material " << m_Params->get_string_param("wagon_wheel_material") << std::endl;
    assert(material);
  }
  const G4double wagon_wheel_sector_phi_offset = m_Params->get_double_param("wagon_wheel_sector_phi_offset_degree") * degree;

  ///////////////////////////////////////////////
  // wagon_wheel_front_frame ring
  ///////////////////////////////////////////////
  if (Verbosity())
  {
    std::cout << __PRETTY_FUNCTION__ << " - wagon_wheel_front_frame z_start = " << z_start << std::endl;
  }

  const G4double wagon_wheel_front_frame_thickness = m_Params->get_double_param("wagon_wheel_front_frame_thickness") * cm;
  const G4double wagon_wheel_front_frame_spoke_width = m_Params->get_double_param("wagon_wheel_front_frame_spoke_width") * cm;

  z_start += wagon_wheel_front_frame_thickness / 2.;
  G4ThreeVector g4vec_wagon_wheel_front_frame(0, 0, z_start);
  z_start += wagon_wheel_front_frame_thickness / 2.;

  const G4double wagon_wheel_front_frame_R_inner = m_Params->get_double_param("wagon_wheel_front_frame_R_inner") * cm;
  const G4double wagon_wheel_front_frame_R_outer = m_Params->get_double_param("wagon_wheel_front_frame_R_outer") * cm;

  for (int ring_id = 0; ring_id <= n_radial_modules; ++ring_id)
  {
    G4double Rin = wagon_wheel_front_frame_R_inner;
    G4double Rout = wagon_wheel_front_frame_R_outer;

    if (ring_id > 0)
    {
      Rin = m_Params->get_double_param(
                boost::str(boost::format("wagon_wheel_front_frame_R_R%1%_outer") % (ring_id))) *
            cm;
    }
    if (ring_id < n_radial_modules)
    {
      Rout = m_Params->get_double_param(
                 boost::str(boost::format("wagon_wheel_front_frame_R_R%1%_inner") % (ring_id + 1))) *
             cm;
    }

    std::string name_base = boost::str(boost::format("%1%_%2%_Ring%3%") % GetName() % "wagon_wheel_front_frame" % ring_id);

    G4VSolid *solid_wagon_wheel_front_frame = new G4Tubs(
        name_base,
        Rin,
        Rout,
        wagon_wheel_front_frame_thickness / 2.,
        0, CLHEP::twopi);

    G4LogicalVolume *log_solid_wagon_wheel_front_frame = new G4LogicalVolume(solid_wagon_wheel_front_frame, material, name_base);
    m_LogicalVolumesSet.insert(log_solid_wagon_wheel_front_frame);

    assmeblyvol->AddPlacedVolume(log_solid_wagon_wheel_front_frame,
                                 g4vec_wagon_wheel_front_frame,
                                 nullptr);
    assert(m_DisplayAction);
    m_DisplayAction->AddVolume(log_solid_wagon_wheel_front_frame, "wagon_wheel");

  }  // for (int ring_id = 0; ring_id <= n_radial_modules; ++ring_id)

  ///////////////////////////////////////////////
  // wagon_wheel_front_frame spoke
  ///////////////////////////////////////////////
  for (int ring_id = 1; ring_id <= n_radial_modules; ++ring_id)
  {
    G4double Rout =
        m_Params->get_double_param(
            boost::str(boost::format("wagon_wheel_front_frame_R_R%1%_outer") % (ring_id))) *
        cm;
    G4double Rin =
        m_Params->get_double_param(
            boost::str(boost::format("wagon_wheel_front_frame_R_R%1%_inner") % (ring_id))) *
        cm;

    const G4double reduced_height = sqrt(Rout * Rout - wagon_wheel_front_frame_spoke_width / 2 * wagon_wheel_front_frame_spoke_width / 2);

    std::vector<G4TwoVector> vertexes;
    vertexes.push_back(G4TwoVector(-wagon_wheel_front_frame_spoke_width / 2, Rin));
    vertexes.push_back(G4TwoVector(+wagon_wheel_front_frame_spoke_width / 2, Rin));
    vertexes.push_back(G4TwoVector(+wagon_wheel_front_frame_spoke_width / 2, reduced_height));
    vertexes.push_back(G4TwoVector(-wagon_wheel_front_frame_spoke_width / 2, reduced_height));

    G4TwoVector zero(0, 0);

    std::string name_base_spoke = boost::str(boost::format("%1%_%2%_Ring%3%_spoke") % GetName() % "wagon_wheel_front_frame" % ring_id);

    G4VSolid *solid_wagon_wheel_front_frame_spoke = new G4ExtrudedSolid(name_base_spoke,
                                                                        vertexes,
                                                                        wagon_wheel_front_frame_thickness / 2.,
                                                                        zero, 1.0,
                                                                        zero, 1.0);
    G4LogicalVolume *log_solid_wagon_wheel_front_frame_spoke = new G4LogicalVolume(solid_wagon_wheel_front_frame_spoke, material, name_base_spoke);
    m_LogicalVolumesSet.insert(log_solid_wagon_wheel_front_frame_spoke);

    const G4double sector_dphi = CLHEP::twopi / n_sectors;
    for (int sector_id = 0; sector_id < n_sectors; ++sector_id)
    {
      G4Transform3D trans_spoke(CLHEP::HepRotationZ(wagon_wheel_sector_phi_offset + sector_dphi * sector_id), g4vec_wagon_wheel_front_frame);

      assmeblyvol->AddPlacedVolume(log_solid_wagon_wheel_front_frame_spoke,
                                   trans_spoke);
      assert(m_DisplayAction);
      m_DisplayAction->AddVolume(log_solid_wagon_wheel_front_frame_spoke, "wagon_wheel");

    }  //     for (int sector_id = 0; sector_id < n_sectors; ++sector_id)

  }  //  for (int ring_id = 0; ring_id < n_radial_modules; ++ring_id)

  ///////////////////////////////////////////////
  // wagon_wheel_rim_outer
  ///////////////////////////////////////////////
  if (Verbosity())
  {
    std::cout << __PRETTY_FUNCTION__ << " - wagon_wheel_rim_outer z_start = " << z_start << std::endl;
  }

  {
    const G4double wagon_wheel_rim_outer_Rin = m_Params->get_double_param("wagon_wheel_rim_outer_Rin") * cm;
    const G4double wagon_wheel_rim_outer_Rout = m_Params->get_double_param("wagon_wheel_rim_outer_Rout") * cm;
    const G4double wagon_wheel_rim_outer_thickness = m_Params->get_double_param("wagon_wheel_rim_outer_thickness") * cm;

    G4ThreeVector g4vec_wagon_wheel_rim_outer(0, 0, z_start + wagon_wheel_rim_outer_thickness / 2.);

    std::string name_base = boost::str(boost::format("%1%_wagon_wheel_rim_outer") % GetName());

    G4VSolid *solid_wagon_wheel = new G4Tubs(
        name_base,
        wagon_wheel_rim_outer_Rin,
        wagon_wheel_rim_outer_Rout,
        wagon_wheel_rim_outer_thickness / 2.,
        0, CLHEP::twopi);

    G4LogicalVolume *log_solid_wagon_wheel = new G4LogicalVolume(solid_wagon_wheel, material, name_base);
    m_LogicalVolumesSet.insert(log_solid_wagon_wheel);

    assmeblyvol->AddPlacedVolume(log_solid_wagon_wheel,
                                 g4vec_wagon_wheel_rim_outer,
                                 nullptr);
    assert(m_DisplayAction);
    m_DisplayAction->AddVolume(log_solid_wagon_wheel, "wagon_wheel");

  }  // wagon_wheel_rim_outer

  ///////////////////////////////////////////////
  // wagon_wheel_spoke
  ///////////////////////////////////////////////
  {
    const G4double wagon_wheel_spoke_width = m_Params->get_double_param("wagon_wheel_spoke_width") * cm;
    const G4double wagon_wheel_spoke_height_inner = m_Params->get_double_param("wagon_wheel_spoke_height_inner") * cm;
    const G4double wagon_wheel_spoke_height_outer = m_Params->get_double_param("wagon_wheel_spoke_height_outer") * cm;
    const G4double wagon_wheel_spoke_R_inner = m_Params->get_double_param("wagon_wheel_spoke_R_inner") * cm;
    const G4double wagon_wheel_spoke_R_outer = m_Params->get_double_param("wagon_wheel_spoke_R_outer") * cm;

    std::string name_base = boost::str(boost::format("%1%_wagon_wheel_spoke") % GetName());

    std::vector<G4TwoVector> vertexes;
    vertexes.push_back(G4TwoVector(0, wagon_wheel_spoke_R_inner));
    vertexes.push_back(G4TwoVector(0, wagon_wheel_spoke_R_outer));
    vertexes.push_back(G4TwoVector(wagon_wheel_spoke_height_outer, wagon_wheel_spoke_R_outer));
    vertexes.push_back(G4TwoVector(wagon_wheel_spoke_height_inner, wagon_wheel_spoke_R_inner));
    G4TwoVector zero(0, 0);

    G4VSolid *solid_wagon_wheel_spoke = new G4ExtrudedSolid(name_base,
                                                            vertexes,
                                                            wagon_wheel_spoke_width / 2.,
                                                            zero, 1.0,
                                                            zero, 1.0);
    G4LogicalVolume *log_solid_wagon_wheel_spoke = new G4LogicalVolume(solid_wagon_wheel_spoke, material, name_base);
    m_LogicalVolumesSet.insert(log_solid_wagon_wheel_spoke);

    G4ThreeVector g4vec_wagon_wheel_spoke(0, 0, z_start + wagon_wheel_spoke_width / 2.);

    const G4double sector_dphi = CLHEP::twopi / n_sectors;
    for (int sector_id = 0; sector_id < n_sectors; ++sector_id)
    {
      G4RotateY3D rotm_spoke(-90 * deg);
      G4Transform3D trans_spoke(CLHEP::HepRotationZ(wagon_wheel_sector_phi_offset + sector_dphi * sector_id),
                                g4vec_wagon_wheel_spoke);
      G4Transform3D trans_spoke_final = trans_spoke * rotm_spoke;

      assmeblyvol->AddPlacedVolume(log_solid_wagon_wheel_spoke,
                                   trans_spoke_final);
      assert(m_DisplayAction);
      m_DisplayAction->AddVolume(log_solid_wagon_wheel_spoke, "wagon_wheel");

    }  //     for (int sector_id = 0; sector_id < n_sectors; ++sector_id)

  }  // wagon_wheel_rim_outer
}

void PHG4TpcEndCapDetector::ConstructElectronics(G4AssemblyVolume *assmeblyvol,
                                                 G4double z_start)
{
  const int n_sectors = m_Params->get_int_param("n_sectors");
  assert(n_sectors >= 1);

  const G4double sector_dphi = CLHEP::twopi / n_sectors;

  const int n_radial_modules = m_Params->get_int_param("n_radial_modules");
  assert(n_radial_modules >= 1);
  const G4double wagon_wheel_sector_phi_offset = m_Params->get_double_param("wagon_wheel_sector_phi_offset_degree") * degree;
  const G4double wagon_wheel_spoke_width = m_Params->get_double_param("wagon_wheel_spoke_width") * cm;

  ///////////////////////////////////////////////
  // electronics_cooling_block_material ring
  ///////////////////////////////////////////////
  const G4double electronics_cooling_block_thickness = m_Params->get_double_param("electronics_cooling_block_thickness") * cm;
  if (electronics_cooling_block_thickness > 0)
  {
    if (Verbosity())
    {
      std::cout << __PRETTY_FUNCTION__ << " - electronics_cooling_block_material z_start = " << z_start << std::endl;
    }

    const std::string electronics_cooling_block_material_name(m_Params->get_string_param("electronics_cooling_block_material"));
    auto material = GetDetectorMaterial(electronics_cooling_block_material_name);
    if (material == nullptr)
    {
      std::cout << __PRETTY_FUNCTION__ << " Fatal Error: missing material " << m_Params->get_string_param("electronics_cooling_block_material_name") << std::endl;
      gSystem->Exit(1);
      exit(1);
    }

    G4ThreeVector g4vec_electronics_cooling_block(0, 0, z_start + electronics_cooling_block_thickness / 2.);

    const G4double electronics_cooling_block_R_inner = m_Params->get_double_param("electronics_cooling_block_R_inner") * cm;
    const G4double electronics_cooling_block_R_outer = m_Params->get_double_param("electronics_cooling_block_R_outer") * cm;

    for (int ring_id = 0; ring_id <= n_radial_modules; ++ring_id)
    {
      G4double Rin = electronics_cooling_block_R_inner;
      G4double Rout = electronics_cooling_block_R_outer;

      if (ring_id > 0)
      {
        Rin = m_Params->get_double_param(
                  boost::str(boost::format("electronics_cooling_block_R_R%1%_outer") % (ring_id))) *
              cm;
      }
      if (ring_id < n_radial_modules)
      {
        Rout = m_Params->get_double_param(
                   boost::str(boost::format("electronics_cooling_block_R_R%1%_inner") % (ring_id + 1))) *
               cm;
      }

      std::string name_base = boost::str(boost::format("%1%_%2%_Ring%3%") % GetName() % "electronics_cooling_block" % ring_id);

      const G4double spoke_phi = atan2(wagon_wheel_spoke_width, Rin);

      //      const G4double sector_dphi = CLHEP::twopi / n_sectors;

      G4VSolid *solid = new G4Tubs(
          name_base,
          Rin,
          Rout,
          electronics_cooling_block_thickness / 2.,
          spoke_phi, sector_dphi - 2 * spoke_phi);

      if (Verbosity())
      {
        std::cout << __PRETTY_FUNCTION__ << " - electronics_cooling_block " << name_base
                  << " Rin = " << Rin << " Rout = " << Rout
                  << " phi = " << spoke_phi << " to " << (sector_dphi - spoke_phi) << std::endl;
      }

      //
      G4LogicalVolume *log_vol = new G4LogicalVolume(solid, material, name_base);
      m_LogicalVolumesSet.insert(log_vol);

      for (int sector_id = 0; sector_id < n_sectors; ++sector_id)
      {
        G4Transform3D trans(
            CLHEP::HepRotationZ(wagon_wheel_sector_phi_offset + sector_dphi * sector_id),
            g4vec_electronics_cooling_block);
        //
        assmeblyvol->AddPlacedVolume(log_vol, trans);
        assert(m_DisplayAction);
        m_DisplayAction->AddVolume(log_vol, "cooling_block");

      }  //     for (int sector_id = 0; sector_id < n_sectors; ++sector_id)
    }    // for (int ring_id = 0; ring_id <= n_radial_modules; ++ring_id)
  }      // electronics_cooling_block_material  if (electronics_cooling_block_thickness>0)

  ///////////////////////////////////////////////
  // electronics
  ///////////////////////////////////////////////
  const G4double electronics_FEE_depth = m_Params->get_double_param("electronics_FEE_depth") * cm;
  const G4double electronics_FEE_Cu_thickness = m_Params->get_double_param("electronics_FEE_Cu_thickness") * cm;
  const G4double electronics_FEE_PCB_thickness = m_Params->get_double_param("electronics_FEE_PCB_thickness") * cm;
  const G4double electronics_FEE_Al_thickness = m_Params->get_double_param("electronics_FEE_Al_thickness") * cm;
  const G4double electronics_assemly_thickness = electronics_FEE_Cu_thickness + electronics_FEE_PCB_thickness + electronics_FEE_Al_thickness;

  if (m_Params->get_int_param("electronics_enable") != 0)
  {
    for (int ring_id = 1; ring_id <= n_radial_modules; ++ring_id)
    {
      const G4double Rout = m_Params->get_double_param(
                                boost::str(boost::format("electronics_cooling_block_R_R%1%_outer") % (ring_id))) *
                                cm -
                            electronics_assemly_thickness;
      const G4double Rin = m_Params->get_double_param(
                               boost::str(boost::format("electronics_cooling_block_R_R%1%_inner") % (ring_id))) *
                               cm +
                           electronics_assemly_thickness;
      const int nFEE = m_Params->get_int_param(boost::str(boost::format("electronics_nFEE_R%1%") % (ring_id)));

      if (nFEE <= 0)
      {
        std::cout << __PRETTY_FUNCTION__ << " warning : ignore FEE construction for module " << ring_id << " as "
                  << boost::str(boost::format("electronics_nFEE_R2%1%") % (ring_id)) << " = " << nFEE << std::endl;

        continue;
      }

      G4AssemblyVolume *assmeblyvol_electronics = new G4AssemblyVolume();
      G4double starting_electronics(0);
      std::string name_base = boost::str(boost::format("%1%_%2%_Ring%3%") % GetName() % "electronics" % ring_id);

      if (Verbosity())
      {
        std::cout << __PRETTY_FUNCTION__ << " - electronics G4_PCB z_start = " << z_start
                  << " starting_electronics = " << starting_electronics << std::endl;
      }
      starting_electronics -= electronics_FEE_PCB_thickness / 2.;
      G4ThreeVector g4vec_electronics;
      g4vec_electronics.set(starting_electronics, (Rout + Rin) * .5, z_start + electronics_FEE_depth / 2.);
      starting_electronics -= electronics_FEE_PCB_thickness / 2.;
      G4VSolid *solid_electronics = nullptr;
      solid_electronics = new G4Box(name_base + "_PCB",
                                    electronics_FEE_PCB_thickness / 2.,
                                    (Rout - Rin) / 2.,
                                    electronics_FEE_depth / 2.);

      G4LogicalVolume *log_electronics = nullptr;
      log_electronics = new G4LogicalVolume(solid_electronics, GetDetectorMaterial("FR4"), name_base + "_PCB");
      m_LogicalVolumesSet.insert(log_electronics);

      assmeblyvol_electronics->AddPlacedVolume(log_electronics,
                                               g4vec_electronics, nullptr);
      m_DisplayAction->AddVolume(log_electronics, "FR4");
      if (Verbosity())
      {
        std::cout << __PRETTY_FUNCTION__ << " - electronics G4_Cu z_start = " << z_start
                  << " starting_electronics = " << starting_electronics << std::endl;
      }
      starting_electronics -= electronics_FEE_Cu_thickness / 2.;
      g4vec_electronics.set(starting_electronics, (Rout + Rin) * .5, z_start + electronics_FEE_depth / 2.);
      starting_electronics -= electronics_FEE_Cu_thickness / 2.;

      solid_electronics = new G4Box(name_base + "_Cu",
                                    electronics_FEE_Cu_thickness / 2.,
                                    (Rout - Rin) / 2.,
                                    electronics_FEE_depth / 2.);

      log_electronics = new G4LogicalVolume(solid_electronics, GetDetectorMaterial("G4_Cu"), name_base + "_Cu");
      m_LogicalVolumesSet.insert(log_electronics);

      assmeblyvol_electronics->AddPlacedVolume(log_electronics,
                                               g4vec_electronics, nullptr);
      m_DisplayAction->AddVolume(log_electronics, "Cu");
      if (Verbosity())
      {
        std::cout << __PRETTY_FUNCTION__ << " - electronics Al z_start = " << z_start
                  << " starting_electronics = " << starting_electronics << std::endl;
      }
      starting_electronics -= electronics_FEE_Al_thickness / 2.;
      g4vec_electronics.set(starting_electronics, (Rout + Rin) * .5, z_start + electronics_FEE_depth / 2.);
      starting_electronics -= electronics_FEE_Al_thickness / 2.;

      solid_electronics = new G4Box(name_base + "_Al",
                                    electronics_FEE_Al_thickness / 2.,
                                    (Rout - Rin) / 2.,
                                    electronics_FEE_depth / 2.);

      log_electronics = new G4LogicalVolume(solid_electronics,
                                            GetDetectorMaterial("G4_Al"),
                                            name_base + "_Al");
      m_LogicalVolumesSet.insert(log_electronics);

      assmeblyvol_electronics->AddPlacedVolume(log_electronics,
                                               g4vec_electronics, nullptr);
      m_DisplayAction->AddVolume(log_electronics, "cooling_block");

      for (int sector_id = 0; sector_id < n_sectors; ++sector_id)
      {
        const G4double sector_phi_shift = wagon_wheel_sector_phi_offset + sector_dphi * sector_id;
        const G4double spoke_phi = atan2(wagon_wheel_spoke_width, Rin);
        const G4double board_dphi = (sector_dphi - 2 * spoke_phi) / (nFEE + 1);
        const G4double board_phi_start = sector_phi_shift + spoke_phi + board_dphi;

        for (int board_id = 0; board_id < nFEE; ++board_id)
        {
          G4Transform3D trans_electronic = G4RotateZ3D(board_phi_start + board_dphi * board_id);

          assmeblyvol->AddPlacedAssembly(assmeblyvol_electronics, trans_electronic);
        }
      }  //     for (int sector_id = 0; sector_id < n_sectors; ++sector_id)

    }  //  for (int ring_id = 0; ring_id < n_radial_modules; ++ring_id)
  }
}

//_______________________________________________________________
void PHG4TpcEndCapDetector::Print(const std::string &what) const
{
  std::cout << "PHG4TpcEndCap Detector:" << std::endl;
  if (what == "ALL" || what == "VOLUME")
  {
    std::cout << "Version 0.1" << std::endl;
    std::cout << "Parameters:" << std::endl;
    m_Params->Print();
  }
  return;
}
