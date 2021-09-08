//____________________________________________________________________________..
//
// This is a working template for the G4 Construct() method which needs to be implemented
// We wedge a method between the G4 Construct() to enable volume hierarchies on the macro
// so here it is called ConstructMe() but there is no functional difference
// Currently this installs a simple G4Box solid, creates a logical volume from it
// and places it. Put your own detector in place (just make sure all active volumes
// get inserted into the m_PhysicalVolumesSet)
//
// Rather than using hardcoded values you should consider using the parameter class
// Parameter names and defaults are set in PHG4TpcEndCapSubsystem::SetDefaultParameters()
// Only parameters defined there can be used (also to override in the macro)
// to avoids typos.
// IMPORTANT: parameters have no inherent units, there is a convention (cm/deg)
// but in any case you need to multiply them here with the correct CLHEP/G4 unit
//
// The place where you put your own detector is marked with
// //begin implement your own here://
// //end implement your own here://
// Do not forget to include the G4 includes for your volumes
//____________________________________________________________________________..

#include "PHG4TpcEndCapDetector.h"
#include "PHG4TpcEndCapDisplayAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>
#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4Color.hh>
#include <Geant4/G4ExtrudedSolid.hh>
#include <Geant4/G4IntersectionSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4MultiUnion.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4String.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4Transform3D.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4TwoVector.hh>
#include <Geant4/G4UserLimits.hh>
#include <Geant4/G4VPhysicalVolume.hh>
#include <Geant4/G4VSolid.hh>
#include <Geant4/G4VisAttributes.hh>

#include <CLHEP/Vector/RotationZ.h>
#include <boost/format.hpp>

#include <cassert>
#include <cmath>
#include <iostream>

class G4VSolid;
class PHCompositeNode;

using namespace std;

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

PHG4TpcEndCapDetector::
    ~PHG4TpcEndCapDetector()
{
  if (m_EndCapAssembly)
  {
    if (Verbosity())
    {
      cout << __PRETTY_FUNCTION__ << " delete m_EndCapAssembly" << endl;
    }

    delete m_EndCapAssembly;
  }
}

//_______________________________________________________________
int PHG4TpcEndCapDetector::IsInDetector(G4VPhysicalVolume *volume) const
{
  set<G4VPhysicalVolume *>::const_iterator iter = m_PhysicalVolumesSet.find(volume);
  if (iter != m_PhysicalVolumesSet.end())
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

  G4TranslateZ3D g4vec_front_z(
      m_Params->get_double_param("envelop_front_surface_z") * cm);

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
  G4AssemblyVolume *assmeblyvol = new G4AssemblyVolume();
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
  const int n_GEM_layers = m_Params->get_int_param("n_GEM_layers");

  for (int gem = 1; gem <= n_GEM_layers; gem++)
  {
    stringstream sid;
    sid << gem;

    //  GEM Copper 1.43 0.0005x6 64 0.134
    AddLayer(assmeblyvol, starting_z, G4String("GEMFrontCu") + G4String(sid.str()), "G4_Cu",
             0.0005 * cm, 64);

    //  GEM Kapton 28.6 0.005x3 64 0.034
    AddLayer(assmeblyvol, starting_z, G4String("GEMKapton") + G4String(sid.str()), "G4_KAPTON",
             0.005 * cm, 64);

    //  GEM Copper 1.43 0.0005x6 64 0.134
    AddLayer(assmeblyvol, starting_z, G4String("GEMBackCu") + G4String(sid.str()), "G4_Cu",
             0.0005 * cm, 64);

    //  GEM frames FR4 17.1 0.15x4 6.5 0.228
    //    AddLayer(assmeblyvol,starting_z,G4String("Frame") + G4String(sid.str()), "G10", 0.15 * cm,
    //             6.5);
    // sPHENIX TPC air gap
    starting_z += 0.2 * cm;
  }

  // 16 layer readout plane by TTM
  // https://indico.bnl.gov/event/8307/contributions/36744/attachments/27646/42337/R3-Review.pptx
  const int n_PCB_layers(16);
  // 35 um / layer Cu
  AddLayer(assmeblyvol, starting_z, G4String("PCBCu"), "G4_Cu", 0.0035 * cm * n_PCB_layers, 80);
  // 7 mil / layer board
  AddLayer(assmeblyvol, starting_z, "PCBBase", "FR4", 0.00254 * cm * 7 * n_PCB_layers, 100);

  ConstructWagonWheel(assmeblyvol, starting_z);
  ConstructElectronics(assmeblyvol, starting_z);

  return assmeblyvol;
}

void PHG4TpcEndCapDetector ::AddLayer(  //
    G4AssemblyVolume *assmeblyvol,
    G4double &z_start,
    std::string _name,         //! name base for this layer
    std::string _material,     //! material name in G4
    G4double _depth,           //! depth in G4 units
    double _percentage_filled  //! percentage filled//
)
{
  z_start += _depth / 2.;
  G4ThreeVector g4vec(0, 0, z_start);
  z_start += _depth / 2.;

  string name_base =
      boost::str(boost::format("%1%_Layer_%2%") % GetName() % _name);

  G4VSolid *solid_layer = new G4Tubs(
      name_base,
      m_Params->get_double_param("envelop_r_min") * cm,
      m_Params->get_double_param("envelop_r_max") * cm,
      _depth * _percentage_filled / 100. / 2.,
      0, CLHEP::twopi);

  auto material = G4Material::GetMaterial(_material);
  if (material == nullptr)
  {
    cout << __PRETTY_FUNCTION__ << " Fatal Error: missing material " << _material << endl;
    assert(material);
  }

  G4LogicalVolume *logical_layer = new G4LogicalVolume(solid_layer, material, name_base);

  assmeblyvol->AddPlacedVolume(logical_layer, g4vec, nullptr);

  assert(m_DisplayAction);
  m_DisplayAction->AddVolume(logical_layer, _material);

  return;
}

void PHG4TpcEndCapDetector::ConstructWagonWheel(G4AssemblyVolume *assmeblyvol,
                                                G4double &z_start)
{
  const int n_sectors = m_Params->get_int_param("n_sectors");
  assert(n_sectors >= 1);
  const int n_radial_modules = m_Params->get_int_param("n_radial_modules");
  assert(n_radial_modules >= 1);

  const string material_name(m_Params->get_string_param("wagon_wheel_material"));
  auto material = G4Material::GetMaterial(material_name);
  if (material == nullptr)
  {
    cout << __PRETTY_FUNCTION__ << " Fatal Error: missing material " << m_Params->get_string_param("wagon_wheel_material") << endl;
    assert(material);
  }
  const G4double wagon_wheel_sector_phi_offset = m_Params->get_double_param("wagon_wheel_sector_phi_offset_degree") * degree;

  ///////////////////////////////////////////////
  // wagon_wheel_front_frame ring
  ///////////////////////////////////////////////
  if (Verbosity())
  {
    cout << __PRETTY_FUNCTION__ << " - wagon_wheel_front_frame z_start = " << z_start << endl;
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

    string name_base = boost::str(boost::format("%1%_%2%_Ring%3%") % GetName() % "wagon_wheel_front_frame" % ring_id);

    G4VSolid *solid_wagon_wheel_front_frame = new G4Tubs(
        name_base,
        Rin,
        Rout,
        wagon_wheel_front_frame_thickness / 2.,
        0, CLHEP::twopi);

    G4LogicalVolume *log_solid_wagon_wheel_front_frame = new G4LogicalVolume(solid_wagon_wheel_front_frame, material, name_base);

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

    string name_base_spoke = boost::str(boost::format("%1%_%2%_Ring%3%_spoke") % GetName() % "wagon_wheel_front_frame" % ring_id);

    G4VSolid *solid_wagon_wheel_front_frame_spoke = new G4ExtrudedSolid(name_base_spoke,
                                                                        vertexes,
                                                                        wagon_wheel_front_frame_thickness / 2.,
                                                                        zero, 1.0,
                                                                        zero, 1.0);
    G4LogicalVolume *log_solid_wagon_wheel_front_frame_spoke = new G4LogicalVolume(solid_wagon_wheel_front_frame_spoke, material, name_base_spoke);

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
    cout << __PRETTY_FUNCTION__ << " - wagon_wheel_rim_outer z_start = " << z_start << endl;
  }

  {
    const G4double wagon_wheel_rim_outer_Rin = m_Params->get_double_param("wagon_wheel_rim_outer_Rin") * cm;
    const G4double wagon_wheel_rim_outer_Rout = m_Params->get_double_param("wagon_wheel_rim_outer_Rout") * cm;
    const G4double wagon_wheel_rim_outer_thickness = m_Params->get_double_param("wagon_wheel_rim_outer_thickness") * cm;

    G4ThreeVector g4vec_wagon_wheel_rim_outer(0, 0, z_start + wagon_wheel_rim_outer_thickness / 2.);

    string name_base = boost::str(boost::format("%1%_wagon_wheel_rim_outer") % GetName());

    G4VSolid *solid_wagon_wheel = new G4Tubs(
        name_base,
        wagon_wheel_rim_outer_Rin,
        wagon_wheel_rim_outer_Rout,
        wagon_wheel_rim_outer_thickness / 2.,
        0, CLHEP::twopi);

    G4LogicalVolume *log_solid_wagon_wheel = new G4LogicalVolume(solid_wagon_wheel, material, name_base);

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

    string name_base = boost::str(boost::format("%1%_wagon_wheel_spoke") % GetName());

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
                                                 G4double &z_start)
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
      cout << __PRETTY_FUNCTION__ << " - electronics_cooling_block_material z_start = " << z_start << endl;
    }

    const string electronics_cooling_block_material_name(m_Params->get_string_param("electronics_cooling_block_material"));
    auto material = G4Material::GetMaterial(electronics_cooling_block_material_name);
    if (material == nullptr)
    {
      cout << __PRETTY_FUNCTION__ << " Fatal Error: missing material " << m_Params->get_string_param("electronics_cooling_block_material_name") << endl;
      assert(material);
    }

    const G4double electronics_cooling_block_thickness = m_Params->get_double_param("electronics_cooling_block_thickness") * cm;

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

      string name_base = boost::str(boost::format("%1%_%2%_Ring%3%") % GetName() % "electronics_cooling_block" % ring_id);

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
        cout << __PRETTY_FUNCTION__ << " - electronics_cooling_block " << name_base
             << " Rin = " << Rin << " Rout = " << Rout
             << " phi = " << spoke_phi << " to " << (sector_dphi - spoke_phi) << endl;
      }

      //
      G4LogicalVolume *log_vol = new G4LogicalVolume(solid, material, name_base);
      for (int sector_id = 0; sector_id < n_sectors; ++sector_id)
      {
        G4Transform3D trans(
            CLHEP::HepRotationZ(wagon_wheel_sector_phi_offset + sector_dphi * sector_id),
            g4vec_electronics_cooling_block);
        //
        assmeblyvol->AddPlacedVolume(log_vol,
                                     trans);
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
    for (int ring_id = 1; ring_id <= n_radial_modules; ++ring_id)
    {
      const G4double Rout =
          m_Params->get_double_param(
              boost::str(boost::format("electronics_cooling_block_R_R%1%_outer") % (ring_id))) *
              cm -
          electronics_assemly_thickness;
      const G4double Rin =
          m_Params->get_double_param(
              boost::str(boost::format("electronics_cooling_block_R_R%1%_inner") % (ring_id))) *
              cm +
          electronics_assemly_thickness;
      const int nFEE =
          m_Params->get_int_param(
              boost::str(boost::format("electronics_nFEE_R%1%") % (ring_id)));

      if (nFEE <= 0)
      {
        cout << __PRETTY_FUNCTION__ << " warning : ignore FEE construction for module " << ring_id << " as "
             << boost::str(boost::format("electronics_nFEE_R2%1%") % (ring_id)) << " = " << nFEE << endl;

        continue;
      }

      G4AssemblyVolume *assmeblyvol_electronics = new G4AssemblyVolume();
      G4double starting_electronics(0);
      string name_base = boost::str(boost::format("%1%_%2%_Ring%3%") % GetName() % "electronics" % ring_id);

      {
        if (Verbosity())
        {
          cout << __PRETTY_FUNCTION__ << " - electronics G4_PCB z_start = " << z_start
               << " starting_electronics = " << starting_electronics << endl;
        }
        starting_electronics -= electronics_FEE_PCB_thickness / 2.;
        G4ThreeVector g4vec_electronics(starting_electronics, (Rout + Rin) * .5, z_start + electronics_FEE_depth / 2.);
        starting_electronics -= electronics_FEE_PCB_thickness / 2.;

        G4VSolid *solid_electronics = new G4Box(name_base + "_PCB",
                                                electronics_FEE_PCB_thickness / 2.,
                                                (Rout - Rin) / 2.,
                                                electronics_FEE_depth / 2.);

        G4LogicalVolume *log_electronics = new G4LogicalVolume(solid_electronics,
                                                               G4Material::GetMaterial("FR4"), name_base + "_PCB");

        assmeblyvol_electronics->AddPlacedVolume(log_electronics,
                                                 g4vec_electronics, nullptr);
        m_DisplayAction->AddVolume(log_electronics, "FR4");
      }
      {
        if (Verbosity())
        {
          cout << __PRETTY_FUNCTION__ << " - electronics G4_Cu z_start = " << z_start
               << " starting_electronics = " << starting_electronics << endl;
        }
        starting_electronics -= electronics_FEE_Cu_thickness / 2.;
        G4ThreeVector g4vec_electronics(starting_electronics, (Rout + Rin) * .5, z_start + electronics_FEE_depth / 2.);
        starting_electronics -= electronics_FEE_Cu_thickness / 2.;

        G4VSolid *solid_electronics = new G4Box(name_base + "_Cu",
                                                electronics_FEE_Cu_thickness / 2.,
                                                (Rout - Rin) / 2.,
                                                electronics_FEE_depth / 2.);

        G4LogicalVolume *log_electronics = new G4LogicalVolume(solid_electronics,
                                                               G4Material::GetMaterial("G4_Cu"), name_base + "_Cu");

        assmeblyvol_electronics->AddPlacedVolume(log_electronics,
                                                 g4vec_electronics, nullptr);
        m_DisplayAction->AddVolume(log_electronics, "Cu");
      }
      {
        if (Verbosity())
        {
          cout << __PRETTY_FUNCTION__ << " - electronics Al z_start = " << z_start
               << " starting_electronics = " << starting_electronics << endl;
        }
        starting_electronics -= electronics_FEE_Al_thickness / 2.;
        G4ThreeVector g4vec_electronics(starting_electronics, (Rout + Rin) * .5, z_start + electronics_FEE_depth / 2.);
        starting_electronics -= electronics_FEE_Al_thickness / 2.;

        G4VSolid *solid_electronics = new G4Box(name_base + "_Al",
                                                electronics_FEE_Al_thickness / 2.,
                                                (Rout - Rin) / 2.,
                                                electronics_FEE_depth / 2.);

        G4LogicalVolume *log_electronics = new G4LogicalVolume(solid_electronics,
                                                               G4Material::GetMaterial("G4_Al"), name_base + "_Al");

        assmeblyvol_electronics->AddPlacedVolume(log_electronics,
                                                 g4vec_electronics, nullptr);
        m_DisplayAction->AddVolume(log_electronics, "cooling_block");
      }

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

//_______________________________________________________________
void PHG4TpcEndCapDetector::Print(const std::string &what) const
{
  cout << "PHG4TpcEndCap Detector:" << endl;
  if (what == "ALL" || what == "VOLUME")
  {
    cout << "Version 0.1" << endl;
    cout << "Parameters:" << endl;
    m_Params->Print();
  }
  return;
}
