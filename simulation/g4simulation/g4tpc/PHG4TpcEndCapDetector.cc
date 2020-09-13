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
  //  //begin implement your own here://
  //  // Do not forget to multiply the parameters with their respective CLHEP/G4 unit !
  //  double xdim = m_Params->get_double_param("size_x") * cm;
  //  double ydim = m_Params->get_double_param("size_y") * cm;
  //  double zdim = m_Params->get_double_param("size_z") * cm;
  //  G4VSolid *solidbox = new G4Box("PHG4TpcEndCapSolid", xdim / 2., ydim / 2., zdim / 2.);
  //  G4LogicalVolume *logical = new G4LogicalVolume(solidbox, G4Material::GetMaterial(m_Params->get_string_param("material")), "PHG4TpcEndCapLogical");
  //
  //  //  G4VisAttributes *vis = new G4VisAttributes(G4Color(G4Colour::Grey()));  // grey is good to see the tracks in the display
  //  //  vis->SetForceSolid(true);
  //  //  logical->SetVisAttributes(vis);
  //  G4RotationMatrix *rotm = new G4RotationMatrix();
  //  rotm->rotateX(m_Params->get_double_param("rot_x") * deg);
  //  rotm->rotateY(m_Params->get_double_param("rot_y") * deg);
  //  rotm->rotateZ(m_Params->get_double_param("rot_z") * deg);
  //
  //  G4VPhysicalVolume *phy = new G4PVPlacement(
  //      rotm,
  //      G4ThreeVector(m_Params->get_double_param("place_x") * cm,
  //                    m_Params->get_double_param("place_y") * cm,
  //                    m_Params->get_double_param("place_z") * cm),
  //      logical, "PHG4TpcEndCap", logicWorld, 0, false, OverlapCheck());
  //  // add it to the list of placed volumes so the IsInDetector method
  //  // picks them up
  //  m_PhysicalVolumesSet.insert(phy);
  //  //end implement your own here://

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

  //  PCB Kapton 28.6 0.005 100 0.017
  AddLayer(assmeblyvol, starting_z, G4String("PCBKapton"), "G4_KAPTON", 0.005 * cm, 100);

  //  PCB Copper 1.43 0.0005 80 0.028
  AddLayer(assmeblyvol, starting_z, G4String("PCBCu"), "G4_Cu", 0.0005 * cm, 80);

  //  Facesheet FR4 17.1 0.025x2 100 0.292
  AddLayer(assmeblyvol, starting_z, "Facesheet", "G10", 0.025 * 2 * cm, 100);

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
  z_start += _depth / 2;
  G4ThreeVector g4vec(0, 0, z_start);
  z_start += _depth / 2;

  string name_base =
      boost::str(boost::format("%1%_Layer_%2%") % GetName() % _name);

  G4VSolid *solid_grid = new G4Tubs(
      name_base,
      m_Params->get_double_param("envelop_r_min") * cm,
      m_Params->get_double_param("envelop_r_max") * cm,
      _depth * _percentage_filled / 2,
      0, CLHEP::twopi);

  auto material = G4Material::GetMaterial(_material);
  if (material == nullptr)
  {
    cout << __PRETTY_FUNCTION__ << " Fatal Error: missing material " << _material << endl;
    assert(material);
  }

  G4LogicalVolume *logical_grid = new G4LogicalVolume(solid_grid, material, name_base);

  assmeblyvol->AddPlacedVolume(logical_grid, g4vec, nullptr);

  assert(m_DisplayAction);
  m_DisplayAction->AddVolume(logical_grid, _material);

  return;
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
