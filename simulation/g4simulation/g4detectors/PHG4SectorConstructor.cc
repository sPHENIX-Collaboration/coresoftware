/*!
 * \file PHG4SectorConstructor.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.6 $
 * \date $Date: 2014/07/31 15:52:37 $
 */

#include "PHG4SectorConstructor.h"
#include "PHG4SectorDisplayAction.h"

#include <g4main/PHG4Detector.h>
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>      // for PHG4Subsystem

#include <Geant4/G4Box.hh>
#include <Geant4/G4DisplacedSolid.hh>     // for G4DisplacedSolid
#include <Geant4/G4Exception.hh>          // for G4Exception
#include <Geant4/G4ExceptionSeverity.hh>  // for FatalException, JustWarning
#include <Geant4/G4IntersectionSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4PhysicalConstants.hh>  // for pi
#include <Geant4/G4Sphere.hh>
#include <Geant4/G4String.hh>
#include <Geant4/G4SystemOfUnits.hh>  // for cm, um, perCent
#include <Geant4/G4ThreeVector.hh>    // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>    // for G4Transform3D, G4RotateX3D
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Types.hh>  // for G4int

#include <algorithm>  // for max
#include <cassert>
#include <climits>
#include <cmath>
#include <iostream>
#include <sstream>

class G4Material;

PHG4Sector::PHG4SectorConstructor::PHG4SectorConstructor(const std::string &name, PHG4Subsystem *subsys)
  : overlapcheck_sector(false)
  , name_base(name)
  , m_DisplayAction(dynamic_cast<PHG4SectorDisplayAction *>(subsys->GetDisplayAction()))
  , m_Verbosity(0)
{
}

void PHG4Sector::PHG4SectorConstructor::Construct_Sectors(G4LogicalVolume *WorldLog)
{
  // geometry checks
  if (geom.get_total_thickness() == 0)
  {
    G4Exception(
        (std::string("PHG4SectorConstructor::Construct_Sectors::") + (name_base)).c_str(),
        __FILE__, FatalException,
        " detector configured with zero thickness!");
  }

  if (geom.get_min_polar_angle() == geom.get_max_polar_angle())
  {
    G4Exception(
        (std::string("PHG4SectorConstructor::Construct_Sectors::") + (name_base)).c_str(),
        __FILE__, FatalException, "min_polar_angle = max_polar_angle!");
  }
  if (geom.get_min_polar_angle() > geom.get_max_polar_angle())
  {
    G4Exception(
        (std::string("PHG4SectorConstructor::Construct_Sectors::") + (name_base)).c_str(),
        __FILE__, JustWarning,
        "min and max polar angle got reversed. Correcting them.");

    const double t = geom.get_max_polar_angle();
    geom.set_max_polar_angle(geom.get_min_polar_angle());
    geom.set_min_polar_angle(t);
  }
  if ((geom.get_min_polar_edge() == Sector_Geometry::kFlatEdge or geom.get_max_polar_edge() == Sector_Geometry::kFlatEdge) and geom.get_N_Sector() <= 2)
  {
    G4Exception(
        (std::string("PHG4SectorConstructor::Construct_Sectors::") + (name_base)).c_str(),
        __FILE__, FatalException,
        "can NOT use flat edge for single or double sector detector!");
  }

  const G4Transform3D transform_Det_to_Hall =
      G4RotateX3D(-geom.get_normal_polar_angle()) * G4TranslateZ3D(
                                                        geom.get_normal_start() + geom.get_total_thickness() / 2);

  const G4Transform3D transform_Hall_to_Det(transform_Det_to_Hall.inverse());

  // during GDML export, numerical value may change at the large digit and go beyond the 0 or pi limit.
  // therefore recess 0/pi by a small amount to avoid such problem
  static const double epsilon = std::numeric_limits<float>::epsilon();
  const double sph_min_polar_angle =
      (geom.get_min_polar_edge() == Sector_Geometry::kConeEdge) ? geom.get_min_polar_angle() : (0 + epsilon);
  const double sph_max_polar_angle =
      (geom.get_max_polar_edge() == Sector_Geometry::kConeEdge) ? geom.get_max_polar_angle() : (pi - epsilon);

  G4VSolid *SecConeBoundary_Hall = new G4Sphere("SecConeBoundary_Hall",                                           //
                                                geom.get_normal_start(), geom.get_max_R(),                        // G4double pRmin, G4double pRmax,
                                                pi / 2 - pi / geom.get_N_Sector(), 2 * pi / geom.get_N_Sector(),  //  G4double pSPhi, G4double pDPhi,
                                                sph_min_polar_angle, sph_max_polar_angle - sph_min_polar_angle    //G4double pSTheta, G4double pDTheta
  );

  G4VSolid *SecConeBoundary_Det = new G4DisplacedSolid("SecConeBoundary_Det",
                                                       SecConeBoundary_Hall, transform_Hall_to_Det);

  G4VSolid *Boundary_Det = SecConeBoundary_Det;

  if (geom.get_min_polar_edge() != Sector_Geometry::kConeEdge or geom.get_max_polar_edge() != Sector_Geometry::kConeEdge)
  {
    // build a flat edge

    const double sph_min_polar_size =
        (geom.get_min_polar_edge() == Sector_Geometry::kConeEdge) ? geom.get_max_R() : geom.get_normal_start() * tan(geom.get_normal_polar_angle() - geom.get_min_polar_angle());
    const double sph_max_polar_size =
        (geom.get_max_polar_edge() == Sector_Geometry::kConeEdge) ? geom.get_max_R() : geom.get_normal_start() * tan(geom.get_max_polar_angle() - geom.get_normal_polar_angle());

    G4VSolid *BoxBoundary_Det = new G4Box("BoxBoundary_Det",
                                          geom.get_max_R(), (sph_min_polar_size + sph_max_polar_size) / 2,
                                          geom.get_total_thickness());
    G4VSolid *BoxBoundary_Det_Place = new G4DisplacedSolid(
        "BoxBoundary_Det_Place", BoxBoundary_Det, nullptr,
        G4ThreeVector(0, (sph_max_polar_size - sph_min_polar_size) / 2, 0));

    Boundary_Det = new G4IntersectionSolid("Boundary_Det",
                                           BoxBoundary_Det_Place, SecConeBoundary_Det);
  }

  G4VSolid *DetectorBox_Det = Construct_Sectors_Plane("DetectorBox_Det",
                                                      -geom.get_total_thickness() / 2, geom.get_total_thickness(),
                                                      Boundary_Det);

  G4Material *p_mat = PHG4Detector::GetDetectorMaterial(geom.get_material());

  G4LogicalVolume *DetectorLog_Det = new G4LogicalVolume(DetectorBox_Det,  //
                                                         p_mat, name_base + "_Log");
  RegisterLogicalVolume(DetectorLog_Det);

  for (G4int sec = 0; sec < geom.get_N_Sector(); sec++)
  {
    RegisterPhysicalVolume(
        new G4PVPlacement(
            G4RotateZ3D(2 * pi / geom.get_N_Sector() * sec) * transform_Det_to_Hall, DetectorLog_Det,
            name_base + "_Physical", WorldLog, false, sec, overlapcheck_sector));
  }

  // construct layers
  double z_start = -geom.get_total_thickness() / 2;

  for (const auto &l : geom.layer_list)
  {
    if (l.percentage_filled > 100. || l.percentage_filled < 0)
    {
      std::ostringstream strstr;
      strstr << name_base << " have invalid layer ";
      strstr << l.name << " with percentage_filled =" << l.percentage_filled;

      G4Exception(
          (std::string("PHG4SectorConstructor::Construct_Sectors::") + (name_base)).c_str(), __FILE__, FatalException,
          strstr.str().c_str());
    }

    const std::string layer_name = name_base + "_" + l.name;

    G4VSolid *LayerSol_Det = Construct_Sectors_Plane(layer_name + "_Sol",
                                                     z_start, l.depth * l.percentage_filled * perCent, Boundary_Det);

    G4LogicalVolume *LayerLog_Det = new G4LogicalVolume(LayerSol_Det,  //
                                                        PHG4Detector::GetDetectorMaterial(l.material), layer_name + "_Log");
    RegisterLogicalVolume(LayerLog_Det);

    RegisterPhysicalVolume(
        new G4PVPlacement(nullptr, G4ThreeVector(), LayerLog_Det,
                          layer_name + "_Physical", DetectorLog_Det, false, 0, overlapcheck_sector),
        l.active);

    z_start += l.depth;
  }

  if (std::abs(z_start - geom.get_total_thickness() / 2) > 1 * um)
  {
    std::ostringstream strstr;
    strstr << name_base << " - accumulated thickness = "
           << (z_start + geom.get_total_thickness() / 2) / um
           << " um expected thickness = " << geom.get_total_thickness() / um
           << " um";
    G4Exception(
        (std::string("PHG4SectorConstructor::Construct_Sectors::") + (name_base)).c_str(),
        __FILE__, FatalException, strstr.str().c_str());
  }

  m_DisplayAction->AddVolume(DetectorLog_Det, "DetectorBox");
  if (Verbosity() > 1)
  {
    std::cout << "PHG4SectorConstructor::Construct_Sectors::" << name_base
              << " - total thickness = " << geom.get_total_thickness() / cm << " cm"
              << std::endl;
    std::cout << "PHG4SectorConstructor::Construct_Sectors::" << name_base << " - "
              << map_log_vol.size() << " logical volume constructed" << std::endl;
    std::cout << "PHG4SectorConstructor::Construct_Sectors::" << name_base << " - "
              << map_phy_vol.size() << " physical volume constructed; "
              << map_active_phy_vol.size() << " is active." << std::endl;
  }
}

G4VSolid *
PHG4Sector::PHG4SectorConstructor::Construct_Sectors_Plane(  //
    const std::string &name,                                 //
    const double start_z,                                    //
    const double thickness,                                  //
    G4VSolid *SecConeBoundary_Det                            //
)
{
  assert(SecConeBoundary_Det);

  G4VSolid *Sol_Raw = new G4Tubs(name + "_Raw",     //const G4String& pName,
                                 0,                 //      G4double pRMin,
                                 geom.get_max_R(),  //      G4double pRMax,
                                 thickness / 2,     //      G4double pDz,
                                 0,                 //      G4double pSPhi,
                                 2 * pi             //      G4double pDPhi
  );

  G4VSolid *Sol_Place = new G4DisplacedSolid(name + "_Place", Sol_Raw, nullptr,
                                             G4ThreeVector(0, 0, start_z + thickness / 2));

  G4VSolid *Sol = new G4IntersectionSolid(name, Sol_Place,
                                          SecConeBoundary_Det);

  return Sol;
  //  return Sol_Place;
}

void PHG4Sector::Sector_Geometry::SetDefault()
{
  N_Sector = 8;
  material = "G4_AIR";

  //! polar angle for the normal vector
  normal_polar_angle = 0;

  //! polar angle for edges
  min_polar_angle = .1;

  //! polar angle for edges
  max_polar_angle = .3;

  //! distance that detector starts from the normal direction
  normal_start = 305 * cm;

  min_polar_edge = kConeEdge;

  max_polar_edge = kConeEdge;
}

G4LogicalVolume *
PHG4Sector::PHG4SectorConstructor::RegisterLogicalVolume(G4LogicalVolume *v)
{
  if (!v)
  {
    std::cout
        << "PHG4SectorConstructor::RegisterVolume - Error - invalid volume!"
        << std::endl;
    return v;
  }
  if (map_log_vol.find(v->GetName()) != map_log_vol.end())
  {
    std::cout << "PHG4SectorConstructor::RegisterVolume - Warning - replacing "
              << v->GetName() << std::endl;
  }

  map_log_vol[v->GetName()] = v;

  return v;
}

G4PVPlacement *
PHG4Sector::PHG4SectorConstructor::RegisterPhysicalVolume(G4PVPlacement *v,
                                                          const bool active)
{
  if (!v)
  {
    std::cout
        << "PHG4SectorConstructor::RegisterPhysicalVolume - Error - invalid volume!"
        << std::endl;
    return v;
  }

  phy_vol_idx_t id(v->GetName(), v->GetCopyNo());

  if (map_phy_vol.find(id) != map_phy_vol.end())
  {
    std::cout
        << "PHG4SectorConstructor::RegisterPhysicalVolume - Warning - replacing "
        << v->GetName() << "[" << v->GetCopyNo() << "]" << std::endl;
  }

  map_phy_vol[id] = v;

  if (active)
    map_active_phy_vol[id] = v;

  return v;
}

double
PHG4Sector::Sector_Geometry::get_total_thickness() const
{
  double sum = 0;
  for (const auto &it : layer_list)
  {
    sum += it.depth;
  }
  return sum;
}

double
PHG4Sector::Sector_Geometry::get_max_R() const
{
  // Geometry check
  assert(std::abs(min_polar_angle - normal_polar_angle) < pi / 2);
  assert(std::abs(max_polar_angle - normal_polar_angle) < pi / 2);

  if (N_Sector <= 2)
  {
    // Geometry check
    if (cos(min_polar_angle + normal_polar_angle) <= 0)
    {
      std::stringstream strstr;
      strstr << "Geometry check failed. " << std::endl;
      strstr << "normal_polar_angle = " << normal_polar_angle << std::endl;
      strstr << "min_polar_angle = " << min_polar_angle << std::endl;
      strstr << "max_polar_angle = " << max_polar_angle << std::endl;
      strstr << "cos(min_polar_angle + normal_polar_angle) = "
             << cos(min_polar_angle + normal_polar_angle) << std::endl;

      G4Exception("Sector_Geometry::get_max_R", __FILE__, FatalException,
                  strstr.str().c_str());
    }
    if (cos(max_polar_angle + normal_polar_angle) <= 0)
    {
      std::stringstream strstr;
      strstr << "Geometry check failed. " << std::endl;
      strstr << "normal_polar_angle = " << normal_polar_angle << std::endl;
      strstr << "min_polar_angle = " << min_polar_angle << std::endl;
      strstr << "max_polar_angle = " << max_polar_angle << std::endl;
      strstr << "cos(max_polar_angle + normal_polar_angle) = "
             << cos(max_polar_angle + normal_polar_angle) << std::endl;

      G4Exception("Sector_Geometry::get_max_R", __FILE__, FatalException,
                  strstr.str().c_str());
    }

    const double max_tan_angle = std::max(
        std::abs(tan(min_polar_angle + normal_polar_angle)),
        std::abs(tan(max_polar_angle + normal_polar_angle)));

    return (get_normal_start() + get_total_thickness()) * sqrt(1 + max_tan_angle * max_tan_angle);
  }
  else
  {
    const double max_angle = std::max(
        std::abs(min_polar_angle - normal_polar_angle),
        std::abs(max_polar_angle - normal_polar_angle));

    return (get_normal_start() + get_total_thickness()) * sqrt(1 + pow(tan(max_angle), 2) + pow(tan(2 * pi / N_Sector), 2)) * 2;
  }
}

//! add Entrace window and drift volume
//! Ref: P. Abbon et al. The COMPASS experiment at CERN. Nucl. Instrum. Meth., A577:455
//! 518, 2007. arXiv:hep-ex/0703049, doi:10.1016/j.nima.2007.03.026. 3
void PHG4Sector::Sector_Geometry::AddLayers_DriftVol_COMPASS(const double drift_vol_thickness)
{
  //  (drift chamber) is enclosed by two Mylarr [52] cathode foils of 25um thickness,
  //  coated with about 10um of graphite,

  AddLayer("EntranceWindow", "G4_MYLAR", 25 * um, false, 100);
  AddLayer("Cathode", "G4_GRAPHITE", 10 * um, false, 100);
  AddLayer("DrftVol", material, drift_vol_thickness, true, 100);
}

//! add HBD GEM to the layer list,
//! Ref: W. Anderson et al. Design, Construction, Operation and Performance of a Hadron
//! Blind Detector for the PHENIX Experiment. Nucl. Instrum. Meth., A646:35 58, 2011.
//! arXiv:1103.4277, doi:10.1016/j.nima.2011.04.015. 3.5.1
void PHG4Sector::Sector_Geometry::AddLayers_HBD_GEM(const int n_GEM_layers)
{
  // Internal HBD structure
  // From doi:10.1016/j.nima.2011.04.015
  // Component Material X0 (cm) Thickness (cm) Area (%) Rad. Length (%)
  //  Mesh SS 1.67 0.003 11.5 0.021  <- not used for GEMs trackers
  //  AddLayer("Mesh", "Steel",
  //          0.003 * cm, false, 11.5);

  //  //  GEM frames FR4 17.1 0.15x4 6.5 0.228 <- not used for GEMs trackers
  //  AddLayer("Frame0", "G10",
  //          0.15 * cm, false, 6.5);

  for (int gem = 1; gem <= n_GEM_layers; gem++)
  {
    std::stringstream sid;
    sid << gem;

    //  GEM Copper 1.43 0.0005x6 64 0.134
    AddLayer(G4String("GEMFrontCu") + G4String(sid.str()), "G4_Cu",
             0.0005 * cm, false, 64);

    //  GEM Kapton 28.6 0.005x3 64 0.034
    AddLayer(G4String("GEMKapton") + G4String(sid.str()), "G4_KAPTON",
             0.005 * cm, false, 64);

    //  GEM Copper 1.43 0.0005x6 64 0.134
    AddLayer(G4String("GEMBackCu") + G4String(sid.str()), "G4_Cu",
             0.0005 * cm, false, 64);

    //  GEM frames FR4 17.1 0.15x4 6.5 0.228
    AddLayer(G4String("Frame") + G4String(sid.str()), "G10", 0.15 * cm, false,
             6.5);
  }

  //  PCB Kapton 28.6 0.005 100 0.017
  AddLayer(G4String("PCBKapton"), "G4_KAPTON", 0.005 * cm, false, 100);

  //  PCB Copper 1.43 0.0005 80 0.028
  AddLayer(G4String("PCBCu"), "G4_Cu", 0.0005 * cm, false, 80);

  //  Facesheet FR4 17.1 0.025x2 100 0.292
  AddLayer("Facesheet", "G10", 0.025 * 2 * cm, false, 100);

  //  Panel core Honeycomb 8170 1.905 100 0.023 <- very thin-X0 stuff, ignore

  //  Total vessel 0.82
  //  Readout
}

//! add HBD readout,
//! Ref: W. Anderson et al. Design, Construction, Operation and Performance of a Hadron
//! Blind Detector for the PHENIX Experiment. Nucl. Instrum. Meth., A646:35 58, 2011.
//! arXiv:1103.4277, doi:10.1016/j.nima.2011.04.015. 3.5.1
void PHG4Sector::Sector_Geometry::AddLayers_HBD_Readout()
{
  //  Readout
  //  Readout board FR4/copper 17.1/1.43 0.05/0.001 100 0.367
  AddLayer(G4String("ReadoutFR4"), "G10", 0.05 * cm, false, 100);
  AddLayer(G4String("ReadoutCu"), "G4_Cu", 0.001 * cm, false, 100);

  //  Preamps + sockets Copper 1.43 0.0005 100 0.66
  //  Total readout 1.03
  AddLayer(G4String("SocketsCu"), "G4_Cu", 0.0005 * cm, false, 100);
}

//! Rough AeroGel detector
//! Ref: T. Iijima et al. A Novel type of proximity focusing RICH counter with multiple
//! refractive index aerogel radiator. Nucl. Instrum. Meth., A548:383 390, 2005. arXiv:
//! physics/0504220, doi:10.1016/j.nima.2005.05.030
void PHG4Sector::Sector_Geometry::AddLayers_AeroGel_ePHENIX(const double radiator_length,
                                                            const double expansion_length, std::string radiator)
{
  if (radiator == "Default")
  {
    radiator = "ePHENIX_AeroGel";
  }

  //  Readout board FR4/copper 17.1/1.43 0.05/0.001 100 0.367
  AddLayer("EntranceWindow", "G10", 0.05 * cm, false, 100);
  AddLayer("AeroGel", radiator, radiator_length, false, 100);
  AddLayer("ExpansionVol", "G4_AIR", expansion_length, false, 100);

  //Some readout
  AddLayer(G4String("ReadoutFR4"), "G10", 0.05 * cm, false, 100);
  AddLayer(G4String("ReadoutCu"), "G4_Cu", 0.001 * cm, false, 100);
  AddLayer(G4String("SocketsCu"), "G4_Cu", 0.0005 * cm, false, 100);
}
