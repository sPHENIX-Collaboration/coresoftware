/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.7 $$
 * \date $$Date: 2015/02/10 15:39:26 $$
 */
#include "PHG4SpacalDetector.h"

#include "PHG4CylinderGeomContainer.h"
#include "PHG4SpacalDisplayAction.h"

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/recoConsts.h>

#include <g4gdml/PHG4GDMLConfig.hh>
#include <g4gdml/PHG4GDMLUtility.hh>

#include <TSystem.h>

#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4String.hh>  // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>  // for G4Transform3D, G4RotateZ3D
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Types.hh>  // for G4double
#include <Geant4/G4UserLimits.hh>

#include <cassert>
#include <cstdlib>   // for exit
#include <iostream>  // for operator<<, basic_ostream
#include <sstream>

class PHG4CylinderGeom;

//_______________________________________________________________
//note this inactive thickness is ~1.5% of a radiation length
PHG4SpacalDetector::PHG4SpacalDetector(PHG4Subsystem *subsys,
                                       PHCompositeNode *Node,
                                       const std::string &dnam,
                                       PHParameters *parameters,
                                       const int lyr,
                                       bool init_geom)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4SpacalDisplayAction *>(subsys->GetDisplayAction()))
  , layer(lyr)
{
  if (init_geom)
  {
    _geom = new SpacalGeom_t();
    if (_geom == nullptr)
    {
      std::cout << "PHG4SpacalDetector::Constructor - Fatal Error - invalid geometry object!" << std::endl;
      gSystem->Exit(1);
      exit(1);
    }
    assert(parameters);
    _geom->ImportParameters(*parameters);

    //  _geom->Print();
  }

  gdml_config = PHG4GDMLUtility::GetOrMakeConfigNode(Node);
  assert(gdml_config);
}

PHG4SpacalDetector::~PHG4SpacalDetector()
{
  // deleting nullptr pointers is allowed (results in NOOP)
  // so checking for not null before deleting is not needed
  delete fiber_core_step_limits;
  delete _geom;
}

//_______________________________________________________________
int PHG4SpacalDetector::IsInCylinderActive(const G4VPhysicalVolume *volume)
{
  //  std::cout << "checking detector" << std::endl;
  if (active && fiber_core_vol.find(volume) != fiber_core_vol.end())
  {
    //      return fiber_core_vol.find(volume)->second;
    return FIBER_CORE;
  }
  if (absorberactive)
  {
    if (fiber_vol.find(volume) != fiber_vol.end())
    {
      return FIBER_CLADING;
    }

    if (block_vol.find(volume) != block_vol.end())
    {
      return ABSORBER;
    }

    if (calo_vol.find(volume) != calo_vol.end())
    {
      return SUPPORT;
    }
  }
  return INACTIVE;
}

//_______________________________________________________________
void PHG4SpacalDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  assert(_geom);

  fiber_core_step_limits = new G4UserLimits(_geom->get_fiber_core_step_size() * cm);

  Verbosity(_geom->get_construction_verbose());

  if ((Verbosity() > 0))
  {
    std::cout << "PHG4SpacalDetector::Construct::" << GetName()
              << " - Start. Print Geometry:" << std::endl;
    Print();
  }

  if ((_geom->get_zmin() * cm + _geom->get_zmax() * cm) / 2 != _geom->get_zpos() * cm)
  {
    std::cout << "PHG4SpacalDetector::Construct - ERROR - not yet support unsymmetric system. Let me know if you need it. - Jin" << std::endl;
    _geom->Print();
    gSystem->Exit(-1);
  }
  if (_geom->get_zmin() * cm >= _geom->get_zmax() * cm)
  {
    std::cout << "PHG4SpacalDetector::Construct - ERROR - zmin >= zmax!" << std::endl;
    _geom->Print();
    gSystem->Exit(-1);
  }

  G4Tubs *cylinder_solid = new G4Tubs(G4String(GetName()),
                                      _geom->get_radius() * cm, _geom->get_max_radius() * cm,
                                      _geom->get_length() * cm / 2.0, 0, twopi);

  recoConsts *rc = recoConsts::instance();
  G4Material *cylinder_mat = GetDetectorMaterial(rc->get_StringFlag("WorldMaterial"));
  assert(cylinder_mat);

  G4LogicalVolume *cylinder_logic = new G4LogicalVolume(cylinder_solid, cylinder_mat, GetName(), nullptr, nullptr, nullptr);
  GetDisplayAction()->AddVolume(cylinder_logic, "SpacalCylinder");
  if (!m_CosmicSetupFlag)
  {
    new G4PVPlacement(nullptr, G4ThreeVector(_geom->get_xpos() * cm, _geom->get_ypos() * cm, _geom->get_zpos() * cm),
                      cylinder_logic, GetName(),
                      logicWorld, false, 0, OverlapCheck());
  }
  // install sectors
  if (_geom->get_sector_map().size() == 0)
  {
    _geom->init_default_sector_map();
  }

  if ((Verbosity() > 0))
  {
    std::cout << "PHG4SpacalDetector::Construct::" << GetName()
              << " - start constructing " << _geom->get_sector_map().size() << " sectors in total. " << std::endl;
    Print();
  }

  std::pair<G4LogicalVolume *, G4Transform3D> psec = Construct_AzimuthalSeg();
  G4LogicalVolume *sec_logic = psec.first;
  const G4Transform3D &sec_trans = psec.second;

  for (const SpacalGeom_t::sector_map_t::value_type &val : _geom->get_sector_map())
  {
    const int sec = val.first;
    const double rot = val.second;

    G4Transform3D sec_place = G4RotateZ3D(rot) * sec_trans;

    std::ostringstream name;
    name << GetName() << "_sec" << sec;
    G4PVPlacement *calo_phys = nullptr;
    if (m_CosmicSetupFlag)
    {
      calo_phys = new G4PVPlacement(nullptr, G4ThreeVector(0, -(_geom->get_radius()) * cm, 0), sec_logic,
                                    G4String(name.str()), logicWorld, false, sec,
                                    OverlapCheck());
    }
    else
    {
      calo_phys = new G4PVPlacement(sec_place, sec_logic,
                                    G4String(name.str()), cylinder_logic, false, sec,
                                    OverlapCheck());
    }
    calo_vol[calo_phys] = sec;

    assert(gdml_config);
    gdml_config->exclude_physical_vol(calo_phys);
  }
  _geom->set_nscint(_geom->get_nscint() * _geom->get_sector_map().size());

  if (active)
  {
    std::ostringstream geonode;
    if (superdetector != "NONE")
    {
      geonode << "CYLINDERGEOM_" << superdetector;
    }
    else
    {
      geonode << "CYLINDERGEOM_" << detector_type << "_" << layer;
    }
    PHG4CylinderGeomContainer *geo = findNode::getClass<PHG4CylinderGeomContainer>(topNode(), geonode.str());
    if (!geo)
    {
      geo = new PHG4CylinderGeomContainer();
      PHNodeIterator iter(topNode());
      PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(geo,
                                                                   geonode.str(), "PHObject");
      runNode->addNode(newNode);
    }
    // here in the detector class we have internal units, convert to cm
    // before putting into the geom object
    PHG4CylinderGeom *mygeom = clone_geom();
    geo->AddLayerGeom(layer, mygeom);
    //    geo->identify();
  }

  if (absorberactive)
  {
    std::ostringstream geonode;
    if (superdetector != "NONE")
    {
      geonode << "CYLINDERGEOM_ABSORBER_" << superdetector;
    }
    else
    {
      geonode << "CYLINDERGEOM_ABSORBER_" << detector_type << "_" << layer;
    }
    PHG4CylinderGeomContainer *geo = findNode::getClass<PHG4CylinderGeomContainer>(topNode(), geonode.str());
    if (!geo)
    {
      geo = new PHG4CylinderGeomContainer();
      PHNodeIterator iter(topNode());
      PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(geo,
                                                                   geonode.str(), "PHObject");
      runNode->addNode(newNode);
    }
    // here in the detector class we have internal units, convert to cm
    // before putting into the geom object
    PHG4CylinderGeom *mygeom = clone_geom();
    geo->AddLayerGeom(layer, mygeom);
    //    geo->identify();
  }

  if ((Verbosity() > 0))
  {
    std::cout << "PHG4SpacalDetector::Construct::" << GetName()
              << " - Completed. Print Geometry:" << std::endl;
    Print();
  }
}

std::pair<G4LogicalVolume *, G4Transform3D>
PHG4SpacalDetector::Construct_AzimuthalSeg()
{
  G4Tubs *sec_solid = new G4Tubs(G4String(GetName() + std::string("_sec")),
                                 _geom->get_radius() * cm, _geom->get_max_radius() * cm,
                                 _geom->get_length() * cm / 2.0, 0, twopi / _geom->get_azimuthal_n_sec());

  G4Material *cylinder_mat = GetDetectorMaterial(_geom->get_absorber_mat());
  assert(cylinder_mat);

  G4LogicalVolume *sec_logic = new G4LogicalVolume(sec_solid, cylinder_mat,
                                                   G4String(G4String(GetName() + std::string("_sec"))), nullptr, nullptr, nullptr);

  GetDisplayAction()->AddVolume(sec_logic, "AzimuthSegment");

  const double fiber_length = _geom->get_thickness() * cm - 2 * _geom->get_fiber_outer_r() * cm;
  G4LogicalVolume *fiber_logic = Construct_Fiber(fiber_length, std::string(""));

  int fiber_count = 0;
  //  double z_step = _geom->get_fiber_distance() * cm * sqrt(3) / 2.;
  double z_step = _geom->get_z_distance() * cm;
  G4double z = _geom->get_zmin() * cm - _geom->get_zpos() * cm + z_step;
  while (z < _geom->get_zmax() * cm - _geom->get_zpos() * cm - z_step)
  {
    const double rot = twopi / _geom->get_azimuthal_n_sec() * ((fiber_count % 2 == 0) ? 1. / 4. : 3. / 4.);

    G4Transform3D fiber_place(G4RotateZ3D(rot) * G4TranslateZ3D(z) * G4TranslateX3D(_geom->get_half_radius() * cm) * G4RotateY3D(halfpi));

    std::ostringstream name;
    name << GetName() << "_fiber_" << fiber_count;

    G4PVPlacement *fiber_physi = new G4PVPlacement(fiber_place, fiber_logic,
                                                   G4String(name.str()), sec_logic, false, fiber_count,
                                                   OverlapCheck());
    fiber_vol[fiber_physi] = fiber_count;
    assert(gdml_config);
    gdml_config->exclude_physical_vol(fiber_physi);

    z += z_step;
    fiber_count++;
  }
  _geom->set_nscint(fiber_count);

  if (Verbosity() > 0)
  {
    std::cout << "PHG4SpacalDetector::Construct_AzimuthalSeg::" << GetName()
              << " - constructed " << fiber_count << " fibers" << std::endl;
    std::cout << "\t"
              << "_geom->get_fiber_distance() = " << _geom->get_fiber_distance()
              << std::endl;
    std::cout << "\t"
              << "fiber_length = " << fiber_length / cm << std::endl;
    std::cout << "\t"
              << "z_step = " << z_step << std::endl;
    std::cout << "\t"
              << "_geom->get_azimuthal_bin() = " << _geom->get_azimuthal_n_sec()
              << std::endl;
    std::cout << "\t"
              << "_geom->get_azimuthal_distance() = "
              << _geom->get_azimuthal_distance() << std::endl;
    std::cout << "\t"
              << "_geom->is_virualize_fiber() = " << _geom->is_virualize_fiber()
              << std::endl;
  }

  return std::make_pair(sec_logic, G4Transform3D::Identity);
}

G4LogicalVolume *
PHG4SpacalDetector::Construct_Fiber(const G4double length, const std::string &id)
{
  G4Tubs *fiber_solid = new G4Tubs(G4String(GetName() + std::string("_fiber") + id),
                                   0, _geom->get_fiber_outer_r() * cm, length / 2.0, 0, twopi);

  G4Material *clading_mat = GetDetectorMaterial(_geom->get_fiber_clading_mat());
  assert(clading_mat);

  G4LogicalVolume *fiber_logic = new G4LogicalVolume(fiber_solid, clading_mat,
                                                     G4String(G4String(GetName() + std::string("_fiber") + id)), nullptr, nullptr,
                                                     nullptr);

  {
    GetDisplayAction()->AddVolume(fiber_logic, "Fiber");
  }

  G4Tubs *core_solid = new G4Tubs(
      G4String(GetName() + std::string("_fiber_core") + id), 0,
      _geom->get_fiber_core_diameter() * cm / 2, length / 2.0, 0, twopi);

  G4Material *core_mat = GetDetectorMaterial(_geom->get_fiber_core_mat());
  assert(core_mat);

  G4LogicalVolume *core_logic = new G4LogicalVolume(core_solid, core_mat,
                                                    G4String(G4String(GetName() + std::string("_fiber_core") + id)), nullptr, nullptr,
                                                    fiber_core_step_limits);

  {
    GetDisplayAction()->AddVolume(core_logic, "FiberCore");
  }

  const bool overlapcheck_fiber = OverlapCheck() and (Verbosity() >= 3);
  G4PVPlacement *core_physi = new G4PVPlacement(nullptr, G4ThreeVector(), core_logic,
                                                G4String(G4String(GetName() + std::string("_fiber_core") + id)), fiber_logic,
                                                false, 0, overlapcheck_fiber);
  fiber_core_vol[core_physi] = 0;

  return fiber_logic;
}

void PHG4SpacalDetector::Print(const std::string & /*what*/) const
{
  std::cout << "PHG4SpacalDetector::Print::" << GetName() << " - Print Geometry:" << std::endl;
  _geom->Print();

  return;
}
