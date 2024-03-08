/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.7 $$
 * \date $$Date: 2015/02/10 15:39:26 $$
 */
#include "PHG4SpacalDetector.h"
#include "PHG4CylinderGeom_Spacalv3.h"

#include "PHG4CylinderCellGeom_Spacalv1.h"
#include "PHG4CylinderGeom_Spacalv1.h"  // for PHG4CylinderGeom_Spaca...

#include "PHG4CellDefs.h"
#include "PHG4CylinderCellGeom.h"
#include "PHG4CylinderCellGeomContainer.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4SpacalDisplayAction.h"

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>
#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/recoConsts.h>

#include <g4gdml/PHG4GDMLConfig.hh>
#include <g4gdml/PHG4GDMLUtility.hh>

#include <calobase/RawTowerDefs.h>           // for convert_name_...
#include <calobase/RawTowerGeom.h>           // for RawTowerGeom
#include <calobase/RawTowerGeomContainer.h>  // for RawTowerGeomC...
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/RawTowerGeomv1.h>

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
// note this inactive thickness is ~1.5% of a radiation length
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
              << "fiber_length = " << fiber_length
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
// This is dulplicated code, we can get rid of it when we have the code to make towergeom for real data reco.
void PHG4SpacalDetector::AddTowerGeometryNode()
{
  PHNodeIterator iter(topNode());
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cout << PHWHERE << "Run Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find Run node in PHG4SpacalDetector::AddGeometryNode");
  }

  PHNodeIterator runIter(runNode);
  PHCompositeNode *RunDetNode = dynamic_cast<PHCompositeNode *>(runIter.findFirst("PHCompositeNode", superdetector));
  if (!RunDetNode)
  {
    RunDetNode = new PHCompositeNode(superdetector);
    runNode->addNode(RunDetNode);
  }

  const RawTowerDefs::CalorimeterId caloid = RawTowerDefs::convert_name_to_caloid(superdetector);

  // get the cell geometry and build up the tower geometry object
  std::string geonodename = "CYLINDERCELLGEOM_" + superdetector;
  PHG4CylinderCellGeomContainer *cellgeos = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode(), geonodename);
  if (!cellgeos)
  {
    std::cout << PHWHERE << " " << geonodename
              << " Node missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find " + geonodename + " node in PHG4SpacalDetector::AddGeometryNode");
  }
  m_TowerGeomNodeName = "TOWERGEOM_" + superdetector;
  m_RawTowerGeomContainer = findNode::getClass<RawTowerGeomContainer>(topNode(), m_TowerGeomNodeName);
  if (!m_RawTowerGeomContainer)
  {
    m_RawTowerGeomContainer = new RawTowerGeomContainer_Cylinderv1(caloid);
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(m_RawTowerGeomContainer, m_TowerGeomNodeName, "PHObject");
    RunDetNode->addNode(newNode);
  }
  // fill the number of layers in the calorimeter
  m_NumLayers = cellgeos->get_NLayers();

  // Create the tower nodes on the tree
  PHG4CylinderCellGeomContainer::ConstIterator miter;
  PHG4CylinderCellGeomContainer::ConstRange begin_end =
      cellgeos->get_begin_end();
  int ifirst = 1;
  int first_layer = -1;
  PHG4CylinderCellGeom *first_cellgeo = nullptr;
  double inner_radius = 0;
  double thickness = 0;
  for (miter = begin_end.first; miter != begin_end.second; ++miter)
  {
    PHG4CylinderCellGeom *cellgeo = miter->second;

    if (Verbosity())
    {
      cellgeo->identify();
    }
    thickness += cellgeo->get_thickness();
    if (ifirst)
    {
      first_cellgeo = miter->second;
      m_CellBinning = cellgeo->get_binning();
      m_NumPhiBins = cellgeo->get_phibins();
      m_PhiMin = cellgeo->get_phimin();
      m_PhiStep = cellgeo->get_phistep();
      if (m_CellBinning == PHG4CellDefs::etaphibinning || m_CellBinning == PHG4CellDefs::etaslatbinning)
      {
        m_NumEtaBins = cellgeo->get_etabins();
        m_EtaMin = cellgeo->get_etamin();
        m_EtaStep = cellgeo->get_etastep();
      }
      else if (m_CellBinning == PHG4CellDefs::sizebinning)
      {
        m_NumEtaBins = cellgeo->get_zbins();  // bin eta in the same number of z bins
      }
      else if (m_CellBinning == PHG4CellDefs::spacalbinning)
      {
        // use eta definiton for each row of towers
        m_NumEtaBins = cellgeo->get_etabins();
      }
      else
      {
        std::cout << "PHG4SpacalDetector::AddGeometryNode::" << GetName()
                  << " - Fatal Error - unsupported cell binning method "
                  << m_CellBinning << std::endl;
      }
      inner_radius = cellgeo->get_radius();
      first_layer = cellgeo->get_layer();
      ifirst = 0;
    }
    else
    {
      if (m_CellBinning != cellgeo->get_binning())
      {
        std::cout << "inconsistent binning method from " << m_CellBinning
                  << " layer " << cellgeo->get_layer() << ": "
                  << cellgeo->get_binning() << std::endl;
        exit(1);
      }
      if (inner_radius > cellgeo->get_radius())
      {
        std::cout << "radius of layer " << cellgeo->get_layer() << " is "
                  << cellgeo->get_radius() << " which smaller than radius "
                  << inner_radius << " of first layer in list " << first_layer
                  << std::endl;
      }
      if (m_NumPhiBins != cellgeo->get_phibins())
      {
        std::cout << "mixing number of phibins, fisrt layer: " << m_NumPhiBins
                  << " layer " << cellgeo->get_layer() << ": "
                  << cellgeo->get_phibins() << std::endl;
        exit(1);
      }
      if (m_PhiMin != cellgeo->get_phimin())
      {
        std::cout << "mixing number of phimin, fisrt layer: " << m_PhiMin
                  << " layer " << cellgeo->get_layer() << ": "
                  << cellgeo->get_phimin() << std::endl;
        exit(1);
      }
      if (m_PhiStep != cellgeo->get_phistep())
      {
        std::cout << "mixing phisteps first layer: " << m_PhiStep << " layer "
                  << cellgeo->get_layer() << ": " << cellgeo->get_phistep()
                  << " diff: " << m_PhiStep - cellgeo->get_phistep() << std::endl;
        exit(1);
      }
      if (m_CellBinning == PHG4CellDefs::etaphibinning || m_CellBinning == PHG4CellDefs::etaslatbinning)
      {
        if (m_NumEtaBins != cellgeo->get_etabins())
        {
          std::cout << "mixing number of EtaBins , first layer: "
                    << m_NumEtaBins << " layer " << cellgeo->get_layer() << ": "
                    << cellgeo->get_etabins() << std::endl;
          exit(1);
        }
        if (fabs(m_EtaMin - cellgeo->get_etamin()) > 1e-9)
        {
          std::cout << "mixing etamin, fisrt layer: " << m_EtaMin << " layer "
                    << cellgeo->get_layer() << ": " << cellgeo->get_etamin()
                    << " diff: " << m_EtaMin - cellgeo->get_etamin() << std::endl;
          exit(1);
        }
        if (fabs(m_EtaStep - cellgeo->get_etastep()) > 1e-9)
        {
          std::cout << "mixing eta steps first layer: " << m_EtaStep
                    << " layer " << cellgeo->get_layer() << ": "
                    << cellgeo->get_etastep() << " diff: "
                    << m_EtaStep - cellgeo->get_etastep() << std::endl;
          exit(1);
        }
      }

      else if (m_CellBinning == PHG4CellDefs::sizebinning)
      {
        if (m_NumEtaBins != cellgeo->get_zbins())
        {
          std::cout << "mixing number of z bins , first layer: " << m_NumEtaBins
                    << " layer " << cellgeo->get_layer() << ": "
                    << cellgeo->get_zbins() << std::endl;
          exit(1);
        }
      }
    }
  }
  m_RawTowerGeomContainer->set_radius(inner_radius);
  m_RawTowerGeomContainer->set_thickness(thickness);
  m_RawTowerGeomContainer->set_phibins(m_NumPhiBins);
  //  m_RawTowerGeomContainer->set_phistep(m_PhiStep);
  //  m_RawTowerGeomContainer->set_phimin(m_PhiMin);
  m_RawTowerGeomContainer->set_etabins(m_NumEtaBins);

  if (!first_cellgeo)
  {
    std::cout << "PHG4SpacalDetector::AddGeometryNode - ERROR - can not find first layer of cells "
              << std::endl;

    exit(1);
  }

  for (int ibin = 0; ibin < first_cellgeo->get_phibins(); ibin++)
  {
    const std::pair<double, double> range = first_cellgeo->get_phibounds(ibin);

    m_RawTowerGeomContainer->set_phibounds(ibin, range);
  }

  if (m_CellBinning == PHG4CellDefs::etaphibinning || m_CellBinning == PHG4CellDefs::etaslatbinning || m_CellBinning == PHG4CellDefs::spacalbinning)
  {
    const double r = inner_radius;

    for (int ibin = 0; ibin < first_cellgeo->get_etabins(); ibin++)
    {
      const std::pair<double, double> range = first_cellgeo->get_etabounds(ibin);

      m_RawTowerGeomContainer->set_etabounds(ibin, range);
    }

    // setup location of all towers
    for (int iphi = 0; iphi < m_RawTowerGeomContainer->get_phibins(); iphi++)
    {
      for (int ieta = 0; ieta < m_RawTowerGeomContainer->get_etabins(); ieta++)
      {
        const RawTowerDefs::keytype key =
            RawTowerDefs::encode_towerid(caloid, ieta, iphi);

        const double x(r * cos(m_RawTowerGeomContainer->get_phicenter(iphi)));
        const double y(r * sin(m_RawTowerGeomContainer->get_phicenter(iphi)));
        const double z(r / tan(PHG4Utils::get_theta(m_RawTowerGeomContainer->get_etacenter(ieta))));

        RawTowerGeom *tg = m_RawTowerGeomContainer->get_tower_geometry(key);
        if (tg)
        {
          if (Verbosity() > 0)
          {
            std::cout << "PHG4SpacalDetector::AddGeometryNode - Tower geometry " << key << " already exists" << std::endl;
          }

          if (fabs(tg->get_center_x() - x) > 1e-4)
          {
            std::cout << "PHG4SpacalDetector::AddGeometryNode - Fatal Error - duplicated Tower geometry " << key << " with existing x = " << tg->get_center_x() << " and expected x = " << x
                      << std::endl;

            exit(1);
          }
          if (fabs(tg->get_center_y() - y) > 1e-4)
          {
            std::cout << "PHG4SpacalDetector::AddGeometryNode - Fatal Error - duplicated Tower geometry " << key << " with existing y = " << tg->get_center_y() << " and expected y = " << y
                      << std::endl;
            exit(1);
          }
          if (fabs(tg->get_center_z() - z) > 1e-4)
          {
            std::cout << "PHG4SpacalDetector::AddGeometryNode - Fatal Error - duplicated Tower geometry " << key << " with existing z= " << tg->get_center_z() << " and expected z = " << z
                      << std::endl;
            exit(1);
          }
        }
        else
        {
          if (Verbosity() > 0)
          {
            std::cout << "PHG4SpacalDetector::AddGeometryNode - building tower geometry " << key << "" << std::endl;
          }

          tg = new RawTowerGeomv1(key);

          tg->set_center_x(x);
          tg->set_center_y(y);
          tg->set_center_z(z);
          m_RawTowerGeomContainer->add_tower_geometry(tg);
        }
      }
    }
  }
  else if (m_CellBinning == PHG4CellDefs::sizebinning)
  {
    const double r = inner_radius;

    for (int ibin = 0; ibin < first_cellgeo->get_zbins(); ibin++)
    {
      const std::pair<double, double> z_range = first_cellgeo->get_zbounds(ibin);
      //          const double r = first_cellgeo->get_radius();
      const double eta1 = -log(tan(atan2(r, z_range.first) / 2.));
      const double eta2 = -log(tan(atan2(r, z_range.second) / 2.));
      m_RawTowerGeomContainer->set_etabounds(ibin, std::make_pair(eta1, eta2));
    }

    // setup location of all towers
    for (int iphi = 0; iphi < m_RawTowerGeomContainer->get_phibins(); iphi++)
    {
      for (int ibin = 0; ibin < first_cellgeo->get_zbins(); ibin++)
      {
        const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(caloid, ibin, iphi);

        const double x(r * cos(m_RawTowerGeomContainer->get_phicenter(iphi)));
        const double y(r * sin(m_RawTowerGeomContainer->get_phicenter(iphi)));
        const double z(first_cellgeo->get_zcenter(ibin));

        RawTowerGeom *tg = m_RawTowerGeomContainer->get_tower_geometry(key);
        if (tg)
        {
          if (Verbosity() > 0)
          {
            std::cout << "PHG4SpacalDetector::AddGeometryNode - Tower geometry " << key << " already exists" << std::endl;
          }

          if (fabs(tg->get_center_x() - x) > 1e-4)
          {
            std::cout << "PHG4SpacalDetector::AddGeometryNode - Fatal Error - duplicated Tower geometry " << key << " with existing x = " << tg->get_center_x() << " and expected x = " << x
                      << std::endl;

            exit(1);
          }
          if (fabs(tg->get_center_y() - y) > 1e-4)
          {
            std::cout << "PHG4SpacalDetector::AddGeometryNode - Fatal Error - duplicated Tower geometry " << key << " with existing y = " << tg->get_center_y() << " and expected y = " << y
                      << std::endl;
            exit(1);
          }
          if (fabs(tg->get_center_z() - z) > 1e-4)
          {
            std::cout << "PHG4SpacalDetector::AddGeometryNode - Fatal Error - duplicated Tower geometry " << key << " with existing z= " << tg->get_center_z() << " and expected z = " << z
                      << std::endl;
            exit(1);
          }
        }
        else
        {
          if (Verbosity() > 0)
          {
            std::cout << "PHG4SpacalDetector::AddGeometryNode - building tower geometry " << key << "" << std::endl;
          }

          tg = new RawTowerGeomv1(key);

          tg->set_center_x(x);
          tg->set_center_y(y);
          tg->set_center_z(z);
          m_RawTowerGeomContainer->add_tower_geometry(tg);
        }
      }
    }
  }
  else
  {
    std::cout << "PHG4SpacalDetector::AddGeometryNode - ERROR - unsupported cell geometry "
              << m_CellBinning << std::endl;
    exit(1);
  }

  if (Verbosity() >= 1)
  {
    m_RawTowerGeomContainer->identify();
  }

  return;
}

void PHG4SpacalDetector::AddCellGeometryNode()
{
  PHNodeIterator iter(topNode());
  std::string detector = superdetector;
  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    exit(1);
  }

  std::string geonodename = "CYLINDERGEOM_" + detector;
  PHG4CylinderGeomContainer *geo = findNode::getClass<PHG4CylinderGeomContainer>(topNode(), geonodename);
  if (!geo)
  {
    std::cout << "PHG4SpacalDetector::AddCellGeometryNode - Could not locate geometry node "
              << geonodename << std::endl;
    topNode()->print();
    exit(1);
  }
  if (Verbosity() > 0)
  {
    std::cout << "PHG4SpacalDetector::AddCellGeometryNode- incoming geometry:"
              << std::endl;
    geo->identify();
  }
  std::string seggeonodename = "CYLINDERCELLGEOM_" + detector;
  PHG4CylinderCellGeomContainer *seggeo = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode(), seggeonodename);
  if (!seggeo)
  {
    seggeo = new PHG4CylinderCellGeomContainer();
    PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(seggeo, seggeonodename, "PHObject");
    runNode->addNode(newNode);
  }

  const PHG4CylinderGeom *layergeom_raw = geo->GetFirstLayerGeom();
  assert(layergeom_raw);

  // a special implimentation of PHG4CylinderGeom is required here.
  const PHG4CylinderGeom_Spacalv3 *layergeom =
      dynamic_cast<const PHG4CylinderGeom_Spacalv3 *>(layergeom_raw);

  if (!layergeom)
  {
    std::cout << "PHG4SpacalDetector::AddCellGeometryNode- Fatal Error -"
              << " require to work with a version of SPACAL geometry of PHG4CylinderGeom_Spacalv3 or higher. "
              << "However the incoming geometry has version "
              << layergeom_raw->ClassName() << std::endl;
    exit(1);
  }
  if (Verbosity() > 1)
  {
    layergeom->identify();
  }

  layergeom->subtower_consistency_check();

  //      int layer = layergeom->get_layer();

  // create geo object and fill with variables common to all binning methods
  PHG4CylinderCellGeom_Spacalv1 *layerseggeo = new PHG4CylinderCellGeom_Spacalv1();
  layerseggeo->set_layer(layergeom->get_layer());
  layerseggeo->set_radius(layergeom->get_radius());
  layerseggeo->set_thickness(layergeom->get_thickness());
  layerseggeo->set_binning(PHG4CellDefs::spacalbinning);

  // construct a map to convert tower_ID into the older eta bins.

  const PHG4CylinderGeom_Spacalv3::tower_map_t &tower_map = layergeom->get_sector_tower_map();
  const PHG4CylinderGeom_Spacalv3::sector_map_t &sector_map = layergeom->get_sector_map();
  const int nphibin = layergeom->get_azimuthal_n_sec()       // sector
                      * layergeom->get_max_phi_bin_in_sec()  // blocks per sector
                      * layergeom->get_n_subtower_phi();     // subtower per block
  const double deltaphi = 2. * M_PI / nphibin;

  using map_z_tower_z_ID_t = std::map<double, int>;
  map_z_tower_z_ID_t map_z_tower_z_ID;
  double phi_min = NAN;

  for (const auto &tower_pair : tower_map)
  {
    const int &tower_ID = tower_pair.first;
    const PHG4CylinderGeom_Spacalv3::geom_tower &tower = tower_pair.second;

    // inspect index in sector 0
    std::pair<int, int> tower_z_phi_ID = layergeom->get_tower_z_phi_ID(tower_ID, 0);

    const int &tower_ID_z = tower_z_phi_ID.first;
    const int &tower_ID_phi = tower_z_phi_ID.second;

    if (tower_ID_phi == 0)
    {
      // assign phi min according phi bin 0
      phi_min = M_PI_2 - deltaphi * (layergeom->get_max_phi_bin_in_sec() * layergeom->get_n_subtower_phi() / 2)  // shift of first tower in sector
                + sector_map.begin()->second;
    }

    if (tower_ID_phi == layergeom->get_max_phi_bin_in_sec() / 2)
    {
      // assign eta min according phi bin 0
      map_z_tower_z_ID[tower.centralZ] = tower_ID_z;
    }
    // ...
  }  //       const auto &tower_pair: tower_map

  assert(!std::isnan(phi_min));
  layerseggeo->set_phimin(phi_min);
  layerseggeo->set_phistep(deltaphi);
  layerseggeo->set_phibins(nphibin);

  PHG4CylinderCellGeom_Spacalv1::tower_z_ID_eta_bin_map_t tower_z_ID_eta_bin_map;
  int eta_bin = 0;
  for (auto &z_tower_z_ID : map_z_tower_z_ID)
  {
    tower_z_ID_eta_bin_map[z_tower_z_ID.second] = eta_bin;
    eta_bin++;
  }
  layerseggeo->set_tower_z_ID_eta_bin_map(tower_z_ID_eta_bin_map);
  layerseggeo->set_etabins(eta_bin * layergeom->get_n_subtower_eta());
  layerseggeo->set_etamin(NAN);
  layerseggeo->set_etastep(NAN);

  // build eta bin maps
  for (const auto &tower_pair : tower_map)
  {
    const int &tower_ID = tower_pair.first;
    const PHG4CylinderGeom_Spacalv3::geom_tower &tower = tower_pair.second;

    // inspect index in sector 0
    std::pair<int, int> tower_z_phi_ID = layergeom->get_tower_z_phi_ID(tower_ID, 0);
    const int &tower_ID_z = tower_z_phi_ID.first;
    const int &tower_ID_phi = tower_z_phi_ID.second;
    const int &etabin = tower_z_ID_eta_bin_map[tower_ID_z];

    if (tower_ID_phi == layergeom->get_max_phi_bin_in_sec() / 2)
    {
      // half z-range
      const double dz = fabs(0.5 * (tower.pDy1 + tower.pDy2) / sin(tower.pRotationAngleX));
      const double tower_radial = layergeom->get_tower_radial_position(tower);

      auto z_to_eta = [&tower_radial](const double &z)
      { return -log(tan(0.5 * atan2(tower_radial, z))); };

      const double eta_central = z_to_eta(tower.centralZ);
      // half eta-range
      const double deta = (z_to_eta(tower.centralZ + dz) - z_to_eta(tower.centralZ - dz)) / 2;
      assert(deta > 0);

      for (int sub_tower_ID_y = 0; sub_tower_ID_y < tower.NSubtowerY;
           ++sub_tower_ID_y)
      {
        assert(tower.NSubtowerY <= layergeom->get_n_subtower_eta());
        // do not overlap to the next bin.
        const int sub_tower_etabin = etabin * layergeom->get_n_subtower_eta() + sub_tower_ID_y;

        const std::pair<double, double> etabounds(eta_central - deta + sub_tower_ID_y * 2 * deta / tower.NSubtowerY,
                                                  eta_central - deta + (sub_tower_ID_y + 1) * 2 * deta / tower.NSubtowerY);

        const std::pair<double, double> zbounds(tower.centralZ - dz + sub_tower_ID_y * 2 * dz / tower.NSubtowerY,
                                                tower.centralZ - dz + (sub_tower_ID_y + 1) * 2 * dz / tower.NSubtowerY);

        layerseggeo->set_etabounds(sub_tower_etabin, etabounds);
        layerseggeo->set_zbounds(sub_tower_etabin, zbounds);
      }
    }
    // ...
  }  //       const auto &tower_pair: tower_map

  // add geo object filled by different binning methods
  seggeo->AddLayerCellGeom(layerseggeo);
  // save this to the run wise tree to store on DST
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  // PHCompositeNode *parNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "PAR"));
  PHNodeIterator runIter(runNode);
  PHCompositeNode *RunDetNode = dynamic_cast<PHCompositeNode *>(runIter.findFirst("PHCompositeNode", detector));
  if (!RunDetNode)
  {
    RunDetNode = new PHCompositeNode(detector);
    runNode->addNode(RunDetNode);
  }
  /*
  std::string paramnodename = "G4CELLPARAM_" + detector;
  SaveToNodeTree(RunDetNode, paramnodename);
  save this to the parNode for use
  PHNodeIterator parIter(parNode);
   PHCompositeNode *ParDetNode = dynamic_cast<PHCompositeNode *>(parIter.findFirst("PHCompositeNode", detector));
   if (!ParDetNode)
   {
     ParDetNode = new PHCompositeNode(detector);
     parNode->addNode(ParDetNode);
   }
   std::string cellgeonodename = "G4CELLGEO_" + detector;
   PutOnParNode(ParDetNode, cellgeonodename);
    */
  return;
}