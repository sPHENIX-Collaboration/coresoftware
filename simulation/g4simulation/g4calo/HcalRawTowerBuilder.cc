#include "HcalRawTowerBuilder.h"
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/RawTowerGeomv1.h>
#include <calobase/RawTowerv1.h>

#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CellDefs.h>
#include <g4detectors/PHG4HcalDefs.h>
#include <phparameter/PHParameters.h>

#include <g4main/PHG4Utils.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <pdbcalbase/PdbParameterMap.h>
#include <pdbcalbase/PdbParameterMapContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

#include <TSystem.h>

#include <cassert>
#include <iostream>
#include <map>
#include <stdexcept>

using namespace std;

HcalRawTowerBuilder::HcalRawTowerBuilder(const std::string &name)
  : SubsysReco(name)
  , PHParameterInterface(name)
  , m_Towers(nullptr)
  , m_RawTowerGeom(nullptr)
  , m_Detector("NONE")
  , m_Emin(NAN)
  , m_ChkEnergyConservationFlag(0)
  , m_TowerEnergySrc(enu_tower_energy_src::unknown)
  , m_NcellToTower(-1)
{
  InitializeParameters();
}

int HcalRawTowerBuilder::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode",
                                                                            "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    gSystem->Exit(1);
  }
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  string paramnodename = "TOWERPARAM_" + m_Detector;

  try
  {
    CreateNodes(topNode);
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    //exit(1);
  }

  // order first default,
  // then parameter from g4detector on node tree
  ReadParamsFromNodeTree(topNode);
  // then macro setting
  UpdateParametersWithMacro();
  PHNodeIterator runIter(runNode);
  PHCompositeNode *RunDetNode = dynamic_cast<PHCompositeNode *>(runIter.findFirst("PHCompositeNode", m_Detector));
  if (!RunDetNode)
  {
    RunDetNode = new PHCompositeNode(m_Detector);
    runNode->addNode(RunDetNode);
  }
  SaveToNodeTree(RunDetNode, paramnodename);
  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "PAR"));
  string geonodename = "TOWERGEO_" + m_Detector;

  PHNodeIterator parIter(parNode);
  PHCompositeNode *ParDetNode = dynamic_cast<PHCompositeNode *>(parIter.findFirst("PHCompositeNode", m_Detector));
  if (!ParDetNode)
  {
    ParDetNode = new PHCompositeNode(m_Detector);
    parNode->addNode(ParDetNode);
  }
  PutOnParNode(ParDetNode, geonodename);
  m_TowerEnergySrc = get_int_param("tower_energy_source");
  m_Emin = get_double_param("emin");
  m_NcellToTower = get_int_param("n_scinti_plates_per_tower");
  if (Verbosity() >= 1)
  {
    cout << "HcalRawTowerBuilder::InitRun :";
    if (m_TowerEnergySrc == kEnergyDeposition)
    {
      cout << "save Geant4 energy deposition in towers" << endl;
    }
    else if (m_TowerEnergySrc == kLightYield)
    {
      cout << "save light yield in towers" << endl;
    }
    else if (m_TowerEnergySrc == kIonizationEnergy)
    {
      cout << "save ionization energy in towers" << endl;
    }
    else
    {
      cout << "unknown energy source" << endl;
    }
  }
  m_TowerGeomNodeName = "TOWERGEOM_" + m_Detector;
  m_RawTowerGeom = findNode::getClass<RawTowerGeomContainer>(topNode,
                                                           m_TowerGeomNodeName);
  if (!m_RawTowerGeom)
  {
    m_RawTowerGeom = new RawTowerGeomContainer_Cylinderv1(RawTowerDefs::convert_name_to_caloid(m_Detector));
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(m_RawTowerGeom,
                                                                 m_TowerGeomNodeName, "PHObject");
    RunDetNode->addNode(newNode);
  }
  double innerrad = get_double_param(PHG4HcalDefs::innerrad);
  double thickness = get_double_param(PHG4HcalDefs::outerrad) - innerrad;
  m_RawTowerGeom->set_radius(innerrad);
  m_RawTowerGeom->set_thickness(thickness);
  m_RawTowerGeom->set_phibins(get_int_param(PHG4HcalDefs::n_towers));
  m_RawTowerGeom->set_etabins(get_int_param("etabins"));
  double geom_ref_radius = innerrad + thickness / 2.;
  double phistart = 0;
  for (int i = 0; i < get_int_param(PHG4HcalDefs::n_towers); i++)
  {
    double phiend = phistart + 2. * M_PI / get_int_param(PHG4HcalDefs::n_towers);
    pair<double, double> range = make_pair(phistart, phiend);
    phistart = phiend;
    m_RawTowerGeom->set_phibounds(i, range);
  }
  double etalowbound = -1.1;
  for (int i = 0; i < get_int_param("etabins"); i++)
  {
    double etahibound = etalowbound + 2.2 / get_int_param("etabins");
    pair<double, double> range = make_pair(etalowbound, etahibound);
    m_RawTowerGeom->set_etabounds(i, range);
    etalowbound = etahibound;
  }
  for (int iphi = 0; iphi < m_RawTowerGeom->get_phibins(); iphi++)
  {
    for (int ieta = 0; ieta < m_RawTowerGeom->get_etabins(); ieta++)
    {
      const RawTowerDefs::keytype key =
          RawTowerDefs::encode_towerid(RawTowerDefs::convert_name_to_caloid(m_Detector), ieta, iphi);

      const double x(geom_ref_radius * cos(m_RawTowerGeom->get_phicenter(iphi)));
      const double y(geom_ref_radius * sin(m_RawTowerGeom->get_phicenter(iphi)));
      const double z(geom_ref_radius / tan(PHG4Utils::get_theta(m_RawTowerGeom->get_etacenter(ieta))));

      RawTowerGeom *tg = m_RawTowerGeom->get_tower_geometry(key);
      if (tg)
      {
        if (Verbosity() > 0)
        {
          cout << "HcalRawTowerBuilder::InitRun - Tower geometry " << key << " already exists" << endl;
        }

        if (fabs(tg->get_center_x() - x) > 1e-4)
        {
          cout << "HcalRawTowerBuilder::InitRun - Fatal Error - duplicated Tower geometry " << key << " with existing x = " << tg->get_center_x() << " and expected x = " << x
               << endl;

          return Fun4AllReturnCodes::ABORTRUN;
        }
        if (fabs(tg->get_center_y() - y) > 1e-4)
        {
          cout << "HcalRawTowerBuilder::InitRun - Fatal Error - duplicated Tower geometry " << key << " with existing y = " << tg->get_center_y() << " and expected y = " << y
               << endl;
          return Fun4AllReturnCodes::ABORTRUN;
        }
        if (fabs(tg->get_center_z() - z) > 1e-4)
        {
          cout << "HcalRawTowerBuilder::InitRun - Fatal Error - duplicated Tower geometry " << key << " with existing z= " << tg->get_center_z() << " and expected z = " << z
               << endl;
          return Fun4AllReturnCodes::ABORTRUN;
        }
      }
      else
      {
        if (Verbosity() > 0)
        {
          cout << "HcalRawTowerBuilder::InitRun - building tower geometry " << key << "" << endl;
        }

        tg = new RawTowerGeomv1(key);

        tg->set_center_x(x);
        tg->set_center_y(y);
        tg->set_center_z(z);
        m_RawTowerGeom->add_tower_geometry(tg);
      }
    }
  }
  if (Verbosity() > 0)
  {
    m_RawTowerGeom->identify();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int HcalRawTowerBuilder::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 3)
  {
    std::cout << PHWHERE << "Process event entered" << std::endl;
  }

  // get cells
  std::string cellnodename = "G4CELL_" + m_Detector;
  PHG4CellContainer *slats = findNode::getClass<PHG4CellContainer>(topNode, cellnodename);
  if (!slats)
  {
    std::cerr << PHWHERE << " " << cellnodename
              << " Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // loop over all slats in an event
  PHG4CellContainer::ConstIterator cell_iter;
  PHG4CellContainer::ConstRange cell_range = slats->getCells();
  for (cell_iter = cell_range.first; cell_iter != cell_range.second;
       ++cell_iter)
  {
    PHG4Cell *cell = cell_iter->second;

    if (Verbosity() > 2)
    {
      std::cout << PHWHERE << " print out the cell:" << std::endl;
      cell->identify();
    }
    short twrrow = get_tower_row(PHG4CellDefs::ScintillatorSlatBinning::get_row(cell->get_cellid()));
    // add the energy to the corresponding tower
    // towers are addressed column/row to make the mapping more intuitive
    RawTower *tower = m_Towers->getTower(PHG4CellDefs::ScintillatorSlatBinning::get_column(cell->get_cellid()), twrrow);
    if (!tower)
    {
      tower = new RawTowerv1();
      tower->set_energy(0);
      m_Towers->AddTower(PHG4CellDefs::ScintillatorSlatBinning::get_column(cell->get_cellid()), twrrow, tower);
    }
    double cell_weight = 0;
    if (m_TowerEnergySrc == kEnergyDeposition)
    {
      cell_weight = cell->get_edep();
    }
    else if (m_TowerEnergySrc == kLightYield)
    {
      cell_weight = cell->get_light_yield();
    }
    else if (m_TowerEnergySrc == kIonizationEnergy)
    {
      cell_weight = cell->get_eion();
    }
    else
    {
      cout << Name() << ": unknown tower energy source "
           << m_TowerEnergySrc << endl;
      gSystem->Exit(1);
    }

    tower->add_ecell(cell->get_cellid(), cell_weight);

    PHG4Cell::ShowerEdepConstRange range = cell->get_g4showers();
    for (PHG4Cell::ShowerEdepConstIterator shower_iter = range.first;
         shower_iter != range.second;
         ++shower_iter)
    {
      tower->add_eshower(shower_iter->first, shower_iter->second);
    }
    tower->set_energy(tower->get_energy() + cell_weight);
  }
  double towerE = 0;
  if (m_ChkEnergyConservationFlag)
  {
    double cellE = slats->getTotalEdep();
    towerE = m_Towers->getTotalEdep();
    if (fabs(cellE - towerE) / cellE > 1e-5)
    {
      cout << "towerE: " << towerE << ", cellE: " << cellE << ", delta: "
           << cellE - towerE << endl;
    }
  }
  if (Verbosity())
  {
    towerE = m_Towers->getTotalEdep();
  }

  m_Towers->compress(m_Emin);
  if (Verbosity())
  {
    cout << "Energy lost by dropping towers with less than " << m_Emin
         << " energy, lost energy: " << towerE - m_Towers->getTotalEdep()
         << endl;
    m_Towers->identify();
    RawTowerContainer::ConstRange begin_end = m_Towers->getTowers();
    RawTowerContainer::ConstIterator iter;
    for (iter = begin_end.first; iter != begin_end.second; ++iter)
    {
      iter->second->identify();
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void HcalRawTowerBuilder::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst(
      "PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cerr << PHWHERE << "Run Node missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find Run node in HcalRawTowerBuilder::CreateNodes");
  }
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst(
      "PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find DST node in HcalRawTowerBuilder::CreateNodes");
  }

  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst(
      "PHCompositeNode", m_Detector));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode(m_Detector);
    dstNode->addNode(DetNode);
  }

  // Create the tower nodes on the tree
  if (m_SimTowerNodePrefix.empty())
  {
    // no prefix, consistent with older convension
    m_TowerNodeName = "TOWER_" + m_Detector;
  }
  else
  {
    m_TowerNodeName = "TOWER_" + m_SimTowerNodePrefix + "_" + m_Detector;
  }
  m_Towers = findNode::getClass<RawTowerContainer>(DetNode, m_TowerNodeName);
  if (!m_Towers)
  {
    m_Towers = new RawTowerContainer(RawTowerDefs::convert_name_to_caloid(m_Detector));

    PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(m_Towers,
                                                                   m_TowerNodeName, "PHObject");
    DetNode->addNode(towerNode);
  }
  return;
}

short HcalRawTowerBuilder::get_tower_row(const short cellrow) const
{
  short twrrow = cellrow / m_NcellToTower;
  return twrrow;
}

void HcalRawTowerBuilder::SetDefaultParameters()
{
  set_default_int_param(PHG4HcalDefs::scipertwr, 5);
  set_default_int_param("tower_energy_source", kLightYield);
  set_default_int_param(PHG4HcalDefs::n_towers, 64);
  set_default_int_param("etabins", 24);

  set_default_double_param("emin", 1.e-6);
  set_default_double_param(PHG4HcalDefs::outerrad, NAN);
  set_default_double_param(PHG4HcalDefs::innerrad, NAN);
}

void HcalRawTowerBuilder::ReadParamsFromNodeTree(PHCompositeNode *topNode)
{
  PHParameters *pars = new PHParameters("temp");
  // we need the number of scintillator plates per tower
  string geonodename = "G4GEOPARAM_" + m_Detector;
  PdbParameterMapContainer *saveparams = findNode::getClass<PdbParameterMapContainer>(topNode, geonodename);
  if (!saveparams)
  {
    cout << "could not find " << geonodename << endl;
    Fun4AllServer *se = Fun4AllServer::instance();
    se->Print("NODETREE");
    return;
  }
  pars->FillFrom(saveparams, 0);
  set_int_param(PHG4HcalDefs::scipertwr, pars->get_int_param(PHG4HcalDefs::scipertwr));
  set_int_param(PHG4HcalDefs::n_towers, pars->get_int_param(PHG4HcalDefs::n_towers));
  set_int_param("etabins", 2 * pars->get_int_param(PHG4HcalDefs::n_scinti_tiles));
  set_double_param(PHG4HcalDefs::innerrad, pars->get_double_param(PHG4HcalDefs::innerrad));
  set_double_param(PHG4HcalDefs::outerrad, pars->get_double_param(PHG4HcalDefs::outerrad));
  delete pars;
  return;
}
