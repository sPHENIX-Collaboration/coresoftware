#include "Prototype2RawTowerBuilder.h"
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/RawTowerGeomv1.h>
#include <calobase/RawTowerv1.h>

#include <g4detectors/PHG4ScintillatorSlat.h>
#include <g4detectors/PHG4ScintillatorSlatContainer.h>
#include <g4detectors/PHG4ScintillatorSlatDefs.h>
#include <g4detectors/PHG4PrototypeHcalDefs.h>

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

#include <iostream>
#include <map>
#include <stdexcept>

using namespace std;

Prototype2RawTowerBuilder::Prototype2RawTowerBuilder(const std::string &name)
  : SubsysReco(name)
  , PHParameterInterface(name)
  , m_Detector("NONE")
  , m_Emin(NAN)
  , m_CheckEnergyConservationFlag(0)
  , m_TowerEnergySrc(kLightYield)
  , m_NumCellToTower(5)

{
  InitializeParameters();
}

int Prototype2RawTowerBuilder::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "PAR"));
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!runNode || !dstNode || !parNode)
  {
    cout << PHWHERE << "Top Nodes missing, quitting after printing node tree" << endl;
    Fun4AllServer *se = Fun4AllServer::instance();
    se->Print("NODETREE");
    gSystem->Exit(1);
  }
  string paramnodename = "G4TOWERPARAM_" + m_Detector;
  string geonodename = "G4TOWERGEO_" + m_Detector;

  if (m_SimTowerNodePrefix.empty())
  {
    // no prefix, consistent with older convension
    m_TowerNodeName = "TOWER_" + m_Detector;
  }
  else
  {
    m_TowerNodeName = "TOWER_" + m_SimTowerNodePrefix + "_" + m_Detector;
  }
  RawTowerContainer *towers = findNode::getClass<RawTowerContainer>(topNode, m_TowerNodeName);
  if (!towers)
  {
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", m_Detector));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode(m_Detector);
      dstNode->addNode(DetNode);
    }
    towers = new RawTowerContainer(RawTowerDefs::convert_name_to_caloid(m_Detector));
    PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(towers, m_TowerNodeName, "PHObject");
    DetNode->addNode(towerNode);
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
  // save this to the parNode for use
  PHNodeIterator parIter(parNode);
  PHCompositeNode *ParDetNode = dynamic_cast<PHCompositeNode *>(parIter.findFirst("PHCompositeNode", m_Detector));
  if (!ParDetNode)
  {
    ParDetNode = new PHCompositeNode(m_Detector);
    parNode->addNode(ParDetNode);
  }
  PutOnParNode(ParDetNode, geonodename);

  m_TowerGeomNodeName = "TOWERGEOM_" + m_Detector;
  RawTowerGeomContainer *rawtowergeom = findNode::getClass<RawTowerGeomContainer>(topNode, m_TowerGeomNodeName);
  if (!rawtowergeom)
  {
    rawtowergeom = new RawTowerGeomContainer_Cylinderv1(RawTowerDefs::convert_name_to_caloid(m_Detector));
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(rawtowergeom, m_TowerGeomNodeName, "PHObject");
    runNode->addNode(newNode);
  }
  rawtowergeom->set_phibins(4);
  rawtowergeom->set_etabins(4);
  for (int irow = 0; irow < rawtowergeom->get_phibins(); irow++)
  {
    for (int icolumn = 0; icolumn < rawtowergeom->get_etabins(); icolumn++)
    {
      RawTowerGeomv1 *tg = new RawTowerGeomv1(RawTowerDefs::encode_towerid(RawTowerDefs::convert_name_to_caloid(m_Detector), icolumn, irow));
      tg->set_center_x(irow * 10 + icolumn);
      tg->set_center_y(irow * 10 + icolumn);
      tg->set_center_z(irow * 10 + icolumn);

      rawtowergeom->add_tower_geometry(tg);
    }
  }

  if (Verbosity() >= 1)
  {
    cout << "Prototype2RawTowerBuilder::InitRun :";
    if (m_TowerEnergySrc == kEnergyDeposition)
    {
      cout << "save Geant4 energy deposition as the weight of the cells"
           << endl;
    }
    else if (m_TowerEnergySrc == kLightYield)
    {
      cout << "save light yield as the weight of the cells" << endl;
    }
  }
  m_NumCellToTower = get_int_param(PHG4PrototypeHcalDefs::scipertwr);
  m_Emin = get_double_param("emin");
  return Fun4AllReturnCodes::EVENT_OK;
}

int Prototype2RawTowerBuilder::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 3)
  {
    std::cout << PHWHERE << "Process event entered" << std::endl;
  }

  // get cells
  string cellnodename = "G4CELL_" + m_Detector;
  PHG4ScintillatorSlatContainer *slats = findNode::getClass<PHG4ScintillatorSlatContainer>(topNode, cellnodename);
  if (!slats)
  {
    cout << PHWHERE << " " << cellnodename << " Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  RawTowerContainer *towers = findNode::getClass<RawTowerContainer>(topNode, m_TowerNodeName);
  if (!towers)
  {
    cout << PHWHERE << " " << m_TowerNodeName << " Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // loop over all slats in an event
  PHG4ScintillatorSlatContainer::ConstIterator cell_iter;
  PHG4ScintillatorSlatContainer::ConstRange cell_range = slats->getScintillatorSlats();
  for (cell_iter = cell_range.first; cell_iter != cell_range.second;
       ++cell_iter)
  {
    PHG4ScintillatorSlat *cell = cell_iter->second;

    if (Verbosity() > 2)
    {
      std::cout << PHWHERE << " print out the cell:" << std::endl;
      cell->identify();
    }
    short twrrow = get_tower_row(cell->get_row());
    // add the energy to the corresponding tower
    // towers are addressed column/row to make the mapping more intuitive
    RawTower *tower = towers->getTower(cell->get_column(), twrrow);
    if (!tower)
    {
      tower = new RawTowerv1();
      tower->set_energy(0);
      towers->AddTower(cell->get_column(), twrrow, tower);
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

    tower->add_ecell(cell->get_key(), cell_weight);

    tower->set_energy(tower->get_energy() + cell_weight);
  }
  double towerE = 0;
  if (m_CheckEnergyConservationFlag)
  {
    double cellE = slats->getTotalEdep();
    towerE = towers->getTotalEdep();
    if (fabs(cellE - towerE) / cellE > 1e-5)
    {
      cout << "towerE: " << towerE << ", cellE: " << cellE << ", delta: "
           << cellE - towerE << endl;
    }
  }
  if (Verbosity())
  {
    towerE = towers->getTotalEdep();
  }

  towers->compress(m_Emin);
  if (Verbosity())
  {
    cout << "Energy lost by dropping towers with less than " << m_Emin
         << " energy, lost energy: " << towerE - towers->getTotalEdep()
         << endl;
    towers->identify();
    RawTowerContainer::ConstRange begin_end = towers->getTowers();
    RawTowerContainer::ConstIterator iter;
    for (iter = begin_end.first; iter != begin_end.second; ++iter)
    {
      iter->second->identify();
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

short Prototype2RawTowerBuilder::get_tower_row(const short cellrow) const
{
  short twrrow = cellrow / m_NumCellToTower;
  return twrrow;
}

void Prototype2RawTowerBuilder::SetDefaultParameters()
{
  set_default_int_param(PHG4PrototypeHcalDefs::scipertwr, 0);
  set_default_double_param("emin", 1.e-6);
  return;
}

void Prototype2RawTowerBuilder::ReadParamsFromNodeTree(PHCompositeNode *topNode)
{
  PHParameters *pars = new PHParameters("temp");
  // we need the number of scintillator plates per tower
  string geonodename = "G4GEOPARAM_" + m_Detector;
  PdbParameterMapContainer *saveparams = findNode::getClass<PdbParameterMapContainer>(topNode,geonodename);
  if (! saveparams)
    {
      cout << "could not find " << geonodename << endl;
      Fun4AllServer *se = Fun4AllServer::instance();
      se->Print("NODETREE");
      return;
    }
  pars->FillFrom(saveparams,0);
  set_int_param(PHG4PrototypeHcalDefs::scipertwr,pars->get_int_param(PHG4PrototypeHcalDefs::scipertwr));
  delete pars;
  return;
}

void Prototype2RawTowerBuilder::Print(const std::string &what) const
{
  cout << "m_NumCellToTower: " << m_NumCellToTower << endl;
  return;
}
