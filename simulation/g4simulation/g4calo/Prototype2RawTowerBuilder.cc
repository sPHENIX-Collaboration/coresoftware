#include "Prototype2RawTowerBuilder.h"
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/RawTowerGeomv1.h>
#include <calobase/RawTowerv1.h>
#include <g4detectors/PHG4ScintillatorSlat.h>
#include <g4detectors/PHG4ScintillatorSlatContainer.h>
#include <g4detectors/PHG4ScintillatorSlatDefs.h>

#include <g4main/PHG4Utils.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>

#include <TSystem.h>

#include <iostream>
#include <stdexcept>
#include <map>

using namespace std;

Prototype2RawTowerBuilder::Prototype2RawTowerBuilder(const std::string& name) :
    SubsysReco(name), 
    PHG4ParameterInterface(name),
    rawtowergeom(nullptr),
    m_Detector("NONE"), 
    emin(NAN),
    chkenergyconservation(0), 
    _tower_energy_src(kLightYield),
    ncell_to_tower(5)

{
  InitializeParameters();
}

int
Prototype2RawTowerBuilder::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
      exit(1);
    }

  try
    {
      CreateNodes(topNode);
    }
  catch (std::exception& e)
    {
      std::cout << e.what() << std::endl;
      //exit(1);
    }

  if (verbosity >= 1)
    {
      cout << "Prototype2RawTowerBuilder::InitRun :";
      if (_tower_energy_src == kEnergyDeposition)
        cout << "save Geant4 energy deposition as the weight of the cells"
            << endl;
      else if (_tower_energy_src == kLightYield)
        cout << "save light yield as the weight of the cells" << endl;
    }
  emin = get_double_param("emin");
  return Fun4AllReturnCodes::EVENT_OK;
}

int
Prototype2RawTowerBuilder::process_event(PHCompositeNode *topNode)
{
  if (verbosity>3)
    {
      std::cout << PHWHERE << "Process event entered" << std::endl;
    }

  // get cells
  string cellnodename = "G4CELL_" + m_Detector;
  PHG4ScintillatorSlatContainer* slats = findNode::getClass<PHG4ScintillatorSlatContainer>(topNode, cellnodename);
  if (!slats)
    {
      cout << PHWHERE << " " << cellnodename << " Node missing, doing nothing." << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  RawTowerContainer *towers = findNode::getClass<RawTowerContainer>(topNode,TowerNodeName);
  if (! towers)
  {
      cout << PHWHERE << " " << TowerNodeName << " Node missing, doing nothing." << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  // loop over all slats in an event
  PHG4ScintillatorSlatContainer::ConstIterator cell_iter;
  PHG4ScintillatorSlatContainer::ConstRange cell_range = slats->getScintillatorSlats();
  for (cell_iter = cell_range.first; cell_iter != cell_range.second;
      ++cell_iter)
    {
      PHG4ScintillatorSlat *cell = cell_iter->second;

      if (verbosity > 2)
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
      if (_tower_energy_src == kEnergyDeposition)
	{
	  cell_weight = cell->get_edep();
	}
      else if (_tower_energy_src == kLightYield)
	{
	  cell_weight = cell->get_light_yield();
	}
      else if (_tower_energy_src == kIonizationEnergy)
	{
	  cell_weight = cell->get_eion();
	}

      tower->add_ecell(cell->get_key(), cell_weight);

      tower->set_energy(tower->get_energy() + cell_weight);

    }
  double towerE = 0;
  if (chkenergyconservation)
    {
      double cellE = slats->getTotalEdep();
      towerE = towers->getTotalEdep();
      if (fabs(cellE - towerE) / cellE > 1e-5)
        {
          cout << "towerE: " << towerE << ", cellE: " << cellE << ", delta: "
              << cellE - towerE << endl;
        }
    }
  if (verbosity)
    {
      towerE = towers->getTotalEdep();
    }

  towers->compress(emin);
  if (verbosity)
    {
      cout << "Energy lost by dropping towers with less than " << emin
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

int
Prototype2RawTowerBuilder::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void
Prototype2RawTowerBuilder::CreateNodes(PHCompositeNode *topNode)
{

  PHNodeIterator iter(topNode);
  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "PAR"));
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
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

  if (_sim_tower_node_prefix.empty())
    {
      // no prefix, consistent with older convension
      TowerNodeName = "TOWER_" + m_Detector;
    }
  else
    {
      TowerNodeName = "TOWER_" + _sim_tower_node_prefix + "_" + m_Detector;
    }
  RawTowerContainer *towers = findNode::getClass<RawTowerContainer>(topNode,TowerNodeName);
  if (! towers)
  {
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", m_Detector));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode(m_Detector);
      dstNode->addNode(DetNode);
    }
 towers = new RawTowerContainer(RawTowerDefs::convert_name_to_caloid(m_Detector));
  PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(towers,TowerNodeName, "PHObject");
  DetNode->addNode(towerNode);
  }

  TowerGeomNodeName = "TOWERGEOM_" + m_Detector;
  rawtowergeom = findNode::getClass<RawTowerGeomContainer>(topNode, TowerGeomNodeName);
  if (!rawtowergeom)
    {

      rawtowergeom = new RawTowerGeomContainer_Cylinderv1(RawTowerDefs::convert_name_to_caloid(m_Detector));
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(rawtowergeom, TowerGeomNodeName.c_str(), "PHObject");
      runNode->addNode(newNode);
    }
  rawtowergeom->set_phibins(4);
  rawtowergeom->set_etabins(4);
  for (int irow = 0; irow < rawtowergeom->get_phibins(); irow++)
    {
      for (int icolumn=0; icolumn<rawtowergeom->get_etabins(); icolumn++)
	{
	  RawTowerGeomv1 * tg = new RawTowerGeomv1(RawTowerDefs::encode_towerid(RawTowerDefs::convert_name_to_caloid(m_Detector), icolumn, irow));
            tg->set_center_x(irow*10+icolumn);
            tg->set_center_y(irow*10+icolumn);
            tg->set_center_z(irow*10+icolumn);

            rawtowergeom->add_tower_geometry(tg);
	}
    }


  return;
}

short
Prototype2RawTowerBuilder::get_tower_row(const short cellrow) const
{
  short twrrow = cellrow/ncell_to_tower;
  return twrrow;
}

void
Prototype2RawTowerBuilder::SetDefaultParameters()
{
set_default_double_param("emin", 1.e-6);
return;
}
