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

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>

#include <iostream>
#include <stdexcept>
#include <map>

using namespace std;

Prototype2RawTowerBuilder::Prototype2RawTowerBuilder(const std::string& name) :
    SubsysReco(name), 
    _towers(NULL),
    rawtowergeom(NULL),
    detector("NONE"), 
    emin(1e-6),
    chkenergyconservation(0), 
    _tower_energy_src(kLightYield),
    ncell_to_tower(5),
    _timer(PHTimeServer::get()->insert_new(name))
{
}

int
Prototype2RawTowerBuilder::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode",
      "DST"));
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
  std::string cellnodename = "G4CELL_" + detector;
  PHG4ScintillatorSlatContainer* slats = findNode::getClass<PHG4ScintillatorSlatContainer>(topNode, cellnodename.c_str());
  if (!slats)
    {
      std::cerr << PHWHERE << " " << cellnodename
          << " Node missing, doing nothing." << std::endl;
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
      RawTower *tower = _towers->getTower(cell->get_column(), twrrow);
      if (!tower)
        {
          tower = new RawTowerv1();
          tower->set_energy(0);
          _towers->AddTower(cell->get_column(), twrrow, tower);
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
      towerE = _towers->getTotalEdep();
      if (fabs(cellE - towerE) / cellE > 1e-5)
        {
          cout << "towerE: " << towerE << ", cellE: " << cellE << ", delta: "
              << cellE - towerE << endl;
        }
    }
  if (verbosity)
    {
      towerE = _towers->getTotalEdep();
    }

  _towers->compress(emin);
  if (verbosity)
    {
      cout << "Energy lost by dropping towers with less than " << emin
          << " energy, lost energy: " << towerE - _towers->getTotalEdep()
          << endl;
      _towers->identify();
      RawTowerContainer::ConstRange begin_end = _towers->getTowers();
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
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst(
      "PHCompositeNode", "RUN"));
  if (!runNode)
    {
      std::cerr << PHWHERE << "Run Node missing, doing nothing." << std::endl;
      throw std::runtime_error(
          "Failed to find Run node in Prototype2RawTowerBuilder::CreateNodes");
    }
  TowerGeomNodeName = "TOWERGEOM_" + detector;
  rawtowergeom = findNode::getClass<RawTowerGeomContainer>(topNode,
      TowerGeomNodeName.c_str());
  if (!rawtowergeom)
    {

      rawtowergeom = new RawTowerGeomContainer_Cylinderv1(RawTowerDefs::convert_name_to_caloid(detector));
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(rawtowergeom,
          TowerGeomNodeName.c_str(), "PHObject");
      runNode->addNode(newNode);
    }
  rawtowergeom->set_phibins(4);
  rawtowergeom->set_etabins(4);
  for (int irow = 0; irow < rawtowergeom->get_phibins(); irow++)
    {
      for (int icolumn=0; icolumn<rawtowergeom->get_etabins(); icolumn++)
	{
	  RawTowerGeomv1 * tg = new RawTowerGeomv1(RawTowerDefs::encode_towerid(RawTowerDefs::convert_name_to_caloid(detector), icolumn, irow));
            tg->set_center_x(irow*10+icolumn);
            tg->set_center_y(irow*10+icolumn);
            tg->set_center_z(irow*10+icolumn);

            rawtowergeom->add_tower_geometry(tg);
	}
    }

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst(
      "PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cerr << PHWHERE << "DST Node missing, doing nothing." << std::endl;
      throw std::runtime_error(
          "Failed to find DST node in Prototype2RawTowerBuilder::CreateNodes");
    }

  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode*>(dstiter.findFirst(
      "PHCompositeNode", detector));
  if (!DetNode)
    {
      DetNode = new PHCompositeNode(detector);
      dstNode->addNode(DetNode);
    }

  // Create the tower nodes on the tree
  _towers = new RawTowerContainer(RawTowerDefs::convert_name_to_caloid(detector));
  if (_sim_tower_node_prefix.length() == 0)
    {
      // no prefix, consistent with older convension
      TowerNodeName = "TOWER_" + detector;
    }
  else
    {
      TowerNodeName = "TOWER_" + _sim_tower_node_prefix + "_" + detector;
    }
  PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(_towers,
      TowerNodeName.c_str(), "PHObject");
  DetNode->addNode(towerNode);

  return;
}

short
Prototype2RawTowerBuilder::get_tower_row(const short cellrow) const
{
  short twrrow = cellrow/ncell_to_tower;
  return twrrow;
}
