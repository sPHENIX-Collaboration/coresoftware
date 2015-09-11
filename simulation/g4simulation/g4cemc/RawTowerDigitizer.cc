#include "RawTowerDigitizer.h"
#include "RawTowerContainer.h"
#include "RawTowerGeomv2.h"
#include "RawTowerv1.h"
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4detectors/PHG4CylinderCell.h>
#include <g4detectors/PHG4CylinderCellDefs.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/getClass.h>

#include <iostream>
#include <stdexcept>
#include <map>

using namespace std;

RawTowerDigitizer::RawTowerDigitizer(const std::string& name) :
    SubsysReco(name), _sim_towers(NULL), _raw_towers(NULL), rawtowergeom(NULL), detector(
        "NONE"), emin(1e-6), _timer(PHTimeServer::get()->insert_new(name))
{
}

int
RawTowerDigitizer::InitRun(PHCompositeNode *topNode)
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
      return Fun4AllReturnCodes::ABORTRUN;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

int
RawTowerDigitizer::process_event(PHCompositeNode *topNode)
{
  if (verbosity)
    {
      std::cout << PHWHERE << "Process event entered" << std::endl;
    }

  double towerE = 0;
  if (verbosity)
    {
      towerE = _raw_towers->getTotalEdep();
    }

  _raw_towers->compress(emin);

  return Fun4AllReturnCodes::EVENT_OK;
}

int
RawTowerDigitizer::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void
RawTowerDigitizer::CreateNodes(PHCompositeNode *topNode)
{

  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = static_cast<PHCompositeNode*>(iter.findFirst(
      "PHCompositeNode", "RUN"));
  if (!runNode)
    {
      std::cerr << PHWHERE << "Run Node missing, doing nothing." << std::endl;
      throw std::runtime_error(
          "Failed to find Run node in RawTowerDigitizer::CreateNodes");
    }

  TowerGeomNodeName = "TOWERGEOM_" + detector;
  rawtowergeom = findNode::getClass<RawTowerGeom>(topNode,
      TowerGeomNodeName.c_str());
  if (!rawtowergeom)
    {
      std::cerr << PHWHERE << " " << TowerGeomNodeName
          << " Node missing, doing bail out!" << std::endl;
      throw std::runtime_error(
          "Failed to find " + TowerGeomNodeName
              + " node in RawTowerDigitizer::CreateNodes");
    }

  if (verbosity >= 1)
    {
      rawtowergeom->identify();
    }

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst(
      "PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cerr << PHWHERE << "DST Node missing, doing nothing." << std::endl;
      throw std::runtime_error(
          "Failed to find DST node in RawTowerDigitizer::CreateNodes");
    }

  SimTowerNodeName = "TOWER_SIM_" + detector;
  _raw_towers = findNode::getClass<RawTowerContainer>(dstNode,
      SimTowerNodeName.c_str());
  if (!_raw_towers)
    {
      std::cerr << PHWHERE << " " << SimTowerNodeName
          << " Node missing, doing bail out!" << std::endl;
      throw std::runtime_error(
          "Failed to find " + SimTowerNodeName
              + " node in RawTowerDigitizer::CreateNodes");
    }

  // Create the tower nodes on the tree
  _raw_towers = new RawTowerContainer();
  RawTowerNodeName = "TOWER_RAW_" + detector;
  PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(_raw_towers,
      RawTowerNodeName.c_str(), "PHObject");
  dstNode->addNode(towerNode);

  return;
}

