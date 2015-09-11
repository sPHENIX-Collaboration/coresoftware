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
#include <phool/PHRandomSeed.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/getClass.h>
#include <fun4all/recoConsts.h>

#include <iostream>
#include <stdexcept>
#include <map>

#include <gsl/gsl_rng.h>

using namespace std;

RawTowerDigitizer::RawTowerDigitizer(const std::string& name) :
    SubsysReco(name), _digi_algorithm(ksimple_photon_digitalization), //
    _sim_towers(NULL), _raw_towers(NULL), rawtowergeom(NULL),//
    detector("NONE"), //
    _photonelec_yield_visible_GeV(1.), //default to apply no digitalization
    _photonelec_ADC(1), //default to apply no digitalization
    _pedstal_central_ADC(0), //default to apply no digitalization
    _pedstal_width_ADC(0), //default to apply no digitalization
    _zero_suppression_ADC(1e-6), //default to apply a moderate zero suppression
    _timer(PHTimeServer::get()->insert_new(name))
{
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  recoConsts *rc = recoConsts::instance();
  if (rc->FlagExist("RANDOMSEED"))
    {
      seed = rc->get_IntFlag("RANDOMSEED");
    }
  else
    {
      seed = PHRandomSeed();
    }
  gsl_rng_set(RandomGenerator, seed);
}

RawTowerDigitizer::~RawTowerDigitizer()
{
  gsl_rng_free(RandomGenerator);
}

void
RawTowerDigitizer::set_seed(const unsigned int iseed)
{
  seed = iseed;
  gsl_rng_set(RandomGenerator, seed);
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

  const int phibins = rawtowergeom->get_phibins();
  const int etabins = rawtowergeom->get_etabins();

  for (int iphi = 0; iphi < phibins; ++iphi)
    for (int ieta = 0; ieta < etabins; ++ieta)
      {
        RawTower *sim_tower = _sim_towers->getTower(ieta, iphi);

        RawTower *digi_tower = NULL;

        if (_digi_algorithm == ksimple_photon_digitalization)
          digi_tower = simple_photon_digitalization(sim_tower);
        else
          {

            std::cout << PHWHERE << " invalid digitalization algorithm #"
                << _digi_algorithm << std::endl;

            if (digi_tower)
              delete digi_tower;

            return Fun4AllReturnCodes::ABORTRUN;
          }

        if (digi_tower)
          _raw_towers->AddTower(ieta, iphi, digi_tower);
      }

  if (verbosity)
    {
      std::cout << PHWHERE
      << "input sum energy = " << _raw_towers->getTotalEdep()
          << ", output sum digitalized value = " << _sim_towers->getTotalEdep()
          << std::endl;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

RawTower *
RawTowerDigitizer::simple_photon_digitalization(RawTower * sim_tower)
{

  return NULL;
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
  _sim_towers = findNode::getClass<RawTowerContainer>(dstNode,
      SimTowerNodeName.c_str());
  if (!_sim_towers)
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

