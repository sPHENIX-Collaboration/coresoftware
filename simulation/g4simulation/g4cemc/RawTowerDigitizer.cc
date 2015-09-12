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
#include <cmath>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

RawTowerDigitizer::RawTowerDigitizer(const std::string& name) :
    SubsysReco(name), _digi_algorithm(kNo_digitalization), //
    _sim_towers(NULL), _raw_towers(NULL), rawtowergeom(NULL), //
    detector("NONE"), //
    _sim_tower_node_prefix("SIM"), _raw_tower_node_prefix("RAW"), //
    _photonelec_yield_visible_GeV(NAN), //default to invalid
    _photonelec_ADC(NAN), //default to invalid
    _pedstal_central_ADC(NAN), //default to invalid
    _pedstal_width_ADC(NAN), //default to invalid
    _zero_suppression_ADC(0), //default to apply no zero suppression
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
      std::cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
          << "DST Node missing, doing nothing." << std::endl;
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
      std::cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
          << "Process event entered. " << "Digitalization method: ";

      if (kNo_digitalization)
        cout << "directly pass the energy of sim tower to digitalized tower";
      else if (kSimple_photon_digitalization)
        cout
            << "simple digitalization with photon statistics, ADC conversion and pedstal";

      std::cout << std::endl;
    }

  const int phibins = rawtowergeom->get_phibins();
  const int etabins = rawtowergeom->get_etabins();

  for (int iphi = 0; iphi < phibins; ++iphi)
    for (int ieta = 0; ieta < etabins; ++ieta)
      {
        RawTower *sim_tower = _sim_towers->getTower(ieta, iphi);

        RawTower *digi_tower = NULL;

        if (_digi_algorithm == kNo_digitalization)
          if (sim_tower)
            digi_tower = new RawTowerv1(*sim_tower);
        if (_digi_algorithm == kSimple_photon_digitalization)
          digi_tower = simple_photon_digitalization(ieta, iphi, sim_tower);
        else
          {

            std::cout << Name() << "::" << detector << "::"
                << __PRETTY_FUNCTION__ << " invalid digitalization algorithm #"
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
      std::cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
          << "input sum energy = " << _sim_towers->getTotalEdep()
          << ", output sum digitalized value = " << _raw_towers->getTotalEdep()
          << std::endl;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

RawTower *
RawTowerDigitizer::simple_photon_digitalization(int ieta, int iphi,
    RawTower * sim_tower)
{
  RawTower *digi_tower = NULL;

  double energy = 0;
  if (sim_tower)
    energy = sim_tower->get_energy();

  const double photon_count_mean = energy * _photonelec_yield_visible_GeV;
  const int photon_count = gsl_ran_poisson(RandomGenerator, photon_count_mean);
  const int signal_ADC = floor(photon_count / _photonelec_ADC);

  const double pedstal = _pedstal_central_ADC
      + ((_pedstal_width_ADC > 0) ?
          gsl_ran_gaussian(RandomGenerator, _pedstal_width_ADC) : 0);
  const int sum_ADC = signal_ADC + (int) pedstal;

  if (sum_ADC > _zero_suppression_ADC)
    {
      // create new digitalizaed tower
      if (sim_tower)
        digi_tower = new RawTowerv1(*sim_tower);
      else
        digi_tower = new RawTowerv1(ieta, iphi);

      digi_tower->set_energy((double) sum_ADC);
    }

  if (verbosity >= 2)
    {
      std::cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
          << " eta " << ieta << " phi " << iphi << std::endl;

      cout << "input: ";
      if (sim_tower)
        sim_tower->identify();
      else
        cout << "None" << endl;
      cout << "output based on " << "sum_ADC = " << sum_ADC << ", zero_sup = "
          << _zero_suppression_ADC << " : ";
      if (digi_tower)
        digi_tower->identify();
      else
        cout << "None" << endl;

    }

  return digi_tower;
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
      std::cerr << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
          << "Run Node missing, doing nothing." << std::endl;
      throw std::runtime_error(
          "Failed to find Run node in RawTowerDigitizer::CreateNodes");
    }

  TowerGeomNodeName = "TOWERGEOM_" + detector;
  rawtowergeom = findNode::getClass<RawTowerGeom>(topNode,
      TowerGeomNodeName.c_str());
  if (!rawtowergeom)
    {
      std::cerr << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
          << " " << TowerGeomNodeName << " Node missing, doing bail out!"
          << std::endl;
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
      std::cerr << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
          << "DST Node missing, doing nothing." << std::endl;
      throw std::runtime_error(
          "Failed to find DST node in RawTowerDigitizer::CreateNodes");
    }

  SimTowerNodeName = "TOWER_" + _sim_tower_node_prefix + "_" + detector;
  _sim_towers = findNode::getClass<RawTowerContainer>(dstNode,
      SimTowerNodeName.c_str());
  if (!_sim_towers)
    {
      std::cerr << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
          << " " << SimTowerNodeName << " Node missing, doing bail out!"
          << std::endl;
      throw std::runtime_error(
          "Failed to find " + SimTowerNodeName
              + " node in RawTowerDigitizer::CreateNodes");
    }

  // Create the tower nodes on the tree
  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode*>(dstiter.findFirst(
      "PHCompositeNode", detector));
  if (!DetNode)
    {
      DetNode = new PHCompositeNode(detector);
      dstNode->addNode(DetNode);
    }

  _raw_towers = new RawTowerContainer();
  RawTowerNodeName = "TOWER_" + _raw_tower_node_prefix + "_" + detector;
  PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(_raw_towers,
      RawTowerNodeName.c_str(), "PHObject");
  DetNode->addNode(towerNode);

  return;
}

