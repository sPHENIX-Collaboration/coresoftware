#include "RawTowerDigitizer.h"

#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDeadMap.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerv1.h>

#include <g4detectors/PHG4CylinderCell.h>
#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4detectors/PHG4CylinderCellDefs.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>

using namespace std;

RawTowerDigitizer::RawTowerDigitizer(const std::string &name)
  : SubsysReco(name)
  , _digi_algorithm(kNo_digitization)
  , _sim_towers(nullptr)
  , _raw_towers(nullptr)
  , rawtowergeom(nullptr)
  , m_deadmap(nullptr)
  , detector("NONE")
  , _sim_tower_node_prefix("SIM")
  , _raw_tower_node_prefix("RAW")
  ,  //
  _photonelec_yield_visible_GeV(NAN)
  ,  //default to invalid
  _photonelec_ADC(NAN)
  ,  //default to invalid
  _pedstal_central_ADC(NAN)
  ,  //default to invalid
  _pedstal_width_ADC(NAN)
  ,  //default to invalid
  _zero_suppression_ADC(0)
  ,  //default to apply no zero suppression
  _tower_type(-1)
  , _timer(PHTimeServer::get()->insert_new(name))
{
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  seed = PHRandomSeed();  // fixed seed handled in PHRandomSeed()
  cout << Name() << " Random Seed: " << seed << endl;
  gsl_rng_set(RandomGenerator, seed);
}

RawTowerDigitizer::~RawTowerDigitizer()
{
  gsl_rng_free(RandomGenerator);
}

void RawTowerDigitizer::set_seed(const unsigned int iseed)
{
  seed = iseed;
  gsl_rng_set(RandomGenerator, seed);
}

int RawTowerDigitizer::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
         << "DST Node missing, doing nothing." << endl;
    exit(1);
  }

  try
  {
    CreateNodes(topNode);
  }
  catch (std::exception &e)
  {
    cout << e.what() << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int RawTowerDigitizer::process_event(PHCompositeNode *topNode)
{
  if (Verbosity())
  {
    cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
         << "Process event entered. "
         << "Digitalization method: ";

    if (kNo_digitization)
    {
      cout << "directly pass the energy of sim tower to digitalized tower";
    }
    else if (kSimple_photon_digitization)
    {
      cout << "simple digitization with photon statistics, ADC conversion and pedstal";
    }
    cout << endl;
  }
  // loop over all possible towers, even empty ones. The digitization can add towers containing
  // pedestals
  RawTowerGeomContainer::ConstRange all_towers = rawtowergeom->get_tower_geometries();

  double deadChanEnergy = 0;

  for (RawTowerGeomContainer::ConstIterator it = all_towers.first;
       it != all_towers.second; ++it)
  {
    const RawTowerDefs::keytype key = it->second->get_id();
    if (_tower_type >= 0)
    {
      // Skip towers that don't match the type we are supposed to digitize
      if (_tower_type != it->second->get_tower_type())
      {
        continue;
      }
    }

    RawTower *sim_tower = _sim_towers->getTower(key);
    if (m_deadmap)
    {
      if (m_deadmap->isDeadTower(key))
      {
        if (sim_tower) deadChanEnergy += sim_tower->get_energy();

        sim_tower = nullptr;

        if (Verbosity() >= VERBOSITY_MORE)
        {
          cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
               << " apply dead tower " << key << endl;
        }
      }
    }

    RawTower *digi_tower = nullptr;

    if (_digi_algorithm == kNo_digitization)
    {
      // for no digitization just copy existing towers
      if (sim_tower)
      {
        digi_tower = new RawTowerv1(*sim_tower);
      }
    }
    else if (_digi_algorithm == kSimple_photon_digitization)
    {
      // for photon digitization towers can be created if sim_tower is null pointer
      digi_tower = simple_photon_digitization(sim_tower);
    }
    else
    {
      cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
           << " invalid digitization algorithm #" << _digi_algorithm
           << endl;

      return Fun4AllReturnCodes::ABORTRUN;
    }

    if (digi_tower)
    {
      _raw_towers->AddTower(key, digi_tower);

      if (Verbosity() >= VERBOSITY_MORE)
      {
        cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
             << " output tower:"
             << endl;
        digi_tower->identify();
      }
    }
  }

  if (Verbosity())
  {
    cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
         << "input sum energy = " << _sim_towers->getTotalEdep() << " GeV"
         << ", dead channel masked energy = " << deadChanEnergy << " GeV"
         << ", output sum digitalized value = " << _raw_towers->getTotalEdep() << " ADC"
         << endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

RawTower *
RawTowerDigitizer::simple_photon_digitization(RawTower *sim_tower)
{
  RawTower *digi_tower = nullptr;

  double energy = 0;
  if (sim_tower)
  {
    energy = sim_tower->get_energy();
  }
  const double photon_count_mean = energy * _photonelec_yield_visible_GeV;
  const int photon_count = gsl_ran_poisson(RandomGenerator, photon_count_mean);
  const int signal_ADC = floor(photon_count / _photonelec_ADC);

  const double pedstal = _pedstal_central_ADC + ((_pedstal_width_ADC > 0) ? gsl_ran_gaussian(RandomGenerator, _pedstal_width_ADC) : 0);
  const int sum_ADC = signal_ADC + (int) pedstal;

  if (sum_ADC > _zero_suppression_ADC)
  {
    // create new digitalizaed tower
    if (sim_tower)
    {
      digi_tower = new RawTowerv1(*sim_tower);
    }
    else
    {
      digi_tower = new RawTowerv1();
    }
    digi_tower->set_energy((double) sum_ADC);
  }

  if (Verbosity() >= 2)
  {
    cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
         << endl;

    cout << "input: ";
    if (sim_tower)
    {
      sim_tower->identify();
    }
    else
    {
      cout << "None" << endl;
    }
    cout << "output based on "
         << "sum_ADC = " << sum_ADC << ", zero_sup = "
         << _zero_suppression_ADC << " : ";
    if (digi_tower)
    {
      digi_tower->identify();
    }
    else
    {
      cout << "None" << endl;
    }
  }

  return digi_tower;
}

void RawTowerDigitizer::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
         << "Run Node missing, doing nothing." << endl;
    throw std::runtime_error(
        "Failed to find Run node in RawTowerDigitizer::CreateNodes");
  }

  TowerGeomNodeName = "TOWERGEOM_" + detector;
  rawtowergeom = findNode::getClass<RawTowerGeomContainer>(topNode, TowerGeomNodeName.c_str());
  if (!rawtowergeom)
  {
    cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
         << " " << TowerGeomNodeName << " Node missing, doing bail out!"
         << endl;
    throw std::runtime_error("Failed to find " + TowerGeomNodeName + " node in RawTowerDigitizer::CreateNodes");
  }

  if (Verbosity() >= 1)
  {
    rawtowergeom->identify();
  }

  const string deadMapName = "DEADMAP_" + detector;
  m_deadmap = findNode::getClass<RawTowerDeadMap>(topNode, deadMapName);
  if (m_deadmap)
  {
    cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
         << " use dead map: ";
    m_deadmap->identify();
  }

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
         << "DST Node missing, doing nothing." << endl;
    throw std::runtime_error("Failed to find DST node in RawTowerDigitizer::CreateNodes");
  }

  SimTowerNodeName = "TOWER_" + _sim_tower_node_prefix + "_" + detector;
  _sim_towers = findNode::getClass<RawTowerContainer>(dstNode, SimTowerNodeName.c_str());
  if (!_sim_towers)
  {
    cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
         << " " << SimTowerNodeName << " Node missing, doing bail out!"
         << endl;
    throw std::runtime_error("Failed to find " + SimTowerNodeName + " node in RawTowerDigitizer::CreateNodes");
  }

  // Create the tower nodes on the tree
  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", detector));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode(detector);
    dstNode->addNode(DetNode);
  }

  // Be careful as a previous digitizer may have been registered for this detector
  RawTowerNodeName = "TOWER_" + _raw_tower_node_prefix + "_" + detector;
  _raw_towers = findNode::getClass<RawTowerContainer>(DetNode, RawTowerNodeName.c_str());
  if (!_raw_towers)
  {
    _raw_towers = new RawTowerContainer(_sim_towers->getCalorimeterID());
    PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(_raw_towers,
                                                                   RawTowerNodeName.c_str(), "PHObject");
    DetNode->addNode(towerNode);
  }

  return;
}
