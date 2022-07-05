#include "RawTowerDigitizer.h"

#include <calobase/RawTower.h>  // for RawTower
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDeadMap.h>
#include <calobase/RawTowerDefs.h>  // for keytype
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerv2.h>

#include <dbfile_calo_calib/CEmcCaloCalibSimpleCorrFilev1.h>
#include <dbfile_calo_calib/HcalCaloCalibSimpleCorrFilev1.h>


#include <fun4all/Fun4AllBase.h>  // for Fun4AllBase::VERBOSITY_MORE
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phparameter/PHParameters.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <cmath>
#include <cstdlib>    // for exit
#include <exception>  // for exception
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>  // for pair

using namespace std;

RawTowerDigitizer::RawTowerDigitizer(const std::string &name)
  : SubsysReco(name)
  , m_DigiAlgorithm(kNo_digitization)
  , m_SimTowers(nullptr)
  , m_RawTowers(nullptr)
  , m_RawTowerGeom(nullptr)
  , m_DeadMap(nullptr)
  , m_Detector("NONE")
  , m_SimTowerNodePrefix("SIM")
  , m_RawTowerNodePrefix("RAW")
  , m_PhotonElecYieldVisibleGeV(NAN)
  , m_PhotonElecADC(NAN)
  , m_PedstalCentralADC(NAN)
  , m_PedstalWidthADC(NAN)
  , m_pedestalFile(false)
  , m_ZeroSuppressionADC(-1000)  //default to apply no zero suppression
  , m_ZeroSuppressionFile(false)
  , m_TowerType(-1)
  , m_SiPMEffectivePixel(40000 * 4)  // sPHENIX EMCal default, 4x Hamamatsu S12572-015P MPPC [sPHENIX TDR]
  , _tower_params(name)
  , m_DoDecal(false)
  , m_DecalInverse(false)
  , m_Decal(true)
  , m_DecalFileName("")
  , m_UseConditionsDB(false)
  , m_CalDBFile(0)
{
  m_RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  m_Seed = PHRandomSeed();  // fixed seed handled in PHRandomSeed()
  // cout << Name() << " Random Seed: " << m_Seed << endl;
  gsl_rng_set(m_RandomGenerator, m_Seed);
}

RawTowerDigitizer::~RawTowerDigitizer()
{
  gsl_rng_free(m_RandomGenerator);
}

void RawTowerDigitizer::set_seed(const unsigned int iseed)
{
  m_Seed = iseed;
  gsl_rng_set(m_RandomGenerator, m_Seed);
}

int RawTowerDigitizer::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
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

  /*
  // this is for getting the file from 
  // the  conditions DB, it's reply is not used right now
  if (m_UseConditionsDB)
    {
      recoConsts *rc = recoConsts::instance();
      uint64_t timestamp = rc->get_IntFlag("RUNNUMBER");
      std::string tag = "example_tag_1";
      std::string cfg = "test";
      xpload::Configurator config(cfg);
      //std::vector<std::string> paths = xpload::fetch(tag, cfg, timestamp, config);
      xpload::Result pathres = xpload::fetch(tag, cfg, timestamp, config);
      if (pathres.paths.empty())
	{
      if (Verbosity())
	{
	  std::cout << "No paths in conditions DB found" << std::endl;
	}
	}
      else
	{
	  if (Verbosity())
	    {
	      std::cout << "Found paths:" << std::endl;
	      
	      for (const std::string &path : pathres.paths)
		{
		  std::cout << path << std::endl;
		}
	    }
	}
    }
*/



  if (m_DoDecal)
    {
      
      if (m_Detector.c_str()[0] == 'H')
	m_CalDBFile = (CaloCalibSimpleCorrFile *)  new  HcalCaloCalibSimpleCorrFilev1();
      else if (m_Detector.c_str()[0] == 'C')
	m_CalDBFile = (CaloCalibSimpleCorrFile *)  new  CEmcCaloCalibSimpleCorrFilev1();
      else 
	{
	  std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
		    << "Calo Decal requested but Detector Name not HCALOUT/IN or CEMC" 
		    << std::endl;
	  return -999;
	}
      
      if (!(m_DecalFileName == ""))
	m_CalDBFile->Open(m_DecalFileName.c_str());
      else
	m_Decal = false;
      //warnings for bad file names handled inside Open
      
    }   
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int RawTowerDigitizer::process_event(PHCompositeNode */*topNode*/)
{
  if (Verbosity())
  {
    cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
         << "Process event entered. "
         << "Digitalization method: ";

    if (m_DigiAlgorithm == kNo_digitization)
    {
      cout << "directly pass the energy of sim tower to digitalized tower";
    }
    else if (m_DigiAlgorithm == kSimple_photon_digitization)
    {
      cout << "simple digitization with photon statistics, ADC conversion and pedstal";
    }
    cout << endl;
  }
  // loop over all possible towers, even empty ones. The digitization can add towers containing
  // pedestals
  RawTowerGeomContainer::ConstRange all_towers = m_RawTowerGeom->get_tower_geometries();

  double deadChanEnergy = 0;

  for (RawTowerGeomContainer::ConstIterator it = all_towers.first;
       it != all_towers.second; ++it)
  {
    const RawTowerDefs::keytype key = it->second->get_id();
    RawTowerDefs::CalorimeterId caloid = RawTowerDefs::decode_caloid(key);
    const int eta = it->second->get_bineta();
    const int phi = it->second->get_binphi();

    if (caloid == RawTowerDefs::LFHCAL) 
    {
      const int l = it->second->get_binl();
      if (m_ZeroSuppressionFile == true)
      {
        const string zsName = "ZS_ADC_eta" + to_string(eta) + "_phi" + to_string(phi) + "_l" + to_string(l);
        m_ZeroSuppressionADC =
          _tower_params.get_double_param(zsName);
      }

      if (m_pedestalFile == true)
      {
        const string pedCentralName = "PedCentral_ADC_eta" + to_string(eta) + "_phi" + to_string(phi) + "_l" + to_string(l);
        m_PedstalCentralADC =
          _tower_params.get_double_param(pedCentralName);
        const string pedWidthName = "PedWidth_ADC_eta" + to_string(eta) + "_phi" + to_string(phi) + "_l" + to_string(l);
        m_PedstalWidthADC =
          _tower_params.get_double_param(pedWidthName);
      }
    }
    else 
    {
      if (m_ZeroSuppressionFile == true)
      {
        const string zsName = "ZS_ADC_eta" + to_string(eta) + "_phi" + to_string(phi);
        m_ZeroSuppressionADC =
          _tower_params.get_double_param(zsName);
      }

      if (m_pedestalFile == true)
      {
        const string pedCentralName = "PedCentral_ADC_eta" + to_string(eta) + "_phi" + to_string(phi);
        m_PedstalCentralADC =
          _tower_params.get_double_param(pedCentralName);
        const string pedWidthName = "PedWidth_ADC_eta" + to_string(eta) + "_phi" + to_string(phi);
        m_PedstalWidthADC =
          _tower_params.get_double_param(pedWidthName);
      }
    } 

    
    
    if (m_TowerType >= 0)
    {
      // Skip towers that don't match the type we are supposed to digitize
      if (m_TowerType != it->second->get_tower_type())
      {
        continue;
      }
    }

    RawTower *sim_tower = m_SimTowers->getTower(key);
    if (m_DeadMap)
    {
      if (m_DeadMap->isDeadTower(key))
      {
        if (sim_tower) deadChanEnergy += sim_tower->get_energy();

        sim_tower = nullptr;

        if (Verbosity() >= VERBOSITY_MORE)
        {
          cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
               << " apply dead tower " << key << endl;
        }
      }
    }

    RawTower *digi_tower = nullptr;

    if (m_DigiAlgorithm == kNo_digitization)
    {
      // for no digitization just copy existing towers
      if (sim_tower)
      {
        digi_tower = new RawTowerv2(*sim_tower);
      }
    }
    else if (m_DigiAlgorithm == kSimple_photon_digitization)
    {
      // for photon digitization towers can be created if sim_tower is null pointer
      digi_tower = simple_photon_digitization(sim_tower);
    }
    else if (m_DigiAlgorithm == kSiPM_photon_digitization)
    {
      // for photon digitization towers can be created if sim_tower is null pointer
      digi_tower = sipm_photon_digitization(sim_tower);
    }
    else
    {
      cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
           << " invalid digitization algorithm #" << m_DigiAlgorithm
           << endl;

      return Fun4AllReturnCodes::ABORTRUN;
    }

    if (digi_tower)
    {

      if (m_DoDecal && m_Decal)
	{
	  float decal_fctr = m_CalDBFile->getCorr(eta,phi);
	  
	  if (m_DecalInverse)
	    decal_fctr = 1.0/decal_fctr;
	  float e_dec = digi_tower->get_energy();
	  digi_tower->set_energy(e_dec*decal_fctr);
	}

      m_RawTowers->AddTower(key, digi_tower);

      if (Verbosity() >= VERBOSITY_MORE)
      {
        cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
             << " output tower:"
             << endl;
        digi_tower->identify();
      }
    }
  }

  if (Verbosity())
  {
    cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
         << "input sum energy = " << m_SimTowers->getTotalEdep() << " GeV"
         << ", dead channel masked energy = " << deadChanEnergy << " GeV"
         << ", output sum digitalized value = " << m_RawTowers->getTotalEdep() << " ADC"
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
  const double photon_count_mean = energy * m_PhotonElecYieldVisibleGeV;
  const int photon_count = gsl_ran_poisson(m_RandomGenerator, photon_count_mean);

  const int signal_ADC = floor(photon_count / m_PhotonElecADC);
  const double pedestal = m_PedstalCentralADC + ((m_PedstalWidthADC > 0) ? gsl_ran_gaussian(m_RandomGenerator, m_PedstalWidthADC) : 0);

  const int sum_ADC = signal_ADC + (int) pedestal;

  double  sum_ADC_d = (double) sum_ADC; 

  // temporary fix to remove digitization
  // for decalibration tests
  // to be replaced permanently for all cases 
  // with sim of pulse extraction/fitting  
  if (m_DoDecal)
    {
  
      double signal_ADC_d = 1.*photon_count;
      signal_ADC_d /= m_PhotonElecADC;
    
      sum_ADC_d  = signal_ADC_d + pedestal;
    }

  if (sum_ADC > m_ZeroSuppressionADC)
  {
    // create new digitalizaed tower
    if (sim_tower)
    {
      digi_tower = new RawTowerv2(*sim_tower);
    }
    else
    {
      digi_tower = new RawTowerv2();
    }
    digi_tower->set_energy(sum_ADC_d);
  }

  if (Verbosity() >= 2)
  {
    cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
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
         << m_ZeroSuppressionADC << " : ";
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

RawTower *
RawTowerDigitizer::sipm_photon_digitization(RawTower *sim_tower)
{
  RawTower *digi_tower = nullptr;

  double energy = 0;
  if (sim_tower)
  {
    energy = sim_tower->get_energy();
  }

  int signal_ADC = 0;

  if (energy > 0)
  {
    const double photon_count_mean = energy * m_PhotonElecYieldVisibleGeV;
    const double poission_param_per_pixel = photon_count_mean / m_SiPMEffectivePixel;
    const double prob_activated_per_pixel = gsl_cdf_poisson_Q(0, poission_param_per_pixel);
    const double active_pixel = gsl_ran_binomial(m_RandomGenerator, prob_activated_per_pixel, m_SiPMEffectivePixel);
    signal_ADC = floor(active_pixel / m_PhotonElecADC);
  }

  const double pedstal = m_PedstalCentralADC + ((m_PedstalWidthADC > 0) ? gsl_ran_gaussian(m_RandomGenerator, m_PedstalWidthADC) : 0);
  const int sum_ADC = signal_ADC + (int) pedstal;

  if (sum_ADC > m_ZeroSuppressionADC)
  {
    // create new digitalizaed tower
    if (sim_tower)
    {
      digi_tower = new RawTowerv2(*sim_tower);
    }
    else
    {
      digi_tower = new RawTowerv2();
    }
    digi_tower->set_energy((double) sum_ADC);
  }

  if (Verbosity() >= 2)
  {
    cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
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
         << m_ZeroSuppressionADC << " : ";
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
    cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
         << "Run Node missing, doing nothing." << endl;
    throw std::runtime_error(
        "Failed to find Run node in RawTowerDigitizer::CreateNodes");
  }

  string TowerGeomNodeName = "TOWERGEOM_" + m_Detector;
  m_RawTowerGeom = findNode::getClass<RawTowerGeomContainer>(topNode, TowerGeomNodeName);
  if (!m_RawTowerGeom)
  {
    cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
         << " " << TowerGeomNodeName << " Node missing, doing bail out!"
         << endl;
    throw std::runtime_error("Failed to find " + TowerGeomNodeName + " node in RawTowerDigitizer::CreateNodes");
  }

  if (Verbosity() >= 1)
  {
    m_RawTowerGeom->identify();
  }

  const string deadMapName = "DEADMAP_" + m_Detector;
  m_DeadMap = findNode::getClass<RawTowerDeadMap>(topNode, deadMapName);
  if (m_DeadMap)
  {
    cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
         << " use dead map: ";
    m_DeadMap->identify();
  }

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
         << "DST Node missing, doing nothing." << endl;
    throw std::runtime_error("Failed to find DST node in RawTowerDigitizer::CreateNodes");
  }

  string SimTowerNodeName = "TOWER_" + m_SimTowerNodePrefix + "_" + m_Detector;
  m_SimTowers = findNode::getClass<RawTowerContainer>(dstNode, SimTowerNodeName);
  if (!m_SimTowers)
  {
    cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
         << " " << SimTowerNodeName << " Node missing, doing bail out!"
         << endl;
    throw std::runtime_error("Failed to find " + SimTowerNodeName + " node in RawTowerDigitizer::CreateNodes");
  }

  // Create the tower nodes on the tree
  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", m_Detector));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode(m_Detector);
    dstNode->addNode(DetNode);
  }

  // Be careful as a previous digitizer may have been registered for this detector
  string RawTowerNodeName = "TOWER_" + m_RawTowerNodePrefix + "_" + m_Detector;
  m_RawTowers = findNode::getClass<RawTowerContainer>(DetNode, RawTowerNodeName);
  if (!m_RawTowers)
  {
    m_RawTowers = new RawTowerContainer(m_SimTowers->getCalorimeterID());
    PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(m_RawTowers,
                                                                   RawTowerNodeName, "PHObject");
    DetNode->addNode(towerNode);
  }

  return;
}
