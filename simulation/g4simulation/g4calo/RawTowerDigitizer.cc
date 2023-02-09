#include "RawTowerDigitizer.h"

#include <calobase/RawTower.h>  // for RawTower
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDeadMap.h>
#include <calobase/RawTowerDefs.h>  // for keytype
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerv2.h>

#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfov1.h>

#include <dbfile_calo_calib/CEmcCaloCalibSimpleCorrFilev1.h>
#include <dbfile_calo_calib/CaloCalibSimpleCorrFile.h>
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

#include <TSystem.h>

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

RawTowerDigitizer::RawTowerDigitizer(const std::string &name)
  : SubsysReco(name)
  , _tower_params(name)
{
  m_RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  m_Seed = PHRandomSeed();  // fixed seed handled in PHRandomSeed()
  // std::cout << Name() << " Random Seed: " << m_Seed << std::endl;
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
    std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
              << "DST Node missing, doing nothing." << std::endl;
    exit(1);
  }

  try
  {
    CreateNodes(topNode);
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (m_DoDecal)
  {
    if (m_Detector.c_str()[0] == 'H')
    {
      m_CalDBFile = new HcalCaloCalibSimpleCorrFilev1();
    }
    else if (m_Detector.c_str()[0] == 'C')
    {
      m_CalDBFile = new CEmcCaloCalibSimpleCorrFilev1();
    }
    else
    {
      std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
                << "Calo Decal requested but Detector Name not HCALOUT/IN or CEMC"
                << std::endl;
      gSystem->Exit(1);
    }

    if (!m_DecalFileName.empty())
    {
      m_CalDBFile->Open(m_DecalFileName.c_str());
    }
    else
    {
      m_Decal = false;
    }
    //warnings for bad file names handled inside Open
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int RawTowerDigitizer::process_event(PHCompositeNode */**topNode*/)
{
  if (Verbosity())
  {
    std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
              << "Process event entered. "
              << "Digitalization method: ";

    if (m_DigiAlgorithm == kNo_digitization)
    {
      std::cout << "directly pass the energy of sim tower to digitalized tower";
    }
    else if (m_DigiAlgorithm == kSimple_photon_digitization)
    {
      std::cout << "simple digitization with photon statistics, ADC conversion and pedstal";
    }
    std::cout << std::endl;
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
	      const std::string zsName = "ZS_ADC_eta" + std::to_string(eta) + "_phi" + std::to_string(phi) + "_l" + std::to_string(l);
	      m_ZeroSuppressionADC =
		_tower_params.get_double_param(zsName);
	    }
	  
	  if (m_pedestalFile == true)
	    {
	      const std::string pedCentralName = "PedCentral_ADC_eta" + std::to_string(eta) + "_phi" + std::to_string(phi) + "_l" + std::to_string(l);
	      m_PedstalCentralADC =
		_tower_params.get_double_param(pedCentralName);
	      const std::string pedWidthName = "PedWidth_ADC_eta" + std::to_string(eta) + "_phi" + std::to_string(phi) + "_l" + std::to_string(l);
	      m_PedstalWidthADC =
		_tower_params.get_double_param(pedWidthName);
	    }
	}
      else
	{
	  if (m_ZeroSuppressionFile == true)
	    {
	      const std::string zsName = "ZS_ADC_eta" + std::to_string(eta) + "_phi" + std::to_string(phi);
	      m_ZeroSuppressionADC = _tower_params.get_double_param(zsName);
	    }
	  
	  if (m_pedestalFile == true)
	    {
	      const std::string pedCentralName = "PedCentral_ADC_eta" + std::to_string(eta) + "_phi" + std::to_string(phi);
	      m_PedstalCentralADC = _tower_params.get_double_param(pedCentralName);
	      const std::string pedWidthName = "PedWidth_ADC_eta" + std::to_string(eta) + "_phi" + std::to_string(phi);
	      m_PedstalWidthADC = _tower_params.get_double_param(pedWidthName);
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
      
      RawTower *digi_tower = nullptr;
      TowerInfo *digi_towerinfo = nullptr;
    
      if(m_UseTowerInfo != 1)
	{

	  RawTower *sim_tower = m_SimTowers->getTower(key);
	  if (m_DeadMap)
	    {
	      if (m_DeadMap->isDeadTower(key))
		{
		  if (sim_tower) deadChanEnergy += sim_tower->get_energy();
		  
		  sim_tower = nullptr;
		  
		  if (Verbosity() >= VERBOSITY_MORE)
		    {
		      std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
				<< " apply dead tower " << key << std::endl;
		    }
		}
	    }
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
	      std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
			<< " invalid digitization algorithm #" << m_DigiAlgorithm
			<< std::endl;
	      
	      return Fun4AllReturnCodes::ABORTRUN;
	    }
	  
	  if (digi_tower)
	    {
	      if (m_DoDecal && m_Decal)
		{
		  float decal_fctr = m_CalDBFile->getCorr(eta, phi);
		  
		  if (m_DecalInverse)
		    {
		      decal_fctr = 1.0 / decal_fctr;
		    }
		  float e_dec = digi_tower->get_energy();
		  digi_tower->set_energy(e_dec * decal_fctr);
		}
	      m_RawTowers->AddTower(key, digi_tower);
	      if (Verbosity() >= VERBOSITY_MORE)
		{
		  std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
			    << " output tower:"
			    << std::endl;
		  digi_tower->identify();
		}
	    }
	}    
      if (m_UseTowerInfo > 0)
	{
	  unsigned int towerkey = (eta << 16U) + phi;
	  unsigned int towerindex = m_SimTowerInfos->decode_key(towerkey);
	  TowerInfo *sim_tower = m_SimTowerInfos->at(towerindex);
	  if (m_DeadMap)
	    {
	      if (m_DeadMap->isDeadTower(key))
		{
		  if (sim_tower) deadChanEnergy += sim_tower->get_energy();
		  sim_tower = nullptr;
		  if (Verbosity() >= VERBOSITY_MORE)
		    {
		      std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
				<< " apply dead tower " << key << std::endl;
		    }
		}
	    }
	  if (m_DigiAlgorithm == kNo_digitization)
	    {
	      // for no digitization just copy existing towers
	      if (sim_tower)
		{
		  digi_towerinfo = new TowerInfov1(*sim_tower);
		}
	    }
	  else if (m_DigiAlgorithm == kSimple_photon_digitization)
	    {
	      // for photon digitization towers can be created if sim_tower is null pointer
	      digi_towerinfo = simple_photon_digitization(sim_tower);
	    }
	  else if (m_DigiAlgorithm == kSiPM_photon_digitization)
	    {
	      // for photon digitization towers can be created if sim_tower is null pointer
	      digi_towerinfo = sipm_photon_digitization(sim_tower);
	    }
	  else
	    {
	      std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
			<< " invalid digitization algorithm #" << m_DigiAlgorithm
			<< std::endl;
	      
	      return Fun4AllReturnCodes::ABORTRUN;
	    }
	
	  if (digi_towerinfo)
	    {
	      if (m_DoDecal && m_Decal)
		{
		  float decal_fctr = m_CalDBFile->getCorr(eta, phi);
		  if (m_DecalInverse)
		    {
		      decal_fctr = 1.0 / decal_fctr;
		    }
		  float e_dec = digi_towerinfo->get_energy();
		  digi_towerinfo->set_energy(e_dec * decal_fctr);
		}
	      TowerInfo *digitized_towerinfo = m_RawTowerInfos->at(towerindex);
	      digitized_towerinfo->set_energy(digi_towerinfo->get_energy());
	    }	
	  delete digi_towerinfo;
	}
    }

  if (Verbosity())
  {
    if (m_UseTowerInfo !=1)
      {
	std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
		  << "input sum energy = " << m_SimTowers->getTotalEdep() << " GeV"
		  << ", dead channel masked energy = " << deadChanEnergy << " GeV"
		  << ", output sum digitalized value = " << m_RawTowers->getTotalEdep() << " ADC"
		  << std::endl;
      }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

RawTower *
RawTowerDigitizer::simple_photon_digitization(RawTower *sim_tower)
{
  RawTower *digi_tower = nullptr;

  double energy = 0.;
  if (sim_tower)
  {
    energy = sim_tower->get_energy();
  }
  const double photon_count_mean = energy * m_PhotonElecYieldVisibleGeV;
  const int photon_count = gsl_ran_poisson(m_RandomGenerator, photon_count_mean);

  const int signal_ADC = floor(photon_count / m_PhotonElecADC);
  const double pedestal = m_PedstalCentralADC + ((m_PedstalWidthADC > 0) ? gsl_ran_gaussian(m_RandomGenerator, m_PedstalWidthADC) : 0);

  const int sum_ADC = signal_ADC + (int) pedestal;

  double sum_ADC_d = (double) sum_ADC;

  // temporary fix to remove digitization
  // for decalibration tests
  // to be replaced permanently for all cases
  // with sim of pulse extraction/fitting
  if (m_DoDecal)
  {
    double signal_ADC_d = 1. * photon_count;
    signal_ADC_d /= m_PhotonElecADC;

    sum_ADC_d = signal_ADC_d + pedestal;
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
    std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
              << std::endl;

    std::cout << "input: ";
    if (sim_tower)
    {
      sim_tower->identify();
    }
    else
    {
      std::cout << "None" << std::endl;
    }
    std::cout << "output based on "
              << "sum_ADC = " << sum_ADC << ", zero_sup = "
              << m_ZeroSuppressionADC << " : ";
    if (digi_tower)
    {
      digi_tower->identify();
    }
    else
    {
      std::cout << "None" << std::endl;
    }
  }
  return digi_tower;
}



TowerInfo *
RawTowerDigitizer::simple_photon_digitization(TowerInfo *sim_tower)
{
  TowerInfo* digi_tower = nullptr;
  // TowerInfo* digi_tower = new TowerInfov1(*sim_tower);
  double energy = 0.;
  if (sim_tower)
  {
    energy = sim_tower->get_energy();
  }
  const double photon_count_mean = energy * m_PhotonElecYieldVisibleGeV;
  const int photon_count = gsl_ran_poisson(m_RandomGenerator, photon_count_mean);

  const int signal_ADC = floor(photon_count / m_PhotonElecADC);
  const double pedestal = m_PedstalCentralADC + ((m_PedstalWidthADC > 0) ? gsl_ran_gaussian(m_RandomGenerator, m_PedstalWidthADC) : 0);

  const int sum_ADC = signal_ADC + (int) pedestal;

  double sum_ADC_d = (double) sum_ADC;

  // temporary fix to remove digitization
  // for decalibration tests
  // to be replaced permanently for all cases
  // with sim of pulse extraction/fitting
  if (m_DoDecal)
  {
    double signal_ADC_d = 1. * photon_count;
    signal_ADC_d /= m_PhotonElecADC;

    sum_ADC_d = signal_ADC_d + pedestal;
  }

  if (sum_ADC > m_ZeroSuppressionADC)
    {  
      if (sim_tower)
	{
	  digi_tower = new TowerInfov1(*sim_tower);
	}
      else
	{
	  digi_tower = new TowerInfov1();
	}

      digi_tower->set_energy(sum_ADC_d);
    }

  if (Verbosity() >= 2)
  {
    std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
              << std::endl;

    std::cout << "input: ";
    if (sim_tower)
    {
      sim_tower->identify();
    }
    else
    {
      std::cout << "None" << std::endl;
    }
    std::cout << "output based on "
              << "sum_ADC = " << sum_ADC << ", zero_sup = "
              << m_ZeroSuppressionADC << " : ";
    if (digi_tower)
    {
      digi_tower->identify();
    }
    else
    {
      std::cout << "None" << std::endl;
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
    std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
              << std::endl;

    std::cout << "input: ";
    if (sim_tower)
    {
      sim_tower->identify();
    }
    else
    {
      std::cout << "None" << std::endl;
    }
    std::cout << "output based on "
              << "sum_ADC = " << sum_ADC << ", zero_sup = "
              << m_ZeroSuppressionADC << " : ";
    if (digi_tower)
    {
      digi_tower->identify();
    }
    else
    {
      std::cout << "None" << std::endl;
    }
  }

  return digi_tower;
}





TowerInfo *
RawTowerDigitizer::sipm_photon_digitization(TowerInfo *sim_tower)
{
  TowerInfo *digi_tower = nullptr;

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
      digi_tower = new TowerInfov1(*sim_tower);
    }
    else
    {
      digi_tower = new TowerInfov1();
    }
    digi_tower->set_energy((double) sum_ADC);
  }

  if (Verbosity() >= 2)
  {
    std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
              << std::endl;

    std::cout << "input: ";
    if (sim_tower)
    {
      sim_tower->identify();
    }
    else
    {
      std::cout << "None" << std::endl;
    }
    std::cout << "output based on "
              << "sum_ADC = " << sum_ADC << ", zero_sup = "
              << m_ZeroSuppressionADC << " : ";
    if (digi_tower)
    {
      digi_tower->identify();
    }
    else
    {
      std::cout << "None" << std::endl;
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
    std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
              << "Run Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find Run node in RawTowerDigitizer::CreateNodes");
  }

  std::string TowerGeomNodeName = "TOWERGEOM_" + m_Detector;
  m_RawTowerGeom = findNode::getClass<RawTowerGeomContainer>(topNode, TowerGeomNodeName);
  if (!m_RawTowerGeom)
  {
    std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
              << " " << TowerGeomNodeName << " Node missing, doing bail out!"
              << std::endl;
    throw std::runtime_error("Failed to find " + TowerGeomNodeName + " node in RawTowerDigitizer::CreateNodes");
  }

  if (Verbosity() >= 1)
  {
    m_RawTowerGeom->identify();
  }

  const std::string deadMapName = "DEADMAP_" + m_Detector;
  m_DeadMap = findNode::getClass<RawTowerDeadMap>(topNode, deadMapName);
  if (m_DeadMap)
  {
    std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
              << " use dead map: ";
    m_DeadMap->identify();
  }

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
              << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find DST node in RawTowerDigitizer::CreateNodes");
  }


  if (m_UseTowerInfo != 1)
    {
      std::string SimTowerNodeName = "TOWER_" + m_SimTowerNodePrefix + "_" + m_Detector;
      m_SimTowers = findNode::getClass<RawTowerContainer>(dstNode, SimTowerNodeName);
      if (!m_SimTowers)
	{
	  std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
		    << " " << SimTowerNodeName << " Node missing, doing bail out!"
		    << std::endl;
	  throw std::runtime_error("Failed to find " + SimTowerNodeName + " node in RawTowerDigitizer::CreateNodes");
	}
    }
  if (m_UseTowerInfo > 0 )
    {
      std::string SimTowerInfoNodeName = "TOWERINFO_" + m_SimTowerNodePrefix + "_" + m_Detector;
      m_SimTowerInfos = findNode::getClass<TowerInfoContainer>(dstNode, SimTowerInfoNodeName);
      if (!m_SimTowerInfos)
	{
	  std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
		    << " " << SimTowerInfoNodeName << " Node missing, doing bail out!"
		    << std::endl;
	  throw std::runtime_error("Failed to find " + SimTowerInfoNodeName + " node in RawTowerDigitizer::CreateNodes");
	}
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

  if (m_UseTowerInfo != 1)
    {
      std::string RawTowerNodeName = "TOWER_" + m_RawTowerNodePrefix + "_" + m_Detector;
      m_RawTowers = findNode::getClass<RawTowerContainer>(DetNode, RawTowerNodeName);
      if (!m_RawTowers)
	{
	  m_RawTowers = new RawTowerContainer(m_SimTowers->getCalorimeterID());
	  PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(m_RawTowers,
									 RawTowerNodeName, "PHObject");
	  DetNode->addNode(towerNode);
	}
    }
  if (m_UseTowerInfo > 0 )
    {
     std::string TowerInfoNodeName = "TOWERINFO_" + m_RawTowerNodePrefix + "_" + m_Detector;
     m_RawTowerInfos = findNode::getClass<TowerInfoContainer>(DetNode,TowerInfoNodeName);
     if (m_RawTowerInfos == nullptr)
       {
	 TowerInfoContainerv1::DETECTOR detec;
	 if (m_Detector == "CEMC")
	   {
	     detec = TowerInfoContainer::DETECTOR::EMCAL;
	   }
	 else if (m_Detector == "HCALIN" || m_Detector == "HCALOUT")
	   {
	     detec = TowerInfoContainer::DETECTOR::HCAL;
	   }
	 else
	   {
	     std::cout << PHWHERE << "Detector not implemented into the TowerInfoContainer object, defaulting to HCal implementation." << std::endl;
	     detec = TowerInfoContainer::DETECTOR::HCAL;
	   }
	 if (!m_RawTowerInfos)
	   {
	     m_RawTowerInfos = new TowerInfoContainerv1(detec);
	     PHIODataNode<PHObject> *towerinfoNode = new PHIODataNode<PHObject>(m_RawTowerInfos, TowerInfoNodeName, "PHObject");
	     DetNode->addNode(towerinfoNode);
	   }
       }
    }
    

  return;
}
