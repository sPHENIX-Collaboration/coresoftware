#include "JetCalib.h"

#include "JetContainer.h"

#include <cdbobjects/CDBTF.h>  // for CDBTF1

#include <ffamodules/CDBInterface.h>

#include <ffaobjects/EventHeader.h>


#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/GlobalVertex.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <TF1.h>
#include <TSystem.h>

#include <boost/format.hpp>

#include <cstdlib>    // for exit
#include <exception>  // for exception
#include <iostream>   // for operator<<, basic_ostream
#include <stdexcept>  // for runtime_error

//____________________________________________________________________________..
JetCalib::JetCalib(const std::string &name)
  : SubsysReco(name)
{
  if (Verbosity() > 0)
  {
    std::cout << "JetCalib::JetCalib(const std::string &name) Calling ctor" << std::endl;
  }
}

//____________________________________________________________________________..
JetCalib::~JetCalib()
{
  
  delete m_etaJesFile;
  delete m_rTrkFile;
  delete m_gammaJetFile;
  
  for(auto & i : m_etaJesFunc) { delete i;}
  delete m_gammaJetFunc;
  delete m_rTrkFunc;
  if (Verbosity() > 0)
  {
    std::cout << "JetCalib::~JetCalib() Calling dtor" << std::endl;
  }
}

//____________________________________________________________________________..
int JetCalib::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator nodeIter(topNode);

  //We most likely don't need the evtHeader right now. Might need some functionality
  //later such as for centrality, but I doubt we'll have the statistics for run-by-run calibrations
  //-a. hodges
  // EventHeader *evtHeader = findNode::getClass<EventHeader>(topNode, "EventHeader");

  // if (evtHeader)
  //   {
  //     m_runNumber = evtHeader->get_RunNumber();
  //   }
  // else
  //   {
  //     m_runNumber = -1;
  //   }
 
  if(m_doEtaJes)
    {
      if(fetchCalibDir("etaJes").empty())
	{
	  std::cout << "JetCalib::process_event - No etaJes calibration available! Will apply calib factor 1" << std::endl;
	}
      else
	{
	  m_etaJesFile = new CDBTF(fetchCalibDir("etaJes"));
	  if(m_etaJesFile)
	    {
	      for(int i = 0; i < m_nEtaBins; i++)
		{
		  m_etaJesFunc[i] = m_etaJesFile->getTF((boost::format("corrFit_eta%d") % i).str().c_str());
		}
	    }
	  else
	    {
	      std::cout << "JetCalib::process_event - Could not open etaJes calibration file!" << std::endl;
	    }
	}
    }

 
    
  if(m_doInsitu)
    {
      //fetch gamma jet calibration if it exists
      if(fetchCalibDir("gammaJet").empty())
	{
	  std::cout << "No gammaJet calibration available! Will apply calib factor 1" << std::endl;
	}
      else
	{
	  m_gammaJetFile = new CDBTF(fetchCalibDir("gammaJet"));
	  if(m_gammaJetFile)
	    {
	      m_gammaJetFunc = m_gammaJetFile->getTF("corrFit");
	    }
	  else
	    {
	      std::cout << "JetCalib::process_event - Could not open rTrk calibration file!" << std::endl;
	    }
	}
      //fetch rtrk calibration if it exists
      if(fetchCalibDir("rTrk").empty())
	{
	  std::cout << "No rTrk calibration available! Will apply calib factor 1" << std::endl;
	}
      else
	{
	  m_rTrkFile = new CDBTF(fetchCalibDir("rTrk"));
	  if(m_rTrkFile)
	    {
	      m_rTrkFunc = m_rTrkFile->getTF("corrFit");
	    }
	  else
	    {
	      std::cout << "JetCalib::process_event - Could not open rTrk calibration file!" << std::endl;
	    }
	}
    }
  
  if(m_doEtaJes)
    {
      m_isEtaDependent = 1;
    }

  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << Name() << "::" << __PRETTY_FUNCTION__
              << "DST Node missing, doing nothing." << std::endl;
    exit(1);
  }
  try
  {
    CreateNodeTree(topNode);
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (Verbosity() > 0)
  {
    topNode->print();
  }

  //initialize eta bins
  for(float i = m_etaStart; i <= m_etaEnd; i += (m_etaStart - m_etaEnd)/m_nEtaBins)
    {
      m_etaBins.push_back(i);
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetCalib::process_event(PHCompositeNode *topNode)
{

  JetContainer *_raw_jets = findNode::getClass<JetContainer>(topNode,(boost::format("AntiKt_TowerInfo_r0%d_sub1") % m_radius).str());
  JetContainer *_calib_jets = findNode::getClass<JetContainer>(topNode,(boost::format("AntiKt_TowerInfo_r0%d_sub1_Calib") % m_radius).str());

  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if(!vertexmap)
    {
      std::cout << "JetCalib::process_event - Error cannot find global vertex node!" << std::endl;
      //return Fun4AllReturnCodes::EVENT_OK;
    }
  else if(vertexmap->empty())
    {
      std::cout << "JetCalib::process_event - global vertex node is empty!" << std::endl;
      //return Fun4AllReturnCodes::EVENT_OK;
    }
  else
    {
      GlobalVertex *vtx = vertexmap -> begin()->second;
      m_zvtx= vtx->get_z();
    }
  
  int ijet = 0; 
  for(auto jet : *_raw_jets)
    {

      float pt = jet->get_pt();      
      float eta  = jet->get_eta();
      int etaBin = getEtaBin(eta);
      auto calib_jet = _calib_jets -> add_jet();

      if(etaBin > -1)
	{//within calibration range
	  float corrFactor = getTotalCorrFactor(m_etaJesFunc[etaBin], m_rTrkFunc, m_gammaJetFunc,pt,m_zvtx);
      
	  calib_jet -> set_px(jet -> get_px()*corrFactor);
	  calib_jet -> set_py(jet -> get_py()*corrFactor);
	  calib_jet -> set_pz(jet -> get_pz()*corrFactor);
	  calib_jet -> set_id(ijet);
	  calib_jet -> set_isCalib(1);
	}
      else
	{
	  if(Verbosity() > 0)
	    {
	      std::cout << "JetCalib::process_event: jet eta outside calibration range, setting eta to -9999" << std::endl;
	    }
	  calib_jet -> set_px(jet -> get_px());
	  calib_jet -> set_py(jet -> get_py());
	  calib_jet -> set_pz(jet -> get_pz());
	  calib_jet -> set_id(ijet);
	  calib_jet -> set_isCalib(0);
	}
      ijet++;
      
    }
  
  if (Verbosity() >0)
    {
      std::cout << "JetCalib::process_event - started process_event with nRawJets: " <<  _raw_jets->size() << std::endl;
      std::cout << "JetCalib::process_event - ended process_event nCalibJets: " <<  _calib_jets->size() << std::endl;
      if(_calib_jets -> size() != _raw_jets->size())
	{
	  std::cout << "JetCalib::process_event - different number of raw jets vs. calib jets! Something is amiss! " << std::endl;
	}
    }


  return Fun4AllReturnCodes::EVENT_OK;
}

int JetCalib::CreateNodeTree(PHCompositeNode *topNode)
{

  PHNodeIterator iter(topNode);
  // Look for the necessary head nodes before continuing. I ABORTRUN if not found because presumably the code 
  // would just crash elsewhere.
  // Looking for the DST node, if this isn't there we're in trouble
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, aborting!" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Looking for the ANTIKT node, which carries all jet stuff
  PHCompositeNode *antiktNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "ANTIKT"));
  if (!antiktNode)
  {
    std::cout << PHWHERE << "ANTIKT node not found, aborting!" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Looking for the TOWER node
  PHCompositeNode *towerNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "TOWER"));
  if (!towerNode)
  {
    std::cout << PHWHERE << "TOWER node not found, aborting!" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // store the new jet collection
  JetContainer *test_jets = findNode::getClass<JetContainer>(topNode, (boost::format("AntiKt_TowerInfo_r0%d_sub1_Calib") %m_radius).str());//test to see if node already exists

  if (!test_jets)
    {
      //if not, we create it
      JetContainer *calib_jets = new JetContainer();
      PHIODataNode<PHObject> *calibjetNode;
     
      if (Verbosity() > 0)
	{
	  std::cout << "JetCalib::CreateNode : creating " <<  (boost::format("AntiKt_TowerInfo_r0%d_sub1_Calib") % m_radius)  << std::endl;
	}
      calibjetNode = new PHIODataNode<PHObject>(calib_jets, (boost::format("AntiKt_TowerInfo_r0%d_sub1_Calib") % m_radius).str(), "PHObject");
	
    towerNode->addNode(calibjetNode);
  }
  else
  {
    std::cout << "JetCalib::CreateNode : " <<  (boost::format("AntiKt_TowerInfo_r0%d_sub1_Calib") % m_radius) << " already exists! " << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int JetCalib::getEtaBin(float jetEta)
{
  int etaBin = -1;
  if(jetEta > m_etaBins.at(m_nEtaBins)) 
    { 
      return -1;
    }
  if(jetEta < m_etaBins.at(0)) 
    { 
      return -1;
    }
  for(int i = 0; i < m_nEtaBins; i++)
    {
      if(jetEta > m_etaBins.at(i) && jetEta <= m_etaBins.at(i+1))
	{
	  etaBin = i;
	}
    }
  return etaBin;
}

std::string JetCalib::fetchCalibDir(const char *calibType)
{
 
  std::string  calibName = (boost::format("%s_%s_r0%d_%s_BGSub%d_MC%d") % calibType % m_jetType % m_radius % m_jetInstrument % m_doBackgroundSub % m_calibyear).str();
  std::cout << "JetCalib::::InitRun Searching for jet correction: " << calibName << "" << std::endl;
  std::string calibdir = CDBInterface::instance()->getUrl(calibName);

  if(calibdir.empty())
    {
      //std::string default_time_independent_calib = "EtaJES_AntiKt_r04_caloJet_BGSub1";
      std::cout << "JetCalib::::InitRun No jet calibration of name: " << calibName << "!" << std::endl;
      return calibdir;
      //if no calibration is selected or user messes up calib options
      //calibdir = CDBInterface::instance()->getUrl(default_time_independent_calib);
    }
  else
    {
      return calibdir;
    }
    
}

 float JetCalib::getTotalCorrFactor(TF1*  etaJesTF, TF1 * rTrkTF, TF1 * gammaJetTF, float jetPt, float zvtx)
{
  float calibEtaJes = 1;
  float calibRtrk = 1;
  float calibGammaJet = 1;
  
  
  if(etaJesTF)
    {
      calibEtaJes = etaJesTF -> Eval(jetPt);
    }
  if(rTrkTF && !std::isnan(zvtx))
    {
      calibRtrk = rTrkTF -> Eval(jetPt);
    }
 if(gammaJetTF  && !std::isnan(zvtx))
    {
      calibGammaJet = gammaJetTF -> Eval(jetPt);
    }
 
 return calibEtaJes*calibRtrk*calibGammaJet;
  
}
