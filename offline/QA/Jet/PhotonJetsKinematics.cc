//____________________________________________________________________________..

#include "PhotonJetsKinematics.h"
#include <fun4all/PHTFileServer.h>

// Fun4all includes
#include <calotrigger/TriggerAnalyzer.h>

#include <fun4all/Fun4AllHistoManager.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

#include <TFile.h>
#include <phool/getClass.h>
#include "TH1.h"
#include "TMath.h"
//#include <jetbase/JetContainer.h>


// Tower includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

// Cluster includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>

// QA includes
#include <qautils/QAHistManagerDef.h>

// c++ includes
#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>
//____________________________________________________________________________..

PhotonJetsKinematics::PhotonJetsKinematics(const std::string &m_modulename, const std::string &m_inputnode) :
  SubsysReco(m_modulename)
  , modulename(m_modulename)
  , inputnode(m_inputnode)
  , histtag("AllTrig")
  , trgToSelect(JetQADefs::GL1::MBDNSJet1)
  , doTrgSelect(false)
{
  if (Verbosity() > 1) std::cout << "PhotonJetsKinematics::PhotonJetsKinematics(const std::string &name, const std::string&outputfilename) Calling ctor" << std::endl;
}
//____________________________________________________________________________..
PhotonJetsKinematics::~PhotonJetsKinematics()
{
  std::cout << "PhotonJetsKinematics::~PhotonJetsKinematics() Calling dtor" << std::endl;
  delete m_analyzer;
}

//____________________________________________________________________________..
int PhotonJetsKinematics::Init(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() > 1) std::cout << "PhotonJetsKinematics::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  // initialize trigger analyzer and hist manager
  delete m_analyzer;
  m_analyzer = new TriggerAnalyzer();
  manager = QAHistManagerDef::getHistoManager();
  if (!manager)
    {
      std::cerr << PHWHERE << "PANIC: couldn't grab histogram manager!" << std::endl;
      assert(manager);
    }

  // make sure module name is lower case
  std::string smallModuleName = modulename;
  std::transform(
		 smallModuleName.begin(),
		 smallModuleName.end(),
		 smallModuleName.begin(),
		 ::tolower);
  // construct histogram names
  std::vector<std::string> vecHistNames = {
    "emcal_cluster_chi2",
    "emcal_cluster_energy",
    "emcal_cluster_eta_phi",
    "emcal_cluster_eta",
    "emcal_cluster_phi"};

  for (auto& histName : vecHistNames)
    {
      histName.insert(0, "h_" + smallModuleName + "_");
      if (!histtag.empty())
	{
	  histName.append("_" + histtag);
	}
    }

  //initializing histograms

  h_emcal_cluster_chi2 = new TH1D(vecHistNames[0].data(), "h_emcal_cluster_chi2", 30, 0, 150);
  h_emcal_cluster_chi2->GetXaxis()->SetTitle("#chi^{2}");
 
  h_emcal_cluster_energy = new TH1D(vecHistNames[1].data(), "h_emcal_cluster_energy", 100, 0, 10);
  h_emcal_cluster_energy->GetXaxis()->SetTitle("E_{T} [GeV]");

  h_emcal_cluster_eta_phi = new TH2F(vecHistNames[2].data(), "h_emcal_cluster_eta_phi", 140, -1.2, 1.2, 64, -M_PI, M_PI);
  h_emcal_cluster_eta_phi->GetXaxis()->SetTitle("#eta");
  h_emcal_cluster_eta_phi->GetYaxis()->SetTitle("#Phi");
  h_emcal_cluster_eta_phi->SetOption("COLZ");

  h_emcal_cluster_eta = new TH1D(vecHistNames[3].data(), "Eta Distribution", 140, -1.2, 1.2);
  h_emcal_cluster_eta->GetXaxis()->SetTitle("#eta");

  h_emcal_cluster_phi = new TH1D(vecHistNames[4].data(), "Phi Distribution", 64, -M_PI, M_PI);
  h_emcal_cluster_phi->GetXaxis()->SetTitle("#phi");

  return Fun4AllReturnCodes::EVENT_OK;
  
}

//____________________________________________________________________________..
int PhotonJetsKinematics::InitRun(PHCompositeNode* /*topNode*/)
{

  if (Verbosity() > 1) std::cout << "PhotonJetsKinematics::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PhotonJetsKinematics::process_event(PHCompositeNode *topNode)
{  

  RawClusterContainer* clusterContainer = findNode::getClass<RawClusterContainer>(topNode, inputnode);
  if (!clusterContainer)
    {
      std::cout << PHWHERE << "funkyCaloStuff::process_event - Fatal Error - CLUSTER_CEMC node is missing. " << std::endl;
      return 0;
    }

  // if needed, check if trigger fired
  if (doTrgSelect)
    {
      m_analyzer->decodeTriggers(topNode);
      bool hasTrigger = JetQADefs::DidTriggerFire(trgToSelect, m_analyzer);
      if (!hasTrigger)
	{
	  return Fun4AllReturnCodes::EVENT_OK;
	}
    }

  RawClusterContainer::ConstIterator clusterIter;
  RawClusterContainer::ConstRange clusterEnd = clusterContainer->getClusters();

  for (clusterIter = clusterEnd.first; clusterIter != clusterEnd.second; clusterIter++)
    {
      RawCluster* recoCluster = clusterIter->second;
      CLHEP::Hep3Vector vertex(0, 0, 0);
      CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*recoCluster, vertex);

      float clusE = E_vec_cluster.mag();
      float clus_eta = E_vec_cluster.pseudoRapidity();
      float clus_phi = E_vec_cluster.phi();
      float clus_chisq = recoCluster->get_chi2();
 
      //Filling Histograms, only taking into account E vec cluster, not reco cluster

       h_emcal_cluster_chi2->Fill(clus_chisq);
       h_emcal_cluster_energy->Fill(clusE);
       h_emcal_cluster_eta_phi->Fill(clus_eta, clus_phi);
       h_emcal_cluster_eta->Fill(clus_eta);  // 1D eta plot
       h_emcal_cluster_phi->Fill(clus_phi);  // 1D phi plot
   }

  // std::cout << "PhotonJetsKinematics::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
   return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PhotonJetsKinematics::ResetEvent(PHCompositeNode* /*topNode*/)
{
  // std::cout << "PhotonJetsKinematics::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PhotonJetsKinematics::EndRun(const int runnumber)
{
  if (Verbosity() > 1) std::cout << "PhotonJetsKinematics::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PhotonJetsKinematics::End(PHCompositeNode* /*topNode*/)
{
  // Commenting out the following line because Jenkins keeps failing the build test :( 
  // if (Verbosity() > 1) std::cout << "PhotonJetsKinematics::End - Output to " << outfilename << std::endl;
  
  //Outputting the histograms
  manager->registerHisto(h_emcal_cluster_eta_phi);
  manager->registerHisto(h_emcal_cluster_energy);
  manager->registerHisto(h_emcal_cluster_chi2);
  manager->registerHisto(h_emcal_cluster_eta);  
  manager->registerHisto(h_emcal_cluster_phi); 

  if (Verbosity() > 1) std::cout << "PhotonJetsKinematics::End(PHCompositeNode *topNode) This is the End..." << std::endl; 
  return Fun4AllReturnCodes::EVENT_OK;
  
}

//____________________________________________________________________________..
int PhotonJetsKinematics::Reset(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() > 1) std::cout << "PhotonJetsKinematics::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void PhotonJetsKinematics::Print(const std::string &what) const
{
  if (Verbosity() > 1) std::cout << "PhotonJetsKinematics::Print(const std::string &what) const Printing info for " << what << std::endl;
}
