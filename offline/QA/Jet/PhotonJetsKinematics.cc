//____________________________________________________________________________..

#include "PhotonJetsKinematics.h"

// calobase includes
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

// calotrigger includes
#include <calotrigger/TriggerAnalyzer.h>

// fun4all includes
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

// phool includes
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

// qautils includes
#include <qautils/QAHistManagerDef.h>

// root includes
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>

// c++ includes
#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>
//____________________________________________________________________________..

PhotonJetsKinematics::PhotonJetsKinematics(const std::string &modulename, const std::string &inputnode, const std::string &histtag) :
  SubsysReco(modulename)
  , m_modulename(modulename)
  , m_inputnode(inputnode)
  , m_histtag(histtag)
  , m_trgToSelect(JetQADefs::GL1::MBDNSJet1)
  , m_doTrgSelect(false)
  , m_doOptHist(false)
{
  if (Verbosity() > 1)
    {
      std::cout << "PhotonJetsKinematics::PhotonJetsKinematics(const std::string &name, const std::string&outputfilename) Calling ctor" << std::endl;
    }
}
//____________________________________________________________________________..
PhotonJetsKinematics::~PhotonJetsKinematics()
{
  if (Verbosity() > 1)
    {
      std::cout << "PhotonJetsKinematics::~PhotonJetsKinematics() Calling dtor" << std::endl;
    }
  delete m_analyzer;
}

//____________________________________________________________________________..
int PhotonJetsKinematics::Init(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() > 1)
    {
      std::cout << "PhotonJetsKinematics::Init(PHCompositeNode *topNode) Initializing" << std::endl;
    }

  // initialize trigger analyzer and hist manager
  delete m_analyzer;
  m_analyzer = new TriggerAnalyzer();

  gStyle->SetOptTitle(0);
  m_manager = QAHistManagerDef::getHistoManager();
  if (!m_manager)
    {
      std::cerr << PHWHERE << "PANIC: couldn't grab histogram manager!" << std::endl;
      assert(m_manager);
    }

  // make sure module name is lower case
  std::string smallModuleName = m_modulename;
  std::transform(
		 smallModuleName.begin(),
		 smallModuleName.end(),
		 smallModuleName.begin(),
		 ::tolower);
  // construct histogram names
  std::vector<std::string> vecHistNames = {
    "emcal_cluster_chi2",
    "emcal_cluster_energy",
    "emcal_cluster_eta_phi", // optional
    "emcal_cluster_eta", // optional
    "emcal_cluster_phi", // optional
    "emcal_cluster_energy_eta", //see what eta distribution looks like with different energy
    "emcal_cluster_chi2_eta", //see what eta distribution looks like with different chi2
    "emcal_cluster_eta_with_cuts",//chi2<3,Et>500 MeV
    "emcal_cluster_energy_phi",
    "emcal_cluster_chi2_phi",
    "emcal_cluster_phi_with_cuts",//chi2 < 3, Et > 500 MeV
    "emcal_cluster_eta_phi_with_cuts",//chi2 < 3, Et > 500 MeV
    "emcal_cluster_eta_with_energy_cut",//optional: Et > 500 MeV only
    "emcal_cluster_chi2_energy"
  };

  for (auto& histName : vecHistNames)
    {
      histName.insert(0, "h_" + smallModuleName + "_");
      if (!m_histtag.empty())
	{
	  histName.append("_" + m_histtag);
	}
    }

  //initializing histograms
  if (m_doOptHist)
  {
    h_emcal_cluster_eta_phi = new TH2F(vecHistNames[2].data(), "", 48, -1.2, 1.2, 64, -M_PI, M_PI);
    h_emcal_cluster_eta_phi->GetXaxis()->SetTitle("#eta");
    h_emcal_cluster_eta_phi->GetYaxis()->SetTitle("#Phi");
    h_emcal_cluster_eta_phi->SetOption("COLZ");

    h_emcal_cluster_eta = new TH1D(vecHistNames[3].data(), "", 48, -1.2, 1.2);
    h_emcal_cluster_eta->GetXaxis()->SetTitle("#eta");

    h_emcal_cluster_phi = new TH1D(vecHistNames[4].data(), "", 64, -M_PI, M_PI);
    h_emcal_cluster_phi->GetXaxis()->SetTitle("#phi");

    h_emcal_cluster_eta_with_energy_cut = new TH1D(vecHistNames[12].data(), "", 48, -1.2, 1.2);
    h_emcal_cluster_eta_with_energy_cut->GetXaxis()->SetTitle("#eta");
  }
  h_emcal_cluster_chi2 = new TH1D(vecHistNames[0].data(), "", 30, 0, 150);
  h_emcal_cluster_chi2->GetXaxis()->SetTitle("#chi^{2}");
 
  h_emcal_cluster_energy = new TH1D(vecHistNames[1].data(), "", 100, 0, 10);
  h_emcal_cluster_energy->GetXaxis()->SetTitle("E_{T} [GeV]");

  h_emcal_cluster_energy_eta = new TH2D(vecHistNames[5].data(), "",48,-1.2,1.2,100,0,10);
  h_emcal_cluster_energy_eta->GetXaxis()->SetTitle("#eta");
  h_emcal_cluster_energy_eta->GetYaxis()->SetTitle("E_{T} [GeV]");

  h_emcal_cluster_chi2_eta = new TH2D(vecHistNames[6].data(), "",48,-1.2,1.2,150,0,150);
  h_emcal_cluster_chi2_eta->GetXaxis()->SetTitle("#eta");
  h_emcal_cluster_chi2_eta->GetYaxis()->SetTitle("#chi^{2}");

  h_emcal_cluster_eta_with_cuts = new TH1D(vecHistNames[7].data(), "", 48, -1.2, 1.2);
  h_emcal_cluster_eta_with_cuts->GetXaxis()->SetTitle("#eta");

  h_emcal_cluster_energy_phi = new TH2D(vecHistNames[8].data(), "",64,-M_PI,M_PI,100,0,10);
  h_emcal_cluster_energy_phi->GetXaxis()->SetTitle("#phi");
  h_emcal_cluster_energy_phi->GetYaxis()->SetTitle("E_{T} [GeV]");

  h_emcal_cluster_chi2_phi = new TH2D(vecHistNames[9].data(), "",64,-M_PI,M_PI,150,0,150);
  h_emcal_cluster_chi2_phi->GetXaxis()->SetTitle("#phi");
  h_emcal_cluster_chi2_phi->GetYaxis()->SetTitle("#chi2^{2}");

  h_emcal_cluster_phi_with_cuts = new TH1D(vecHistNames[10].data(), "", 64,-M_PI, M_PI);
  h_emcal_cluster_phi_with_cuts->GetXaxis()->SetTitle("#phi");

  h_emcal_cluster_eta_phi_with_cuts = new TH2F(vecHistNames[11].data(), "", 48, -1.2, 1.2, 64, -M_PI, M_PI);
  h_emcal_cluster_eta_phi_with_cuts->GetXaxis()->SetTitle("#eta");
  h_emcal_cluster_eta_phi_with_cuts->GetYaxis()->SetTitle("#phi");
  h_emcal_cluster_eta_phi_with_cuts->SetOption("COLZ");

  h_emcal_cluster_chi2_energy = new TH2F(vecHistNames[13].data(),"",100,0,10,150,0,150);
  h_emcal_cluster_chi2_energy->GetXaxis()->SetTitle("E_{T} [GeV]");
  h_emcal_cluster_chi2_energy->GetYaxis()->SetTitle("#chi2^{2}");

  return Fun4AllReturnCodes::EVENT_OK;
  
}

//____________________________________________________________________________..
int PhotonJetsKinematics::InitRun(PHCompositeNode* /*topNode*/)
{

  if (Verbosity() > 1)
    {
      std::cout << "PhotonJetsKinematics::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PhotonJetsKinematics::process_event(PHCompositeNode *topNode)
{  

  RawClusterContainer* clusterContainer = findNode::getClass<RawClusterContainer>(topNode, m_inputnode);
  if (!clusterContainer)
    {
      std::cout << PHWHERE << "PhotonJetsKinematics::process_event - Fatal Error - CLUSTER_CEMC node is missing. " << std::endl;
      return 0;
    }

  // if needed, check if trigger fired
  if (m_doTrgSelect)
    {
      m_analyzer->decodeTriggers(topNode);
      bool hasTrigger = JetQADefs::DidTriggerFire(m_trgToSelect, m_analyzer);
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
      float clusEt = clusE/std::cosh(clus_eta);//this is what really wanted to be plotted on histogram related to energy
 
      //Filling Histograms, only taking into account E vec cluster, not reco cluster

       h_emcal_cluster_chi2->Fill(clus_chisq);
       h_emcal_cluster_energy->Fill(clusEt);
       if (m_doOptHist)
       {
         h_emcal_cluster_eta_phi->Fill(clus_eta, clus_phi);
         h_emcal_cluster_eta->Fill(clus_eta);  // 1D eta plot
         h_emcal_cluster_phi->Fill(clus_phi);  // 1D phi plot
       }
       h_emcal_cluster_energy_eta->Fill(clus_eta, clusEt);//2D energy eta plot
       h_emcal_cluster_chi2_eta->Fill(clus_eta,clus_chisq);//2D chi2 eta plot
       h_emcal_cluster_energy_phi->Fill(clus_phi, clusEt);//2D energy phi plot
       h_emcal_cluster_chi2_phi->Fill(clus_phi,clus_chisq);//2D chi2 phi plot
       h_emcal_cluster_chi2_energy->Fill(clusEt,clus_chisq);
       if(clus_chisq<3 && clusEt> 0.5){
       h_emcal_cluster_eta_with_cuts->Fill(clus_eta);
       h_emcal_cluster_phi_with_cuts->Fill(clus_phi);
       h_emcal_cluster_eta_phi_with_cuts->Fill(clus_eta, clus_phi);
       }
       if (m_doOptHist)
       {
         if (clusEt>0.5)
         {
           h_emcal_cluster_eta_with_energy_cut->Fill(clus_eta);
         }
       }
     }

  if (Verbosity() > 1)
    {
      std::cout << "PhotonJetsKinematics::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PhotonJetsKinematics::ResetEvent(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() > 1)
    {
      std::cout << "PhotonJetsKinematics::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PhotonJetsKinematics::EndRun(const int runnumber)
{
  if (Verbosity() > 1) 
    {
      std::cout << "PhotonJetsKinematics::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PhotonJetsKinematics::End(PHCompositeNode* /*topNode*/)
{
  // Commenting out the following line because Jenkins keeps failing the build test :( 
  // if (Verbosity() > 1) std::cout << "PhotonJetsKinematics::End - Output to " << outfilename << std::endl;
  
  //Outputting the histograms
  if (m_doOptHist)
  {
    m_manager->registerHisto(h_emcal_cluster_eta_phi);
    m_manager->registerHisto(h_emcal_cluster_eta);  
    m_manager->registerHisto(h_emcal_cluster_phi);
    m_manager->registerHisto(h_emcal_cluster_eta_with_energy_cut);
  }
  m_manager->registerHisto(h_emcal_cluster_energy);
  m_manager->registerHisto(h_emcal_cluster_energy_eta);
  m_manager->registerHisto(h_emcal_cluster_chi2);
  m_manager->registerHisto(h_emcal_cluster_chi2_eta);  
  m_manager->registerHisto(h_emcal_cluster_eta_with_cuts);
  m_manager->registerHisto(h_emcal_cluster_phi_with_cuts);
  m_manager->registerHisto(h_emcal_cluster_energy_phi);
  m_manager->registerHisto(h_emcal_cluster_chi2_phi);
  m_manager->registerHisto(h_emcal_cluster_eta_phi_with_cuts);
  m_manager->registerHisto(h_emcal_cluster_chi2_energy);

  if (Verbosity() > 1)
    {
      std::cout << "PhotonJetsKinematics::End(PHCompositeNode *topNode) This is the End..." << std::endl; 
    }
  return Fun4AllReturnCodes::EVENT_OK;
  
}

//____________________________________________________________________________..
int PhotonJetsKinematics::Reset(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() > 1)
    {
      std::cout << "PhotonJetsKinematics::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void PhotonJetsKinematics::Print(const std::string &what) const
{
  if (Verbosity() > 1)
    {
      std::cout << "PhotonJetsKinematics::Print(const std::string &what) const Printing info for " << what << std::endl;
    }
}
