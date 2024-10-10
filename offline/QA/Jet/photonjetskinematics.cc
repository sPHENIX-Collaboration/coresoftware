//____________________________________________________________________________..

#include "photonjetskinematics.h"
#include <fun4all/PHTFileServer.h>

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


//____________________________________________________________________________..
photonjetskinematics::photonjetskinematics(const std::string &outputfilename):
SubsysReco("photonjetskinematics")
, outfilename(outputfilename)
{
 std::cout << "photonjetskinematics::photonjetskinematics(const std::string &name, const std::string&outputfilename) Calling ctor" << std::endl;
}
//____________________________________________________________________________..
photonjetskinematics::~photonjetskinematics()
{
  std::cout << "photonjetskinematics::~photonjetskinematics() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int photonjetskinematics::Init(PHCompositeNode *topNode)
{
  std::cout << "photonjetskinematics::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  //outfile = new TFile(outfilename.c_str(), "RECREATE");
  PHTFileServer::get().open(outfilename, "RECREATE");

  //initializing histograms

  h_emcal_cluster_chi2 = new TH1D("h_emcal_cluster_chi2", "h_emcal_cluster_chi2", 30, 0, 150); //1000,0,5e6
  h_emcal_cluster_chi2->GetXaxis()->SetTitle("#chi^{2}");
 

  h_emcal_cluster_energy = new TH1D("h_emcal_cluster_energy", "h_emcal_cluster_energy", 100, 0, 10);
  h_emcal_cluster_energy->GetXaxis()->SetTitle("E_{T} [GeV]");


  h_emcal_cluster_eta_phi = new TH2F("h_emcal_cluster_eta_phi", "h_emcal_cluster_eta_phi", 140, -1.2, 1.2, 64, -M_PI, M_PI);
  //eta used to be 24 increase to 100
  //TH2D changed to TH2F 
  h_emcal_cluster_eta_phi->GetXaxis()->SetTitle("#eta");
  h_emcal_cluster_eta_phi->GetYaxis()->SetTitle("#Phi");
  h_emcal_cluster_eta_phi->SetOption("COLZ");


  h_emcal_cluster_eta = new TH1D("h_emcal_cluster_eta", "Eta Distribution", 140, -1.2, 1.2);
  h_emcal_cluster_eta->GetXaxis()->SetTitle("#eta");

  h_emcal_cluster_phi = new TH1D("h_emcal_cluster_phi", "Phi Distribution", 64, -M_PI, M_PI);
  h_emcal_cluster_phi->GetXaxis()->SetTitle("#phi");

  return Fun4AllReturnCodes::EVENT_OK;
  
}

//____________________________________________________________________________..
int photonjetskinematics::InitRun(PHCompositeNode *topNode)
{
  std::cout << "photonjetskinematics::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int photonjetskinematics::process_event(PHCompositeNode *topNode)
{  

  RawClusterContainer* clusterContainer = findNode::getClass<RawClusterContainer>(topNode, "CLUSTERINFO_CEMC");
  if (!clusterContainer)
    {
      std::cout << PHWHERE << "funkyCaloStuff::process_event - Fatal Error - CLUSTER_CEMC node is missing. " << std::endl;
      return 0;
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

       // std::cout << "E vec phi = " << clus_phi << std::endl;

       // std::cout << "E vec energy = " << clusE << std::endl;

       // std::cout << "E vec chi2  = " << clus_chisq << std::endl;

       // std::cout << "E vec eta = " << clus_eta << std::endl;
         
  }

  // std::cout << "photonjetskinematics::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int photonjetskinematics::ResetEvent(PHCompositeNode *topNode)
{
  // std::cout << "photonjetskinematics::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int photonjetskinematics::EndRun(const int runnumber)
{
  std::cout << "photonjetskinematics::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int photonjetskinematics::End(PHCompositeNode *topNode)
{
  std::cout << "photonjetskinematics::End - Output to " << outfilename << std::endl;
  PHTFileServer::get().cd(outfilename);

  //outfile->cd("photonjetskinematics");

  //Outputting the histograms
  h_emcal_cluster_eta_phi->Write("h_emcal_cluster_eta_phi");
  h_emcal_cluster_energy->Write("h_emcal_cluster_energy");
  h_emcal_cluster_chi2->Write("h_emcal_cluster_chi2");
  h_emcal_cluster_eta->Write("h_emcal_cluster_eta");  // Save 1D eta plot
  h_emcal_cluster_phi->Write("h_emcal_cluster_phi");  // Save 1D phi plot

  std::cout << "photonjetskinematics::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  //  outfile->Write();                                                                                                                                              
  // outfile->Close(); 
  return Fun4AllReturnCodes::EVENT_OK;
  
}

//____________________________________________________________________________..
int photonjetskinematics::Reset(PHCompositeNode *topNode)
{
 std::cout << "photonjetskinematics::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void photonjetskinematics::Print(const std::string &what) const
{
  std::cout << "photonjetskinematics::Print(const std::string &what) const Printing info for " << what << std::endl;
}
