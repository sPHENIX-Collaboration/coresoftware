///===========================================================================
/*!Migrated from analysis/JS-Jet/JetUESize-v2Psi2Centrality 
 * Original Author: Micah Meskowitz
 * Author: Jinglin Liu
 * First Commit Date: 08.17.2025
 *
 *
 * Underlying event QA module
 */ 
///===========================================================================

#include "UEvsEtaCentrality.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>

#include <phool/PHCompositeNode.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

//#include <g4jets/JetMap.h>
//#include <g4jets/Jetv1.h>

#include <centrality/CentralityInfo.h>

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>

#include <calotrigger/TriggerAnalyzer.h>

#include <jetbackground/TowerBackground.h>

#include <qautils/QAHistManagerDef.h>

#include <TTree.h>
#include <TStyle.h>
#include <TH1F.h>
#include <TH2F.h>

//____________________________________________________________________________..
UEvsEtaCentrality::UEvsEtaCentrality(const std::string& modulename)
  : SubsysReco(modulename)
  , m_moduleName(modulename)
{
  if (m_config.debug && (Verbosity() > 1)){
    std::cout << "UEvsEtaCentrality::UEvsEtaCentrality(const std::string &moduleName) Calling ctor" << std::endl;
  }
}

//____________________________________________________________________________..
// Module constructor accepting a configuration
UEvsEtaCentrality::UEvsEtaCentrality(const Config& config)
  : SubsysReco(config.moduleName)
  , m_config(config)
{
  // print debug message
  if (m_config.debug && (Verbosity() > 1))
  {
    std::cout << "UEvsEtaCentrality::UEvsEtaCentrality(Config&) Calling ctor" << std::endl;
  }
}  // end ctor(Config&)

//____________________________________________________________________________..
UEvsEtaCentrality::~UEvsEtaCentrality()
{
  // print debug message
  if (m_config.debug && (Verbosity() > 1)){
    std::cout << "UEvsEtaCentrality::~UEvsEtaCentrality() Calling dtor" << std::endl;
  }
}

//____________________________________________________________________________..
int UEvsEtaCentrality::Init(PHCompositeNode* /*topNode*/)
{
  if (m_config.debug && (Verbosity() > 1)){
    std::cout << "UEvsEtaCentrality::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  }

  // initialize trigger analyzer
  delete m_analyzer;
  m_analyzer = new TriggerAnalyzer();

  InitHistManager();

  // make sure module name is lower case
  std::string smallModuleName = m_moduleName;
  std::transform(
      smallModuleName.begin(),
      smallModuleName.end(),
      smallModuleName.begin(),
      ::tolower);

  // construct histogram names
  std::vector<std::string> vecHistNames = {
      "v2_cent",
      "psi2_cent",
      "ihcaleta_cent0_20",
      "ohcaleta_cent0_20",
      "emcaleta_cent0_20",
      "ihcaleta_cent20_50",
      "ohcaleta_cent20_50",
      "emcaleta_cent20_50",
      "ihcaleta_cent50_100",
      "ohcaleta_cent50_100",
      "emcaleta_cent50_100",
      "ihcaleta",
      "ohcaleta",
      "emcaleta",
  };
  for (auto &vecHistName : vecHistNames)
  {
    vecHistName.insert(0, "h_" + smallModuleName + "_");
    if (!m_config.histTag.empty())
    {
      vecHistName.append("_" + m_config.histTag);
    }
  }

  // Initiate histograms
  hv2_cent = new TH2F(vecHistNames[0].data(),"",10,0,100,50,0,0.5);
  hPsi2_cent = new TH2F(vecHistNames[1].data(),"",10,0,100,50,-1.57,1.57);
  hUEiHcalEta_Cent0_20 = new TH2F(vecHistNames[2].data(), "", 24, -1.1,1.1, 48, 0, 0.25);
  hUEoHcalEta_Cent0_20 = new TH2F(vecHistNames[3].data(), "", 24, -1.1,1.1, 48, 0, 0.5);
  hUEemcalEta_Cent0_20 = new TH2F(vecHistNames[4].data(), "", 24, -1.1,1.1, 48, 0, 1.5);
  
  hUEiHcalEta_Cent20_50 = new TH2F(vecHistNames[5].data(), "", 24, -1.1,1.1, 48, 0, 0.25);
  hUEoHcalEta_Cent20_50 = new TH2F(vecHistNames[6].data(), "", 24, -1.1,1.1, 48, 0, 0.5);
  hUEemcalEta_Cent20_50 = new TH2F(vecHistNames[7].data(), "", 24, -1.1,1.1, 48, 0, 1.5);

  hUEiHcalEta_Cent50_100 = new TH2F(vecHistNames[8].data(), "", 24, -1.1,1.1, 48, 0, 0.25);
  hUEoHcalEta_Cent50_100 = new TH2F(vecHistNames[9].data(), "", 24, -1.1,1.1, 48, 0, 0.5);
  hUEemcalEta_Cent50_100 = new TH2F(vecHistNames[10].data(), "", 24, -1.1,1.1, 48, 0, 1.5);

  hUEiHcalEta = new TH2F(vecHistNames[11].data(), "", 24, -1.1,1.1, 48, 0, 0.25);
  hUEoHcalEta = new TH2F(vecHistNames[12].data(), "", 24, -1.1,1.1, 48, 0, 0.5);
  hUEemcalEta = new TH2F(vecHistNames[13].data(), "", 24, -1.1,1.1, 48, 0, 1.5);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int UEvsEtaCentrality::InitRun(PHCompositeNode* /*topNode*/)
{
  if (m_config.debug && (Verbosity() > 1))
  {
    std::cout << "UEvsEtaCentrality::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int UEvsEtaCentrality::process_event(PHCompositeNode *topNode)
{ 

  if (m_config.debug && (Verbosity() > 1))
  {
    std::cout << "UEvsEtaCentrality::process_event(PHCompositeNode* topNode) Processing Event" << std::endl;
  }

  // if needed, check if selected trigger fired
  if (m_config.doTrgSelect)
  {
    m_analyzer->decodeTriggers(topNode);
    bool hasTrigger = JetQADefs::DidTriggerFire(m_config.trgToSelect, m_analyzer);
    if (!hasTrigger)
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }

  //centrality
  CentralityInfo* cent_node = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
  if (!cent_node)
  {
    std::cout
      << "UEvsEtaCentrality::process_event - Error can not find centrality node "
      << std::endl;
    exit(-1);
  }  

  //calorimeter towers
  TowerInfoContainer *towersEM3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
  TowerInfoContainer *towersIH3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
  TowerInfoContainer *towersOH3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
  RawTowerGeomContainer *tower_geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  RawTowerGeomContainer *tower_geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  if(!towersEM3 || !towersIH3 || !towersOH3){
    std::cout
      <<"UEvsEtaCentrality::process_event - Error can not find raw tower node "
      << std::endl;
    exit(-1);
  }

  if(!tower_geom || !tower_geomOH){
    std::cout
      <<"UEvsEtaCentrality::process_event - Error can not find raw tower geometry "
      << std::endl;
    exit(-1);
  }

  //underlying event
  TowerBackground *background = findNode::getClass<TowerBackground>(topNode, "TowerInfoBackground_Sub2");
  if(!background){
    std::cout<<"Can't get background. Exiting"<<std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }
   
  float background_v2 = 0;
  float background_Psi2 = 0;
  float m_centrality = 0;

  m_centrality =  cent_node->get_centile(CentralityInfo::PROP::mbd_NS); 
  background_v2 = background->get_v2();
  background_Psi2 = background->get_Psi2();

  // print centrality debug message
  if (m_config.debug && (Verbosity() > 1))
  {
    std::cout << "UEvsEtaCentrality::process_event, centrality = " << m_centrality << std::endl;
  }

  hv2_cent->Fill(m_centrality, background_v2);
  hPsi2_cent->Fill(m_centrality, background_Psi2);
    
  int i=0;    
  for (i=0;i<=23;i++){  
    float UEi = background->get_UE(1).at(i); //24 eta bins for each detector Inner Hcal
    float UEo = background->get_UE(2).at(i); //24 eta bins for each detector, Outer Hcal
    float UEe = background->get_UE(0).at(i); //24 eta bins for each detector, EMCal
   
    double eta = tower_geom->get_etacenter(i); //eta comes from IHCal

    // all centrality histograms
    hUEiHcalEta->Fill(eta, UEi);
    hUEoHcalEta->Fill(eta, UEo);
    hUEemcalEta->Fill(eta, UEe);

    if (m_centrality > 0 && m_centrality <= 20){
      hUEiHcalEta_Cent0_20->Fill(eta, UEi);                                                                                  
      hUEoHcalEta_Cent0_20->Fill(eta, UEo);                                                                                  
      hUEemcalEta_Cent0_20->Fill(eta, UEe);
    }
    if (m_centrality > 20 && m_centrality <= 50){
      hUEiHcalEta_Cent20_50->Fill(eta, UEi);                                                                                 
      hUEoHcalEta_Cent20_50->Fill(eta, UEo);                                                                                
      hUEemcalEta_Cent20_50->Fill(eta, UEe);
    }
    if (m_centrality > 50 && m_centrality <= 100){
      hUEiHcalEta_Cent50_100->Fill(eta, UEi);
      hUEoHcalEta_Cent50_100->Fill(eta, UEo);
      hUEemcalEta_Cent50_100->Fill(eta, UEe);
    }
  } // end of eta bin loop

  // std::cout << "UEvsEtaCentrality::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
} // end process_event

//____________________________________________________________________________..
int UEvsEtaCentrality::ResetEvent(PHCompositeNode* /*topNode*/)
{
  // std::cout << "UEvsEtaCentrality::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;

  // IMPORTANT!! CLEAR YOUR VECTORS AND RESET YOUR TREE VARIABLES HERE!!

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int UEvsEtaCentrality::EndRun(const int runnumber)
{
  if (m_config.debug && (Verbosity() > 1))
  {
    std::cout << "UEvsEtaCentrality::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int UEvsEtaCentrality::End(PHCompositeNode* /*topNode*/)
{

  m_manager->registerHisto(hv2_cent);
  m_manager->registerHisto(hPsi2_cent);
  m_manager->registerHisto(hUEiHcalEta_Cent0_20);
  m_manager->registerHisto(hUEoHcalEta_Cent0_20);
  m_manager->registerHisto(hUEemcalEta_Cent0_20);
  m_manager->registerHisto(hUEiHcalEta_Cent20_50);
  m_manager->registerHisto(hUEoHcalEta_Cent20_50);
  m_manager->registerHisto(hUEemcalEta_Cent20_50);
  m_manager->registerHisto(hUEiHcalEta_Cent50_100);
  m_manager->registerHisto(hUEoHcalEta_Cent50_100);
  m_manager->registerHisto(hUEemcalEta_Cent50_100);
  m_manager->registerHisto(hUEiHcalEta);
  m_manager->registerHisto(hUEoHcalEta);
  m_manager->registerHisto(hUEemcalEta);

  if (m_config.debug && (Verbosity() > 1))
  { 
    std::cout << "UEvsEtaCentrality::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  }
 
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int UEvsEtaCentrality::Reset(PHCompositeNode* /*topNode*/)
{
  if (m_config.debug && (Verbosity() > 1))
  { 
    std::cout << "UEvsEtaCentrality::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void UEvsEtaCentrality::Print(const std::string &what) const
{
  std::cout << "UEvsEtaCentrality::Print(const std::string &what) const Printing info for " << what << std::endl;
}

// Private methods =============================================================

//____________________________________________________________________________..
// Initialize histogram manager
void UEvsEtaCentrality::InitHistManager()
{

  // print debug message
  if (m_config.debug && (Verbosity() > 1))
  {
    std::cout << "UEvsEtaCentrality::InitHistManager() Initializing histogram manager" << std::endl;
  }
  
  gStyle->SetOptTitle(0);
  m_manager = QAHistManagerDef::getHistoManager();
  if (!m_manager)
  {
    std::cerr << PHWHERE << ": PANIC! Couldn't grab histogram manager!" << std::endl;
    assert(m_manager);
  }
  return;

}  // end 'InitHistManager()'
