#include "JetSeedCount.h"

#include <jetbase/JetContainer.h>

#include <centrality/CentralityInfo.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TH1.h>
#include <TH2.h>

#include <algorithm>
#include <iostream>

JetSeedCount::JetSeedCount(const std::string &moduleName, const std::string &recojetname, const std::string &rawSeedName, const std::string &subSeedName, const std::string &truthjetname, const std::string &outputfilename)
  : SubsysReco(moduleName)
  , m_moduleName(moduleName)
  , m_recoJetName(recojetname)
  , m_rawSeedName(rawSeedName)
  , m_subSeedName(subSeedName)
  , m_truthJetName(truthjetname)
  , m_outputFileName(outputfilename)
  , m_histTag("AllTrig_AntiKt_Tower_r04_Sub1")
  , m_etaRange(-1, 1)
  , m_ptRange(5, 100)
  , m_doTrgSelect(false)
  , m_trgToSelect(JetQADefs::GL1::MBDNSJet1)
{
  // std::cout << "JetSeedCount::JetSeedCount(const std::string &name) Calling ctor" << std::endl;
}

int JetSeedCount::Init(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 1)
  {
    std::cout << "JetSeedCount::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  }
  if (m_writeToOutputFile)
  {
    std::cout << "Opening output file named " << m_outputFileName << std::endl;
    PHTFileServer::get().open(m_outputFileName, "RECREATE");
  }
  m_manager = QAHistManagerDef::getHistoManager();
  if (!m_manager)
  {
    std::cout << "No m_manager" << std::endl;
    assert(m_manager);
  }
  
  // make sure module name is lower case
  std::string smallModuleName = m_moduleName;
  std::transform(
      smallModuleName.begin(),
      smallModuleName.end(),
      smallModuleName.begin(),
      ::tolower);

  // construct histogram names
  std::vector<std::string> vecHistNames = {
      "rawseedcount",
      "rawpt",
      "rawptall",
      "rawetavsphi",
      "subseedcount",
      "subpt",
      "subptall",
      "subetavsphi",
      "rawseedenergyvscent",
      "subseedenergyvscent",
      "centmbd",
      "rawseedvscent",
      "subseedvscent"};
  for (auto &vecHistName : vecHistNames)
  {
    vecHistName.insert(0, "h_" + smallModuleName + "_");
    if (!m_histTag.empty())
    {
      vecHistName.append("_" + m_histTag);
    }
  }

  // make histograms
  m_hRawSeedCount = new TH1F(vecHistNames[0].data(), "Raw Seed Count per Event", 100, 0.00, 50.00);
  m_hRawSeedCount->GetXaxis()->SetTitle("Raw Seed Count per Event");
  m_hRawSeedCount->GetYaxis()->SetTitle("Number of Entries");
  m_manager->registerHisto(m_hRawSeedCount);

  m_hRawPt = new TH1F(vecHistNames[1].data(), "Raw p_{T}", 1000, 0.00, 50.00);
  m_hRawPt->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  m_hRawPt->GetYaxis()->SetTitle("Number of Entries");
  m_manager->registerHisto(m_hRawPt);

  m_hRawPt_All = new TH1F(vecHistNames[2].data(), "Raw p_{T} (all jet seeds)", 1000, 0.00, 50.00);
  m_hRawPt_All->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  m_hRawPt_All->GetYaxis()->SetTitle("Number of Entries");
  m_manager->registerHisto(m_hRawPt_All);

  m_hRawEtaVsPhi = new TH2F(vecHistNames[3].data(), "Raw Seed Eta Vs Phi", 220, -1.1, 1.1, 628, -3.14, 3.14);
  m_hRawEtaVsPhi->GetXaxis()->SetTitle("Jet #eta [Rads.]");
  m_hRawEtaVsPhi->GetYaxis()->SetTitle("Jet #phi [Rads.]");
  m_manager->registerHisto(m_hRawEtaVsPhi);

  m_hSubSeedCount = new TH1F(vecHistNames[4].data(), "Sub Seed Count per Event", 100, 0.00, 50.00);
  m_hSubSeedCount->GetXaxis()->SetTitle("Sub Seed Count per Event");
  m_hSubSeedCount->GetYaxis()->SetTitle("Number of Entries");
  m_manager->registerHisto(m_hSubSeedCount);

  m_hSubPt = new TH1F(vecHistNames[5].data(), "Sub. p_{T}", 1000, 0.00, 50.00);
  m_hSubPt->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  m_hSubPt->GetYaxis()->SetTitle("Number of Entries");
  m_manager->registerHisto(m_hSubPt);

  m_hSubPt_All = new TH1F(vecHistNames[6].data(), "Sub. p_{T} (all jet seeds)", 1000, 0.00, 50.00);
  m_hSubPt_All->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  m_hSubPt_All->GetYaxis()->SetTitle("Number of Entries");
  m_manager->registerHisto(m_hSubPt_All);

  m_hSubEtaVsPhi = new TH2F(vecHistNames[7].data(), "Sub. Seed Eta Vs Phi", 220, -1.1, 1.1, 628, -3.14, 3.14);
  m_hSubEtaVsPhi->GetXaxis()->SetTitle("Jet #eta [Rads.]");
  m_hSubEtaVsPhi->GetYaxis()->SetTitle("Jet #phi [Rads.]");
  m_manager->registerHisto(m_hSubEtaVsPhi);

  // If not in pp mode, plot quantities vs. centrality
  if (!m_inPPMode)
  {
    m_hRawSeedEnergyVsCent = new TH2F(vecHistNames[8].data(), "Raw Seed Energy Vs Centrality", 10.00, 0.00, 100.00, 100, 0.00, 50.00);
    m_hRawSeedEnergyVsCent->GetXaxis()->SetTitle("Centrality");
    m_hRawSeedEnergyVsCent->GetYaxis()->SetTitle("RawSeedEnergy");
    m_manager->registerHisto(m_hRawSeedEnergyVsCent);

    m_hSubSeedEnergyVsCent = new TH2F(vecHistNames[9].data(), "Sub Seed Energy Vs Centrality", 10.00, 0.00, 100.00, 100, 0.00, 50.00);
    m_hSubSeedEnergyVsCent->GetXaxis()->SetTitle("Centrality");
    m_hSubSeedEnergyVsCent->GetYaxis()->SetTitle("SubSeedEnergy");
    m_manager->registerHisto(m_hSubSeedEnergyVsCent);

    m_hCentMbd = new TH1F(vecHistNames[10].data(), "hCentMbd", 10, 0.00, 100.00);
    m_hCentMbd->GetXaxis()->SetTitle("Centrality (Mbd)");
    m_hCentMbd->GetYaxis()->SetTitle("Number of Entries");
    m_manager->registerHisto(m_hCentMbd);

    m_hRawSeedVsCent = new TH2F(vecHistNames[11].data(), "Raw Seed Vs Centrality", 10, 0.00, 100.00, 101, -0.5, 100.5);
    m_hRawSeedVsCent->GetXaxis()->SetTitle("Centrality");
    m_hRawSeedVsCent->GetYaxis()->SetTitle("Raw Seed Count");
    m_manager->registerHisto(m_hRawSeedVsCent);

    m_hSubSeedVsCent = new TH2F(vecHistNames[12].data(), "Sub Seed Vs Centrality", 10, 0.00, 100.00, 101, -0.5, 100.5);
    m_hSubSeedVsCent->GetXaxis()->SetTitle("Centrality");
    m_hSubSeedVsCent->GetYaxis()->SetTitle("Sub Seed Count");
    m_manager->registerHisto(m_hSubSeedVsCent);
  }  // end if not in pp mode
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetSeedCount::InitRun(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 1)
  {
    std::cout << "JetSeedCount::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetSeedCount::process_event(PHCompositeNode *topNode)
{
  // if needed, check if selected trigger fired
  if (m_doTrgSelect)
  {
    bool hasTrigger = JetQADefs::DidTriggerFire(m_trgToSelect, topNode);
    if (!hasTrigger)
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }

  // std::cout << "JetSeedCount::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  ++m_event;
  // Calling Raw Jet Seeds
  JetContainer *seedjetsraw = findNode::getClass<JetContainer>(topNode, m_rawSeedName);
  if (!seedjetsraw)
  {
    std::cout << "JetSeedCount::process_event - Error can not find DST raw seed jets" << std::endl;
    exit(-1);
  }

  // Calling Sub jet seeds
  JetContainer *seedjetssub = findNode::getClass<JetContainer>(topNode, m_subSeedName);
  if (!seedjetssub)
  {
    std::cout << "JetSeedCount::process_event - Error can not find DST sub seed jets" << std::endl;
    exit(-1);
  }

  // If not in pp mode, call Centrality Info
  CentralityInfo *cent_node = nullptr;
  if (!m_inPPMode)
  {
    cent_node = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
    if (!cent_node)
    {
      std::cout << "JetSeedCount::process_event - Error can not find CentralityInfo" << std::endl;
      exit(-1);
    }
  }

  // Z vertex info
  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexmap)
  {
    std::cout
        << "JetSeedCount::process_event - Error can not find global vertex  node "
        << std::endl;
    exit(-1);
  }
  if (vertexmap->empty())
  {
    if (Verbosity() > 1)
    {
      std::cout
          << "JetSeedCount::process_event - global vertex node is empty "
          << std::endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // If not in pp mode, save centrailty info
  float cent_mbd = std::numeric_limits<float>::max();
  if (!m_inPPMode)
  {
    cent_mbd = cent_node->get_centile(CentralityInfo::PROP::bimp);
    m_hCentMbd->Fill(cent_mbd);
  }

  // Raw Seed Count
  uint64_t n_seed_raw = 0;
  //  float Counter = 0;
  // for (JetMap::Iter iter = seedjetsraw->begin(); iter != seedjetsraw->end(); ++iter){
  for (auto jet : *seedjetsraw)
  {
    // Jet* jet = iter->second;
    int passesCut = jet->get_property(seedjetsraw->property_index(Jet::PROPERTY::prop_SeedItr));
    // Counter += 1;
    m_hRawPt_All->Fill(jet->get_pt());
    if (passesCut == 1)
    {
      m_hRawPt->Fill(jet->get_pt());
      m_hRawEtaVsPhi->Fill(jet->get_eta(), jet->get_phi());
      if (!m_inPPMode)
      {
        m_hRawSeedEnergyVsCent->Fill(cent_mbd, jet->get_e());
      }
      n_seed_raw++;
    }
  }
  m_hRawSeedCount->Fill(n_seed_raw);
  if (!m_inPPMode)
  {
    m_hRawSeedVsCent->Fill(cent_mbd, n_seed_raw);
  }

  // Sub Seed Count
  uint64_t n_seed_sub = 0;
  //  Counter = 0;
  // for (unsigned int iter = 0; iter < seedjetssub->size(); ++iter){
  // Jet* jet = seedjetsub->get_jet(iter);
  for (auto jet : *seedjetssub)
  {
    // Jet* jet = iter->second;
    int passesCut = jet->get_property(seedjetssub->property_index(Jet::PROPERTY::prop_SeedItr));
    //  Counter += 1;
    m_hSubPt_All->Fill(jet->get_pt());
    if (passesCut == 2)
    {
      m_hSubPt->Fill(jet->get_pt());
      m_hSubEtaVsPhi->Fill(jet->get_eta(), jet->get_phi());
      if (!m_inPPMode)
      {
        m_hSubSeedEnergyVsCent->Fill(cent_mbd, jet->get_e());
      }
      n_seed_sub++;
    }
  }
  m_hSubSeedCount->Fill(n_seed_sub);
  if (!m_inPPMode)
  {
    m_hSubSeedVsCent->Fill(cent_mbd, n_seed_sub);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetSeedCount::End(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 1)
  {
    std::cout << "JetSeedCount::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  }
  if (m_writeToOutputFile)
  {
    PHTFileServer::get().cd(m_outputFileName);
    m_hRawSeedCount->Write();
    m_hRawPt->Write();
    m_hRawPt_All->Write();
    m_hRawEtaVsPhi->Write();
    m_hSubSeedCount->Write();
    m_hSubPt->Write();
    m_hSubPt_All->Write();
    m_hSubEtaVsPhi->Write();
    if (!m_inPPMode)
    {
      m_hRawSeedEnergyVsCent->Write();
      m_hSubSeedEnergyVsCent->Write();
      m_hCentMbd->Write();
      m_hRawSeedVsCent->Write();
      m_hSubSeedVsCent->Write();
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void JetSeedCount::Print(const std::string &what) const
{
  std::cout << "JetSeedCount::Print(const std::string &what) const Printing info for " << what << std::endl;
}
