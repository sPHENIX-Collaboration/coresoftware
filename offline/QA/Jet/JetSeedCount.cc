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

  GlobalVertex *vtx = vertexmap->begin()->second;
  z_vtx = vtx->get_z();

  // If not in pp mode, save centrailty info
  float cent_mbd = std::numeric_limits<float>::max();
  if (!m_inPPMode)
  {
    cent_mbd = cent_node->get_centile(CentralityInfo::PROP::bimp);
    m_centrality.push_back(cent_mbd);
  }

  // Raw Seed Count
  m_seed_raw = 0;
  //  float Counter = 0;
  // for (JetMap::Iter iter = seedjetsraw->begin(); iter != seedjetsraw->end(); ++iter){
  for (auto jet : *seedjetsraw)
  {
    // Jet* jet = iter->second;
    int passesCut = jet->get_property(seedjetsraw->property_index(Jet::PROPERTY::prop_SeedItr));
    // Counter += 1;
    m_rawpt_all.push_back(jet->get_pt());
    if (passesCut == 1)
    {
      m_rawpt.push_back(jet->get_pt());
      m_rawenergy.push_back(jet->get_e());
      m_RawEta.push_back(jet->get_eta());
      m_RawPhi.push_back(jet->get_phi());
      if (!m_inPPMode)
      {
        m_rawcent.push_back(cent_mbd);
      }
      m_seed_raw++;
    }
  }
  m_raw_counts.push_back(m_seed_raw);

  // Sub Seed Count
  m_seed_sub = 0;
  //  Counter = 0;
  // for (unsigned int iter = 0; iter < seedjetssub->size(); ++iter){
  // Jet* jet = seedjetsub->get_jet(iter);
  for (auto jet : *seedjetssub)
  {
    // Jet* jet = iter->second;
    int passesCut = jet->get_property(seedjetssub->property_index(Jet::PROPERTY::prop_SeedItr));
    //  Counter += 1;
    m_subpt_all.push_back(jet->get_pt());
    if (passesCut == 2)
    {
      m_subpt.push_back(jet->get_pt());
      m_subenergy.push_back(jet->get_e());
      m_SubEta.push_back(jet->get_eta());
      m_SubPhi.push_back(jet->get_phi());
      if (!m_inPPMode)
      {
        m_subcent.push_back(cent_mbd);
      }
      m_seed_sub++;
    }
  }
  m_sub_counts.push_back(m_seed_sub);
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

  TH1 *hRawSeedCount = new TH1F(vecHistNames[0].data(), "Raw Seed Count per Event", 100, 0.00, 50.00);
  hRawSeedCount->GetXaxis()->SetTitle("Raw Seed Count per Event");
  hRawSeedCount->GetYaxis()->SetTitle("Number of Entries");
  for (int m_raw_count : m_raw_counts)
  {
    hRawSeedCount->Fill(m_raw_count);
  }
  m_manager->registerHisto(hRawSeedCount);

  TH1 *hRawPt = new TH1F(vecHistNames[1].data(), "Raw p_{T}", 1000, 0.00, 50.00);
  hRawPt->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  hRawPt->GetYaxis()->SetTitle("Number of Entries");
  for (double j : m_rawpt)
  {
    hRawPt->Fill(j);
  }
  m_manager->registerHisto(hRawPt);

  TH1 *hRawPt_All = new TH1F(vecHistNames[2].data(), "Raw p_{T} (all jet seeds)", 1000, 0.00, 50.00);
  hRawPt_All->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  hRawPt_All->GetYaxis()->SetTitle("Number of Entries");
  for (double j : m_rawpt_all)
  {
    hRawPt_All->Fill(j);
  }
  m_manager->registerHisto(hRawPt_All);

  TH2 *hRawEtaVsPhi = new TH2F(vecHistNames[3].data(), "Raw Seed Eta Vs Phi", 220, -1.1, 1.1, 628, -3.14, 3.14);
  hRawEtaVsPhi->GetXaxis()->SetTitle("Jet #eta [Rads.]");
  hRawEtaVsPhi->GetYaxis()->SetTitle("Jet #phi [Rads.]");
  for (int j = 0; j < (int) m_RawEta.size(); j++)
  {
    hRawEtaVsPhi->Fill(m_RawEta.at(j), m_RawPhi.at(j));
  }
  m_manager->registerHisto(hRawEtaVsPhi);

  TH1 *hSubSeedCount = new TH1F(vecHistNames[4].data(), "Sub Seed Count per Event", 100, 0.00, 50.00);
  hSubSeedCount->GetXaxis()->SetTitle("Sub Seed Count per Event");
  hSubSeedCount->GetYaxis()->SetTitle("Number of Entries");
  for (int m_sub_count : m_sub_counts)
  {
    hSubSeedCount->Fill(m_sub_count);
  }
  m_manager->registerHisto(hSubSeedCount);

  TH1 *hSubPt = new TH1F(vecHistNames[5].data(), "Sub. p_{T}", 1000, 0.00, 50.00);
  hSubPt->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  hSubPt->GetYaxis()->SetTitle("Number of Entries");
  for (double j : m_subpt)
  {
    hSubPt->Fill(j);
  }
  m_manager->registerHisto(hSubPt);

  TH1 *hSubPt_All = new TH1F(vecHistNames[6].data(), "Sub. p_{T} (all jet seeds)", 1000, 0.00, 50.00);
  hSubPt_All->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  hSubPt_All->GetYaxis()->SetTitle("Number of Entries");
  for (double j : m_subpt_all)
  {
    hSubPt_All->Fill(j);
  }
  m_manager->registerHisto(hSubPt_All);

  TH2 *hSubEtaVsPhi = new TH2F(vecHistNames[7].data(), "Sub. Seed Eta Vs Phi", 220, -1.1, 1.1, 628, -3.14, 3.14);
  hSubEtaVsPhi->GetXaxis()->SetTitle("Jet #eta [Rads.]");
  hSubEtaVsPhi->GetYaxis()->SetTitle("Jet #phi [Rads.]");
  for (int j = 0; j < (int) m_SubEta.size(); j++)
  {
    hSubEtaVsPhi->Fill(m_SubEta.at(j), m_SubPhi.at(j));
  }
  m_manager->registerHisto(hSubEtaVsPhi);

  // If not in pp mode, plot quantities vs. centrality
  if (!m_inPPMode)
  {
    TH2 *hRawSeedEnergyVsCent = new TH2F(vecHistNames[8].data(), "Raw Seed Energy Vs Centrality", 10.00, 0.00, 100.00, 100, 0.00, 50.00);
    hRawSeedEnergyVsCent->GetXaxis()->SetTitle("Centrality");
    hRawSeedEnergyVsCent->GetYaxis()->SetTitle("RawSeedEnergy");
    for (int j = 0; j < (int) m_rawenergy.size(); j++)
    {
      hRawSeedEnergyVsCent->Fill(m_rawcent.at(j), m_rawenergy.at(j));
    }
    m_manager->registerHisto(hRawSeedEnergyVsCent);

    TH2 *hSubSeedEnergyVsCent = new TH2F(vecHistNames[9].data(), "Sub Seed Energy Vs Centrality", 10.00, 0.00, 100.00, 100, 0.00, 50.00);
    hSubSeedEnergyVsCent->GetXaxis()->SetTitle("Centrality");
    hSubSeedEnergyVsCent->GetYaxis()->SetTitle("SubSeedEnergy");
    for (int j = 0; j < (int) m_subenergy.size(); j++)
    {
      hSubSeedEnergyVsCent->Fill(m_subcent.at(j), m_subenergy.at(j));
    }
    m_manager->registerHisto(hSubSeedEnergyVsCent);

    TH1 *hCentMbd = new TH1F(vecHistNames[10].data(), "hCentMbd", 10, 0.00, 100.00);
    hCentMbd->GetXaxis()->SetTitle("Centrality (Mbd)");
    hCentMbd->GetYaxis()->SetTitle("Number of Entries");
    for (int j : m_centrality)
    {
      hCentMbd->Fill(j);
    }
    m_manager->registerHisto(hCentMbd);

    TH2 *hRawSeedVsCent = new TH2F(vecHistNames[10].data(), "Raw Seed Vs Centrality", 10, 0.00, 100.00, 101, -0.5, 100.5);
    hRawSeedVsCent->GetXaxis()->SetTitle("Centrality");
    hRawSeedVsCent->GetYaxis()->SetTitle("Raw Seed Count");
    for (int j = 0; j < (int) m_raw_counts.size(); j++)
    {
      hRawSeedVsCent->Fill(m_centrality.at(j), m_raw_counts.at(j));
    }
    m_manager->registerHisto(hRawSeedVsCent);

    TH2 *hSubSeedVsCent = new TH2F(vecHistNames[12].data(), "Sub Seed Vs Centrality", 10, 0.00, 100.00, 101, -0.5, 100.5);
    hSubSeedVsCent->GetXaxis()->SetTitle("Centrality");
    hSubSeedVsCent->GetYaxis()->SetTitle("Sub Seed Count");
    for (int j = 0; j < (int) m_sub_counts.size(); j++)
    {
      hSubSeedVsCent->Fill(m_centrality.at(j), m_sub_counts.at(j));
    }
    m_manager->registerHisto(hSubSeedVsCent);
  }  // end if not in pp mode
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void JetSeedCount::Print(const std::string &what) const
{
  std::cout << "JetSeedCount::Print(const std::string &what) const Printing info for " << what << std::endl;
}
