#include "StructureinJets.h"

#include <calotrigger/TriggerAnalyzer.h>

#include <centrality/CentralityInfo.h>

#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>
#include <jetbase/JetInput.h>

#include <qautils/QAHistManagerDef.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeed.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TVector3.h>
#include <TStyle.h>//for gStyle 

#include <boost/format.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <utility>

//____________________________________________________________________________..
StructureinJets::StructureinJets(const std::string& moduleName, const std::string& recoJetName, const std::string& trkNodeName, const std::string& histTag, const std::string& outputfilename)
  : SubsysReco(moduleName)
  , m_moduleName(moduleName)
  , m_recoJetName(recoJetName)
  , m_trkNodeName(trkNodeName)
  , m_histTag(histTag)
  , m_outputFileName(outputfilename)
{
  if(Verbosity() > 1 )
  {
    std::cout << "StructureinJets::StructureinJets(const std::string& x 5) Calling ctor" << std::endl;
  }
}

//____________________________________________________________________________..
StructureinJets::~StructureinJets()
{
  if (Verbosity() > 1)
  {
    std::cout << "StructureinJets::~StructureinJets() Calling dtor" << std::endl;
  }
}

//____________________________________________________________________________..
int StructureinJets::Init(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() > 0)
  {
    std::cout << "StructureinJets::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  }

  if (m_writeToOutputFileFlag)
  {
    PHTFileServer::open(m_outputFileName, "RECREATE");
  }

  delete m_analyzer;
  m_analyzer = new TriggerAnalyzer();
  m_manager = QAHistManagerDef::getHistoManager();

  // make sure module name is lower case
  std::string smallModuleName = m_moduleName;
  std::transform(
      smallModuleName.begin(),
      smallModuleName.end(),
      smallModuleName.begin(),
      ::tolower);

  // construct histogram names
  std::vector<std::string> vecHistNames = {
      "sumtrkvsjetptvscent",
      "sumtrkvsjetpt",
      "sumtrkoverjtpt",
      "sumtrkpt"};
  for (auto& vecHistName : vecHistNames)
  {
    vecHistName.insert(0, "h_" + smallModuleName + "_");
    if (!m_histTag.empty())
    {
      vecHistName.append("_" + m_histTag);
    }
  }

  m_hSumTrkVsJetPtVsCent = new TH3F(vecHistNames[0].data(), "", 100, 0, 100, 200, 0, 100, 10, 0, 100);
  m_hSumTrkVsJetPtVsCent->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  m_hSumTrkVsJetPtVsCent->GetYaxis()->SetTitle("Sum track p_{T} [GeV]");
  m_hSumTrkVsJetPt = new TH2F(vecHistNames[1].data(), "", 100, 0, 100, 200, 0, 100);
  m_hSumTrkVsJetPt->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  m_hSumTrkVsJetPt->GetYaxis()->SetTitle("Sum track p_{T} [GeV]");
  m_hSumTrkOverJetPt = new TH1F(vecHistNames[2].data(), "", 50, 0, 5);
  m_hSumTrkOverJetPt->GetXaxis()->SetTitle("Sum track p_{T} / Jet p_{T}");
  m_hSumTrkOverJetPt->GetYaxis()->SetTitle("Counts");
  m_hSumTrkPt = new TH1F(vecHistNames[3].data(), "", 200, 0, 100);
  m_hSumTrkPt->GetXaxis()->SetTitle("Sum track p_{T}");
  m_hSumTrkPt->GetYaxis()->SetTitle("Counts");

  if (m_isAAFlag)
  {
    m_manager->registerHisto(m_hSumTrkVsJetPtVsCent);
  }
  m_manager->registerHisto(m_hSumTrkVsJetPt);
  m_manager->registerHisto(m_hSumTrkOverJetPt);
  m_manager->registerHisto(m_hSumTrkPt);
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int StructureinJets::InitRun(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() > 0)
  {
    std::cout << "StructureinJets::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int StructureinJets::process_event(PHCompositeNode* topNode)
{

  if (Verbosity() > 1)
  {
    std::cout << "StructureinJets::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  }

  // if needed, check if selected trigger fired
  if (m_doTrgSelect)
  {
    m_analyzer->decodeTriggers(topNode);
    bool hasTrigger = JetQADefs::DidTriggerFire(m_trgToSelect, m_analyzer);
    if (!hasTrigger)
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }

  // interface to reco jets
  JetContainer* jets = findNode::getClass<JetContainer>(topNode, m_recoJetName);
  if (!jets)
  {
    std::cout
        << "StructureInJets::process_event - Error can not find DST Reco JetContainer node "
        << m_recoJetName << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // get reco tracks
  SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trkNodeName);
  if (!trackmap)
  {
    // if specified node name not found, try looking for this node;
    // and if that's not found, exit
    trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
    if (!trackmap)
    {
      std::cout
          << "StructureinJets::process_event - Error can not find DST trackmap node SvtxTrackMap" << std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }

  // get event centrality
  int cent = -1;
  if (m_isAAFlag)
  {
    CentralityInfo* cent_node = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
    if (!cent_node)
    {
      std::cout
          << "StructureinJets::process_event - Error can not find centrality node "
          << std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }
    cent = cent_node->get_centile(CentralityInfo::PROP::mbd_NS);
  }

  // Loop through jets
  for (auto *jet : *jets)
  {
    if (!jet)
    {
      if (Verbosity() > 2)
      {
        std::cout << "WARNING!!! Jet not found" << std::endl;
      }
      continue;
    }

    // remove noise
    if (jet->get_pt() < m_ptJetRange.first)
    {
      continue;
    }

    // apply jet kinematic cuts
    const bool inJetEtaCut = (jet->get_eta() >= m_etaJetRange.first) && (jet->get_eta() <= m_etaJetRange.second);
    const bool inJetPtCut = (jet->get_pt() >= m_ptJetRange.first) && (jet->get_pt() <= m_ptJetRange.second);
    if (!inJetEtaCut || !inJetPtCut)
    {
      continue;
    }

    // sum up tracks in jet
    TVector3 sumtrk(0, 0, 0);
    for (auto& iter : *trackmap)
    {
      SvtxTrack* track = iter.second;
      float quality = track->get_quality();
      auto *silicon_seed = track->get_silicon_seed();
      auto *tpc_seed = track->get_tpc_seed();

      // get no. of clusters in silicon seed
      int nsiliconclusts = 0;
      if (silicon_seed)
      {
        nsiliconclusts = silicon_seed->size_cluster_keys();
      }

      // get no. of clusters in tpc seed
      int ntpcclusts = 0;
      if (tpc_seed)
      {
        ntpcclusts = tpc_seed->size_cluster_keys();
      }

      // do some basic quality selection on tracks
      const bool inTrkPtCut = (track->get_pt() >= m_trk_pt_cut);
      const bool inTrkQualCut = (quality <= m_trk_qual_cut);
      const bool inTrkNSilCut = (nsiliconclusts >= m_trk_nsil_cut);
      const bool inTrkNTPCCut = (ntpcclusts >= m_trk_ntpc_cut);
      if (!inTrkPtCut || !inTrkQualCut || !inTrkNSilCut || !inTrkNTPCCut)
      {
        continue;
      }

      // Calculate delta eta, delta phi, and dR
      TVector3 v(track->get_px(), track->get_py(), track->get_pz());
      double dEta = v.Eta() - jet->get_eta();
      double dPhi = v.Phi() - jet->get_phi();
      while (dPhi > M_PI)
      {
        dPhi -= 2 * M_PI;
      }
      while (dPhi < -M_PI)
      {
        dPhi += 2 * M_PI;
      }
      double dR = sqrt((dEta * dEta) + (dPhi * dPhi));

      // Check if track is within jet radius
      if (dR < m_jetRadius)
      {
        sumtrk += v;
      }
    }

    // Fill histogram for the current jet
    assert(m_hSumTrkVsJetPtVsCent);
    assert(m_hSumTrkVsJetPt);
    assert(m_hSumTrkOverJetPt);
    assert(m_hSumTrkPt);

    // Fill TH3 histogram for Au+Au collisions
    if (m_isAAFlag)
    {
      m_hSumTrkVsJetPtVsCent->Fill(jet->get_pt(), sumtrk.Perp(), cent);
    }

    // Always fill others
    m_hSumTrkVsJetPt->Fill(jet->get_pt(), sumtrk.Perp());
    m_hSumTrkOverJetPt->Fill(sumtrk.Perp()/jet->get_pt());
    m_hSumTrkPt->Fill(sumtrk.Perp());
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int StructureinJets::ResetEvent(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() > 1)
  {
    std::cout << "StructureinJets::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int StructureinJets::EndRun(const int runnumber)
{
  if (Verbosity() > 0)
  {
  std::cout << "StructureinJets::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int StructureinJets::End(PHCompositeNode* /*topNode*/)
{
  // if flag is true, write to output file
  // otherwise rely on histogram manager
  if (m_writeToOutputFileFlag)
  {
    if (Verbosity() > 1)
    {
      std::cout << "StructureinJets::End - Output to " << m_outputFileName << std::endl;
    }
    PHTFileServer::cd(m_outputFileName);
  }
  else
  {
    if (Verbosity() > 1)
    {
      std::cout << "StructureinJets::End - Output to histogram manager" << std::endl;
    }
  }

  if (m_isAAFlag)
  {
    TH2* h_proj;
    for (int i = 0; i < m_hSumTrkVsJetPtVsCent->GetNbinsZ(); i++)
    {
      // construct histogram name for projection
      std::string name = m_hSumTrkVsJetPtVsCent->GetName();
      name.append("_centbin" + std::to_string(i + 1));

      m_hSumTrkVsJetPtVsCent->GetZaxis()->SetRange(i + 1, i + 1);
      h_proj = (TH2*) m_hSumTrkVsJetPtVsCent->Project3D("yx");
      h_proj->SetName(name.data());
      if (m_writeToOutputFileFlag)
      {
        h_proj->Write();
      }
      else
      {
        m_manager->registerHisto(h_proj);
      }
    }
  }
  else
  {
    if (m_writeToOutputFileFlag)
    {
      m_hSumTrkVsJetPt->Write();  // if pp, do not project onto centrality bins
    }
  }
  if (Verbosity() > 0)
  {
    std::cout << "StructureinJets::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int StructureinJets::Reset(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() > 0)
  {
    std::cout << "StructureinJets::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void StructureinJets::Print(const std::string& what) const
{
  std::cout << "StructureinJets::Print(const std::string &what) const Printing info for " << what << std::endl;
}
