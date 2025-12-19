#include "QAKFParticle.h"
#include "QAKFParticleTrackPtAsymmetry.h"

#include <qautils/QAHistManagerDef.h>  // for getHistoManager

#include <kfparticle_sphenix/KFParticle_Container.h>
#include <kfparticle_sphenix/KFParticle_Tools.h>

#include <calotrigger/TriggerRunInfo.h>

#include <ffarawobjects/Gl1Packet.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <trackbase_historic/SvtxTrack.h>  // for SvtxTrack
#include <trackbase_historic/SvtxTrackMap.h>

#include <trackbase/TrkrDefs.h>  // for cluskey

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>

#include <KFParticle.h>

#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>
#include <HepMC/GenVertex.h>
#include <HepMC/IteratorRange.h>

#include <TH1.h>
#include <TH2.h>

#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Vector/ThreeVector.h>

#include <cassert>
#include <cstdlib>  // for abs, NULL
#include <iostream>
#include <map>
#include <string>
#include <utility>  // for pair
#include <vector>

KFParticle_Tools kfpTools;

namespace
{
  inline void fillBunchCrossingRanges(std::vector<std::pair<double, double>> &ranges)
  {
    ranges.clear();

    int start = -100;
    int stop = 800;
    int step = 100;

    for (int lower = start; lower < stop;)
    {
      int upper = lower + step;

      if (lower == 0)
      {
        lower += 1;
      }

      ranges.emplace_back(static_cast<double>(lower), static_cast<double>(upper));

      if (lower == 1)
      {
        lower = 0;
      }

      lower = upper;
    }
  }
}  // namespace

QAKFParticle::QAKFParticle(const std::string &name, const std::string &mother_name, double min_m, double max_m)
  : SubsysReco(name)
  , m_mother_id(kfpTools.getParticleID(mother_name))
  , m_min_mass(min_m)
  , m_max_mass(max_m)
{
}

int QAKFParticle::InitRun(PHCompositeNode *topNode)
{
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  if (!m_trackMap)
  {
    std::cout << __PRETTY_FUNCTION__ << " Fatal Error : missing " << m_trackMapName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAKFParticle::Init(PHCompositeNode * /*topNode*/)
{
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  fillBunchCrossingRanges(bunchCrossingRanges);

  std::cout << "Bunch crossing ranges for mass histograms: " << std::endl;
  if (Verbosity() > 0)
  {
    std::cout << "Total number of ranges: " << bunchCrossingRanges.size() << std::endl;
    for (const auto &range : bunchCrossingRanges)
    {
      std::cout << "  [" << range.first << ", " << range.second << ")" << std::endl;
    }
  }

  TH1 *h(nullptr);

  TH2 *h2(nullptr);

  h = new TH1F(TString(get_histo_prefix()) + "InvMass_KFP",  //
               ";mass [GeV];Entries", 40, m_min_mass, m_max_mass);
  hm->registerHisto(h);

  float eta_min = -1.3;
  float eta_max = 1.3;
  float phi_min = -3.14;
  float phi_max = 3.14;
  float pt_min = 0;
  float pt_max = 5;

  h2 = new TH2F(TString(get_histo_prefix()) + "InvMass_KFP_Eta",  //
                ";mass [GeV]; #eta", 100, m_min_mass, m_max_mass, 100, eta_min, eta_max);
  hm->registerHisto(h2);

  h2 = new TH2F(TString(get_histo_prefix()) + "InvMass_KFP_Phi",  //
                ";mass [GeV]; #phi [rad.]", 100, m_min_mass, m_max_mass, 100, phi_min, phi_max);
  hm->registerHisto(h2);

  h2 = new TH2F(TString(get_histo_prefix()) + "InvMass_KFP_pT",  //
                ";mass [GeV]; p_{T} [GeV]", 100, m_min_mass, m_max_mass, 100, pt_min, pt_max);
  hm->registerHisto(h2);

  // mass v.s bunch crossing
  h2 = new TH2F(TString(get_histo_prefix()) + "BunchCrossing_InvMass_KFP",                                         //
                ";Beam crossing (106 ns/crossing);mass [GeV];", 899, -149.5, 749.5, 100, m_min_mass, m_max_mass);  // the axis title is the same as the approved plot https://www.sphenix.bnl.gov/sites/default/files/2025-04/sphenix-perf-4-25-Kscrossingcut.pdf
  hm->registerHisto(h2);

  h = new TH1F(TString(get_histo_prefix()) + "InvMass_KFP_crossing0",  //
               ";mass [GeV];Entries", 100, m_min_mass, m_max_mass);
  hm->registerHisto(h);

  h = new TH1F(TString(get_histo_prefix()) + "InvMass_KFP_non_crossing0",  //
               ";mass [GeV];Entries", 100, m_min_mass, m_max_mass);
  hm->registerHisto(h);

  h = new TH1F(TString(get_histo_prefix()) + "InvMass_KFP_ZDC_Coincidence",  //
               ";mass [GeV];Entries", 100, m_min_mass, m_max_mass);
  hm->registerHisto(h);

  h = new TH1F(TString(get_histo_prefix()) + "InvMass_KFP_MBD_NandS_geq_1_vtx_l_10_cm",  //
               ";mass [GeV];Entries", 100, m_min_mass, m_max_mass);
  hm->registerHisto(h);

  h = new TH1F(TString(get_histo_prefix()) + "InvMass_KFP_Jet_6_GeV_MBD_NandS_geq_1_vtx_l_10_cm",  //
               ";mass [GeV];Entries", 100, m_min_mass, m_max_mass);
  hm->registerHisto(h);

  // create the histograms for different bunch crossing ranges
  for (const auto &range : bunchCrossingRanges)
  {
    std::string histname = get_histo_prefix() + "InvMass_KFP_crossing_" + std::to_string(static_cast<int>(range.first)) + "_to_" + std::to_string(static_cast<int>(range.second));
    h = new TH1F(histname.c_str(), ";mass [GeV];Entries", 100, m_min_mass, m_max_mass);
    hm->registerHisto(h);
  }

  h_mass_KFP = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "InvMass_KFP"));
  assert(h_mass_KFP);

  h_mass_KFP_eta = dynamic_cast<TH2F *>(hm->getHisto(get_histo_prefix() + "InvMass_KFP_Eta"));
  assert(h_mass_KFP_eta);

  h_mass_KFP_phi = dynamic_cast<TH2F *>(hm->getHisto(get_histo_prefix() + "InvMass_KFP_Phi"));
  assert(h_mass_KFP_phi);

  h_mass_KFP_pt = dynamic_cast<TH2F *>(hm->getHisto(get_histo_prefix() + "InvMass_KFP_pT"));
  assert(h_mass_KFP_pt);

  h_bunchcrossing_mass_KFP = dynamic_cast<TH2F *>(hm->getHisto(get_histo_prefix() + "BunchCrossing_InvMass_KFP"));
  assert(h_bunchcrossing_mass_KFP);

  h_mass_KFP_crossing0 = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "InvMass_KFP_crossing0"));
  assert(h_mass_KFP_crossing0);

  h_mass_KFP_non_crossing0 = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "InvMass_KFP_non_crossing0"));
  assert(h_mass_KFP_non_crossing0);

  h_mass_KFP_ZDC_Coincidence = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "InvMass_KFP_ZDC_Coincidence"));
  assert(h_mass_KFP_ZDC_Coincidence);

  h_mass_KFP_MBD_NandS_geq_1_vtx_l_10_cm = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "InvMass_KFP_MBD_NandS_geq_1_vtx_l_10_cm"));
  assert(h_mass_KFP_MBD_NandS_geq_1_vtx_l_10_cm);

  h_mass_KFP_Jet_6_GeV_MBD_NandS_geq_1_vtx_l_10_cm = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "InvMass_KFP_Jet_6_GeV_MBD_NandS_geq_1_vtx_l_10_cm"));
  assert(h_mass_KFP_Jet_6_GeV_MBD_NandS_geq_1_vtx_l_10_cm);

  // histograms for different bunch crossing ranges
for (const auto &range : bunchCrossingRanges)
  {
    h_mass_KFP_crossingrange.push_back(dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "InvMass_KFP_crossing_" + std::to_string(static_cast<int>(range.first)) + "_to_" + std::to_string(static_cast<int>(range.second)))));
    assert(h_mass_KFP_crossingrange.back());
  }

  // pt asymmetry analyzer
  if (m_doTrackPtAsymmetry)
  {
    m_trackPtAsymmetryAnalyzer = std::make_unique<QAKFParticleTrackPtAsymmetry>(get_histo_prefix(),     //
                                                                                m_min_mass,             //
                                                                                m_max_mass,             //
                                                                                m_trackPtAsymEtaBins,   //
                                                                                m_trackPtAsymPhiBins);  //
    m_trackPtAsymmetryAnalyzer->bookHistograms(hm);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAKFParticle::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 2)
  {
    std::cout << "QAKFParticle::process_event() entered" << std::endl;
  }

  // load relevant nodes from NodeTree
  load_nodes(topNode);

  if (counter == 0)
  {
    initializeTriggerInfo(topNode);
  }
  // histogram manager

  if (hasTriggerInfo)
  {
    triggeranalyzer->decodeTriggers(topNode);
  }

  for (auto &iter : *m_kfpContainer)
  {
    if (iter.second->GetPDG() == m_mother_id)
    {
      // filling mother histogram information
      float eta = 0;
      float etaErr = 0;
      float phi = 0;
      float phiErr = 0;

      iter.second->GetEta(eta, etaErr);
      iter.second->GetPhi(phi, phiErr);

      h_mass_KFP->Fill(iter.second->GetMass());

      h_mass_KFP_eta->Fill(iter.second->GetMass(), eta);

      h_mass_KFP_phi->Fill(iter.second->GetMass(), phi);

      h_mass_KFP_pt->Fill(iter.second->GetMass(), iter.second->GetPt());

      const std::vector<int> track_ids = iter.second->DaughterIds();
      SvtxTrack *kfpTrack = nullptr;
      for (auto &iter2 : *m_trackMap)
      {
        if (iter2.first == (unsigned int) track_ids[0])
        {
          kfpTrack = iter2.second;
          break;
        }
      }

      if (kfpTrack)
      {
        if (kfpTrack->get_crossing() == 0)
        {
          h_mass_KFP_crossing0->Fill(iter.second->GetMass());
        }
        else
        {
          h_mass_KFP_non_crossing0->Fill(iter.second->GetMass());
        }

        h_bunchcrossing_mass_KFP->Fill(kfpTrack->get_crossing(), iter.second->GetMass());

        // fill the appropriate crossing range histogram
        size_t range_index = 0;
        for (const auto &range : bunchCrossingRanges)
        {
          if (kfpTrack->get_crossing() >= range.first && kfpTrack->get_crossing() < range.second)
          {
            h_mass_KFP_crossingrange[range_index]->Fill(iter.second->GetMass());
            break;
          }
          ++range_index;
        }
      }

      if (hasTriggerInfo)
      {
        for (int i = 0; i < nTriggerBits; ++i)
        {
          if (m_ZDC_Coincidence_bit == i && triggeranalyzer->didTriggerFire(i) == 1)
          {
            h_mass_KFP_ZDC_Coincidence->Fill(iter.second->GetMass());
          }
          else if (m_MBD_NandS_geq_1_vtx_l_10_cm_bit == i && triggeranalyzer->didTriggerFire(i) == 1)
          {
            h_mass_KFP_MBD_NandS_geq_1_vtx_l_10_cm->Fill(iter.second->GetMass());
          }
          else if (m_Jet_6_GeV_MBD_NandS_geq_1_vtx_l_10_cm_bit == i && triggeranalyzer->didTriggerFire(i) == 1)
          {
            h_mass_KFP_Jet_6_GeV_MBD_NandS_geq_1_vtx_l_10_cm->Fill(iter.second->GetMass());
          }
        }
      }

      // hook up the pt asymmetry analysis (can be used for any particle with 2-body decay)
      if (m_doTrackPtAsymmetry)
      {
        m_trackPtAsymmetryAnalyzer->setVerbosity(Verbosity());
        m_trackPtAsymmetryAnalyzer->analyzeTrackPtAsymmetry(m_trackMap, iter.second);
      }
    }
  }

  ++counter;

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAKFParticle::load_nodes(PHCompositeNode *topNode)
{
  std::string kfpContainerNodeName = m_KFParticleNodeName + "_KFParticle_Container";
  m_kfpContainer = findNode::getClass<KFParticle_Container>(topNode, kfpContainerNodeName);
  if (!m_kfpContainer)
  {
    std::cout << kfpContainerNodeName << " - Fatal Error - "
              << "unable to find DST node "
              << "KFParticle_QA" << std::endl;
    assert(m_kfpContainer);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

std::string QAKFParticle::get_histo_prefix() { return std::string("h_") + Name() + std::string("_") + m_trackMapName + std::string("_"); }

void QAKFParticle::initializeTriggerInfo(PHCompositeNode *topNode)
{
  triggeranalyzer = new TriggerAnalyzer();

  // Check whether we actually have the right information
  auto *gl1packet = findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");
  if (!gl1packet)
  {
    // std::cout << "No GL1RAWHIT" << std::endl;
    gl1packet = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    if (!gl1packet)
    {
      // std::cout << "No GL1Packet" << std::endl;
      return;
    }
  }

  auto *triggerruninfo = findNode::getClass<TriggerRunInfo>(topNode, "TriggerRunInfo");
  if (!triggerruninfo)
  {
    hasTriggerInfo = false;
    // std::cout << "No triggerRunInfo" << std::endl;
    return;
  }

  size_t pos;
  std::string undrscr = "_";
  std::string nothing;
  std::map<std::string, std::string> forbiddenStrings;
  forbiddenStrings[" "] = undrscr;
  forbiddenStrings[","] = nothing;
  forbiddenStrings["/"] = undrscr;
  forbiddenStrings["&"] = "and";
  forbiddenStrings["="] = "eq";
  forbiddenStrings["<"] = "l";
  forbiddenStrings[">"] = "g";
  forbiddenStrings["+"] = "plus";
  forbiddenStrings["-"] = "minus";
  forbiddenStrings["*"] = "star";

  triggeranalyzer->decodeTriggers(topNode);

  for (int i = 0; i < nTriggerBits; ++i)
  {
    std::string triggerName = triggeranalyzer->getTriggerName(i);

    if (triggerName.find("unknown") != std::string::npos)
    {
      continue;
    }

    for (auto const &[badString, goodString] : forbiddenStrings)
    {
      while ((pos = triggerName.find(badString)) != std::string::npos)
      {
        triggerName.replace(pos, 1, goodString);
      }
    }

    if (triggerName == "ZDC_Coincidence")
    {
      m_ZDC_Coincidence_bit = i;
    }
    else if (triggerName == "MBD_NandS_geq_1_vtx_l_10_cm")
    {
      m_MBD_NandS_geq_1_vtx_l_10_cm_bit = i;
    }
    else if (triggerName == "Jet_6_GeV_MBD_NandS_geq_1_vtx_l_10_cm")
    {
      m_Jet_6_GeV_MBD_NandS_geq_1_vtx_l_10_cm_bit = i;
    }
  }
}
