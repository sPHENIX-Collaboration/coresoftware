//____________________________________________________________________________..

#include "JetKinematicCheck.h"

#include <calotrigger/TriggerAnalyzer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <jetbase/JetContainer.h>
#include <jetbase/Jetv2.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLegend.h>
#include <TPad.h>
#include <TStyle.h>

#include <boost/format.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
//____________________________________________________________________________..

JetKinematicCheck::JetKinematicCheck(const std::string &moduleName,
                                     const std::string &recojetnameR02,
                                     const std::string &recojetnameR03,
                                     const std::string &recojetnameR04,
                                     const std::string &recojetnameR05)
  : SubsysReco(moduleName)
  , m_moduleName(moduleName)
  , m_recoJetNameR02(recojetnameR02)
  , m_recoJetNameR03(recojetnameR03)
  , m_recoJetNameR04(recojetnameR04)
  , m_recoJetNameR05(recojetnameR05)
{
  if (Verbosity() > 1)
  {
    std::cout << "JetKinematicCheck::JetKinematicCheck(const std::string &name) Calling ctor" << std::endl;
  }
}

//____________________________________________________________________________..
JetKinematicCheck::~JetKinematicCheck()
{
  if (Verbosity() > 1)
  {
    std::cout << "JetKinematicCheck::~JetKinematicCheck() Calling dtor" << std::endl;
  }
  delete m_analyzer;
}

//____________________________________________________________________________..
int JetKinematicCheck::Init(PHCompositeNode * /*unused*/)
{
  delete m_analyzer;
  m_analyzer = new TriggerAnalyzer();
  hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // make sure module name is lower case
  std::string smallModuleName = m_moduleName;
  std::transform(
      smallModuleName.begin(),
      smallModuleName.end(),
      smallModuleName.begin(),
      ::tolower);

  // construct histogram names
  std::vector<std::string> vecHistNames = {
      "r02_spectra",
      "r03_spectra",
      "r04_spectra",
      "r05_spectra",
      "r02_etavsphi",
      "r03_etavsphi",
      "r04_etavsphi",
      "r05_etavsphi",
      "r02_jetmassvspt",
      "r03_jetmassvspt",
      "r04_jetmassvspt",
      "r05_jetmassvspt",
      "r02_jetmassvseta",
      "r03_jetmassvseta",
      "r04_jetmassvseta",
      "r05_jetmassvseta",
      "r02_jetmassvsptprofile",
      "r03_jetmassvsptprofile",
      "r04_jetmassvsptprofile",
      "r05_jetmassvsptprofile",
      "r02_jetmassvsetaprofile",
      "r03_jetmassvsetaprofile",
      "r04_jetmassvsetaprofile",
      "r05_jetmassvsetaprofile"};
  for (auto &vecHistName : vecHistNames)
  {
    vecHistName.insert(0, "h_" + smallModuleName + "_");
    if (!m_histTag.empty())
    {
      vecHistName.append("_" + m_histTag);
    }
  }

  // initialize histograms

  h_jet_pt_spectra_r02 = new TH1D(vecHistNames[0].data(), "", 19, 5, 100);
  h_jet_pt_spectra_r02->GetXaxis()->SetTitle("p_{T} [GeV/c]");

  h_jet_pt_spectra_r03 = new TH1D(vecHistNames[1].data(), "", 19, 5, 100);
  h_jet_pt_spectra_r03->GetXaxis()->SetTitle("p_{T} [GeV/c]");

  h_jet_pt_spectra_r04 = new TH1D(vecHistNames[2].data(), "", 19, 5, 100);
  h_jet_pt_spectra_r04->GetXaxis()->SetTitle("p_{T} [GeV/c]");

  h_jet_pt_spectra_r05 = new TH1D(vecHistNames[3].data(), "", 19, 5, 100);
  h_jet_pt_spectra_r05->GetXaxis()->SetTitle("p_{T} [GeV/c]");

  h_jet_eta_phi_r02 = new TH2D(vecHistNames[4].data(), "", 24, -1.1, 1.1, 64, -M_PI, M_PI);
  h_jet_eta_phi_r02->GetXaxis()->SetTitle("#eta");
  h_jet_eta_phi_r02->GetYaxis()->SetTitle("#Phi");

  h_jet_eta_phi_r03 = new TH2D(vecHistNames[5].data(), "", 24, -1.1, 1.1, 64, -M_PI, M_PI);
  h_jet_eta_phi_r03->GetXaxis()->SetTitle("#eta");
  h_jet_eta_phi_r03->GetYaxis()->SetTitle("#Phi");

  h_jet_eta_phi_r04 = new TH2D(vecHistNames[6].data(), "", 24, -1.1, 1.1, 64, -M_PI, M_PI);
  h_jet_eta_phi_r04->GetXaxis()->SetTitle("#eta");
  h_jet_eta_phi_r04->GetYaxis()->SetTitle("#Phi");

  h_jet_eta_phi_r05 = new TH2D(vecHistNames[7].data(), "", 24, -1.1, 1.1, 64, -M_PI, M_PI);
  h_jet_eta_phi_r05->GetXaxis()->SetTitle("#eta");
  h_jet_eta_phi_r05->GetYaxis()->SetTitle("#Phi");

  h_jet_mass_pt_r02 = new TH2D(vecHistNames[8].data(), "", 19, 5, 100, 15, 0, 15);
  h_jet_mass_pt_r02->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  h_jet_mass_pt_r02->GetYaxis()->SetTitle("Jet Mass [GeV/c^{2}]");

  h_jet_mass_pt_r03 = new TH2D(vecHistNames[9].data(), "", 19, 5, 100, 15, 0, 15);
  h_jet_mass_pt_r03->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  h_jet_mass_pt_r03->GetYaxis()->SetTitle("Jet Mass [GeV/c^{2}]");

  h_jet_mass_pt_r04 = new TH2D(vecHistNames[10].data(), "", 19, 5, 100, 15, 0, 15);
  h_jet_mass_pt_r04->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  h_jet_mass_pt_r04->GetYaxis()->SetTitle("Jet Mass [GeV/c^{2}]");

  h_jet_mass_pt_r05 = new TH2D(vecHistNames[11].data(), "", 19, 5, 100, 15, 0, 15);
  h_jet_mass_pt_r05->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  h_jet_mass_pt_r05->GetYaxis()->SetTitle("Jet Mass [GeV/c^{2}]");

  h_jet_mass_eta_r02 = new TH2D(vecHistNames[12].data(), "", 24, -1.1, 1.1, 15, 0, 15);
  h_jet_mass_eta_r02->GetXaxis()->SetTitle("#eta");
  h_jet_mass_eta_r02->GetYaxis()->SetTitle("Jet Mass [GeV/c^{2}]");

  h_jet_mass_eta_r03 = new TH2D(vecHistNames[13].data(), "", 24, -1.1, 1.1, 15, 0, 15);
  h_jet_mass_eta_r03->GetXaxis()->SetTitle("#eta");
  h_jet_mass_eta_r03->GetYaxis()->SetTitle("Jet Mass [GeV/c^{2}]");

  h_jet_mass_eta_r04 = new TH2D(vecHistNames[14].data(), "", 24, -1.1, 1.1, 15, 0, 15);
  h_jet_mass_eta_r04->GetXaxis()->SetTitle("#eta");
  h_jet_mass_eta_r04->GetYaxis()->SetTitle("Jet Mass [GeV/c^{2}]");

  h_jet_mass_eta_r05 = new TH2D(vecHistNames[15].data(), "", 24, -1.1, 1.1, 15, 0, 15);
  h_jet_mass_eta_r05->GetXaxis()->SetTitle("#eta");
  h_jet_mass_eta_r05->GetYaxis()->SetTitle("Jet Mass [GeV/c^{2}]");

  h_jet_average_mass_pt_1D_r02 = new TH1D(vecHistNames[16].data(), "", 19, 5, 100);
  h_jet_average_mass_pt_1D_r02->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  h_jet_average_mass_pt_1D_r02->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");

  h_jet_average_mass_pt_1D_r03 = new TH1D(vecHistNames[17].data(), "", 19, 5, 100);
  h_jet_average_mass_pt_1D_r03->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  h_jet_average_mass_pt_1D_r03->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");

  h_jet_average_mass_pt_1D_r04 = new TH1D(vecHistNames[18].data(), "", 19, 5, 100);
  h_jet_average_mass_pt_1D_r04->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  h_jet_average_mass_pt_1D_r04->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");

  h_jet_average_mass_pt_1D_r05 = new TH1D(vecHistNames[19].data(), "", 19, 5, 100);
  h_jet_average_mass_pt_1D_r05->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  h_jet_average_mass_pt_1D_r05->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");

  h_jet_average_mass_eta_1D_r02 = new TH2D(vecHistNames[20].data(), "", 24, -1.1, 1.1, 15, 0, 15);
  h_jet_average_mass_eta_1D_r02->GetXaxis()->SetTitle("#eta");
  h_jet_average_mass_eta_1D_r02->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");

  h_jet_average_mass_eta_1D_r03 = new TH2D(vecHistNames[21].data(), "", 24, -1.1, 1.1, 15, 0, 15);
  h_jet_average_mass_eta_1D_r03->GetXaxis()->SetTitle("#eta");
  h_jet_average_mass_eta_1D_r03->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");

  h_jet_average_mass_eta_1D_r04 = new TH2D(vecHistNames[22].data(), "", 24, -1.1, 1.1, 15, 0, 15);
  h_jet_average_mass_eta_1D_r04->GetXaxis()->SetTitle("#eta");
  h_jet_average_mass_eta_1D_r04->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");

  h_jet_average_mass_eta_1D_r05 = new TH2D(vecHistNames[23].data(), "", 24, -1.1, 1.1, 15, 0, 15);
  h_jet_average_mass_eta_1D_r05->GetXaxis()->SetTitle("#eta");
  h_jet_average_mass_eta_1D_r05->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");

  if (Verbosity() > 1)
  {
    std::cout << "JetKinematicCheck::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetKinematicCheck::InitRun(PHCompositeNode * /*unused*/)
{
  if (Verbosity() > 1)
  {
    std::cout << "JetKinematicCheck::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetKinematicCheck::process_event(PHCompositeNode *topNode)
{
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

  std::vector<std::string> m_recoJetName_array = {m_recoJetNameR02, m_recoJetNameR03, m_recoJetNameR04, m_recoJetNameR05};
  m_radii = {0.2, 0.3, 0.4, 0.5};
  int n_radii = m_radii.size();

  // Loop over each reco jet radii from array
  for (int i = 0; i < n_radii; i++)
  {

    // update eta range based on resolution parameter
    std::pair<float, float> etaRangeUse;
    if (m_restrictEtaRange)
    {
      etaRangeUse = {m_etaRange.first + m_radii[i], m_etaRange.second - m_radii[i]};
    }
    else
    {
      etaRangeUse = {m_etaRange.first, m_etaRange.second};
    }

    const std::string& recoJetName = m_recoJetName_array[i];

    JetContainer *jets = findNode::getClass<JetContainer>(topNode, recoJetName);
    if (!jets)
    {
      std::cout
          << "JetKinematicCheck::process_event - Error can not find DST Reco JetContainer node "
          << recoJetName << std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }

    // loop over jets
    for (auto *jet : *jets)
    {
      bool eta_cut = (jet->get_eta() >= etaRangeUse.first) && (jet->get_eta() <= etaRangeUse.second);
      bool pt_cut = (jet->get_pt() >= m_ptRange.first) && (jet->get_pt() <= m_ptRange.second);
      if ((! eta_cut) || (! pt_cut))
      {
        continue;
      }
      if (jet->get_pt() < m_ptRange.first)
      {
        continue;  // to remove noise jets
      }

      if (i == 0)
      {
        h_jet_pt_spectra_r02->Fill(jet->get_pt());
        h_jet_eta_phi_r02->Fill(jet->get_eta(), jet->get_phi());
        h_jet_mass_pt_r02->Fill(jet->get_pt(), jet->get_mass());
        h_jet_mass_eta_r02->Fill(jet->get_eta(), jet->get_mass());
      }

      else if (i == 1)
      {
        h_jet_pt_spectra_r03->Fill(jet->get_pt());
        h_jet_eta_phi_r03->Fill(jet->get_eta(), jet->get_phi());
        h_jet_mass_pt_r03->Fill(jet->get_pt(), jet->get_mass());
        h_jet_mass_eta_r03->Fill(jet->get_eta(), jet->get_mass());
      }

      else if (i == 2)
      {
        h_jet_pt_spectra_r04->Fill(jet->get_pt());
        h_jet_eta_phi_r04->Fill(jet->get_eta(), jet->get_phi());
        h_jet_mass_pt_r04->Fill(jet->get_pt(), jet->get_mass());
        h_jet_mass_eta_r04->Fill(jet->get_eta(), jet->get_mass());
      }

      else if (i == 3)
      {
        h_jet_pt_spectra_r05->Fill(jet->get_pt());
        h_jet_eta_phi_r05->Fill(jet->get_eta(), jet->get_phi());
        h_jet_mass_pt_r05->Fill(jet->get_pt(), jet->get_mass());
        h_jet_mass_eta_r05->Fill(jet->get_eta(), jet->get_mass());
      }
    }
  }

  if (Verbosity() > 1)
  {
    std::cout << "JetKinematicCheck::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetKinematicCheck::EndRun(const int runnumber)
{
  if (Verbosity() > 1)
  {
    std::cout << "JetKinematicCheck::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int JetKinematicCheck::End(PHCompositeNode * /*unused*/)
{
  if (Verbosity() > 1)
  {
    std::cout << "JetKinematicCheck::End(PHCompositeNode *topNode) Entering the end" << std::endl;
  }

  // for jet spectra [R02]
  h_jet_pt_spectra_r02->SetMarkerStyle(8);
  h_jet_pt_spectra_r02->SetMarkerColor(1);
  h_jet_pt_spectra_r02->SetLineColor(1);
  h_jet_pt_spectra_r02->SetTitle("");
  h_jet_pt_spectra_r02->GetYaxis()->SetTitle("Counts");
  h_jet_pt_spectra_r02->SetStats(false);

  // for jet spectra [R03]
  h_jet_pt_spectra_r03->SetMarkerStyle(8);
  h_jet_pt_spectra_r03->SetMarkerColor(1);
  h_jet_pt_spectra_r03->SetLineColor(1);
  h_jet_pt_spectra_r03->SetTitle("");
  h_jet_pt_spectra_r03->GetYaxis()->SetTitle("Counts");
  h_jet_pt_spectra_r03->SetStats(false);

  // for jet spectra [R04]
  h_jet_pt_spectra_r04->SetMarkerStyle(8);
  h_jet_pt_spectra_r04->SetMarkerColor(1);
  h_jet_pt_spectra_r04->SetLineColor(1);
  h_jet_pt_spectra_r04->SetTitle("");
  h_jet_pt_spectra_r04->GetYaxis()->SetTitle("Counts");
  h_jet_pt_spectra_r04->SetStats(false);

  // for jet spectra [R05]
  h_jet_pt_spectra_r05->SetMarkerStyle(8);
  h_jet_pt_spectra_r05->SetMarkerColor(1);
  h_jet_pt_spectra_r05->SetLineColor(1);
  h_jet_pt_spectra_r05->SetTitle("");
  h_jet_pt_spectra_r05->GetYaxis()->SetTitle("Counts");
  h_jet_pt_spectra_r05->SetStats(false);

  // for jet eta-phi [R02]
  h_jet_eta_phi_r02->SetStats(false);
  h_jet_eta_phi_r02->SetTitle("");

  // for jet eta-phi [R03]
  h_jet_eta_phi_r03->SetStats(false);
  h_jet_eta_phi_r03->SetTitle("");

  // for jet eta-phi [R04]
  h_jet_eta_phi_r04->SetStats(false);
  h_jet_eta_phi_r04->SetTitle("");

  // for jet eta-phi [R05]
  h_jet_eta_phi_r05->SetStats(false);
  h_jet_eta_phi_r05->SetTitle("");

  // for jet mass vs pt [R02]
  h_jet_mass_pt_r02->SetStats(false);
  h_jet_mass_pt_r02->SetTitle("");

  // for average jet mass vs pt [R02]
  h_jet_average_mass_pt_1D_r02 = (TH1D *) h_jet_mass_pt_r02->ProfileX();
  h_jet_average_mass_pt_1D_r02->SetStats(false);
  h_jet_average_mass_pt_1D_r02->SetTitle("");
  h_jet_average_mass_pt_1D_r02->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  h_jet_average_mass_pt_1D_r02->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");
  h_jet_average_mass_pt_1D_r02->SetMarkerStyle(8);
  h_jet_average_mass_pt_1D_r02->SetMarkerColor(1);
  h_jet_average_mass_pt_1D_r02->SetLineColor(1);

  // for jet mass vs pt [R03]
  h_jet_mass_pt_r03->SetStats(false);
  h_jet_mass_pt_r03->SetTitle("");

  // for average jet mass vs pt [R03]
  h_jet_average_mass_pt_1D_r03 = (TH1D *) h_jet_mass_pt_r03->ProfileX();
  h_jet_average_mass_pt_1D_r03->SetStats(false);
  h_jet_average_mass_pt_1D_r03->SetTitle("");
  h_jet_average_mass_pt_1D_r03->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  h_jet_average_mass_pt_1D_r03->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");
  h_jet_average_mass_pt_1D_r03->SetMarkerStyle(8);
  h_jet_average_mass_pt_1D_r03->SetMarkerColor(1);
  h_jet_average_mass_pt_1D_r03->SetLineColor(1);

  // for jet mass vs pt [R04]
  h_jet_mass_pt_r04->SetStats(false);
  h_jet_mass_pt_r04->SetTitle("");

  // for average jet mass vs pt [R04]
  h_jet_average_mass_pt_1D_r04 = (TH1D *) h_jet_mass_pt_r04->ProfileX();
  h_jet_average_mass_pt_1D_r04->SetStats(false);
  h_jet_average_mass_pt_1D_r04->SetTitle("");
  h_jet_average_mass_pt_1D_r04->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  h_jet_average_mass_pt_1D_r04->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");
  h_jet_average_mass_pt_1D_r04->SetMarkerStyle(8);
  h_jet_average_mass_pt_1D_r04->SetMarkerColor(1);
  h_jet_average_mass_pt_1D_r04->SetLineColor(1);

  // for jet mass vs pt [R05]
  h_jet_mass_pt_r05->SetStats(false);
  h_jet_mass_pt_r05->SetTitle("");

  // for average jet mass vs pt [R05]
  h_jet_average_mass_pt_1D_r05 = (TH1D *) h_jet_mass_pt_r05->ProfileX();
  h_jet_average_mass_pt_1D_r05->SetStats(false);
  h_jet_average_mass_pt_1D_r05->SetTitle("");
  h_jet_average_mass_pt_1D_r05->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  h_jet_average_mass_pt_1D_r05->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");
  h_jet_average_mass_pt_1D_r05->SetMarkerStyle(8);
  h_jet_average_mass_pt_1D_r05->SetMarkerColor(1);
  h_jet_average_mass_pt_1D_r05->SetLineColor(1);

  // for jet mass vs eta [R02]
  h_jet_mass_eta_r02->SetStats(false);
  h_jet_mass_eta_r02->SetTitle("");

  // for average jet mass vs eta [R02]
  h_jet_average_mass_eta_1D_r02 = (TH1D *) h_jet_mass_eta_r02->ProfileX();
  h_jet_average_mass_eta_1D_r02->SetStats(false);
  h_jet_average_mass_eta_1D_r02->SetTitle("");
  h_jet_average_mass_eta_1D_r02->GetXaxis()->SetTitle("#eta");
  h_jet_average_mass_eta_1D_r02->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");
  h_jet_average_mass_eta_1D_r02->SetMarkerStyle(8);
  h_jet_average_mass_eta_1D_r02->SetMarkerColor(1);
  h_jet_average_mass_eta_1D_r02->SetLineColor(1);

  // for jet mass vs eta [R03]
  h_jet_mass_eta_r03->SetStats(false);
  h_jet_mass_eta_r03->SetTitle("");

  // for average jet mass vs eta [R03]
  h_jet_average_mass_eta_1D_r03 = (TH1D *) h_jet_mass_eta_r03->ProfileX();
  h_jet_average_mass_eta_1D_r03->SetStats(false);
  h_jet_average_mass_eta_1D_r03->SetTitle("");
  h_jet_average_mass_eta_1D_r03->GetXaxis()->SetTitle("#eta");
  h_jet_average_mass_eta_1D_r03->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");
  h_jet_average_mass_eta_1D_r03->SetMarkerStyle(8);
  h_jet_average_mass_eta_1D_r03->SetMarkerColor(1);
  h_jet_average_mass_eta_1D_r03->SetLineColor(1);

  // for jet mass vs eta [R04]
  h_jet_mass_eta_r04->SetStats(false);
  h_jet_mass_eta_r04->SetTitle("");

  // for average jet mass vs eta [R04]
  h_jet_average_mass_eta_1D_r04 = (TH1D *) h_jet_mass_eta_r04->ProfileX();
  h_jet_average_mass_eta_1D_r04->SetStats(false);
  h_jet_average_mass_eta_1D_r04->SetTitle("");
  h_jet_average_mass_eta_1D_r04->GetXaxis()->SetTitle("#eta");
  h_jet_average_mass_eta_1D_r04->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");
  h_jet_average_mass_eta_1D_r04->SetMarkerStyle(8);
  h_jet_average_mass_eta_1D_r04->SetMarkerColor(1);
  h_jet_average_mass_eta_1D_r04->SetLineColor(1);

  // for jet mass vs eta [R05]
  h_jet_mass_eta_r05->SetStats(false);
  h_jet_mass_eta_r05->SetTitle("");

  // for average jet mass vs eta [R05]
  h_jet_average_mass_eta_1D_r05 = (TH1D *) h_jet_mass_eta_r05->ProfileX();
  h_jet_average_mass_eta_1D_r05->SetStats(false);
  h_jet_average_mass_eta_1D_r05->SetTitle("");
  h_jet_average_mass_eta_1D_r05->GetXaxis()->SetTitle("#eta");
  h_jet_average_mass_eta_1D_r05->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");
  h_jet_average_mass_eta_1D_r05->SetMarkerStyle(8);
  h_jet_average_mass_eta_1D_r05->SetMarkerColor(1);
  h_jet_average_mass_eta_1D_r05->SetLineColor(1);

  hm->registerHisto(h_jet_pt_spectra_r02);
  hm->registerHisto(h_jet_pt_spectra_r03);
  hm->registerHisto(h_jet_pt_spectra_r04);
  hm->registerHisto(h_jet_pt_spectra_r05);
  hm->registerHisto(h_jet_eta_phi_r02);
  hm->registerHisto(h_jet_eta_phi_r03);
  hm->registerHisto(h_jet_eta_phi_r04);
  hm->registerHisto(h_jet_eta_phi_r05);
  hm->registerHisto(h_jet_mass_pt_r02);
  hm->registerHisto(h_jet_mass_pt_r03);
  hm->registerHisto(h_jet_mass_pt_r04);
  hm->registerHisto(h_jet_mass_pt_r05);
  hm->registerHisto(h_jet_average_mass_pt_1D_r02);
  hm->registerHisto(h_jet_average_mass_pt_1D_r03);
  hm->registerHisto(h_jet_average_mass_pt_1D_r04);
  hm->registerHisto(h_jet_average_mass_pt_1D_r05);
  hm->registerHisto(h_jet_mass_eta_r02);
  hm->registerHisto(h_jet_mass_eta_r03);
  hm->registerHisto(h_jet_mass_eta_r04);
  hm->registerHisto(h_jet_mass_eta_r05);
  hm->registerHisto(h_jet_average_mass_eta_1D_r02);
  hm->registerHisto(h_jet_average_mass_eta_1D_r03);
  hm->registerHisto(h_jet_average_mass_eta_1D_r04);
  hm->registerHisto(h_jet_average_mass_eta_1D_r05);

  if (Verbosity() > 1)
  {
    std::cout << "JetKinematicCheck::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
