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
      "spectra_r02",
      "spectra_r03",
      "spectra_r04",
      "spectra_r05",
      "etavsphi_r02",
      "etavsphi_r03",
      "etavsphi_r04",
      "etavsphi_r05",
      "jetmassvspt_r02",
      "jetmassvspt_r03",
      "jetmassvspt_r04",
      "jetmassvspt_r05",
      "jetmassvseta_r02",
      "jetmassvseta_r03",
      "jetmassvseta_r04",
      "jetmassvseta_r05",
      "jetmassvsptprofile_r02",
      "jetmassvsptprofile_r03",
      "jetmassvsptprofile_r04",
      "jetmassvsptprofile_r05",
      "jetmassvsetaprofile_r02",
      "jetmassvsetaprofile_r03",
      "jetmassvsetaprofile_r04",
      "jetmassvsetaprofile_r05"};
  for (auto &vecHistName : vecHistNames)
  {
    vecHistName.insert(0, "h_" + smallModuleName + "_");
    if (!m_histTag.empty())
    {
      vecHistName.append("_" + m_histTag);
    }
  }

  // initialize histograms

  jet_spectra_r02 = new TH1D(vecHistNames[0].data(), "", 19, 5, 100);
  jet_spectra_r02->GetXaxis()->SetTitle("p_{T} [GeV/c]");

  jet_spectra_r03 = new TH1D(vecHistNames[1].data(), "", 19, 5, 100);
  jet_spectra_r03->GetXaxis()->SetTitle("p_{T} [GeV/c]");

  jet_spectra_r04 = new TH1D(vecHistNames[2].data(), "", 19, 5, 100);
  jet_spectra_r04->GetXaxis()->SetTitle("p_{T} [GeV/c]");

  jet_spectra_r05 = new TH1D(vecHistNames[3].data(), "", 19, 5, 100);
  jet_spectra_r05->GetXaxis()->SetTitle("p_{T} [GeV/c]");

  jet_eta_phi_r02 = new TH2D(vecHistNames[4].data(), "", 24, -1.1, 1.1, 64, -M_PI, M_PI);
  jet_eta_phi_r02->GetXaxis()->SetTitle("#eta");
  jet_eta_phi_r02->GetYaxis()->SetTitle("#Phi");

  jet_eta_phi_r03 = new TH2D(vecHistNames[5].data(), "", 24, -1.1, 1.1, 64, -M_PI, M_PI);
  jet_eta_phi_r03->GetXaxis()->SetTitle("#eta");
  jet_eta_phi_r03->GetYaxis()->SetTitle("#Phi");

  jet_eta_phi_r04 = new TH2D(vecHistNames[6].data(), "", 24, -1.1, 1.1, 64, -M_PI, M_PI);
  jet_eta_phi_r04->GetXaxis()->SetTitle("#eta");
  jet_eta_phi_r04->GetYaxis()->SetTitle("#Phi");

  jet_eta_phi_r05 = new TH2D(vecHistNames[7].data(), "", 24, -1.1, 1.1, 64, -M_PI, M_PI);
  jet_eta_phi_r05->GetXaxis()->SetTitle("#eta");
  jet_eta_phi_r05->GetYaxis()->SetTitle("#Phi");

  jet_mass_pt_r02 = new TH2D(vecHistNames[8].data(), "", 19, 5, 100, 15, 0, 15);
  jet_mass_pt_r02->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  jet_mass_pt_r02->GetYaxis()->SetTitle("Jet Mass [GeV/c^{2}]");

  jet_mass_pt_r03 = new TH2D(vecHistNames[9].data(), "", 19, 5, 100, 15, 0, 15);
  jet_mass_pt_r03->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  jet_mass_pt_r03->GetYaxis()->SetTitle("Jet Mass [GeV/c^{2}]");

  jet_mass_pt_r04 = new TH2D(vecHistNames[10].data(), "", 19, 5, 100, 15, 0, 15);
  jet_mass_pt_r04->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  jet_mass_pt_r04->GetYaxis()->SetTitle("Jet Mass [GeV/c^{2}]");

  jet_mass_pt_r05 = new TH2D(vecHistNames[11].data(), "", 19, 5, 100, 15, 0, 15);
  jet_mass_pt_r05->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  jet_mass_pt_r05->GetYaxis()->SetTitle("Jet Mass [GeV/c^{2}]");

  jet_mass_eta_r02 = new TH2D(vecHistNames[12].data(), "", 24, -1.1, 1.1, 15, 0, 15);
  jet_mass_eta_r02->GetXaxis()->SetTitle("#eta");
  jet_mass_eta_r02->GetYaxis()->SetTitle("Jet Mass [GeV/c^{2}]");

  jet_mass_eta_r03 = new TH2D(vecHistNames[13].data(), "", 24, -1.1, 1.1, 15, 0, 15);
  jet_mass_eta_r03->GetXaxis()->SetTitle("#eta");
  jet_mass_eta_r03->GetYaxis()->SetTitle("Jet Mass [GeV/c^{2}]");

  jet_mass_eta_r04 = new TH2D(vecHistNames[14].data(), "", 24, -1.1, 1.1, 15, 0, 15);
  jet_mass_eta_r04->GetXaxis()->SetTitle("#eta");
  jet_mass_eta_r04->GetYaxis()->SetTitle("Jet Mass [GeV/c^{2}]");

  jet_mass_eta_r05 = new TH2D(vecHistNames[15].data(), "", 24, -1.1, 1.1, 15, 0, 15);
  jet_mass_eta_r05->GetXaxis()->SetTitle("#eta");
  jet_mass_eta_r05->GetYaxis()->SetTitle("Jet Mass [GeV/c^{2}]");

  jet_mass_pt_1D_r02 = new TH1D(vecHistNames[16].data(), "", 19, 5, 100);
  jet_mass_pt_1D_r02->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  jet_mass_pt_1D_r02->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");

  jet_mass_pt_1D_r03 = new TH1D(vecHistNames[17].data(), "", 19, 5, 100);
  jet_mass_pt_1D_r03->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  jet_mass_pt_1D_r03->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");

  jet_mass_pt_1D_r04 = new TH1D(vecHistNames[18].data(), "", 19, 5, 100);
  jet_mass_pt_1D_r04->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  jet_mass_pt_1D_r04->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");

  jet_mass_pt_1D_r05 = new TH1D(vecHistNames[19].data(), "", 19, 5, 100);
  jet_mass_pt_1D_r05->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  jet_mass_pt_1D_r05->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");

  jet_mass_eta_1D_r02 = new TH2D(vecHistNames[20].data(), "", 24, -1.1, 1.1, 15, 0, 15);
  jet_mass_eta_1D_r02->GetXaxis()->SetTitle("#eta");
  jet_mass_eta_1D_r02->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");

  jet_mass_eta_1D_r03 = new TH2D(vecHistNames[21].data(), "", 24, -1.1, 1.1, 15, 0, 15);
  jet_mass_eta_1D_r03->GetXaxis()->SetTitle("#eta");
  jet_mass_eta_1D_r03->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");

  jet_mass_eta_1D_r04 = new TH2D(vecHistNames[22].data(), "", 24, -1.1, 1.1, 15, 0, 15);
  jet_mass_eta_1D_r04->GetXaxis()->SetTitle("#eta");
  jet_mass_eta_1D_r04->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");

  jet_mass_eta_1D_r05 = new TH2D(vecHistNames[23].data(), "", 24, -1.1, 1.1, 15, 0, 15);
  jet_mass_eta_1D_r05->GetXaxis()->SetTitle("#eta");
  jet_mass_eta_1D_r05->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");

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

    std::string recoJetName = m_recoJetName_array[i];

    JetContainer *jets = findNode::getClass<JetContainer>(topNode, recoJetName);
    if (!jets)
    {
      std::cout
          << "JetKinematicCheck::process_event - Error can not find DST Reco JetContainer node "
          << recoJetName << std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }

    // loop over jets
    for (auto jet : *jets)
    {
      bool eta_cut = (jet->get_eta() >= etaRangeUse.first) and (jet->get_eta() <= etaRangeUse.second);
      bool pt_cut = (jet->get_pt() >= m_ptRange.first) and (jet->get_pt() <= m_ptRange.second);
      if ((not eta_cut) or (not pt_cut))
      {
        continue;
      }
      if (jet->get_pt() < 1)
      {
        continue;  // to remove noise jets
      }

      if (i == 0)
      {
        jet_spectra_r02->Fill(jet->get_pt());
        jet_eta_phi_r02->Fill(jet->get_eta(), jet->get_phi());
        jet_mass_pt_r02->Fill(jet->get_pt(), jet->get_mass());
        jet_mass_eta_r02->Fill(jet->get_eta(), jet->get_mass());
      }

      else if (i == 1)
      {
        jet_spectra_r03->Fill(jet->get_pt());
        jet_eta_phi_r03->Fill(jet->get_eta(), jet->get_phi());
        jet_mass_pt_r03->Fill(jet->get_pt(), jet->get_mass());
        jet_mass_eta_r03->Fill(jet->get_eta(), jet->get_mass());
      }

      else if (i == 2)
      {
        jet_spectra_r04->Fill(jet->get_pt());
        jet_eta_phi_r04->Fill(jet->get_eta(), jet->get_phi());
        jet_mass_pt_r04->Fill(jet->get_pt(), jet->get_mass());
        jet_mass_eta_r04->Fill(jet->get_eta(), jet->get_mass());
      }

      else if (i == 3)
      {
        jet_spectra_r05->Fill(jet->get_pt());
        jet_eta_phi_r05->Fill(jet->get_eta(), jet->get_phi());
        jet_mass_pt_r05->Fill(jet->get_pt(), jet->get_mass());
        jet_mass_eta_r05->Fill(jet->get_eta(), jet->get_mass());
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
  TLegend *leg1 = new TLegend(.7, .9, .9, 1);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.06);
  leg1->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < p_{T} < %2% [GeV/c]") % m_ptRange.first % m_ptRange.second).c_str(), "");
  leg1->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < #eta < %2%") % m_etaRange.first % m_etaRange.second).c_str(), "");
  jet_spectra_r02->SetMarkerStyle(8);
  jet_spectra_r02->SetMarkerColor(1);
  jet_spectra_r02->SetLineColor(1);
  jet_spectra_r02->SetTitle("Jet Spectra [R = 0.2]");
  jet_spectra_r02->GetYaxis()->SetTitle("Counts");
  jet_spectra_r02->GetListOfFunctions()->Add(leg1);
  jet_spectra_r02->SetStats(false);

  // for jet spectra [R03]
  TLegend *leg2 = new TLegend(.7, .9, .9, 1);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.06);
  leg2->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < p_{T} < %2% [GeV/c]") % m_ptRange.first % m_ptRange.second).c_str(), "");
  leg2->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < #eta < %2%") % m_etaRange.first % m_etaRange.second).c_str(), "");
  jet_spectra_r03->SetMarkerStyle(8);
  jet_spectra_r03->SetMarkerColor(1);
  jet_spectra_r03->SetLineColor(1);
  jet_spectra_r03->SetTitle("Jet Spectra [R = 0.3]");
  jet_spectra_r03->GetYaxis()->SetTitle("Counts");
  jet_spectra_r03->GetListOfFunctions()->Add(leg2);
  jet_spectra_r03->SetStats(false);

  // for jet spectra [R04]
  TLegend *leg3 = new TLegend(.7, .9, .9, 1);
  leg3->SetFillStyle(0);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.06);
  leg3->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < p_{T} < %2% [GeV/c]") % m_ptRange.first % m_ptRange.second).c_str(), "");
  leg3->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < #eta < %2%") % m_etaRange.first % m_etaRange.second).c_str(), "");
  jet_spectra_r04->SetMarkerStyle(8);
  jet_spectra_r04->SetMarkerColor(1);
  jet_spectra_r04->SetLineColor(1);
  jet_spectra_r04->SetTitle("Jet Spectra [R = 0.4]");
  jet_spectra_r04->GetYaxis()->SetTitle("Counts");
  jet_spectra_r04->GetListOfFunctions()->Add(leg3);
  jet_spectra_r04->SetStats(false);

  // for jet spectra [R05]
  TLegend *leg4 = new TLegend(.7, .9, .9, 1);
  leg4->SetFillStyle(0);
  leg4->SetBorderSize(0);
  leg4->SetTextSize(0.06);
  leg4->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < p_{T} < %2% [GeV/c]") % m_ptRange.first % m_ptRange.second).c_str(), "");
  leg4->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < #eta < %2%") % m_etaRange.first % m_etaRange.second).c_str(), "");
  jet_spectra_r05->SetMarkerStyle(8);
  jet_spectra_r05->SetMarkerColor(1);
  jet_spectra_r05->SetLineColor(1);
  jet_spectra_r05->SetTitle("Jet Spectra [R = 0.5]");
  jet_spectra_r05->GetYaxis()->SetTitle("Counts");
  jet_spectra_r05->GetListOfFunctions()->Add(leg4);
  jet_spectra_r05->SetStats(false);

  // for jet eta-phi [R02]
  TLegend *leg5 = new TLegend(.7, .9, .9, 1);
  leg5->SetFillStyle(0);
  leg5->SetBorderSize(0);
  leg5->SetTextSize(0.06);
  leg5->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < p_{T} < %2% [GeV/c]") % m_ptRange.first % m_ptRange.second).c_str(), "");
  leg5->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < #eta < %2%") % m_etaRange.first % m_etaRange.second).c_str(), "");
  jet_eta_phi_r02->SetStats(false);
  jet_eta_phi_r02->SetTitle("Jet Eta-Phi [R = 0.2]");
  jet_eta_phi_r02->GetListOfFunctions()->Add(leg5);

  // for jet eta-phi [R03]
  TLegend *leg6 = new TLegend(.7, .9, .9, 1);
  leg6->SetFillStyle(0);
  leg6->SetBorderSize(0);
  leg6->SetTextSize(0.06);
  leg6->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < p_{T} < %2% [GeV/c]") % m_ptRange.first % m_ptRange.second).c_str(), "");
  leg6->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < #eta < %2%") % m_etaRange.first % m_etaRange.second).c_str(), "");
  jet_eta_phi_r03->SetStats(false);
  jet_eta_phi_r03->SetTitle("Jet Eta-Phi [R = 0.3]");
  jet_eta_phi_r03->GetListOfFunctions()->Add(leg6);

  // for jet eta-phi [R04]
  TLegend *leg7 = new TLegend(.7, .9, .9, 1);
  leg7->SetFillStyle(0);
  leg7->SetBorderSize(0);
  leg7->SetTextSize(0.06);
  leg7->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < p_{T} < %2% [GeV/c]") % m_ptRange.first % m_ptRange.second).c_str(), "");
  leg7->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < #eta < %2%") % m_etaRange.first % m_etaRange.second).c_str(), "");
  jet_eta_phi_r04->SetStats(false);
  jet_eta_phi_r04->SetTitle("Jet Eta-Phi [R = 0.4]");
  jet_eta_phi_r04->GetListOfFunctions()->Add(leg7);

  // for jet eta-phi [R05]
  TLegend *leg8 = new TLegend(.7, .9, .9, 1);
  leg8->SetFillStyle(0);
  leg8->SetBorderSize(0);
  leg8->SetTextSize(0.06);
  leg8->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < p_{T} < %2% [GeV/c]") % m_ptRange.first % m_ptRange.second).c_str(), "");
  leg8->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < #eta < %2%") % m_etaRange.first % m_etaRange.second).c_str(), "");
  jet_eta_phi_r05->SetStats(false);
  jet_eta_phi_r05->SetTitle("Jet Eta-Phi [R = 0.5]");
  jet_eta_phi_r05->GetListOfFunctions()->Add(leg8);

  // for jet mass vs pt [R02]
  TLegend *leg9 = new TLegend(.7, .9, .9, 1);
  leg9->SetFillStyle(0);
  leg9->SetBorderSize(0);
  leg9->SetTextSize(0.06);
  leg9->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < p_{T} < %2% [GeV/c]") % m_ptRange.first % m_ptRange.second).c_str(), "");
  leg9->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < #eta < %2%") % m_etaRange.first % m_etaRange.second).c_str(), "");
  jet_mass_pt_r02->SetStats(false);
  jet_mass_pt_r02->SetTitle("Jet Mass vs p_{T} [R = 0.2]");
  jet_mass_pt_r02->GetListOfFunctions()->Add(leg9);

  // for average jet mass vs pt [R02]
  TLegend *leg10 = new TLegend(.7, .9, .9, 1);
  leg10->SetFillStyle(0);
  leg10->SetBorderSize(0);
  leg10->SetTextSize(0.06);
  leg10->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < p_{T} < %2% [GeV/c]") % m_ptRange.first % m_ptRange.second).c_str(), "");
  leg10->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < #eta < %2%") % m_etaRange.first % m_etaRange.second).c_str(), "");
  jet_mass_pt_1D_r02 = (TH1D *) jet_mass_pt_r02->ProfileX();
  jet_mass_pt_1D_r02->SetStats(false);
  jet_mass_pt_1D_r02->SetTitle("Average Jet Mass vs p_{T} [R = 0.2]");
  jet_mass_pt_1D_r02->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  jet_mass_pt_1D_r02->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");
  jet_mass_pt_1D_r02->SetMarkerStyle(8);
  jet_mass_pt_1D_r02->SetMarkerColor(1);
  jet_mass_pt_1D_r02->SetLineColor(1);
  jet_mass_pt_1D_r02->GetListOfFunctions()->Add(leg10);

  // for jet mass vs pt [R03]
  TLegend *leg11 = new TLegend(.7, .9, .9, 1);
  leg11->SetFillStyle(0);
  leg11->SetBorderSize(0);
  leg11->SetTextSize(0.06);
  leg11->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < p_{T} < %2% [GeV/c]") % m_ptRange.first % m_ptRange.second).c_str(), "");
  leg11->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < #eta < %2%") % m_etaRange.first % m_etaRange.second).c_str(), "");
  jet_mass_pt_r03->SetStats(false);
  jet_mass_pt_r03->SetTitle("Jet Mass vs p_{T} [R = 0.3]");
  jet_mass_pt_r03->GetListOfFunctions()->Add(leg11);

  // for average jet mass vs pt [R03]
  TLegend *leg12 = new TLegend(.7, .9, .9, 1);
  leg12->SetFillStyle(0);
  leg12->SetBorderSize(0);
  leg12->SetTextSize(0.06);
  leg12->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < p_{T} < %2% [GeV/c]") % m_ptRange.first % m_ptRange.second).c_str(), "");
  leg12->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < #eta < %2%") % m_etaRange.first % m_etaRange.second).c_str(), "");
  jet_mass_pt_1D_r03 = (TH1D *) jet_mass_pt_r03->ProfileX();
  jet_mass_pt_1D_r03->SetStats(false);
  jet_mass_pt_1D_r03->SetTitle("Average Jet Mass vs p_{T} [R = 0.3]");
  jet_mass_pt_1D_r03->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  jet_mass_pt_1D_r03->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");
  jet_mass_pt_1D_r03->SetMarkerStyle(8);
  jet_mass_pt_1D_r03->SetMarkerColor(1);
  jet_mass_pt_1D_r03->SetLineColor(1);
  jet_mass_pt_1D_r03->GetListOfFunctions()->Add(leg12);

  // for jet mass vs pt [R04]
  TLegend *leg13 = new TLegend(.7, .9, .9, 1);
  leg13->SetFillStyle(0);
  leg13->SetBorderSize(0);
  leg13->SetTextSize(0.06);
  leg13->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < p_{T} < %2% [GeV/c]") % m_ptRange.first % m_ptRange.second).c_str(), "");
  leg13->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < #eta < %2%") % m_etaRange.first % m_etaRange.second).c_str(), "");
  jet_mass_pt_r04->SetStats(false);
  jet_mass_pt_r04->SetTitle("Jet Mass vs p_{T} [R = 0.4]");
  jet_mass_pt_r04->GetListOfFunctions()->Add(leg13);

  // for average jet mass vs pt [R04]
  TLegend *leg14 = new TLegend(.7, .9, .9, 1);
  leg14->SetFillStyle(0);
  leg14->SetBorderSize(0);
  leg14->SetTextSize(0.06);
  leg14->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < p_{T} < %2% [GeV/c]") % m_ptRange.first % m_ptRange.second).c_str(), "");
  leg14->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < #eta < %2%") % m_etaRange.first % m_etaRange.second).c_str(), "");
  jet_mass_pt_1D_r04 = (TH1D *) jet_mass_pt_r04->ProfileX();
  jet_mass_pt_1D_r04->SetStats(false);
  jet_mass_pt_1D_r04->SetTitle("Average Jet Mass vs p_{T} [R = 0.4]");
  jet_mass_pt_1D_r04->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  jet_mass_pt_1D_r04->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");
  jet_mass_pt_1D_r04->SetMarkerStyle(8);
  jet_mass_pt_1D_r04->SetMarkerColor(1);
  jet_mass_pt_1D_r04->SetLineColor(1);
  jet_mass_pt_1D_r04->GetListOfFunctions()->Add(leg14);

  // for jet mass vs pt [R05]
  TLegend *leg15 = new TLegend(.7, .9, .9, 1);
  leg15->SetFillStyle(0);
  leg15->SetBorderSize(0);
  leg15->SetTextSize(0.06);
  leg15->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < p_{T} < %2% [GeV/c]") % m_ptRange.first % m_ptRange.second).c_str(), "");
  leg15->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < #eta < %2%") % m_etaRange.first % m_etaRange.second).c_str(), "");
  jet_mass_pt_r05->SetStats(false);
  jet_mass_pt_r05->SetTitle("Jet Mass vs p_{T} [R = 0.5]");
  jet_mass_pt_r05->GetListOfFunctions()->Add(leg15);

  // for average jet mass vs pt [R05]
  TLegend *leg16 = new TLegend(.7, .9, .9, 1);
  leg16->SetFillStyle(0);
  leg16->SetBorderSize(0);
  leg16->SetTextSize(0.06);
  leg16->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < p_{T} < %2% [GeV/c]") % m_ptRange.first % m_ptRange.second).c_str(), "");
  leg16->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < #eta < %2%") % m_etaRange.first % m_etaRange.second).c_str(), "");
  jet_mass_pt_1D_r05 = (TH1D *) jet_mass_pt_r05->ProfileX();
  jet_mass_pt_1D_r05->SetStats(false);
  jet_mass_pt_1D_r05->SetTitle("Average Jet Mass vs p_{T} [R = 0.5]");
  jet_mass_pt_1D_r05->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  jet_mass_pt_1D_r05->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");
  jet_mass_pt_1D_r05->SetMarkerStyle(8);
  jet_mass_pt_1D_r05->SetMarkerColor(1);
  jet_mass_pt_1D_r05->SetLineColor(1);
  jet_mass_pt_1D_r05->GetListOfFunctions()->Add(leg16);

  // for jet mass vs eta [R02]
  TLegend *leg17 = new TLegend(.7, .9, .9, 1);
  leg17->SetFillStyle(0);
  leg17->SetBorderSize(0);
  leg17->SetTextSize(0.06);
  leg17->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < p_{T} < %2% [GeV/c]") % m_ptRange.first % m_ptRange.second).c_str(), "");
  leg17->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < #eta < %2%") % m_etaRange.first % m_etaRange.second).c_str(), "");
  jet_mass_eta_r02->SetStats(false);
  jet_mass_eta_r02->SetTitle("Jet Mass vs #eta [R = 0.2]");
  jet_mass_eta_r02->GetListOfFunctions()->Add(leg17);

  // for average jet mass vs eta [R02]
  TLegend *leg18 = new TLegend(.7, .9, .9, 1);
  leg18->SetFillStyle(0);
  leg18->SetBorderSize(0);
  leg18->SetTextSize(0.06);
  leg18->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < p_{T} < %2% [GeV/c]") % m_ptRange.first % m_ptRange.second).c_str(), "");
  leg18->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < #eta < %2%") % m_etaRange.first % m_etaRange.second).c_str(), "");
  jet_mass_eta_1D_r02 = (TH1D *) jet_mass_eta_r02->ProfileX();
  jet_mass_eta_1D_r02->SetStats(false);
  jet_mass_eta_1D_r02->SetTitle("Average Jet Mass vs #eta [R = 0.2]");
  jet_mass_eta_1D_r02->GetXaxis()->SetTitle("#eta");
  jet_mass_eta_1D_r02->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");
  jet_mass_eta_1D_r02->SetMarkerStyle(8);
  jet_mass_eta_1D_r02->SetMarkerColor(1);
  jet_mass_eta_1D_r02->SetLineColor(1);
  jet_mass_eta_1D_r02->GetListOfFunctions()->Add(leg18);

  // for jet mass vs eta [R03]
  TLegend *leg19 = new TLegend(.7, .9, .9, 1);
  leg19->SetFillStyle(0);
  leg19->SetBorderSize(0);
  leg19->SetTextSize(0.06);
  leg19->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < p_{T} < %2% [GeV/c]") % m_ptRange.first % m_ptRange.second).c_str(), "");
  leg19->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < #eta < %2%") % m_etaRange.first % m_etaRange.second).c_str(), "");
  jet_mass_eta_r03->SetStats(false);
  jet_mass_eta_r03->SetTitle("Jet Mass vs #eta [R = 0.3]");
  jet_mass_eta_r03->GetListOfFunctions()->Add(leg19);

  // for average jet mass vs eta [R03]
  TLegend *leg20 = new TLegend(.7, .9, .9, 1);
  leg20->SetFillStyle(0);
  leg20->SetBorderSize(0);
  leg20->SetTextSize(0.06);
  leg20->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < p_{T} < %2% [GeV/c]") % m_ptRange.first % m_ptRange.second).c_str(), "");
  leg20->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < #eta < %2%") % m_etaRange.first % m_etaRange.second).c_str(), "");
  jet_mass_eta_1D_r03 = (TH1D *) jet_mass_eta_r03->ProfileX();
  jet_mass_eta_1D_r03->SetStats(false);
  jet_mass_eta_1D_r03->SetTitle("Average Jet Mass vs #eta [R = 0.3]");
  jet_mass_eta_1D_r03->GetXaxis()->SetTitle("#eta");
  jet_mass_eta_1D_r03->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");
  jet_mass_eta_1D_r03->SetMarkerStyle(8);
  jet_mass_eta_1D_r03->SetMarkerColor(1);
  jet_mass_eta_1D_r03->SetLineColor(1);
  jet_mass_eta_1D_r03->GetListOfFunctions()->Add(leg20);

  // for jet mass vs eta [R04]
  TLegend *leg21 = new TLegend(.7, .9, .9, 1);
  leg21->SetFillStyle(0);
  leg21->SetBorderSize(0);
  leg21->SetTextSize(0.06);
  leg21->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < p_{T} < %2% [GeV/c]") % m_ptRange.first % m_ptRange.second).c_str(), "");
  leg21->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < #eta < %2%") % m_etaRange.first % m_etaRange.second).c_str(), "");
  jet_mass_eta_r04->SetStats(false);
  jet_mass_eta_r04->SetTitle("Jet Mass vs #eta [R = 0.4]");
  jet_mass_eta_r04->GetListOfFunctions()->Add(leg21);

  // for average jet mass vs eta [R04]
  TLegend *leg22 = new TLegend(.7, .9, .9, 1);
  leg22->SetFillStyle(0);
  leg22->SetBorderSize(0);
  leg22->SetTextSize(0.06);
  leg22->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < p_{T} < %2% [GeV/c]") % m_ptRange.first % m_ptRange.second).c_str(), "");
  leg22->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < #eta < %2%") % m_etaRange.first % m_etaRange.second).c_str(), "");
  jet_mass_eta_1D_r04 = (TH1D *) jet_mass_eta_r04->ProfileX();
  jet_mass_eta_1D_r04->SetStats(false);
  jet_mass_eta_1D_r04->SetTitle("Average Jet Mass vs #eta [R = 0.4]");
  jet_mass_eta_1D_r04->GetXaxis()->SetTitle("#eta");
  jet_mass_eta_1D_r04->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");
  jet_mass_eta_1D_r04->SetMarkerStyle(8);
  jet_mass_eta_1D_r04->SetMarkerColor(1);
  jet_mass_eta_1D_r04->SetLineColor(1);
  jet_mass_eta_1D_r04->GetListOfFunctions()->Add(leg22);

  // for jet mass vs eta [R05]
  TLegend *leg23 = new TLegend(.7, .9, .9, 1);
  leg23->SetFillStyle(0);
  leg23->SetBorderSize(0);
  leg23->SetTextSize(0.06);
  leg23->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < p_{T} < %2% [GeV/c]") % m_ptRange.first % m_ptRange.second).c_str(), "");
  leg23->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < #eta < %2%") % m_etaRange.first % m_etaRange.second).c_str(), "");
  jet_mass_eta_r05->SetStats(false);
  jet_mass_eta_r05->SetTitle("Jet Mass vs #eta [R = 0.5]");
  jet_mass_eta_r05->GetListOfFunctions()->Add(leg23);

  // for average jet mass vs eta [R05]
  TLegend *leg24 = new TLegend(.7, .9, .9, 1);
  leg24->SetFillStyle(0);
  leg24->SetBorderSize(0);
  leg24->SetTextSize(0.06);
  leg24->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < p_{T} < %2% [GeV/c]") % m_ptRange.first % m_ptRange.second).c_str(), "");
  leg24->AddEntry((TObject *) nullptr, boost::str(boost::format("%1% < #eta < %2%") % m_etaRange.first % m_etaRange.second).c_str(), "");
  jet_mass_eta_1D_r05 = (TH1D *) jet_mass_eta_r05->ProfileX();
  jet_mass_eta_1D_r05->SetStats(false);
  jet_mass_eta_1D_r05->SetTitle("Average Jet Mass vs #eta [R = 0.5]");
  jet_mass_eta_1D_r05->GetXaxis()->SetTitle("#eta");
  jet_mass_eta_1D_r05->GetYaxis()->SetTitle("Average Jet Mass [GeV/c^{2}]");
  jet_mass_eta_1D_r05->SetMarkerStyle(8);
  jet_mass_eta_1D_r05->SetMarkerColor(1);
  jet_mass_eta_1D_r05->SetLineColor(1);
  jet_mass_eta_1D_r05->GetListOfFunctions()->Add(leg24);

  hm->registerHisto(jet_spectra_r02);
  hm->registerHisto(jet_spectra_r03);
  hm->registerHisto(jet_spectra_r04);
  hm->registerHisto(jet_spectra_r05);
  hm->registerHisto(jet_eta_phi_r02);
  hm->registerHisto(jet_eta_phi_r03);
  hm->registerHisto(jet_eta_phi_r04);
  hm->registerHisto(jet_eta_phi_r05);
  hm->registerHisto(jet_mass_pt_r02);
  hm->registerHisto(jet_mass_pt_r03);
  hm->registerHisto(jet_mass_pt_r04);
  hm->registerHisto(jet_mass_pt_r05);
  hm->registerHisto(jet_mass_pt_1D_r02);
  hm->registerHisto(jet_mass_pt_1D_r03);
  hm->registerHisto(jet_mass_pt_1D_r04);
  hm->registerHisto(jet_mass_pt_1D_r05);
  hm->registerHisto(jet_mass_eta_r02);
  hm->registerHisto(jet_mass_eta_r03);
  hm->registerHisto(jet_mass_eta_r04);
  hm->registerHisto(jet_mass_eta_r05);
  hm->registerHisto(jet_mass_eta_1D_r02);
  hm->registerHisto(jet_mass_eta_1D_r03);
  hm->registerHisto(jet_mass_eta_1D_r04);
  hm->registerHisto(jet_mass_eta_1D_r05);

  if (Verbosity() > 1)
  {
    std::cout << "JetKinematicCheck::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
