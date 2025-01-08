#include "RhosinEvent.h"

// fun4all includes
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
// #include <fun4all/PHTFileServer.h>

#include <qautils/QAHistManagerDef.h>

// phool includes
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

// jetbackground includes
#include <jetbackground/TowerRho.h>
#include <jetbackground/TowerRhov1.h>

#include <TH1.h>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <string>
#include <vector>

RhosinEvent::RhosinEvent(const std::string& moduleName, const std::string& tag)
  : SubsysReco(moduleName)
  , m_moduleName(moduleName)
  , m_histTag(tag)
  // , m_name(outputfilename)
  , m_do_mult_rho(true)
  , m_do_area_rho(true)
  , m_mult_rho_node("TowerRho_MULT")
  , m_area_rho_node("TowerRho_AREA")
  , m_doTrgSelect(false)
  , m_trgToSelect(JetQADefs::GL1::MBDNSJet1)
  , h1_mult_rho(nullptr)
  , h1_mult_rho_sigma(nullptr)
  , h1_area_rho(nullptr)
  , h1_area_rho_sigma(nullptr)
{
}

int RhosinEvent::Init(PHCompositeNode* /*topNode*/)
{
  // if (Verbosity() > 0)
  // {
  // std::cout << "RhosinEvent::Init - Output to " << m_outputFileName << std::endl;
  // }

  // create output file
  // PHTFileServer::get().open(m_outputFileName, "RECREATE");

  m_manager = QAHistManagerDef::getHistoManager();
  if (!m_manager)
  {
    std::cerr << PHWHERE << ": PANIC: couldn't grab histogram manager!" << std::endl;
    assert(m_manager);
  }

  // Initialize histograms
  const int N_rho_mult = 200;
  Double_t N_rho_mult_bins[N_rho_mult + 1];
  for (int i = 0; i <= N_rho_mult; i++)
  {
    N_rho_mult_bins[i] = (0.3 / 200.0) * i;
  }

  const int N_rho_area = 200;
  Double_t N_rho_area_bins[N_rho_area + 1];
  for (int i = 0; i <= N_rho_area; i++)
  {
    N_rho_area_bins[i] = (10.0 / 200.0) * i;
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
      "rhomult",
      "sigmamult",
      "rhoarea",
      "sigmaarea"};
  for (auto& vecHistName : vecHistNames)
  {
    vecHistName.insert(0, "h_" + smallModuleName + "_");
    if (!m_histTag.empty())
    {
      vecHistName.append("_" + m_histTag);
    }
  }

  h1_mult_rho = new TH1D(vecHistNames[0].data(), "h1_mult_rho", N_rho_mult, N_rho_mult_bins);
  h1_mult_rho->GetXaxis()->SetTitle("rho_M");
  h1_mult_rho->GetYaxis()->SetTitle("Counts");

  h1_mult_rho_sigma = new TH1D(vecHistNames[1].data(), "h1_mult_rho_sigma", N_rho_mult, N_rho_mult_bins);
  h1_mult_rho_sigma->GetXaxis()->SetTitle("sigma_M");
  h1_mult_rho_sigma->GetYaxis()->SetTitle("Counts");

  h1_area_rho = new TH1D(vecHistNames[2].data(), "h1_area_rho", N_rho_area, N_rho_area_bins);
  h1_area_rho->GetXaxis()->SetTitle("rho_A");
  h1_area_rho->GetYaxis()->SetTitle("Counts");

  h1_area_rho_sigma = new TH1D(vecHistNames[3].data(), "h1_area_rho_sigma", N_rho_area, N_rho_area_bins);
  h1_area_rho_sigma->GetXaxis()->SetTitle("sigma_A");
  h1_area_rho_sigma->GetYaxis()->SetTitle("Counts");

  if (Verbosity() > 0)
  {
    std::cout << "RhosinEvent::Init - Histograms created" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int RhosinEvent::process_event(PHCompositeNode* topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "RhosinEvent::process_event - Process event..." << std::endl;
  }

  // if needed, check if selected trigger fired
  if (m_doTrgSelect)
  {
    bool hasTrigger = JetQADefs::DidTriggerFire(m_trgToSelect, topNode);
    if (!hasTrigger)
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }

  if (m_do_area_rho)
  {
    TowerRho* towerrho = findNode::getClass<TowerRhov1>(topNode, m_area_rho_node);
    if (!towerrho)
    {
      std::cout << "RhosinEvent::process_event - Error can not find towerrho " << m_area_rho_node << std::endl;
      exit(-1);
    }
    h1_area_rho->Fill(towerrho->get_rho());
    h1_area_rho_sigma->Fill(towerrho->get_sigma());
  }

  if (m_do_mult_rho)
  {
    TowerRho* towerrho = findNode::getClass<TowerRhov1>(topNode, m_mult_rho_node);
    if (!towerrho)
    {
      std::cout << "RhosinEvent::process_event - Error can not find towerrho " << m_mult_rho_node << std::endl;
      exit(-1);
    }
    h1_mult_rho->Fill(towerrho->get_rho());
    h1_mult_rho_sigma->Fill(towerrho->get_sigma());
  }

  if (Verbosity() > 1)
  {
    std::cout << "RhosinEvent::process_event - Event processed" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int RhosinEvent::End(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() > 0)
  {
    std::cout << "RhosinEvent::EndRun - End run " << std::endl;
    // std::cout << "RhosinEvent::EndRun - Writing to " << m_outputFileName << std::endl;
  }

  // PHTFileServer::get().cd(m_outputFileName);

  // h1_mult_rho->Write();
  // h1_mult_rho_sigma->Write();
  // h1_area_rho->Write();
  // h1_area_rho_sigma->Write();

  // register histograms
  m_manager->registerHisto(h1_mult_rho);
  m_manager->registerHisto(h1_mult_rho_sigma);
  m_manager->registerHisto(h1_area_rho);
  m_manager->registerHisto(h1_area_rho_sigma);

  if (Verbosity() > 0)
  {
    std::cout << "RhosinEvent::EndRun - Done" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
