#include "TpcSiliconQA.h"

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <qautils/QAHistManagerDef.h>
#include <qautils/QAUtil.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TH1.h>
#include <TH2.h>

#include <boost/format.hpp>

//____________________________________________________________________________..
TpcSiliconQA::TpcSiliconQA(const std::string& name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int TpcSiliconQA::InitRun(PHCompositeNode* /*topNode*/)
{
  createHistos();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcSiliconQA::process_event(PHCompositeNode* topNode)
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  auto silseedmap = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  if (!silseedmap)
  {
    std::cout << "Silicon seed map not found, aborting event" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  auto tpcseedmap = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  if (!tpcseedmap)
  {
    std::cout << "TPC seed map not found, aborting event" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  for (const auto& silseed : *silseedmap)
  {
    if (!silseed)
    {
      continue;
    }

    m_crossing = (float) silseed->get_crossing();
    h_crossing->Fill(m_crossing);

    m_silseedx = silseed->get_x();
    m_silseedy = silseed->get_y();
    m_silseedz = silseed->get_z();
    m_silseedphi = silseed->get_phi();
    m_silseedeta = silseed->get_eta();

    for (const auto& tpcseed : *tpcseedmap)
    {
      if (!tpcseed)
      {
        continue;
      }

      m_tpcseedx = tpcseed->get_x();
      m_tpcseedy = tpcseed->get_y();
      m_tpcseedz = tpcseed->get_z();
      m_tpcseedphi = tpcseed->get_phi();
      m_tpcseedeta = tpcseed->get_eta();

      h_phiDiff[0]->Fill(m_tpcseedphi - m_silseedphi);
      h_etaDiff[0]->Fill(m_tpcseedeta - m_silseedeta);
      h_xDiff[0]->Fill(m_tpcseedx - m_silseedx);
      h_yDiff[0]->Fill(m_tpcseedy - m_silseedy);
      h_zDiff[0]->Fill(m_tpcseedz - m_silseedz);

      if (abs(m_tpcseedx - m_silseedx) > m_xcut || abs(m_tpcseedy - m_silseedy) > m_ycut)
      {
        continue;
      }

      h_phiDiff[1]->Fill(m_tpcseedphi - m_silseedphi);
      h_etaDiff[1]->Fill(m_tpcseedeta - m_silseedeta);
      h_xDiff[1]->Fill(m_tpcseedx - m_silseedx);
      h_yDiff[1]->Fill(m_tpcseedy - m_silseedy);
      h_zDiff[1]->Fill(m_tpcseedz - m_silseedz);

      if (abs(m_tpcseedeta - m_silseedeta) > m_etacut)
      {
        continue;
      }

      h_phiDiff[2]->Fill(m_tpcseedphi - m_silseedphi);
      h_etaDiff[2]->Fill(m_tpcseedeta - m_silseedeta);
      h_xDiff[2]->Fill(m_tpcseedx - m_silseedx);
      h_yDiff[2]->Fill(m_tpcseedy - m_silseedy);
      h_zDiff[2]->Fill(m_tpcseedz - m_silseedz);

      if (abs(m_tpcseedphi - m_silseedphi) > m_phicut)
      {
        continue;
      }

      h_phiDiff[3]->Fill(m_tpcseedphi - m_silseedphi);
      h_etaDiff[3]->Fill(m_tpcseedeta - m_silseedeta);
      h_xDiff[3]->Fill(m_tpcseedx - m_silseedx);
      h_yDiff[3]->Fill(m_tpcseedy - m_silseedy);
      h_zDiff[3]->Fill(m_tpcseedz - m_silseedz);
    }
  }

  m_event++;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcSiliconQA::EndRun(const int /*runnumber*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

std::string TpcSiliconQA::getHistoPrefix() const
{
  return std::string("h_") + Name() + std::string("_");
}

void TpcSiliconQA::createHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  std::vector<std::string> cutNames = {"", "_xyCut", "_etaCut", "_phiCut"};

  {
    h_crossing = new TH1F(std::string(getHistoPrefix() + "crossing").c_str(),
                          "Track Crossing Value", 1000, -500, 500);
    h_crossing->GetXaxis()->SetTitle("Track Crossing");
    h_crossing->GetYaxis()->SetTitle("Entries");
    hm->registerHisto(h_crossing);
  }
  /*
  {
    h_trackMatch = new TH1F(std::string(getHistoPrefix() + "trackMatch").c_str(),
                        "TPC and Silicon Seed Exist", 2, -0.5, 1.5);
    h_trackMatch->GetXaxis()->SetTitle("1 - TPC+Sil Seed, 0 - Missing TPC and/or Sil");
    h_trackMatch->GetYaxis()->SetTitle("Entries");
    hm->registerHisto(h_trackMatch);
  }
  */
  int i = 0;
  for (const std::string& name : cutNames)
  {
    h_phiDiff[i] = new TH1F(std::string(getHistoPrefix() + "phiDiff" + name).c_str(),
                            "TPC-Silicon #phi Difference", 100, -0.5, 0.5);
    h_phiDiff[i]->GetXaxis()->SetTitle("TPC Seed #phi - Silicon Seed #phi");
    h_phiDiff[i]->GetYaxis()->SetTitle("Entries");
    hm->registerHisto(h_phiDiff[i]);
    i++;
  }
  i = 0;
  for (const std::string& name : cutNames)
  {
    h_etaDiff[i] = new TH1F(std::string(getHistoPrefix() + "etaDiff" + name).c_str(),
                            "TPC-Silicon #eta Difference", 100, -0.1, 0.1);
    h_etaDiff[i]->GetXaxis()->SetTitle("TPC Seed #eta - Silicon Seed #eta");
    h_etaDiff[i]->GetYaxis()->SetTitle("Entries");
    hm->registerHisto(h_etaDiff[i]);
    i++;
  }
  i = 0;
  for (const std::string& name : cutNames)
  {
    h_xDiff[i] = new TH1F(std::string(getHistoPrefix() + "xDiff" + name).c_str(),
                          "TPC-Silicon x Difference", 100, -2, 2);
    h_xDiff[i]->GetXaxis()->SetTitle("TPC Seed x - Silicon Seed x [cm]");
    h_xDiff[i]->GetYaxis()->SetTitle("Entries");
    hm->registerHisto(h_xDiff[i]);
    i++;
  }
  i = 0;
  for (const std::string& name : cutNames)
  {
    h_yDiff[i] = new TH1F(std::string(getHistoPrefix() + "yDiff" + name).c_str(),
                          "TPC-Silicon y Difference", 100, -2, 2);
    h_yDiff[i]->GetXaxis()->SetTitle("TPC Seed y - Silicon Seed y [cm]");
    h_yDiff[i]->GetYaxis()->SetTitle("Entries");
    hm->registerHisto(h_yDiff[i]);
    i++;
  }
  i = 0;
  for (const std::string& name : cutNames)
  {
    h_zDiff[i] = new TH1F(std::string(getHistoPrefix() + "zDiff" + name).c_str(),
                          "TPC-Silicon z Difference", 500, -100, 100);
    h_zDiff[i]->GetXaxis()->SetTitle("TPC Seed z - Silicon Seed z [cm]");
    h_zDiff[i]->GetYaxis()->SetTitle("Entries");
    hm->registerHisto(h_zDiff[i]);
    i++;
  }

  return;
}
