#include "TpcSiliconQA.h"

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeedHelper.h>

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

    {
      const auto position = TrackSeedHelper::get_xyz(silseed);
      m_silseedx = position.x();
      m_silseedy = position.y();
      m_silseedz = position.z();
      m_silseedphi = silseed->get_phi();
      m_silseedeta = silseed->get_eta();
    }

    for (const auto& tpcseed : *tpcseedmap)
    {
      if (!tpcseed)
      {
        continue;
      }

      const auto position = TrackSeedHelper::get_xyz(tpcseed);
      m_tpcseedx = position.x();
      m_tpcseedy = position.y();
      m_tpcseedz = position.z();
      m_tpcseedphi = tpcseed->get_phi();
      m_tpcseedeta = tpcseed->get_eta();

      h_phiDiff[0]->Fill(m_tpcseedphi - m_silseedphi);
      h_etaDiff[0]->Fill(m_tpcseedeta - m_silseedeta);
      h_xDiff[0]->Fill(m_tpcseedx - m_silseedx);
      h_yDiff[0]->Fill(m_tpcseedy - m_silseedy);
      h_zDiff[0]->Fill(m_tpcseedz - m_silseedz);

      if (m_tpcseedeta > 0 && m_silseedeta > 0)
      {
        h_phiDiff[4]->Fill(m_tpcseedphi - m_silseedphi);
        h_etaDiff[4]->Fill(m_tpcseedeta - m_silseedeta);
        h_xDiff[4]->Fill(m_tpcseedx - m_silseedx);
        h_yDiff[4]->Fill(m_tpcseedy - m_silseedy);
        h_zDiff[4]->Fill(m_tpcseedz - m_silseedz);
      }
      else if (m_tpcseedeta < 0 && m_silseedeta < 0)
      {
        h_phiDiff[5]->Fill(m_tpcseedphi - m_silseedphi);
        h_etaDiff[5]->Fill(m_tpcseedeta - m_silseedeta);
        h_xDiff[5]->Fill(m_tpcseedx - m_silseedx);
        h_yDiff[5]->Fill(m_tpcseedy - m_silseedy);
        h_zDiff[5]->Fill(m_tpcseedz - m_silseedz);
      }

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

      if (m_tpcseedeta > 0 && m_silseedeta > 0)
      {
        h_phiDiff[6]->Fill(m_tpcseedphi - m_silseedphi);
        h_etaDiff[6]->Fill(m_tpcseedeta - m_silseedeta);
        h_xDiff[6]->Fill(m_tpcseedx - m_silseedx);
        h_yDiff[6]->Fill(m_tpcseedy - m_silseedy);
        h_zDiff[6]->Fill(m_tpcseedz - m_silseedz);
      }
      else if (m_tpcseedeta < 0 && m_silseedeta < 0)
      {
        h_phiDiff[7]->Fill(m_tpcseedphi - m_silseedphi);
        h_etaDiff[7]->Fill(m_tpcseedeta - m_silseedeta);
        h_xDiff[7]->Fill(m_tpcseedx - m_silseedx);
        h_yDiff[7]->Fill(m_tpcseedy - m_silseedy);
        h_zDiff[7]->Fill(m_tpcseedz - m_silseedz);
      }
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

  std::stringstream stream1, stream2, stream3, stream4;
  stream1 << std::fixed << std::setprecision(2) << m_xcut;
  stream2 << std::fixed << std::setprecision(2) << m_ycut;
  stream3 << std::fixed << std::setprecision(2) << m_etacut;
  stream4 << std::fixed << std::setprecision(2) << m_phicut;

  std::vector<std::string> cutNames = {"", "_xyCut", "_etaCut", "_phiCut", "North", "South", "NorthAllCuts", "SouthAllCuts"};
  std::vector<std::string> cutVals = {"All Track Seeds",
                         std::string("|xdiff| < " + stream1.str() + "cm , |ydiff| < " + stream2.str() + "cm"),
                         std::string("xy cuts and |etadiff| < " + stream3.str()),
                         std::string("xy, eta cuts and |phidiff| < " + stream4.str()),
                         "All Track Seeds (North Only)",
                         "All Track Seeds (South Only)",
                         "North All Cuts (x,y,eta,phi)",
                         "South All Cuts (x,y,eta,phi)"};

  {
    h_crossing = new TH1F(std::string(getHistoPrefix() + "crossing").c_str(),
                          "Track Crossing Value", 1000, -500, 500);
    h_crossing->GetXaxis()->SetTitle("Track Crossing");
    h_crossing->GetYaxis()->SetTitle("Entries");
    hm->registerHisto(h_crossing);
  }
  int i = 0;
  for (const std::string& name : cutNames)
  {
    std::string histoTitle = std::string(cutVals[i]);
    h_phiDiff[i] = new TH1F(std::string(getHistoPrefix() + "phiDiff" + name).c_str(),
                            histoTitle.c_str(), 100, -0.5, 0.5);
    h_phiDiff[i]->GetXaxis()->SetTitle("TPC Seed #phi - Silicon Seed #phi");
    h_phiDiff[i]->GetYaxis()->SetTitle("Entries");
    hm->registerHisto(h_phiDiff[i]);
    i++;
  }
  i = 0;
  for (const std::string& name : cutNames)
  {
    std::string histoTitle = std::string(cutVals[i]);
    h_etaDiff[i] = new TH1F(std::string(getHistoPrefix() + "etaDiff" + name).c_str(),
                            histoTitle.c_str(), 100, -0.1, 0.1);
    h_etaDiff[i]->GetXaxis()->SetTitle("TPC Seed #eta - Silicon Seed #eta");
    h_etaDiff[i]->GetYaxis()->SetTitle("Entries");
    hm->registerHisto(h_etaDiff[i]);
    i++;
  }
  i = 0;
  for (const std::string& name : cutNames)
  {
    std::string histoTitle = std::string(cutVals[i]);
    h_xDiff[i] = new TH1F(std::string(getHistoPrefix() + "xDiff" + name).c_str(),
                          histoTitle.c_str(), 100, -2, 2);
    h_xDiff[i]->GetXaxis()->SetTitle("TPC Seed x - Silicon Seed x [cm]");
    h_xDiff[i]->GetYaxis()->SetTitle("Entries");
    hm->registerHisto(h_xDiff[i]);
    i++;
  }
  i = 0;
  for (const std::string& name : cutNames)
  {
    std::string histoTitle = std::string(cutVals[i]);
    h_yDiff[i] = new TH1F(std::string(getHistoPrefix() + "yDiff" + name).c_str(),
                          histoTitle.c_str(), 100, -2, 2);
    h_yDiff[i]->GetXaxis()->SetTitle("TPC Seed y - Silicon Seed y [cm]");
    h_yDiff[i]->GetYaxis()->SetTitle("Entries");
    hm->registerHisto(h_yDiff[i]);
    i++;
  }
  i = 0;
  for (const std::string& name : cutNames)
  {
    std::string histoTitle = std::string(cutVals[i]);
    h_zDiff[i] = new TH1F(std::string(getHistoPrefix() + "zDiff" + name).c_str(),
                          histoTitle.c_str(), 500, -100, 100);
    h_zDiff[i]->GetXaxis()->SetTitle("TPC Seed z - Silicon Seed z [cm]");
    h_zDiff[i]->GetYaxis()->SetTitle("Entries");
    hm->registerHisto(h_zDiff[i]);
    i++;
  }

  return;
}
