#include "TpcSiliconQA.h"

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
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
TpcSiliconQA::TpcSiliconQA(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int TpcSiliconQA::InitRun(PHCompositeNode * /*topNode*/)
{
  createHistos();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcSiliconQA::process_event(PHCompositeNode *topNode)
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
  auto trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  if (!trackmap)
  {
    std::cout << "Track map not found, aborting event" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT; 
  }

  for (const auto& [key, track] : *trackmap)
  {
    if (!track)
    {
      continue;
    }

    m_crossing = (float) track->get_crossing();
    h_crossing->Fill(m_crossing);

    m_silseedx = std::numeric_limits<float>::quiet_NaN();
    m_silseedy = std::numeric_limits<float>::quiet_NaN();
    m_silseedz = std::numeric_limits<float>::quiet_NaN();
    m_silseedphi = std::numeric_limits<float>::quiet_NaN();
    m_silseedeta = std::numeric_limits<float>::quiet_NaN();
    m_tpcseedx = std::numeric_limits<float>::quiet_NaN();
    m_tpcseedy = std::numeric_limits<float>::quiet_NaN();
    m_tpcseedz = std::numeric_limits<float>::quiet_NaN();
    m_tpcseedphi = std::numeric_limits<float>::quiet_NaN();
    m_tpcseedeta = std::numeric_limits<float>::quiet_NaN();
 
    auto tpcseed = track->get_tpc_seed();
    auto silseed = track->get_silicon_seed();
    if (silseed && tpcseed)
    {
      h_trackMatch->Fill(1);
      m_silseedx = silseed->get_x();
      m_silseedy = silseed->get_y();
      m_silseedz = silseed->get_z();
      m_silseedphi = silseed->get_phi();
      m_silseedeta = silseed->get_eta();
      m_tpcseedx = tpcseed->get_x();
      m_tpcseedy = tpcseed->get_y();
      m_tpcseedz = tpcseed->get_z();
      m_tpcseedphi = tpcseed->get_phi();
      m_tpcseedeta = tpcseed->get_eta();
  
      h_phiDiff->Fill(m_tpcseedphi - m_silseedphi);
      h_etaDiff->Fill(m_tpcseedeta - m_silseedeta);
      h_xDiff->Fill(m_tpcseedx - m_silseedx);
      h_yDiff->Fill(m_tpcseedy - m_silseedy);
      h_zDiff->Fill(m_tpcseedz - m_silseedz);
    }
    else
    {
      h_trackMatch->Fill(0);
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
 
  {
  h_crossing = new TH1F(std::string(getHistoPrefix() + "crossing").c_str(),
                      "Track Crossing Value", 1000, -500, 500);
  h_crossing->GetXaxis()->SetTitle("Track Crossing");
  h_crossing->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_crossing);
  }
  { 
  h_trackMatch = new TH1F(std::string(getHistoPrefix() + "trackMatch").c_str(),
                      "TPC and Silicon Seed Exist", 2, -0.5, 1.5);
  h_trackMatch->GetXaxis()->SetTitle("1 - TPC+Sil Seed, 0 - Missing TPC and/or Sil");
  h_trackMatch->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_trackMatch);
  }
  {
  h_phiDiff = new TH1F(std::string(getHistoPrefix() + "phiDiff").c_str(),
                      "TPC-Silicon #phi Difference", 100, -0.5, 0.5);
  h_phiDiff->GetXaxis()->SetTitle("TPC Seed #phi - Silicon Seed #phi");
  h_phiDiff->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_phiDiff);
  }
  {
  h_etaDiff = new TH1F(std::string(getHistoPrefix() + "etaDiff").c_str(),
                      "TPC-Silicon #eta Difference", 100, -0.1, 0.1);
  h_etaDiff->GetXaxis()->SetTitle("TPC Seed #eta - Silicon Seed #eta");
  h_etaDiff->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_etaDiff);
  }
  {
  h_xDiff = new TH1F(std::string(getHistoPrefix() + "xDiff").c_str(),
                      "TPC-Silicon x Difference", 100, -2, 2);
  h_xDiff->GetXaxis()->SetTitle("TPC Seed x - Silicon Seed x [cm]");
  h_xDiff->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_xDiff);
  }
  {
  h_yDiff = new TH1F(std::string(getHistoPrefix() + "yDiff").c_str(),
                      "TPC-Silicon y Difference", 100, -2, 2);
  h_yDiff->GetXaxis()->SetTitle("TPC Seed y - Silicon Seed y [cm]");
  h_yDiff->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_yDiff);
  }
  {
  h_zDiff = new TH1F(std::string(getHistoPrefix() + "zDiff").c_str(),
                      "TPC-Silicon z Difference", 500, -100, 100);
  h_zDiff->GetXaxis()->SetTitle("TPC Seed z - Silicon Seed z [cm]");
  h_zDiff->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_zDiff);
  } 

  return;
}
