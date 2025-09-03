#include "TpcLaserQA.h"

#include <qautils/QAHistManagerDef.h>

#include <tpc/LaserEventInfo.h>
#include <trackbase/LaserCluster.h>
#include <trackbase/LaserClusterContainer.h>
#include <trackbase/TpcDefs.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TH1.h>
#include <TH2.h>

#include <cassert>
#include <cmath>
#include <format>
#include <iostream>  // for basic_ostream, operator<<
#include <map>       // for _Rb_tree_const_iterator
#include <string>
#include <utility>  // for get

TpcLaserQA::TpcLaserQA(const std::string &name)
  : SubsysReco(name)
{
}

int TpcLaserQA::InitRun(PHCompositeNode * /*topNode*/)
{
  createHistos();

  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  m_nLaserEvents = dynamic_cast<TH1 *>(hm->getHisto(std::format("{}nLaserEvents", getHistoPrefix())));

  for (int side = 0; side < 2; side++)
  {
    m_TPCWheel[side] = dynamic_cast<TH2 *>(hm->getHisto(std::format("{}TPCWheel_{}", getHistoPrefix(), (side == 1 ? "North" : "South"))));

    m_nLaserClusters[side] = dynamic_cast<TH1 *>(hm->getHisto(std::format("{}nLaserClusters_{}", getHistoPrefix(), (side == 1 ? "North" : "South"))));
    m_saturation[side] = dynamic_cast<TH2 *>(hm->getHisto(std::format("{}saturation_{}", getHistoPrefix(), (side == 1 ? "North" : "South"))));

    for (int sec = 0; sec < 12; sec++)
    {
      m_sample_R1[side][sec] = dynamic_cast<TH1 *>(hm->getHisto(std::format("{}sample_R1_{}_{}", getHistoPrefix(), (side == 1 ? "North" : "South"), sec)));
      m_sample_R2[side][sec] = dynamic_cast<TH1 *>(hm->getHisto(std::format("{}sample_R2_{}_{}", getHistoPrefix(), (side == 1 ? "North" : "South"), sec)));
      m_sample_R3[side][sec] = dynamic_cast<TH1 *>(hm->getHisto(std::format("{}sample_R3_{}_{}", getHistoPrefix(), (side == 1 ? "North" : "South"), sec)));
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcLaserQA::process_event(PHCompositeNode *topNode)
{
  LaserEventInfo *lei = findNode::getClass<LaserEventInfo>(topNode, "LaserEventInfo");
  if (!lei)
  {
    std::cout << PHWHERE << "LaserEventInfo Node missing, abort." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  LaserClusterContainer *lcc = findNode::getClass<LaserClusterContainer>(topNode, "LASER_CLUSTER");
  if (!lcc)
  {
    std::cout << PHWHERE << "LASER_CLUSTER Node missing, abort." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (!lei->isLaserEvent())
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  /*
  if(lcc->size() < 1000)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  */

  m_nLaserEvents->Fill(1.0);

  int nN = 0;
  int nS = 0;

  auto clusrange = lcc->getClusters();
  for (auto cmitr = clusrange.first;
       cmitr != clusrange.second;
       ++cmitr)
  {
    const auto &[cmkey, cmclus_orig] = *cmitr;
    LaserCluster *cmclus = cmclus_orig;
    int side = TpcDefs::getSide(cmkey);

    if (side)
    {
      nN++;
    }
    else
    {
      nS++;
    }

    const unsigned int nhits = cmclus->getNhits();
    for (unsigned int i = 0; i < nhits; i++)
    {
      float layer = cmclus->getHitLayer(i);
      float hitAdc = cmclus->getHitAdc(i);
      float hitIT = cmclus->getHitIT(i);

      double phi = atan2(cmclus->getHitY(i), cmclus->getHitX(i));
      if (phi < -M_PI / 12.) phi += 2 * M_PI;

      int mod = -1;
      double RValue = -999;
      if (layer >= 7 && layer <= 22)
      {
        mod = 0;
        RValue = 0.4;
      }
      else if (layer >= 23 && layer <= 38)
      {
        mod = 1;
        RValue = 0.6;
      }
      else if (layer >= 39 && layer <= 54)
      {
        mod = 2;
        RValue = 0.8;
      }

      int sector = -1;
      for (int sec = 0; sec < 12; sec++)
      {
        if (phi >= sec * M_PI / 6.0 - M_PI / 12 && phi < sec * M_PI / 6.0 + M_PI / 12)
        {
          sector = sec;
          break;
        }
      }

      if (mod == -1) continue;

      m_TPCWheel[side]->Fill(phi, RValue, hitAdc);
      if (hitAdc > 900) m_saturation[side]->Fill(phi, RValue);

      if (hitAdc > 100)
      {
        if (mod == 0)
        {
          m_sample_R1[side][sector]->Fill(hitIT);
        }
        else if (mod == 1)
        {
          m_sample_R2[side][sector]->Fill(hitIT);
        }
        if (mod == 2)
        {
          m_sample_R3[side][sector]->Fill(hitIT);
        }
      }
    }
  }

  m_nLaserClusters[0]->Fill(nS);
  m_nLaserClusters[1]->Fill(nN);

  return Fun4AllReturnCodes::EVENT_OK;
}

std::string TpcLaserQA::getHistoPrefix() const { return std::string("h_") + Name() + std::string("_"); }  // Define prefix to all histos in HistoManager

void TpcLaserQA::createHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  auto h1 = new TH1F(std::format("{}nLaserEvents", getHistoPrefix()).c_str(), "Number of Laser Events", 1, 0.5, 1.5);
  hm->registerHisto(h1);

  for (int side = 0; side < 2; side++)
  {
    auto h = new TH2F(std::format("{}TPCWheel_{}", getHistoPrefix(), (side == 1 ? "North" : "South")).c_str(), std::format("Laser Hit ADC per event {}", (side == 1 ? "North" : "South")).c_str(), 12, -M_PI / 12., 23. * M_PI / 12., 4, rBinEdges);
    hm->registerHisto(h);

    auto h2 = new TH1F(std::format("{}nLaserClusters_{}", getHistoPrefix(), (side == 1 ? "North" : "South")).c_str(), std::format("Number of Laser Clusters per Event {}", (side == 1 ? "North" : "South")).c_str(), 81, -50, 8050);
    if (side == 1)
    {
      h2->SetLineColor(kRed);
    }
    else
    {
      h2->SetLineColor(kBlue);
    }
    hm->registerHisto(h2);

    auto h3 = new TH2F(std::format("{}saturation_{}", getHistoPrefix(), (side == 1 ? "North" : "South")).c_str(), std::format("Number of Saturated Hits per event {}", (side == 1 ? "North" : "South")).c_str(), 12, -M_PI / 12., 23. * M_PI / 12., 4, rBinEdges);
    hm->registerHisto(h3);

    for (int sec = 0; sec < 12; sec++)
    {
      auto h4 = new TH1F(std::format("{}sample_R1_{}_{}", getHistoPrefix(), (side == 1 ? "North" : "South"), sec).c_str(), std::format("Diffuse Laser Time Sample Sector {} R1 {}", sec, (side == 1 ? "North" : "South")).c_str(), 51, 305.5, 356.5);
      h4->SetLineColor(kRed);
      auto h5 = new TH1F(std::format("{}sample_R2_{}_{}", getHistoPrefix(), (side == 1 ? "North" : "South"), sec).c_str(), std::format("Diffuse Laser Time Sample Sector {} R2 {}", sec, (side == 1 ? "North" : "South")).c_str(), 51, 305.5, 356.5);
      h5->SetLineColor(kBlue);
      auto h6 = new TH1F(std::format("{}sample_R3_{}_{}", getHistoPrefix(), (side == 1 ? "North" : "South"), sec).c_str(), std::format("Diffuse Laser Time Sample Sector {} R3 {}", sec, (side == 1 ? "North" : "South")).c_str(), 51, 305.5, 356.5);
      h6->SetLineColor(kGreen + 2);

      hm->registerHisto(h4);
      hm->registerHisto(h5);
      hm->registerHisto(h6);
    }
  }
}
