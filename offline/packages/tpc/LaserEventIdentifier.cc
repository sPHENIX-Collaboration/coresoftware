#include "LaserEventIdentifier.h"

#include "LaserEventInfo.h"
#include "LaserEventInfov2.h"

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>  // for hitkey, getLayer
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <ffarawobjects/Gl1Packet.h>

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>  // for PHIODataNode
#include <phool/PHNode.h>        // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>

#include <Acts/Definitions/Units.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <TF1.h>
#include <TFile.h>

LaserEventIdentifier::LaserEventIdentifier(const std::string &name)
  : SubsysReco(name)
{
}

int LaserEventIdentifier::InitRun(PHCompositeNode *topNode)
{
  // get node containing the digitized hits
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!m_hits)
  {
    std::cout << PHWHERE << "ERROR: Can't find node TRKR_HITSET" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_geom_container =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!m_geom_container)
  {
    std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_tGeometry = findNode::getClass<ActsGeometry>(topNode,
                                                 "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << PHWHERE
              << "ActsGeometry not found on node tree. Exiting"
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode =
      dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode("TRKR");
    dstNode->addNode(DetNode);
  }

  LaserEventInfo *laserEventInfo = new LaserEventInfov2();

  PHIODataNode<PHObject> *laserEventInfoNode = new PHIODataNode<PHObject>(laserEventInfo, "LaserEventInfo", "PHObject");
  DetNode->addNode(laserEventInfoNode);

  if (m_debug)
  {
    m_debugFile = new TFile(m_debugFileName.c_str(), "RECREATE");
  }
  float timeHistMax = m_time_samples_max;
  timeHistMax -= 0.5;
  m_itHist_0 = new TH1I("m_itHist_0", "side 0;it", m_time_samples_max, -0.5, timeHistMax);
  m_itHist_1 = new TH1I("m_itHist_1", "side 1;it", m_time_samples_max, -0.5, timeHistMax);

  if (m_debug)
  {
    m_hitTree = new TTree("hitTree", "hitTree");
    m_hitTree->Branch("itHist_0", &m_itHist_0);
    m_hitTree->Branch("itHist_1", &m_itHist_1);
    m_hitTree->Branch("isLaserEvent", &isLaserEvent);
    m_hitTree->Branch("isGl1LaserEvent", &isGl1LaserEvent);
    m_hitTree->Branch("isGl1LaserPileupEvent", &isGl1LaserPileupEvent);
    m_hitTree->Branch("peakSample_0", &peakSample0);
    m_hitTree->Branch("peakSample_1", &peakSample1);
    m_hitTree->Branch("peakWidth_0", &peakWidth0);
    m_hitTree->Branch("peakWidth_1", &peakWidth1);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int LaserEventIdentifier::process_event(PHCompositeNode *topNode)
{
  m_laserEventInfo = findNode::getClass<LaserEventInfo>(topNode, "LaserEventInfo");
  if (!m_laserEventInfo)
  {
    std::cout << "no laser event info node" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  Gl1Packet *gl1pkt = findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");
  if (!gl1pkt)
  {
    std::cout << "no GL1RAWHIT node" << std::endl;
    m_laserEventInfo->setIsGl1LaserEvent(false);
    m_laserEventInfo->setIsGl1LaserPileupEvent(false);
    //return Fun4AllReturnCodes::ABORTRUN;
  }
  else if(m_runnumber > 66153)
  {
    if ((gl1pkt->getGTMAllBusyVector() & (1<<14)) == 0)
    {
      m_laserEventInfo->setIsGl1LaserEvent(true);
      m_laserEventInfo->setIsGl1LaserPileupEvent(false);
      isGl1LaserEvent = true;
      isGl1LaserPileupEvent = false;
      prev_BCO = gl1pkt->getBCO();
    }
    else if ((gl1pkt->getBCO() - prev_BCO) < 350.0/30*16)
    {
      m_laserEventInfo->setIsGl1LaserEvent(false);
      m_laserEventInfo->setIsGl1LaserPileupEvent(true);
      isGl1LaserEvent = false;
      isGl1LaserPileupEvent = true;
      prev_BCO = 0;
    }
    else
    {
      m_laserEventInfo->setIsGl1LaserEvent(false);
      m_laserEventInfo->setIsGl1LaserPileupEvent(false);
      isGl1LaserEvent = false;
      isGl1LaserPileupEvent = false;
      prev_BCO = 0;
    }
  }
  else
  {
    m_laserEventInfo->setIsGl1LaserEvent(false);
    m_laserEventInfo->setIsGl1LaserPileupEvent(false);
  }

  TrkrHitSetContainer::ConstRange hitsetrange = m_hits->getHitSets(TrkrDefs::TrkrId::tpcId);
  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr)
  {
    TrkrHitSet *hitset = hitsetitr->second;
    int side = TpcDefs::getSide(hitsetitr->first);

    TrkrHitSet::ConstRange hitrangei = hitset->getHits();
    for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
         hitr != hitrangei.second;
         ++hitr)
    {
      int it = TpcDefs::getTBin(hitr->first);

      if (side == 0)
      {
        m_itHist_0->Fill(it);
      }
      else
      {
        m_itHist_1->Fill(it);
      }
    }
  }

  double itMeanContent_0 = 0.0;
  double itMeanContent_1 = 0.0;

  for (int i = 1; i <= m_time_samples_max; i++)
  {
    itMeanContent_0 += m_itHist_0->GetBinContent(i);
    itMeanContent_1 += m_itHist_1->GetBinContent(i);
  }

  itMeanContent_0 = 1.0 * itMeanContent_0 / m_time_samples_max;
  itMeanContent_1 = 1.0 * itMeanContent_1 / m_time_samples_max;

  m_itHist_0->GetXaxis()->SetRange(320, m_time_samples_max - 0.5);
  double itMax_0 = m_itHist_0->GetBinCenter(m_itHist_0->GetMaximumBin());
  double itMaxContent_0 = m_itHist_0->GetMaximum();
  m_itHist_0->GetXaxis()->SetRange(0, 0);

  m_itHist_1->GetXaxis()->SetRange(320, m_time_samples_max - 0.5);
  double itMax_1 = m_itHist_1->GetBinCenter(m_itHist_1->GetMaximumBin());
  double itMaxContent_1 = m_itHist_1->GetMaximum();
  m_itHist_1->GetXaxis()->SetRange(0, 0);

  auto f0 = std::make_unique<TF1>("f0", "gausn(0)");
  f0->SetParameters(itMaxContent_0, itMax_0, 1);
  f0->SetParLimits(1, itMax_0 - 2, itMax_0 + 2);
  m_itHist_0->Fit(f0.get(), "Bq0");

  auto f1 = std::make_unique<TF1>("f1", "gausn(0)");
  f1->SetParameters(itMaxContent_1, itMax_1, 1);
  f1->SetParLimits(1, itMax_1 - 2, itMax_1 + 2);
  m_itHist_1->Fit(f1.get(), "Bq0");


  if ((itMaxContent_0 / itMeanContent_0 >= 7 && itMaxContent_0 > 1000) || (itMaxContent_1 / itMeanContent_1 >= 7 && itMaxContent_1 > 1000))
  {
    m_laserEventInfo->setIsLaserEvent(true);
    m_laserEventInfo->setPeakSample(false, (int) itMax_0);
    m_laserEventInfo->setPeakWidth(false, f0->GetParameter(2));
    m_laserEventInfo->setPeakSample(true, (int) itMax_1);
    m_laserEventInfo->setPeakWidth(true, f1->GetParameter(2));

    isLaserEvent = true;
    peakSample0 = (int) itMax_0;
    peakSample1 = (int) itMax_1;
    peakWidth0 = f0->GetParameter(2);
    peakWidth1 = f1->GetParameter(2);
  }
  else
  {
    m_laserEventInfo->setIsLaserEvent(false);

    isLaserEvent = false;
    peakSample0 = -999;
    peakSample1 = -999;
    peakWidth0 = -999;
    peakWidth1 = -999;
  }

  if (m_debug)
  {
    m_hitTree->Fill();
  }
  

  return Fun4AllReturnCodes::EVENT_OK;
}

int LaserEventIdentifier::ResetEvent(PHCompositeNode * /*topNode*/)
{
  m_itHist_0->Reset();
  m_itHist_1->Reset();

  return Fun4AllReturnCodes::EVENT_OK;
}

int LaserEventIdentifier::End(PHCompositeNode * /*topNode*/)
{
  if (m_debug)
  {
    m_debugFile->cd();
    m_hitTree->Write();
    m_debugFile->Close();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
