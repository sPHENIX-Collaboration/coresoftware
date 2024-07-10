#include "GlobalQA.h"

// Calo includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

#include <zdcinfo/ZdcReco.h>
#include <zdcinfo/Zdcinfo.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMapv1.h>

#include <qautils/QAHistManagerDef.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <ffarawobjects/Gl1Packet.h>
#include <ffaobjects/EventHeader.h>

#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TSystem.h>

#include <boost/format.hpp>

#include <cassert>
#include <cmath>  // for log10, pow, sqrt, abs, M_PI
#include <cstdint>
#include <iostream>  // for operator<<, endl, basic_...
#include <limits>
#include <map>  // for operator!=, _Rb_tree_con...
#include <string>
#include <utility>  // for pair

GlobalQA::GlobalQA(const std::string &name)
  : SubsysReco(name)
  , detector("HCALIN")
{
  evtcount = 0;
}

GlobalQA::~GlobalQA() = default;

int GlobalQA::Init(PHCompositeNode * /*unused*/)
{
  if (m_debug)
  {
    std::cout << "In GlobalQA::Init" << std::endl;
  }

  createHistos();

  if (m_debug)
  {
    std::cout << "Leaving GlobalQA::Init" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int GlobalQA::process_event(PHCompositeNode *topNode)
{
  _eventcounter++;
  process_towers(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int GlobalQA::process_towers(PHCompositeNode *topNode)
{
  if (m_debug)
  {
    std::cout << _eventcounter << std::endl;
  }

  //---------------------------Event header--------------------------------//
  EventHeader *eventheader = findNode::getClass<EventHeader>(topNode, "EventHeader");
  int event_number = 0;
  if (eventheader)
  {
    if (eventheader->isValid())
    {
      event_number = eventheader->get_EvtSequence();
    }
  }
  else{
    std::cout << "GlobalQA::process_event()  No event header" << std::endl;
  }

  //--------------------------- MBD vertex------------------------------//
  MbdVertexMap *mbdmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");
  MbdVertex *bvertex = nullptr;
  float mbd_zvtx = -999;
  if (mbdmap)
  {
    for (MbdVertexMap::ConstIter mbditer = mbdmap->begin(); mbditer != mbdmap->end(); ++mbditer)
    {
      bvertex = mbditer->second;
    }
    if (bvertex)
    {
      mbd_zvtx = bvertex->get_z();
    }
  }
  h_GlobalQA_mbd_zvtx->Fill(mbd_zvtx);
  h_GlobalQA_mbd_zvtx_wide->Fill(mbd_zvtx);
  if (mbd_zvtx == -999)
  {
    h_GlobalQA_mbd_zvtxq->Fill(0);
  }
  else
  {
    h_GlobalQA_mbd_zvtxq->Fill(1);
  }

  //--------------------------- trigger and GL1-------------------------------//
  bool scaledBits[64] = {false};
  Gl1Packet *gl1PacketInfo = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
  if (!gl1PacketInfo)
  {
    std::cout << PHWHERE << "GlobalQA::process_event: GL1Packet node is missing" << std::endl;
  }
  if (gl1PacketInfo)
  {
    uint64_t triggervec = gl1PacketInfo->getScaledVector();
    for (int i = 0; i < 64; i++)
    {
      bool trig_decision = ((triggervec & 0x1U) == 0x1U);
      scaledBits[i] = trig_decision;
      if (trig_decision)
      {
        h_GlobalQA_triggerVec->Fill(i);
      }
      triggervec = (triggervec >> 1U) & 0xffffffffU;
    }
  }

  // ------------------------------------- ZDC -----------------------------------------//
  {
    Zdcinfo *_zdcinfo = findNode::getClass<Zdcinfo>(topNode, "Zdcinfo");
    float totalzdcsouthcalib = 0.;
    float totalzdcnorthcalib = 0.;
    if (_zdcinfo)
    {
      totalzdcsouthcalib = _zdcinfo->get_zdc_energy(0);
      totalzdcnorthcalib = _zdcinfo->get_zdc_energy(1);
    }
    h_GlobalQA_zdc_zvtx->Fill(999);
    h_GlobalQA_zdc_energy_s->Fill(totalzdcsouthcalib);
    h_GlobalQA_zdc_energy_n->Fill(totalzdcnorthcalib);
  }

  //--------------------------- MBD ----------------------------------------//
  MbdPmtContainer *bbcpmts = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
  if (!bbcpmts)
  {
    std::cout << "GlobalQA::process_event: Could not find MbdPmtContainer," << std::endl;
    // return Fun4AllReturnCodes::ABORTEVENT;
  }

  int hits = 0;
  int hits_n = 0;
  int hits_s = 0;
  int hits_n_t = 0;
  int hits_s_t = 0;
  std::vector<float> time_sum_s = {};
  std::vector<float> time_sum_n = {};
  float sum_s = 0.;
  float sum_n = 0.;
  float sum_s2 = 0.;
  float sum_n2 = 0.;
  float tot_charge_s = 0.;
  float tot_charge_n = 0.;

  float charge_thresh = 0.4;
  if (bbcpmts)
  {
    int nPMTs = bbcpmts->get_npmt();
    for (int i = 0; i < nPMTs; i++)
    {
      MbdPmtHit *mbdpmt = bbcpmts->get_pmt(i);
      float q = mbdpmt->get_q();
      float t = mbdpmt->get_time();
      if (i < 64)
      {
        tot_charge_s += q;
        if (q > charge_thresh)
        {
          hits_s++;
        }
        if (i == 56 || isnan(t))
        {
          continue;
        }
        hits_s_t++;
        time_sum_s.push_back(t);
        sum_s += t;
        sum_s2 += t * t;
      }
      else if (i >= 64)
      {
        tot_charge_n += q;
        if (q > charge_thresh)
        {
          hits_n++;
        }
        if (i == 120 || isnan(t))
        {
          continue;
        }
        hits_n_t++;
        time_sum_n.push_back(t);
        sum_n += t;
        sum_n2 += t * t;
      }

      // float pmtadc = mbdpmt->get_q();
      if (q > 0.4)
      {
        hits++;
      }
    }
  }

  // Calculating the zvtx
  std::sort(time_sum_n.begin(), time_sum_n.end());
  std::sort(time_sum_s.begin(), time_sum_s.end());
  unsigned length_s = time_sum_s.size();
  unsigned length_n = time_sum_n.size();
  float mean_north = 999;
  float mean_south = 999;
  int central_cut = 4;
  float sigma_cut = 1.5;

  if (hits_s_t >= central_cut)
  {
    mean_south = sum_s / static_cast<float>(hits_s_t);
    float rms_s = sqrt(sum_s2 / static_cast<float>(hits_s_t) - TMath::Power(mean_south, 2));
    int nhit_s_center = 0;
    float sum_s_center = 0.;

    for (unsigned int is = 0; is < length_s; is++)
    {
      if (fabs(time_sum_s.at(is) - mean_south) < sigma_cut * rms_s)
      {
        sum_s_center += time_sum_s.at(is);
        nhit_s_center++;
      }
    }

    if (nhit_s_center > 0)
    {
      float mean_south_center = sum_s_center / static_cast<float>(nhit_s_center);
      mean_south = mean_south_center;
    }
  }
  else if (hits_s >= 2 && (hits_s_t >= 1))
  {
    mean_south = sum_s / static_cast<float>(hits_s_t);
  }

  if (hits_n_t >= central_cut)
  {
    mean_north = sum_n / static_cast<float>(hits_n_t);
    float rms_n = sqrt(sum_n2 / static_cast<float>(hits_n_t) - TMath::Power(mean_north, 2));
    int nhit_n_center = 0;
    float sum_n_center = 0.;

    for (unsigned int ino = 0; ino < length_n; ino++)
    {
      if (fabs(time_sum_n.at(ino) - mean_north) < sigma_cut * rms_n)
      {
        sum_n_center += time_sum_n.at(ino);
        nhit_n_center++;
      }
    }

    if (nhit_n_center > 0)
    {
      float mean_north_center = sum_n_center / static_cast<float>(nhit_n_center);
      mean_north = mean_north_center;
    }
  }
  else if (hits_n >= 2 && hits_n_t >= 1)
  {
    mean_north = sum_n / static_cast<float>(hits_n_t);
  }
  float calc_zvtx;
  if (mean_north != 999 && mean_south != 999)
  {
    calc_zvtx = 15 * (mean_south - mean_north);
  }
  else
  {
    calc_zvtx = 999;
  }

  h_GlobalQA_calc_zvtx->Fill(calc_zvtx);
  h_GlobalQA_calc_zvtx_wide->Fill(calc_zvtx);
  h_GlobalQA_mbd_charge_s->Fill(tot_charge_s);
  h_GlobalQA_mbd_charge_n->Fill(tot_charge_n);
  h_GlobalQA_mbd_nhit_s->Fill(hits_s);
  h_GlobalQA_mbd_nhit_n->Fill(hits_n);


  //---------------------------- Trigger / alignment -------------------------------------//
  float leading_cluster_ecore = 0;
  int evtNum_overK = event_number / 1000;

  RawClusterContainer *clusterContainer = findNode::getClass<RawClusterContainer>(topNode, "CLUSTERINFO_CEMC");
  if (clusterContainer)
  {
    RawClusterContainer::ConstRange clusterEnd = clusterContainer->getClusters();
    RawClusterContainer::ConstIterator clusterIter;
    RawClusterContainer::ConstIterator clusterIter2;

    for (clusterIter = clusterEnd.first; clusterIter != clusterEnd.second; clusterIter++)
    {
      RawCluster *recoCluster = clusterIter->second;

      CLHEP::Hep3Vector vertex(0, 0, 0);
      CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*recoCluster, vertex);

      float clusE = E_vec_cluster.mag();
      if (clusE > leading_cluster_ecore)
      {
        leading_cluster_ecore = clusE;
      }
    }
    for (int i = 0; i < 64; i++)
    {
      if (scaledBits[i])
      {
        h_ldClus_trig[i]->Fill(leading_cluster_ecore);
        pr_evtNum_ldClus_trig[i]->Fill(evtNum_overK, leading_cluster_ecore);
        pr_ldClus_trig->Fill(i, leading_cluster_ecore);
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int GlobalQA::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void GlobalQA::createHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // MBD QA
  h_GlobalQA_mbd_zvtxq = new TH1D("h_GlobalQA_mbd_zvtxq", ";Has zvtx?;percentage", 2, -0.5, 1.5);
  h_GlobalQA_mbd_zvtx = new TH1D("h_GlobalQA_mbd_zvtx", ";zvtx [cm]", 100, -50, 50);
  h_GlobalQA_mbd_zvtx_wide = new TH1D("h_GlobalQA_mbd_zvtx_wide", ";zvtx [cm]", 100, -300, 300);
  h_GlobalQA_calc_zvtx = new TH1D("h_GlobalQA_calc_zvtx", ";zvtx [cm]", 100, -50, 50);
  h_GlobalQA_calc_zvtx_wide = new TH1D("h_GlobalQA_calc_zvtx_wide", ";zvtx [cm]", 100, -300, 300);
  h_GlobalQA_mbd_charge_s = new TH1D("h_GlobalQA_mbd_charge_s", ";charge", 100, 0, 10);
  h_GlobalQA_mbd_charge_n = new TH1D("h_GlobalQA_mbd_charge_n", ";charge", 100, 0, 10);
  h_GlobalQA_mbd_nhit_s = new TH1D("h_GlobalQA_mbd_nhit_s", ";nhit", 30, -0.5, 29.5);
  h_GlobalQA_mbd_nhit_n = new TH1D("h_GlobalQA_mbd_nhit_n", ";nhit", 30, -0.5, 29.5);
  hm->registerHisto(h_GlobalQA_mbd_zvtx);
  hm->registerHisto(h_GlobalQA_mbd_zvtxq);
  hm->registerHisto(h_GlobalQA_mbd_zvtx_wide);
  hm->registerHisto(h_GlobalQA_calc_zvtx);
  hm->registerHisto(h_GlobalQA_calc_zvtx_wide);
  hm->registerHisto(h_GlobalQA_mbd_charge_s);
  hm->registerHisto(h_GlobalQA_mbd_charge_n);
  hm->registerHisto(h_GlobalQA_mbd_nhit_s);
  hm->registerHisto(h_GlobalQA_mbd_nhit_n);

  // ZDC QA
  h_GlobalQA_zdc_zvtx = new TH1D("h_GlobalQA_zdc_zvtx", ";zvtx [cm]", 100, -300, 300);
  h_GlobalQA_zdc_energy_s = new TH1D("h_GlobalQA_zdc_energy_s", ";Energy [Gev]", 100, 10, 340);
  h_GlobalQA_zdc_energy_n = new TH1D("h_GlobalQA_zdc_energy_n", ";Energy [Gev]", 100, 10, 340);
  hm->registerHisto(h_GlobalQA_zdc_zvtx);
  hm->registerHisto(h_GlobalQA_zdc_energy_s);
  hm->registerHisto(h_GlobalQA_zdc_energy_n);

  h_GlobalQA_triggerVec = new TH1F("h_GlobalQA_triggerVec", "", 64, 0, 64);
  hm->registerHisto(h_GlobalQA_triggerVec);
  pr_ldClus_trig = new TProfile("pr_GlobalQA_ldClus_trig", "", 64, 0, 64, 0, 10);
  hm->registerHisto(pr_ldClus_trig);

  for (int i = 0; i < 64; i++)
  {
    h_ldClus_trig[i] = new TH1F(boost::str(boost::format("h_GlobalQA_ldClus_trig%d") % i).c_str(), "", 100, 0, 10);
    hm->registerHisto(h_ldClus_trig[i]);
    pr_evtNum_ldClus_trig[i] = new TProfile(boost::str(boost::format("pr_GlobalQA_evtNum_ldClus_trig%d") % i).c_str(), "", 100000, 0, 100000, 0, 10);
    hm->registerHisto(pr_evtNum_ldClus_trig[i]);
  }
}
