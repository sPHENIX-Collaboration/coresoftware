#include "GlobalQA.h"

// Calo includes
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

#include <calotrigger/MinimumBiasInfo.h>
#include <calotrigger/TriggerAnalyzer.h>

#include <centrality/CentralityInfo.h>

#include <mbd/MbdOut.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

#include <zdcinfo/Zdcinfo.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>

#include <qautils/QAHistManagerDef.h>

#include <ffaobjects/RunHeader.h>

#include <ffarawobjects/Gl1Packet.h>

#include <ffamodules/CDBInterface.h>

#include <cdbobjects/CDBTTree.h>  // for CDBTTree

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/RunnumberRange.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE
#include <phool/recoConsts.h>

#include <TH1.h>
#include <TH2.h>
#include <TProfile2D.h>
#include <TSystem.h>

#include <boost/format.hpp>

#include <algorithm>  // for sort
#include <cassert>
#include <cmath>  // for log10, pow, sqrt, abs, M_PI
#include <cstdint>
#include <format>
#include <iostream>  // for operator<<, endl, basic_...
#include <limits>
#include <map>  // for operator!=, _Rb_tree_con...
#include <string>
#include <utility>  // for pair

GlobalQA::GlobalQA(const std::string &name)
  : SubsysReco(name)
  , detector("HCALIN")
{
}

GlobalQA::~GlobalQA()
{
  delete m_trigana;
  delete m_zdccalib;
}

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

int GlobalQA::InitRun(PHCompositeNode *topNode)
{
  RunHeader *runhdr = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if (runhdr)
  {
    m_runnumber = runhdr->get_RunNumber();
  }
  if ((m_runnumber >= RunnumberRange::RUN2AUAU_FIRST && m_runnumber <= RunnumberRange::RUN2AUAU_LAST) ||
      (m_runnumber >= RunnumberRange::RUN3AUAU_FIRST && m_runnumber <= RunnumberRange::RUN3AUAU_LAST))
  {
    m_species = "AuAu";
  }
  else if (m_runnumber >= RunnumberRange::RUN3OO_FIRST && m_runnumber <= RunnumberRange::RUN3OO_LAST)
  {
    m_species = "OO";
  }
  else
  {
    m_species = "pp";
  }
  if (!m_mbbooked)
  {
    m_trigana = new TriggerAnalyzer();
    createMBHistos();
    m_mbbooked = true;
  }
  m_trigdecided = false;
  // ZDC compartment calibration of this run, only when the CDB is set up
  recoConsts *rc = recoConsts::instance();
  if (rc->FlagExist("CDB_GLOBALTAG") && rc->FlagExist("TIMESTAMP"))
  {
    const std::string url = CDBInterface::instance()->getUrl("data_driven_zdc_calib");
    if (url != m_zdccalib_url)
    {
      delete m_zdccalib;
      m_zdccalib = nullptr;
      m_zdccalib_url = url;
      if (!url.empty())
      {
        m_zdccalib = new CDBTTree(url);
        m_zdccalib->LoadCalibrations();
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void GlobalQA::decideTriggers()
{
  m_mbbit = (m_runnumber >= RunnumberRange::RUN2PP_FIRST && m_runnumber <= RunnumberRange::RUN2PP_LAST) ? 10 : 14;
  if (m_species == "AuAu" && m_trigana->getTriggerPrescale(14) < 0)
  {
    m_mbbit = 10;
  }
  m_mbbit_or = (m_trigana->getTriggerPrescale(m_mbbit) > 1000 && m_trigana->getTriggerPrescale(12) >= 0) ? 12 : -1;
  m_zdcbits.clear();
  m_mbdbits.clear();
  for (int b = 0; b < 64; b++)
  {
    const std::string name = m_trigana->getTriggerName(b);
    if (name.find("ZDC") != std::string::npos && m_trigana->getTriggerPrescale(b) >= 0)
    {
      m_zdcbits.push_back(b);
    }
    if (m_species != "pp" && name.find("MBD N&S") != std::string::npos)
    {
      m_mbdbits.push_back(b);
    }
  }
  m_trigdecided = true;
}

int GlobalQA::process_event(PHCompositeNode *topNode)
{
  _eventcounter++;
  if (m_debug)
  {
    std::cout << _eventcounter << std::endl;
  }

  //--------------------------- trigger and GL1-------------------------------//
  Gl1Packet *gl1PacketInfo = findNode::getClass<Gl1Packet>(topNode, 14001);
  if (!gl1PacketInfo)
  {
    gl1PacketInfo = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    if (!gl1PacketInfo)
    {
      std::cout << PHWHERE << "GlobalQA::process_event: GL1Packet node is missing" << std::endl;
    }
  }

  triggervec = 0;
  if (gl1PacketInfo)
  {
    triggervec = gl1PacketInfo->getScaledVector();
    for (int i = 0; i < 64; i++)
    {
      bool trig_decision = ((triggervec & 0x1U) == 0x1U);

      if (trig_decision)
      {
        h_GlobalQA_triggerVec->Fill(i);
      }
      triggervec = (triggervec >> 1U) & 0xffffffffU;
    }
    triggervec = gl1PacketInfo->getScaledVector();
  }

  process_towers(topNode);

  if ( triggervec & mbdtrig )
  {
    process_mbd(topNode);
  }

  if (gl1PacketInfo && m_trigana && m_trigana->decodeTriggers(topNode) == 0)
  {
    if (!m_trigdecided)
    {
      decideTriggers();
    }
    if (m_trigana->didTriggerFire(m_mbbit) || (m_mbbit_or >= 0 && m_trigana->didTriggerFire(m_mbbit_or)))
    {
      process_mb(topNode);
    }
    bool zdctrig = false;
    for (int b : m_zdcbits)
    {
      zdctrig = zdctrig || m_trigana->didTriggerFire(b);
    }
    for (int b : m_mbdbits)
    {
      zdctrig = zdctrig && !m_trigana->checkRawTrigger(b);
    }
    if (zdctrig)
    {
      process_zdctrig(topNode);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void GlobalQA::process_mb(PHCompositeNode *topNode)
{
  MbdOut *mbdout = findNode::getClass<MbdOut>(topNode, "MbdOut");
  m_mbdt0 = std::numeric_limits<float>::quiet_NaN();
  if (mbdout)
  {
    m_mbdt0 = mbdout->get_t0();
    h_GlobalQA_mb_mbd_zvtx->Fill(mbdout->get_zvtx());
    h_GlobalQA_mb_mbd_t0->Fill(m_mbdt0);
  }

  MbdPmtContainer *pmts = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
  if (pmts)
  {
    std::array<float, 2> qsum{0, 0};
    std::array<int, 2> nhit{0, 0};
    const int npmt = pmts->get_npmt();
    for (int i = 0; i < npmt; i++)
    {
      MbdPmtHit *pmt = pmts->get_pmt(i);
      const float q = pmt->get_q();
      const int arm = i / 64;
      if (std::isfinite(pmt->get_tt()))
      {
        h2_GlobalQA_mb_mbd_q_pmt->Fill(i, q);
        if (q > 0)
        {
          h2_GlobalQA_mb_mbd_t_pmt->Fill(i, pmt->get_time());
        }
      }
      qsum[arm] += q;
      if (q > 0.4)
      {
        nhit[arm]++;
      }
    }
    h_GlobalQA_mb_mbd_qavg_s->Fill(qsum[0] / 64.);
    h_GlobalQA_mb_mbd_qavg_n->Fill(qsum[1] / 64.);
    h_GlobalQA_mb_mbd_npmt_s->Fill(nhit[0]);
    h_GlobalQA_mb_mbd_npmt_n->Fill(nhit[1]);
  }

  // centrality: MB events passing the MinimumBiasClassifier selection (O+O: |MBD z| < 60 cm)
  MinimumBiasInfo *mbinfo = findNode::getClass<MinimumBiasInfo>(topNode, "MinimumBiasInfo");
  CentralityInfo *cent = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
  const float centvz = (m_species == "OO") ? 60. : std::numeric_limits<float>::max();
  if (mbinfo && cent && mbdout && mbinfo->isAuAuMinimumBias() && std::abs(mbdout->get_zvtx()) < centvz)
  {
    h_GlobalQA_mb_cent->Fill(100. * cent->get_centile(CentralityInfo::PROP::mbd_NS));
    h_GlobalQA_mb_cent_mbd_qsum->Fill(mbdout->get_q(0) + mbdout->get_q(1));
  }

  // ZDC as in its calibration: heavy ions |z| < 10 cm
  Zdcinfo *zdcinfo = findNode::getClass<Zdcinfo>(topNode, "Zdcinfo");
  if (zdcinfo)
  {
    const std::array<float, 2> ezdc{zdcinfo->get_zdc_energy(0), zdcinfo->get_zdc_energy(1)};
    h_GlobalQA_mb_zdc_energy_s->Fill(ezdc[0]);
    h_GlobalQA_mb_zdc_energy_n->Fill(ezdc[1]);
    h2_GlobalQA_mb_zdc_energy_ns->Fill(ezdc[0], ezdc[1]);
    float vtxz = std::numeric_limits<float>::quiet_NaN();
    GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
    if (vertexmap && !vertexmap->empty() && vertexmap->begin()->second)
    {
      vtxz = vertexmap->begin()->second->get_z();
    }
    TowerInfoContainer *ohcal = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
    if (ohcal && (m_species == "pp" || std::abs(vtxz) < 10))
    {
      float esum = 0;
      for (unsigned int ch = 0; ch < ohcal->size(); ch++)
      {
        TowerInfo *tower = ohcal->get_tower_at_channel(ch);
        if (tower->get_isGood())
        {
          esum += tower->get_energy();
        }
      }
      for (int a = 0; a < 2; a++)
      {
        h2_GlobalQA_mb_zdc_energy_ohcal[a]->Fill(esum, ezdc[a]);
        if (std::isfinite(m_mbdt0))
        {
          h2_GlobalQA_mb_zdc_energy_t0[a]->Fill(m_mbdt0, ezdc[a]);
          h2Profile_GlobalQA_mb_zdc_t0_ohcal[a]->Fill(m_mbdt0, esum, ezdc[a]);
        }
      }
    }
  }

  // sEPD calibrated tiles (nMIP)
  TowerInfoContainer *sepd = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_SEPD");
  if (!sepd)
  {
    sepd = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_EPD");
  }
  if (sepd)
  {
    std::array<float, 2> qsum{0, 0};
    std::array<int, 2> nhit{0, 0};
    const int size = sepd->size();
    for (int ch = 0; ch < size; ch++)
    {
      const float q = sepd->get_tower_at_channel(ch)->get_energy();
      const int arm = TowerInfoDefs::get_epd_arm(sepd->encode_key(ch));
      h2_GlobalQA_mb_sEPD_q_channel->Fill(ch, q);
      qsum[arm] += q;
      if (q > 0.5)
      {
        nhit[arm]++;
      }
    }
    h_GlobalQA_mb_sEPD_qsum_s->Fill(qsum[0]);
    h_GlobalQA_mb_sEPD_qsum_n->Fill(qsum[1]);
    h2_GlobalQA_mb_sEPD_qsum_ns->Fill(qsum[0], qsum[1]);
    h_GlobalQA_mb_sEPD_nhit_s->Fill(nhit[0]);
    h_GlobalQA_mb_sEPD_nhit_n->Fill(nhit[1]);
  }
}

void GlobalQA::process_zdctrig(PHCompositeNode *topNode)
{
  Zdcinfo *zdcinfo = findNode::getClass<Zdcinfo>(topNode, "Zdcinfo");
  if (zdcinfo)
  {
    const float es = zdcinfo->get_zdc_energy(0);
    const float en = zdcinfo->get_zdc_energy(1);
    if (es > 0)
    {
      h_GlobalQA_zdctrig_zdc_energy_s->Fill(es);
    }
    if (en > 0)
    {
      h_GlobalQA_zdctrig_zdc_energy_n->Fill(en);
    }
    h2_GlobalQA_zdctrig_zdc_energy_ns->Fill(es, en);
  }
  // compartment energies: TOWERS_ZDC x data_driven_zdc_calib, in-time pulses
  TowerInfoContainer *zdc = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_ZDC");
  if (!zdc || !m_zdccalib || zdc->size() < 16)
  {
    return;
  }
  for (int a = 0; a < 2; a++)
  {
    for (int m = 0; m < 3; m++)
    {
      const int ch = 8 * a + 2 * m;
      TowerInfo *tower = zdc->get_tower_at_channel(ch);
      const float t = tower->get_time();
      if (t > 5. && t < 9.)
      {
        const float e = tower->get_energy() * m_zdccalib->GetFloatValue(ch, "zdc_calib");
        if (e > 0)
        {
          h_GlobalQA_zdctrig_zdc_module_energy[3 * a + m]->Fill(e);
        }
      }
    }
  }
}

int GlobalQA::process_towers(PHCompositeNode *topNode)
{

  //--------------------------- sEPD ------------------------------//

  TowerInfoContainer *_sepd_towerinfo = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_SEPD");
  unsigned int ntowers = 0;
  if (_sepd_towerinfo)
  {
    ntowers = _sepd_towerinfo->size();
    if (ntowers != 744)
    {
      std::cout << "sEPD container has unexpected size - skipping sEPD!" << std::endl;
      //gSystem->Exit(1);     // commenting out exit so that ZDC and MBD aren't prevented
    }
  }
  if ( (triggervec & mbdtrig) && (ntowers == 744) )  // Any MBD trigger (bits 10-15)
  {
    //--------------------------- sEPD ------------------------------//

    float sepdsouthadcsum = 0.;
    float sepdnorthadcsum = 0.;
    if (_sepd_towerinfo)
    {
      for (unsigned int i = 0; i < ntowers; i++)
      {
        TowerInfo *_tower = _sepd_towerinfo->get_tower_at_channel(i);
        float _e = _tower->get_energy();
        bool isZS = _tower->get_isZS();
        unsigned int key = TowerInfoDefs::encode_epd(i);
        int arm = TowerInfoDefs::get_epd_arm(key);
        float rbin = (float) (TowerInfoDefs::get_epd_rbin(key));
        float phibin = (float) (TowerInfoDefs::get_epd_phibin(key));
        if (!isZS)
        {
          if (arm == 0)
          {
            sepdsouthadcsum += _e;
            h2Profile_GlobalQA_sEPD_tiles_south->Fill(rbin, phibin, _e);
          }
          else if (arm == 1)
          {
            sepdnorthadcsum += _e;
            h2Profile_GlobalQA_sEPD_tiles_north->Fill(rbin, phibin, _e);
          }
        }
      }

      h_GlobalQA_sEPD_adcsum_s->Fill(sepdsouthadcsum);
      h_GlobalQA_sEPD_adcsum_n->Fill(sepdnorthadcsum);
      h2_GlobalQA_sEPD_adcsum_ns->Fill(sepdsouthadcsum, sepdnorthadcsum);
    }
  }

  // ------------------------------------- ZDC
  // -----------------------------------------//

  Zdcinfo *zdcinfo = findNode::getClass<Zdcinfo>(topNode, "Zdcinfo");
  float totalzdcsouthcalib = 0.;
  float totalzdcnorthcalib = 0.;
  if (zdcinfo)
  {
    totalzdcsouthcalib = zdcinfo->get_zdc_energy(0);
    totalzdcnorthcalib = zdcinfo->get_zdc_energy(1);
    zdc_zvtx = zdcinfo->get_zvertex();
    h_GlobalQA_zdc_zvtx->Fill(zdc_zvtx);
    h_GlobalQA_zdc_zvtx_wide->Fill(zdc_zvtx);
    h_GlobalQA_zdc_energy_s->Fill(totalzdcsouthcalib);
    h_GlobalQA_zdc_energy_n->Fill(totalzdcnorthcalib);
  }


  return Fun4AllReturnCodes::EVENT_OK;
}


int GlobalQA::process_mbd(PHCompositeNode *topNode)
{
  //--------------------------- MBD ------------------------------//
  MbdVertexMap *mbdmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");
  MbdVertex *bvertex = nullptr;
  float mbd_zvtx = std::numeric_limits<float>::quiet_NaN();
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
  if (!std::isfinite(mbd_zvtx))
  {
    h_GlobalQA_mbd_zvtxq->SetBinContent(1, h_GlobalQA_mbd_zvtxq->GetBinContent(1) + 1);
  }
  else
  {
    h_GlobalQA_mbd_zvtxq->SetBinContent(2, h_GlobalQA_mbd_zvtxq->GetBinContent(2) + 1);
  }

  MbdPmtContainer *bbcpmts = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
  if (!bbcpmts)
  {
    std::cout << "GlobalQA::process_event: Could not find MbdPmtContainer,"
              << std::endl;
    // return Fun4AllReturnCodes::ABORTEVENT;
  }

  //    int hits = 0;
  int hits_n = 0;
  int hits_s = 0;
  int hits_n_t = 0;
  int hits_s_t = 0;
  std::vector<float> time_sum_s;
  std::vector<float> time_sum_n;
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
        if (i == 56 || std::isnan(t))
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
        if (i == 120 || std::isnan(t))
        {
          continue;
        }
        hits_n_t++;
        time_sum_n.push_back(t);
        sum_n += t;
        sum_n2 += t * t;
      }
      // if (q > charge_thresh) {
      //   hits++;
      // }
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
    float rms_s = std::sqrt((sum_s2 / static_cast<float>(hits_s_t)) - (mean_south * mean_south));
    int nhit_s_center = 0;
    float sum_s_center = 0.;

    for (unsigned int is = 0; is < length_s; is++)
    {
      if (std::fabs(time_sum_s.at(is) - mean_south) < sigma_cut * rms_s)
      {
        sum_s_center += time_sum_s.at(is);
        nhit_s_center++;
      }
    }

    if (nhit_s_center > 0)
    {
      float mean_south_center =
          sum_s_center / static_cast<float>(nhit_s_center);
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
    float rms_n = std::sqrt((sum_n2 / static_cast<float>(hits_n_t)) - (mean_north * mean_north));
    int nhit_n_center = 0;
    float sum_n_center = 0.;

    for (unsigned int ino = 0; ino < length_n; ino++)
    {
      if (std::abs(time_sum_n.at(ino) - mean_north) < sigma_cut * rms_n)
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
  h_GlobalQA_mbd_charge_sum->Fill(tot_charge_s + tot_charge_n);
  h2_GlobalQA_mbd_charge_NS_correlation->Fill(tot_charge_s, tot_charge_n);
  h2_GlobalQA_mbd_nhits_NS_correlation->Fill(hits_s, hits_n);

  return Fun4AllReturnCodes::EVENT_OK;
}


void GlobalQA::createHistos()
{
  auto *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // MBD QA
  h_GlobalQA_mbd_zvtxq = new TH1D("h_GlobalQA_mbd_zvtxq", ";Scaled Trigger 10: MBD Coincidence    Has zvtx?;percentage", 2, -0.5, 1.5);

  h_GlobalQA_mbd_zvtx = new TH1D( "h_GlobalQA_mbd_zvtx", ";Scaled Trigger 10: MBD Coincidence    zvtx [cm]", 100, -50, 50);

  h_GlobalQA_mbd_zvtx_wide = new TH1D( "h_GlobalQA_mbd_zvtx_wide", ";Scaled Trigger 10: MBD Coincidence    zvtx [cm]", 100, -300, 300);

  h_GlobalQA_calc_zvtx = new TH1D("h_GlobalQA_calc_zvtx", ";Scaled Trigger 10: MBD Coincidence    zvtx [cm]", 100, -50, 50);

  h_GlobalQA_calc_zvtx_wide = new TH1D("h_GlobalQA_calc_zvtx_wide", ";Scaled Trigger 10: MBD Coincidence    zvtx [cm]", 100, -300, 300);

  h_GlobalQA_mbd_charge_s = new TH1D("h_GlobalQA_mbd_charge_s", ";Scaled Trigger 10: MBD Coincidence    charge", 1500, 0, 1500);

  h_GlobalQA_mbd_charge_n = new TH1D("h_GlobalQA_mbd_charge_n", ";Scaled Trigger 10: MBD Coincidence    charge", 1500, 0, 1500);

  h_GlobalQA_mbd_nhit_s = new TH1D("h_GlobalQA_mbd_nhit_s", ";Scaled Trigger 10: MBD Coincidence    nhit", 65, -0.5, 64.5);

  h_GlobalQA_mbd_nhit_n = new TH1D("h_GlobalQA_mbd_nhit_n", ";Scaled Trigger 10: MBD Coincidence    nhit", 65, -0.5, 64.5);

  h_GlobalQA_mbd_charge_sum = new TH1F("h_GlobalQA_mbd_charge_sum", " ; MBD Total Charge ; Counts", 2500, 0., 2500);

  h2_GlobalQA_mbd_charge_NS_correlation = new TH2F("h2_GlobalQA_mbd_charge_NS_correlation",
      "MBD Charge Correlation ; Total Charge (South); Total Charge (North)",
      500, 0, 1500, 500, 0, 1500);

  h2_GlobalQA_mbd_nhits_NS_correlation = new TH2F("h2_GlobalQA_mbd_nhits_NS_correlation",
      "MBD Number Of Hits Correlation ; Number Of Hits (South); " "Number Of Hits (North)",
      65, -0.5, 64.5, 65, -0.5, 64.5);

  hm->registerHisto(h_GlobalQA_mbd_zvtx);
  hm->registerHisto(h_GlobalQA_mbd_zvtxq);
  hm->registerHisto(h_GlobalQA_mbd_zvtx_wide);
  hm->registerHisto(h_GlobalQA_calc_zvtx);
  hm->registerHisto(h_GlobalQA_calc_zvtx_wide);
  hm->registerHisto(h_GlobalQA_mbd_charge_s);
  hm->registerHisto(h_GlobalQA_mbd_charge_n);
  hm->registerHisto(h_GlobalQA_mbd_nhit_s);
  hm->registerHisto(h_GlobalQA_mbd_nhit_n);
  hm->registerHisto(h_GlobalQA_mbd_charge_sum);
  hm->registerHisto(h2_GlobalQA_mbd_charge_NS_correlation);
  hm->registerHisto(h2_GlobalQA_mbd_nhits_NS_correlation);

  // --- 6e6 for AuAu, 2e4 for pp
  // sEPD QA
  // h_GlobalQA_sEPD_adcsum_s = new TH1D(
  //     "h_GlobalQA_sEPD_adcsum_s", " ; sEPD ADC sum ; Counts", 200, 0, 6e6);

  // h_GlobalQA_sEPD_adcsum_n = new TH1D(
  //     "h_GlobalQA_sEPD_adcsum_n", " ; sEPD ADC sum ; Counts", 200, 0, 6e6);

  // h2_GlobalQA_sEPD_adcsum_ns =
  //     new TH2D("h2_GlobalQA_sEPD_adcsum_ns",
  //              "sEPD NS ADC sum correlation ; ADC sum (south); ADC sum (north)",
  //              500, 0, 6e6, 500, 0, 6e6);

  h_GlobalQA_sEPD_adcsum_s = new TH1D(
      "h_GlobalQA_sEPD_adcsum_s", " ; sEPD ADC sum ; Counts", 200, 0, 2e4);

  h_GlobalQA_sEPD_adcsum_n = new TH1D(
      "h_GlobalQA_sEPD_adcsum_n", " ; sEPD ADC sum ; Counts", 200, 0, 2e4);

  h2_GlobalQA_sEPD_adcsum_ns =
      new TH2D("h2_GlobalQA_sEPD_adcsum_ns",
               "sEPD NS ADC sum correlation ; ADC sum (south); ADC sum (north)",
               500, 0, 2e4, 500, 0, 2e4);

  h2Profile_GlobalQA_sEPD_tiles_north =
      new TProfile2D("h2Profile_GlobalQA_sEPD_tiles_north",
                     "sEPD Tile Mean Energy (north); #eta; #phi", 16, -0.5,
                     15.5, 24, -0.5, 23.5);

  h2Profile_GlobalQA_sEPD_tiles_south =
      new TProfile2D("h2Profile_GlobalQA_sEPD_tiles_south",
                     "sEPD Tile Mean Energy (south); #eta; #phi", 16, -0.5,
                     15.5, 24, -0.5, 23.5);

  h2_GlobalQA_sEPD_ADC_channel_south =
      new TH2D("h2_GlobalQA_sEPD_ADC_channel_south",
               "h2_GlobalQA_sEPD_ADC_channel_south ; #eta; #phi", 16, -0.5,
               15.5, 24, -0.5, 23.5);
  h2_GlobalQA_sEPD_ADC_channel_north =
      new TH2D("h2_GlobalQA_sEPD_ADC_channel_north",
               "h2_GlobalQA_sEPD_ADC_channel_north ; #eta; #phi", 16, -0.5,
               15.5, 24, -0.5, 23.5);

  hm->registerHisto(h_GlobalQA_sEPD_adcsum_s);
  hm->registerHisto(h_GlobalQA_sEPD_adcsum_n);
  hm->registerHisto(h2_GlobalQA_sEPD_adcsum_ns);
  hm->registerHisto(h2Profile_GlobalQA_sEPD_tiles_south);
  hm->registerHisto(h2Profile_GlobalQA_sEPD_tiles_north);
  hm->registerHisto(h2_GlobalQA_sEPD_ADC_channel_south);
  hm->registerHisto(h2_GlobalQA_sEPD_ADC_channel_north);

  // ZDC QA
  h_GlobalQA_zdc_zvtx = new TH1D(
      "h_GlobalQA_zdc_zvtx", ";Scaled Trigger 3: ZDC Coincidence    zvtx [cm]",
      100, -1000, 1000);
  h_GlobalQA_zdc_zvtx_wide = new TH1D(
      "h_GlobalQA_zdc_zvtx_wide",
      ";Scaled Trigger 3: ZDC Coincidence    zvtx [cm]", 100, -2000, 2000);
  h_GlobalQA_zdc_energy_s = new TH1D(
      "h_GlobalQA_zdc_energy_s",
      ";Scaled Trigger 3: ZDC Coincidence    Energy [Gev]", 100, 10, 510);
  h_GlobalQA_zdc_energy_n = new TH1D(
      "h_GlobalQA_zdc_energy_n",
      ";Scaled Trigger 3: ZDC Coincidence    Energy [Gev]", 100, 10, 510);
  hm->registerHisto(h_GlobalQA_zdc_zvtx);
  hm->registerHisto(h_GlobalQA_zdc_zvtx_wide);
  hm->registerHisto(h_GlobalQA_zdc_energy_s);
  hm->registerHisto(h_GlobalQA_zdc_energy_n);

  h_GlobalQA_triggerVec = new TH1F("h_GlobalQA_triggerVec", "", 65, -0.5, 64.5);
  h_GlobalQA_triggerVec->SetDirectory(nullptr);
  hm->registerHisto(h_GlobalQA_triggerVec);
}

void GlobalQA::createMBHistos()
{
  auto *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  const bool pp = (m_species == "pp");
  const bool oo = (m_species == "OO");
  float qmax = 200.;
  float zdcmax = 6000.;
  float sepdqmax = 50.;
  float sepdsummax = 5000.;
  float ohmax = 400.;
  if (pp)
  {
    qmax = 25.;
    zdcmax = 1000.;
    sepdqmax = 15.;
    sepdsummax = 300.;
    ohmax = 20.;
  }
  else if (oo)
  {
    qmax = 50.;
    zdcmax = 2500.;
    sepdqmax = 30.;
    sepdsummax = 1000.;
    ohmax = 40.;
  }
  const std::string zsel = pp ? "MB" : "MB |z| < 10 cm";

  h2_GlobalQA_mb_mbd_q_pmt = new TH2F("h2_GlobalQA_mb_mbd_q_pmt", "MBD PMT charge, valid TDC, MB;PMT;q", 128, 0, 128, std::lround((qmax + 1) / 0.05), -1, qmax);
  h2_GlobalQA_mb_mbd_t_pmt = new TH2F("h2_GlobalQA_mb_mbd_t_pmt", "MBD PMT time q>0, MB;PMT;t (ns)", 128, 0, 128, 500, -25, 25);
  h_GlobalQA_mb_mbd_qavg_s = new TH1D("h_GlobalQA_mb_mbd_qavg_s", "MBD south average PMT charge, MB;<q>;events", 500, 0, qmax);
  h_GlobalQA_mb_mbd_qavg_n = new TH1D("h_GlobalQA_mb_mbd_qavg_n", "MBD north average PMT charge, MB;<q>;events", 500, 0, qmax);
  h_GlobalQA_mb_mbd_npmt_s = new TH1D("h_GlobalQA_mb_mbd_npmt_s", "MBD south PMTs q>0.4, MB;N_{PMT};events", 65, -0.5, 64.5);
  h_GlobalQA_mb_mbd_npmt_n = new TH1D("h_GlobalQA_mb_mbd_npmt_n", "MBD north PMTs q>0.4, MB;N_{PMT};events", 65, -0.5, 64.5);
  h_GlobalQA_mb_mbd_zvtx = new TH1D("h_GlobalQA_mb_mbd_zvtx", "MBD z-vertex, MB;z (cm);events", 300, -150, 150);
  h_GlobalQA_mb_mbd_t0 = new TH1D("h_GlobalQA_mb_mbd_t0", "MBD t0, MB;t0 (ns);events", 200, -20, 20);
  h_GlobalQA_mb_cent = new TH1D("h_GlobalQA_mb_cent", "centrality, MB trigger and MB event selection;centrality (%);events", 100, 0, 100);
  h_GlobalQA_mb_cent_mbd_qsum = new TH1D("h_GlobalQA_mb_cent_mbd_qsum", "MBD charge sum, MB trigger and MB event selection;MBD charge sum (N+S);events", 1000, 0, oo ? 500. : 3000.);
  for (TH1 *h : {static_cast<TH1 *>(h2_GlobalQA_mb_mbd_q_pmt), static_cast<TH1 *>(h2_GlobalQA_mb_mbd_t_pmt), h_GlobalQA_mb_mbd_qavg_s, h_GlobalQA_mb_mbd_qavg_n,
                 h_GlobalQA_mb_mbd_npmt_s, h_GlobalQA_mb_mbd_npmt_n, h_GlobalQA_mb_mbd_zvtx, h_GlobalQA_mb_mbd_t0, h_GlobalQA_mb_cent, h_GlobalQA_mb_cent_mbd_qsum})
  {
    hm->registerHisto(h);
  }

  h_GlobalQA_mb_zdc_energy_s = new TH1D("h_GlobalQA_mb_zdc_energy_s", "ZDC south energy, MB;E (GeV);events", std::lround(zdcmax), 0, zdcmax);
  h_GlobalQA_mb_zdc_energy_n = new TH1D("h_GlobalQA_mb_zdc_energy_n", "ZDC north energy, MB;E (GeV);events", std::lround(zdcmax), 0, zdcmax);
  h2_GlobalQA_mb_zdc_energy_ns = new TH2F("h2_GlobalQA_mb_zdc_energy_ns", "ZDC north vs south energy, MB;E_{south} (GeV);E_{north} (GeV)", 200, 0, zdcmax, 200, 0, zdcmax);
  hm->registerHisto(h_GlobalQA_mb_zdc_energy_s);
  hm->registerHisto(h_GlobalQA_mb_zdc_energy_n);
  hm->registerHisto(h2_GlobalQA_mb_zdc_energy_ns);
  for (int a = 0; a < 2; a++)
  {
    const std::string arm = a ? "north" : "south";
    const std::string sn = a ? "n" : "s";
    h2_GlobalQA_mb_zdc_energy_ohcal[a] = new TH2F(std::format("h2_GlobalQA_mb_zdc_energy_ohcal_{}", sn).c_str(),
                                                  std::format("ZDC {} energy vs OHCal energy sum, {};#SigmaE_{{OHCal}} (GeV);E (GeV)", arm, zsel).c_str(), 200, 0, ohmax, 250, 0, zdcmax);
    h2_GlobalQA_mb_zdc_energy_t0[a] = new TH2F(std::format("h2_GlobalQA_mb_zdc_energy_t0_{}", sn).c_str(),
                                               std::format("ZDC {} energy vs MBD t0, {};t0 (ns);E (GeV)", arm, zsel).c_str(), 40, -10, 10, 250, 0, zdcmax);
    h2Profile_GlobalQA_mb_zdc_t0_ohcal[a] = new TProfile2D(std::format("h2Profile_GlobalQA_mb_zdc_t0_ohcal_{}", sn).c_str(),
                                                           std::format("ZDC {} mean energy, {};t0 (ns);#SigmaE_{{OHCal}} (GeV)", arm, zsel).c_str(), 40, -10, 10, 40, 0, ohmax);
    hm->registerHisto(h2_GlobalQA_mb_zdc_energy_ohcal[a]);
    hm->registerHisto(h2_GlobalQA_mb_zdc_energy_t0[a]);
    hm->registerHisto(h2Profile_GlobalQA_mb_zdc_t0_ohcal[a]);
  }

  const float sepdqbinw = (m_species == "AuAu") ? 0.1 : 0.05;
  h2_GlobalQA_mb_sEPD_q_channel = new TH2F("h2_GlobalQA_mb_sEPD_q_channel", "sEPD tile charge, MB;channel;q (nMIP)", 744, -0.5, 743.5, std::lround((sepdqmax + 1) / sepdqbinw), -1, sepdqmax);
  h_GlobalQA_mb_sEPD_qsum_s = new TH1D("h_GlobalQA_mb_sEPD_qsum_s", "sEPD south charge sum per event, MB;#Sigmaq (nMIP);events", 500, 0, sepdsummax);
  h_GlobalQA_mb_sEPD_qsum_n = new TH1D("h_GlobalQA_mb_sEPD_qsum_n", "sEPD north charge sum per event, MB;#Sigmaq (nMIP);events", 500, 0, sepdsummax);
  h2_GlobalQA_mb_sEPD_qsum_ns = new TH2F("h2_GlobalQA_mb_sEPD_qsum_ns", "sEPD north vs south charge sum, MB;#Sigmaq_{south} (nMIP);#Sigmaq_{north} (nMIP)", 200, 0, sepdsummax, 200, 0, sepdsummax);
  h_GlobalQA_mb_sEPD_nhit_s = new TH1D("h_GlobalQA_mb_sEPD_nhit_s", "sEPD south tiles with q > 0.5 nMIP, MB;N_{tile};events", 373, -0.5, 372.5);
  h_GlobalQA_mb_sEPD_nhit_n = new TH1D("h_GlobalQA_mb_sEPD_nhit_n", "sEPD north tiles with q > 0.5 nMIP, MB;N_{tile};events", 373, -0.5, 372.5);
  for (TH1 *h : {static_cast<TH1 *>(h2_GlobalQA_mb_sEPD_q_channel), h_GlobalQA_mb_sEPD_qsum_s, h_GlobalQA_mb_sEPD_qsum_n, static_cast<TH1 *>(h2_GlobalQA_mb_sEPD_qsum_ns),
                 h_GlobalQA_mb_sEPD_nhit_s, h_GlobalQA_mb_sEPD_nhit_n})
  {
    hm->registerHisto(h);
  }

  const int nzdc = std::lround((zdcmax + 50.) / 2.5);
  h_GlobalQA_zdctrig_zdc_energy_s = new TH1D("h_GlobalQA_zdctrig_zdc_energy_s", "ZDC south energy, ZDC trigger;E (GeV);events", nzdc, -50., zdcmax);
  h_GlobalQA_zdctrig_zdc_energy_n = new TH1D("h_GlobalQA_zdctrig_zdc_energy_n", "ZDC north energy, ZDC trigger;E (GeV);events", nzdc, -50., zdcmax);
  h2_GlobalQA_zdctrig_zdc_energy_ns = new TH2F("h2_GlobalQA_zdctrig_zdc_energy_ns", "ZDC north vs south energy, ZDC trigger;E_{south} (GeV);E_{north} (GeV)", 200, 0, zdcmax, 200, 0, zdcmax);
  hm->registerHisto(h_GlobalQA_zdctrig_zdc_energy_s);
  hm->registerHisto(h_GlobalQA_zdctrig_zdc_energy_n);
  hm->registerHisto(h2_GlobalQA_zdctrig_zdc_energy_ns);
  for (int a = 0; a < 2; a++)
  {
    const std::string arm = a ? "north" : "south";
    for (int m = 0; m < 3; m++)
    {
      TH1 *h = new TH1D(std::format("h_GlobalQA_zdctrig_zdc_module_energy_{}_{}", a ? "n" : "s", m + 1).c_str(),
                        std::format("ZDC {} compartment {} energy, ZDC trigger;E (GeV);events", arm, m + 1).c_str(), nzdc, -50., zdcmax);
      h_GlobalQA_zdctrig_zdc_module_energy[3 * a + m] = h;
      hm->registerHisto(h);
    }
  }
}
