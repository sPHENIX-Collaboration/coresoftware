#include "TpcRawHitQA.h"

#include <qautils/QAHistManagerDef.h>
#include <qautils/QAUtil.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TH2.h>
#include <TMath.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <cassert>
#include <cmath>

//____________________________________________________________________________..
TpcRawHitQA::TpcRawHitQA(const std::string &name)
  : SubsysReco(name)
{
  M.setMapNames("AutoPad-R1-RevA.sch.ChannelMapping.csv", "AutoPad-R2-RevA-Pads.sch.ChannelMapping.csv", "AutoPad-R3-RevA.sch.ChannelMapping.csv");
}

//____________________________________________________________________________..
int TpcRawHitQA::InitRun(PHCompositeNode *topNode)
{
  createHistos();

  rawhitcont = findNode::getClass<TpcRawHitContainer>(topNode, "TPCRAWHIT");

  if (!rawhitcont)
  {
    std::cout << PHWHERE << "Missing TpcRawHitContainer node!!!" << std::endl;
  }

  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  for (int s = 0; s < 24; s++)
  {
    h_nhits_sectors[s] = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nhits_sec" + std::to_string(s)).c_str()));
    h_nhits_sectors_fees[s] = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "nhits_sec" + std::to_string(s) + "_fees").c_str()));
    for (int r = 0; r < 3; r++)
    {
      h_nhits_sam[s][r] = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nhits_sample_sec" + std::to_string(s) + "_R" + std::to_string(r)).c_str()));
      h_adc[s][r] = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "adc_sec" + std::to_string(s) + "_R" + std::to_string(r)).c_str()));
    }
  }

  h_bco = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "bco").c_str()));
  h_xy_N = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "xyPos_North").c_str()));
  h_xy_S = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "xyPos_South").c_str()));

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcRawHitQA::process_event(PHCompositeNode * /*unused*/)
{
  // std::vector < TpcRawHit* > hits;
  std::vector<uint64_t> bcos;
  std::vector<int> sectors;
  std::vector<uint16_t> fees;

  // hits.clear();
  bcos.clear();
  sectors.clear();
  fees.clear();

  unsigned int raw_hit_num = 0;
  if (rawhitcont)
  {
    raw_hit_num = rawhitcont->get_nhits();
    for (unsigned int i = 0; i < raw_hit_num; i++)
    {
      auto hit = rawhitcont->get_hit(i);
      uint64_t bco = hit->get_gtm_bco();
      uint16_t sam = hit->get_samples();
      int32_t packet_id = hit->get_packetid();
      int ep = (packet_id - 4000) % 10;
      int sector = (packet_id - 4000 - ep) / 10;
      uint16_t fee = hit->get_fee();
      int channel = hit->get_channel();
      int region = FEE_R[fee];
      int feeM = FEE_map[fee];
      if (FEE_R[fee] == 2)
      {
        feeM += 6;
      }
      if (FEE_R[fee] == 3)
      {
        feeM += 14;
      }
      double R = M.getR(feeM, channel);
      double phi = 0;
      if (sector < 12)  // NS
      {
        phi = M.getPhi(feeM, channel) + (sector) *M_PI / 6;
      }
      else if (sector >= 12)  // SS
      {
        phi = M.getPhi(feeM, channel) + (18 - sector) * M_PI / 6;
      }

      std::vector<int> values;
      values.reserve(sam);
      for (int sampleN = 0; sampleN < sam; sampleN++)
      {
        values.push_back((int) hit->get_adc(sampleN));
      }
      std::sort(values.begin(), values.end());
      size_t size = values.size();
      float median;
      if (size % 2 == 0)
      {
        median = (values[size / 2 - 1] + values[size / 2]) / 2.0;
      }
      else
      {
        median = values[size / 2];
      }
      std::vector<int> selectedValues;
      for (int value : values)
      {
        if (value >= median - 40 && value <= median + 40)
        {
          selectedValues.push_back(value);
        }
      }

      float stdDev = 3;
      if (selectedValues.size() > 0)
      {
        float mean = 0.0;
        for (int value : selectedValues)
        {
          mean += value;
        }
        mean /= selectedValues.size();
        float sumSquares = 0.0;
        for (int value : selectedValues)
        {
          float diff = value - mean;
          sumSquares += std::pow(diff, 2);
        }
        float variance = sumSquares / selectedValues.size();
        stdDev = std::sqrt(variance);
      }

      for (int sampleN = 0; sampleN < sam; sampleN++)
      {
        float adc = hit->get_adc(sampleN);
        if (adc - median <= (std::max(5 * stdDev, (float) 20.)))
        {
          continue;
        }

        if (sector < 12)
        {
          h_xy_N->Fill(R * cos(phi), R * sin(phi));
        }
        if (sector >= 12)
        {
          h_xy_S->Fill(R * cos(phi), R * sin(phi));
        }
        h_nhits_sam[sector][region - 1]->Fill(sampleN);
        h_adc[sector][region - 1]->Fill(adc);
      }

      // hits.push_back( hit );
      bcos.push_back(bco);
      sectors.push_back(sector);
      fees.push_back(fee);
    }
  }

  // if no raw hit is found, skip this event
  if (raw_hit_num == 0)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  float nhit_sectors[24] = {0};
  float nhit_sectors_fees[24][26] = {{0}};

  for (int i = 0; i < (int) raw_hit_num; i++)
  {
    h_bco->Fill(bcos[i]);
    nhit_sectors[sectors[i]]++;
    nhit_sectors_fees[sectors[i]][fees[i]]++;
  }
  for (int s = 0; s < 24; s++)
  {
    h_nhits_sectors[s]->Fill(nhit_sectors[s]);
    for (int f = 0; f < 26; f++)
    {
      h_nhits_sectors_fees[s]->Fill(f, nhit_sectors_fees[s][f]);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcRawHitQA::EndRun(const int /*runnumber*/)
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcRawHitQA::End(PHCompositeNode * /*unused*/) { return Fun4AllReturnCodes::EVENT_OK; }

std::string TpcRawHitQA::getHistoPrefix() const { return std::string("h_") + Name() + std::string("_"); }

void TpcRawHitQA::createHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  for (int s = 0; s < 24; s++)
  {
    {
      auto h = new TH1F(std::string(getHistoPrefix() + "nhits_sec" + std::to_string(s)).c_str(),
                        std::string("Number of Hits in Sector " + std::to_string(s) + ";Number of Hits;Entries").c_str(), 100, 0, 10000);
      hm->registerHisto(h);
    }
    {
      auto h = new TH2F(std::string(getHistoPrefix() + "nhits_sec" + std::to_string(s) + "_fees").c_str(),
                        std::string("Sector " + std::to_string(s) + " Fee Hit Distribution;FEE;Number of Hits").c_str(), 26, -0.5, 25.5, 100, 0, 1000);
      hm->registerHisto(h);
    }
    for (int r = 0; r < 3; r++)
    {
      auto h = new TH1F(std::string(getHistoPrefix() + "nhits_sample_sec" + std::to_string(s) + "_R" + std::to_string(r)).c_str(),
                        std::string("Sector " + std::to_string(s) + " Sample Time Distribution;Time Bin [1/20 MHz];Entries").c_str(), 400, -0.5, 399.5);
      hm->registerHisto(h);
      auto h2 = new TH1F(std::string(getHistoPrefix() + "adc_sec" + std::to_string(s) + "_R" + std::to_string(r)).c_str(),
                         std::string("Sector " + std::to_string(s) + " ADC Distribution;ADC-pedestal [ADU];Entries").c_str(), 281, -100, 1024);
      hm->registerHisto(h2);
    }
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "bco").c_str(), "BCO distribution;BCO;Entries", 100, 0, TMath::Power(2, 40));
    hm->registerHisto(h);
  }
  {
    auto h = new TH2F(std::string(getHistoPrefix() + "xyPos_North").c_str(), "Hit XY distribution (North);X [mm];Y [mm]", 400, -800, 800, 400, -800, 800);
    hm->registerHisto(h);
  }
  {
    auto h = new TH2F(std::string(getHistoPrefix() + "xyPos_South").c_str(), "Hit XY distribution (South);X [mm];Y [mm]", 400, -800, 800, 400, -800, 800);
    hm->registerHisto(h);
  }
}
