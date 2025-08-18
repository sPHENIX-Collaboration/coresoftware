#include "TpcRawHitQA.h"

#include <qautils/QAHistManagerDef.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <ffarawobjects/TpcRawHit.h>
#include <ffarawobjects/TpcRawHitContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHPointerListIterator.h>

#include <TH1.h>
#include <TH2.h>

#include <stddef.h>   // for size_t
#include <stdint.h>   // for uint16_t, int32_t
#include <algorithm>  // for max, sort
#include <cassert>
#include <cmath>
#include <iostream>  // for basic_ostream, operator<<
#include <memory>    // for unique_ptr

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

  PHNodeIterator trkr_itr(topNode);
  PHCompositeNode *tpc_node = dynamic_cast<PHCompositeNode *>(trkr_itr.findFirst("PHCompositeNode", "TPC"));
  if (!tpc_node)
  {
    std::cout << __PRETTY_FUNCTION__ << " : ERROR : "
              << "No TPC node found, unregistering module" << std::endl;
    Fun4AllServer::instance()->unregisterSubsystem(this);
    return Fun4AllReturnCodes::EVENT_OK;
  }

  PHNodeIterator tpc_itr(tpc_node);
  {
    PHPointerListIterator<PHNode> iter(tpc_itr.ls());

    PHNode *thisNode_raw;
    while ((thisNode_raw = iter()))
    {
      if (thisNode_raw->getType() != "PHIODataNode")
      {
        continue;
      }

      PHIODataNode<TpcRawHitContainer> *thisNode = static_cast<PHIODataNode<TpcRawHitContainer> *>(thisNode_raw);
      if (thisNode)
      {
        std::cout << __PRETTY_FUNCTION__ << " : Found TpcRawHitContainer Node "
                  << thisNode->getName() << std::endl;

        TpcRawHitContainer *rawhitcont = (TpcRawHitContainer *) thisNode->getData();
        if (rawhitcont)
        {
          rawhitcont_vec.push_back(rawhitcont);
        }
      }
    }
  }

  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  for (int s = 0; s < 24; s++)
  {
    h_nhits_sectors[s] = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nhits_sec" + std::to_string(s)).c_str()));
    h_nhits_sectors_fees[s] = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "nhits_sec" + std::to_string(s) + "_fees").c_str()));
    h_nhits_sectors_laser[s] = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nhits_sec" + std::to_string(s) + "_laser").c_str()));
    h_nhits_sectors_fees_laser[s] = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "nhits_sec" + std::to_string(s) + "_fees_laser").c_str()));
    for (int r = 0; r < 3; r++)
    {
      h_nhits_sam[s][r] = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nhits_sample_sec" + std::to_string(s) + "_R" + std::to_string(r)).c_str()));
      h_adc[s][r] = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "adc_sec" + std::to_string(s) + "_R" + std::to_string(r)).c_str()));
    }
  }

  h_xy_N = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "xyPos_North").c_str()));
  h_xy_S = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "xyPos_South").c_str()));

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcRawHitQA::process_event(PHCompositeNode * /*unused*/)
{
  // std::vector < TpcRawHit* > hits;
  std::vector<int> sectors;
  std::vector<uint16_t> fees;

  float nhit_sectors[24] = {0};
  float nhit_sectors_fees[24][26] = {{0}};
  float nhit_sectors_laser[24] = {0};
  float nhit_sectors_fees_laser[24][26] = {{0}};

  unsigned int raw_hit_num = 0;
  for (TpcRawHitContainer *&rawhitcont : rawhitcont_vec)
  {
    raw_hit_num = rawhitcont->get_nhits();
    for (unsigned int i = 0; i < raw_hit_num; i++)
    {
      auto hit = rawhitcont->get_hit(i);
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

      float median = 60;
      float stdDev = 3;
      if (longPresample)
      {
        std::vector<int> values;
        values.reserve(sam);
        // for (int sampleN = 0; sampleN < sam; sampleN++)
        // {
        for (std::unique_ptr<TpcRawHit::AdcIterator> adc_iterator(hit->CreateAdcIterator());
             !adc_iterator->IsDone();
             adc_iterator->Next())
        {
          values.push_back((int) adc_iterator->CurrentAdc());
        }
        std::sort(values.begin(), values.end());
        size_t size = values.size();
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
      }

      // for (int sampleN = 0; sampleN < sam; sampleN++)
      // {
      //   float adc = hit->get_adc(sampleN);
      for (std::unique_ptr<TpcRawHit::AdcIterator> adc_iterator(hit->CreateAdcIterator());
           !adc_iterator->IsDone();
           adc_iterator->Next())
      {
        const uint16_t sampleN = adc_iterator->CurrentTimeBin();
        const uint16_t adc = adc_iterator->CurrentAdc();
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
        h_adc[sector][region - 1]->Fill(adc - median);
        if (sampleN >= 400 && sampleN <= 430)
        {
          nhit_sectors_laser[sector]++;
          nhit_sectors_fees_laser[sector][fee]++;
          continue;
        }
        nhit_sectors[sector]++;
        nhit_sectors_fees[sector][fee]++;
      }
    }
  }

  // if no raw hit is found, skip this event
  if (raw_hit_num == 0)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  for (int s = 0; s < 24; s++)
  {
    h_nhits_sectors[s]->Fill(nhit_sectors[s]);
    h_nhits_sectors_laser[s]->Fill(nhit_sectors_laser[s]);
    for (int f = 0; f < 26; f++)
    {
      h_nhits_sectors_fees[s]->Fill(f, nhit_sectors_fees[s][f]);
      h_nhits_sectors_fees_laser[s]->Fill(f, nhit_sectors_fees_laser[s][f]);
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

std::string TpcRawHitQA::getHistoPrefix() const { return std::string("h_") + Name() + std::string("_"); }

void TpcRawHitQA::createHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  for (int s = 0; s < 24; s++)
  {
      hm->registerHisto(new TH1F(std::string(getHistoPrefix() + "nhits_sec" + std::to_string(s)).c_str(),
					  std::string("Number of Hits in Sector " + std::to_string(s) + ";Number of Hits/Event;Entries").c_str(), 100, 0, 30000));
      hm->registerHisto(new TH2F(std::string(getHistoPrefix() + "nhits_sec" + std::to_string(s) + "_fees").c_str(),
				 std::string("Sector " + std::to_string(s) + " Fee Hit Distribution;FEE;Number of Hits/Event").c_str(), 26, -0.5, 25.5, 100, 0, 3000));
      hm->registerHisto(new TH1F(std::string(getHistoPrefix() + "nhits_sec" + std::to_string(s) + "_laser").c_str(),
				 std::string("Laser Hits in Sector " + std::to_string(s) + ";Number of Hits/Event;Entries").c_str(), 100, 0, 1000));
      hm->registerHisto(new TH2F(std::string(getHistoPrefix() + "nhits_sec" + std::to_string(s) + "_fees_laser").c_str(),
				 std::string("Sector " + std::to_string(s) + " Fee Laser Hits;FEE;Number of Hits/Event").c_str(), 26, -0.5, 25.5, 100, 0, 500));
    for (int r = 0; r < 3; r++)
    {
      hm->registerHisto(new TH1F(std::string(getHistoPrefix() + "nhits_sample_sec" + std::to_string(s) + "_R" + std::to_string(r)).c_str(),
				 std::string("Sector " + std::to_string(s) + " Sample Time Distribution;Time Bin [1/20 MHz];Number of Total Hits").c_str(), 1051, -0.5, 1050.5));
      hm->registerHisto(new TH1F(std::string(getHistoPrefix() + "adc_sec" + std::to_string(s) + "_R" + std::to_string(r)).c_str(),
				 std::string("Sector " + std::to_string(s) + " ADC Distribution;ADC-pedestal [ADU];Entries").c_str(), 281, -100, 1024));
    }
  }

    hm->registerHisto(new TH2F(std::string(getHistoPrefix() + "xyPos_North").c_str(), "Hit XY distribution (North);X [mm];Y [mm]", 400, -800, 800, 400, -800, 800));
    hm->registerHisto(new TH2F(std::string(getHistoPrefix() + "xyPos_South").c_str(), "Hit XY distribution (South);X [mm];Y [mm]", 400, -800, 800, 400, -800, 800));
}
