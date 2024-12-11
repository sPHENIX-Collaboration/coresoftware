#include "InttRawHitQA.h"
#include "InttQaCommon.h"  // for HistConfig, kFee_num

#include <qautils/QAHistManagerDef.h>

#include <ffarawobjects/InttRawHit.h>           // for InttRawHit
#include <ffarawobjects/InttRawHitContainer.h>  // for InttRawHitContainer

#include <phool/PHCompositeNode.h>
#include <phool/PHPointerListIterator.h>
#include <fun4all/Fun4AllHistoManager.h>  // for Fun4AllHistoManager
#include <fun4all/Fun4AllReturnCodes.h>   // for EVENT_OK, ABORTEVENT

#include <phool/getClass.h>  // for getClass
#include <phool/phool.h>     // for PHWHERE

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>

//____________________________________________________________________________..
InttRawHitQA::InttRawHitQA(const std::string &name)
  : SubsysReco(name)
{
}

std::vector<InttRawHit *> InttRawHitQA::GetHits(InttRawHitContainer* container)
{
  std::vector<InttRawHit *> hits;
  if(container == nullptr)
  {
    return hits;
  }
  auto raw_hit_num = container->get_nhits();
  for (unsigned int i = 0; i < raw_hit_num; i++)
  {
    auto hit = container->get_hit(i);
    hits.push_back(hit);
  }

  return hits;
}

int InttRawHitQA::Init(PHCompositeNode * /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int InttRawHitQA::InitRun(PHCompositeNode *topNode)
{
  createHistos();

  /////////////////////////////////////////////////////////////////////////
  // INTT raw hit
  /////////////////////////////////////////////////////////////////////////
  PHNodeIterator trkr_itr(topNode);
  PHCompositeNode *intt_node = dynamic_cast<PHCompositeNode *>(
      trkr_itr.findFirst("PHCompositeNode", "INTT"));  
  if(!intt_node)
  {
    std::cout << PHWHERE << " No INTT node found, exit" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHNodeIterator intt_itr(intt_node);
  PHPointerListIterator<PHNode> iter(intt_itr.ls());
  PHNode *thisnode;
  while((thisnode = iter()))
  {
    if(thisnode->getType() !="PHIODataNode")
    {
      continue;
    }
    PHIODataNode<InttRawHitContainer> *theNode = static_cast<PHIODataNode<InttRawHitContainer> *>(thisnode);
    if(theNode)
    {
      std::cout << PHWHERE << " Found INTT Raw hit container node " << theNode->getName() << std::endl;
      auto cont = (InttRawHitContainer*)theNode->getData();
      if(cont)
      {
        m_rawhit_containers.push_back(cont);
      }
    }
  }
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  for (int felix = 0; felix < InttQa::kFelix_num; felix++)
  {
    std::string name = getHistoPrefix() + "intt" + std::to_string(felix);
    hist_fee_chip_chan_[felix] = dynamic_cast<TH3D *>(hm->getHisto(name.c_str()));

    std::string name_bco_event = name + "_ladder_bco_full_event_counter";
    hist_fee_bco_full_event_counter_[felix] = dynamic_cast<TH3D *>(hm->getHisto(name_bco_event.c_str()));

    std::string name_bco_event_diff = name + "_ladder_bco_full_event_counter_diff";
    hist_fee_bco_full_event_counter_diff_[felix] = dynamic_cast<TH3D *>(hm->getHisto(name_bco_event_diff.c_str()));

    std::string name_event_counter = name + "_event_counter";
    hist_event_counter_[felix] = dynamic_cast<TH1D *>(hm->getHisto(name_event_counter.c_str()));

    std::string name_event_counter_diff = name + "_event_counter_diff";
    hist_event_counter_diff_[felix] = dynamic_cast<TH1D *>(hm->getHisto(name_event_counter_diff.c_str()));
  }

  for (int felix = 0; felix < InttQa::kFelix_num; felix++)
  {
    for (int ladder = 0; ladder < InttQa::kFee_num; ladder++)
    {
      std::string name = getHistoPrefix() + "intt" + std::to_string(felix) + "_" + std::to_string(ladder);
      hist_hitmap_[felix][ladder] = dynamic_cast<TH2I *>(hm->getHisto(name.c_str()));
    }
  }

  hist_nhit_ = dynamic_cast<TH1D *>(hm->getHisto(std::string(getHistoPrefix() + "nhit").c_str()));

  hist_nhit_south_ = dynamic_cast<TH1D *>(hm->getHisto(std::string(getHistoPrefix() + "nhit_south").c_str()));
  hist_nhit_north_ = dynamic_cast<TH1D *>(hm->getHisto(std::string(getHistoPrefix() + "nhit_north").c_str()));
  hist_pid_ = dynamic_cast<TH1D *>(hm->getHisto(std::string(getHistoPrefix() + "pid").c_str()));
  hist_adc_ = dynamic_cast<TH1D *>(hm->getHisto(std::string(getHistoPrefix() + "adc").c_str()));
  hist_bco_ = dynamic_cast<TH1D *>(hm->getHisto(std::string(getHistoPrefix() + "bco").c_str()));
  hist_bco_full_ = dynamic_cast<TH1D *>(hm->getHisto(std::string(getHistoPrefix() + "bco_full").c_str()));

  return Fun4AllReturnCodes::EVENT_OK;
}

void InttRawHitQA::createHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  for (int felix = 0; felix < InttQa::kFelix_num; felix++)
  {
    std::string name = getHistoPrefix() + "intt" + std::to_string(felix);
    std::string title = name + ";FELIX CH;Chip;Channel;Entries";
    {
      auto h = new TH3D(name.c_str(), title.c_str(),
                        InttQa::kFee_num, 0, InttQa::kFee_num,
                        InttQa::kChip_num, 1, InttQa::kChip_num + 1,
                        InttQa::kChan_num, 0, InttQa::kChan_num);
      hm->registerHisto(h);
    }

    std::string name_bco_event = name + "_ladder_bco_full_event_counter";
    std::string title_bco_event = name + ";FELIX_CH;BCO full;Event Counter;Entries";
    {
      auto h = new TH3D(name_bco_event.c_str(), title_bco_event.c_str(),
                        InttQa::kFee_num, 0, InttQa::kFee_num,
                        100, 0, std::pow(2, 40),
                        1e4, 0, 1e7);
      hm->registerHisto(h);
    }

    std::string name_bco_event_diff = name + "_ladder_bco_full_event_counter_diff";
    std::string title_bco_event_diff = name + ";FELIX_CH;#Delta BCO full;#Delta Event Counter;Entries";
    int max = 1000;
    {
      auto h = new TH3D(name_bco_event_diff.c_str(), title_bco_event_diff.c_str(),
                        InttQa::kFee_num, 0, InttQa::kFee_num,
                        2 * max / 100, -max, max,
                        2 * max / 100, -max, max);
      hm->registerHisto(h);
    }

    std::string name_event_counter = name + "_event_counter";
    std::string title_event_counter = name_event_counter + ";Event Counter;Entries";
    {
      auto h = new TH1D(name_event_counter.c_str(), title_event_counter.c_str(), 1e4, 0, 1e7);
      hm->registerHisto(h);
    }

    std::string name_event_counter_diff = name + "_event_counter_diff";
    std::string title_event_counter_diff = name_event_counter_diff + ";Event Counter;Entries";
    {
      auto h = new TH1D(name_event_counter_diff.c_str(), title_event_counter_diff.c_str(), 2 * max / 100, -max, max);
      hm->registerHisto(h);
    }
  }

  for (int felix = 0; felix < InttQa::kFelix_num; felix++)
  {
    for (int ladder = 0; ladder < InttQa::kFee_num; ladder++)
    {
      std::string name = getHistoPrefix() + "intt" + std::to_string(felix) + "_" + std::to_string(ladder);
      std::string title = name + ";Chip;Channel;Entries";
      auto h = new TH2I(name.c_str(), title.c_str(),
                        InttQa::kChip_num, 1, InttQa::kChip_num,
                        InttQa::kChan_num, 0, InttQa::kChan_num);
      hm->registerHisto(h);
    }
  }

  {
    auto h = new TH1D(std::string(getHistoPrefix() + "nhit").c_str(), "#INTTRAWHIT per event;#hit;Entries", 1e4, 0, 1e4);
    InttQa::HistConfig(h);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1D(std::string(getHistoPrefix() + "nhit_south").c_str(), "#INTTRAWHIT South;event;#hit", 1e4, 0, 1e7);
    InttQa::HistConfig(h);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1D(std::string(getHistoPrefix() + "nhit_north").c_str(), "#INTTRAWHIT North;event;#hit", 1e4, 0, 1e7);
    InttQa::HistConfig(h);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1D(std::string(getHistoPrefix() + "pid").c_str(), "Packet ID distribution;pid;Entries", InttQa::kFelix_num, InttQa::kFirst_pid, InttQa::kFirst_pid + InttQa::kFelix_num);
    InttQa::HistConfig(h);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1D(std::string(getHistoPrefix() + "adc").c_str(), "ADC distribution;ADC;Entries", 8, 0, 8);
    InttQa::HistConfig(h);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1D(std::string(getHistoPrefix() + "bco").c_str(), "BCO distribution;BCO;Entries", InttQa::kBco_max + 10, -5, InttQa::kBco_max + 5);
    InttQa::HistConfig(h);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1D(std::string(getHistoPrefix() + "bco_full").c_str(), "BCO full distribution;BCO full;Entries", 100, 0, std::pow(2, 40));
    InttQa::HistConfig(h);
    hm->registerHisto(h);
  }
}

int InttRawHitQA::process_event(PHCompositeNode * /*unused*/)
{
  for(auto& cont : m_rawhit_containers)
  {
   
  auto hits = GetHits(cont);

  auto raw_hit_num = hits.size();
  hist_nhit_->Fill(raw_hit_num);

  // if no raw hit is found, skip this event
  if (raw_hit_num == 0 || cont == nullptr)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  event_counter_by_myself_++;

  //////////////////////////////////////////////////////////////////
  // processes for each event                                     //
  //////////////////////////////////////////////////////////////////
  uint64_t bco_full = (cont->get_hit(0)->get_bco());
  hist_bco_full_->Fill(bco_full);

  //////////////////////////////////////////////////////////////////
  // primary raw hit sweep to get some reference values           //
  //////////////////////////////////////////////////////////////////
  uint32_t event_counter_ref = cont->get_hit(0)->get_event_counter();

  //////////////////////////////////////////////////////////////////
  // processes for each raw hit                                   //
  //////////////////////////////////////////////////////////////////
  // loop over all raw hits
  for (auto hit : hits)
  {
    int felix = hit->get_packetid() - InttQa::kFirst_pid;
    int felix_ch = hit->get_fee();

    // uint16_t InttRawHit::get_chip_id
    int chip = hit->get_chip_id();
    if (chip > InttQa::kChip_num)
    {
      chip = chip - InttQa::kChip_num;
    }

    int chan = hit->get_channel_id();
    auto adc = hit->get_adc();
    auto bco = hit->get_FPHX_BCO();
    int event_counter = hit->get_event_counter();  // uint32_t IttRawHit::get_event_counter()
    if (is_first_event_ == true)
    {
      event_counter = 0;
    }
    else if (event_counter - previous_event_counter_ > 1000)
    {
      event_counter = -1;  // it means bad
    }
    else
    {
      last_event_counter_ = event_counter;
    }

    /*
        int bco_diff = 0;
        if( (bco_full & 0x7f )  > bco )
          bco_diff = int(bco_full & 0x7f ) - bco;
        else
          bco_diff = int(bco_full & 0x7f ) + (128 - bco);
    */

    //////////////////////////////////////////////////////////////////
    // Filling hists                                                //
    //////////////////////////////////////////////////////////////////

    hist_fee_chip_chan_[felix]->Fill(felix_ch, chip, chan);

    hist_hitmap_[felix][felix_ch]->Fill(chip, chan);

    hist_pid_->AddBinContent(felix + 1);

    hist_adc_->Fill(adc);
    hist_bco_->Fill(bco);

    hist_fee_bco_full_event_counter_[felix]
        ->Fill(felix_ch,  // chip, chan );
               hit->get_bco(),
               event_counter);

    hist_fee_bco_full_event_counter_diff_[felix]
        ->Fill(felix_ch,
               hit->get_bco() - bco_full,
               event_counter - event_counter_ref);

    hist_event_counter_[felix]
        ->Fill(event_counter);

    hist_event_counter_diff_[felix]
        ->Fill(event_counter - event_counter_ref);

    if (felix < 4)
    {
      hist_nhit_south_->Fill(event_counter);
    }
    else
    {
      hist_nhit_north_->Fill(event_counter);
    }
  }

  is_first_event_ = false;

  if (last_event_counter_ - previous_event_counter_ < 1000)
  {
    previous_event_counter_ = last_event_counter_;  // in the case of reasonable event counter
  }
  else
  {
    previous_event_counter_ = -1;  // in the case of a crazy event counter
  }
  }
  // cout << "-------------------------------------------------" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int InttRawHitQA::ResetEvent(PHCompositeNode * /*unused*/)
{
  // Intitialize for Clone hit counter
  return Fun4AllReturnCodes::EVENT_OK;
}

int InttRawHitQA::EndRun(const int /*unused*/)
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttRawHitQA::End(PHCompositeNode * /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

std::string InttRawHitQA::getHistoPrefix() const { return std::string("h_") + Name() + std::string("_"); }

int InttRawHitQA::Reset(PHCompositeNode * /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
