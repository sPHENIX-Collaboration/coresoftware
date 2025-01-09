#include "InttStreamQA.h"

/// Fun4All includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <ffarawobjects/Gl1Packet.h>
#include <ffarawobjects/InttRawHit.h>
#include <ffarawobjects/InttRawHitContainer.h>
#include <phool/PHPointerListIterator.h>
#include <qautils/QAHistManagerDef.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include <map>
#include <set>
#include <sstream>

/**
 * Constructor of module
 */
InttStreamQA::InttStreamQA(const std::string& name)
  : SubsysReco(name)
{
}

int InttStreamQA::InitRun(PHCompositeNode* topNode)
{
  if (Verbosity() > 5)
  {
    std::cout << "Beginning InitRun in InttStreamQA" << std::endl;
  }

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

  createHistos();

  return Fun4AllReturnCodes::EVENT_OK;
}

/**
 * Main workhorse function where each event is looped over and
 * data from each event is collected from the node tree for analysis
 */
int InttStreamQA::process_event(PHCompositeNode* topNode)
{
  static int event = 0;
  static int nskip = 0;

  Gl1Packet* gl1 = findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");

  uint64_t bco_gl1 = (gl1 != nullptr) ? gl1->getBCO() : std::numeric_limits<uint64_t>::max();
  // int      evt_gl1  = (gl1 !=nullptr) ? gl1->getEvtSequence() : -1;
  int bunch_gl1 = (gl1 != nullptr) ? gl1->getBunchNumber() : -1;

  // uint64_t trig =  (gl1!=nullptr) ? gl1->getLiveVector() : 0;
  // uint64_t trig =  (gl1!=nullptr) ? gl1->getScaledVector() : 0;

for(auto& rawhitmap : m_rawhit_containers)
{
  if (rawhitmap == nullptr)
  {
    std::cout << PHWHERE << "rawhit  is null" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  uint64_t bcointt = (rawhitmap->get_nhits() > 0)
                         ? rawhitmap->get_hit(0)->get_bco()
                         : std::numeric_limits<uint64_t>::max();

  static uint64_t prebcointt = 0;

  if (bco_gl1 == std::numeric_limits<uint64_t>::max())
  {
    std::cout << "StreamQA bco is max. not valid" << std::endl;
  }
  if (bcointt == std::numeric_limits<uint64_t>::max())
  {
    std::cout << "StreamQA inttbco is max. no intt data valid. skip nth: " << nskip << std::endl;
    nskip++;
  }

  // to 40bit
  bco_gl1 &= 0xFFFFFFFFFFULL;
  bcointt &= 0xFFFFFFFFFFULL;

  uint64_t bcodiff = bcointt - prebcointt;

  int evtcnt = (rawhitmap->get_nhits() > 0)
                   ? rawhitmap->get_hit(0)->get_event_counter()
                   : std::numeric_limits<int>::max();

  //////////////////////////////

  int64_t diff_inttgl1 = bcointt - bco_gl1;
  h_bcointtgl1_diff->Fill(diff_inttgl1);

  if (Verbosity() > 5)
  {
    std::cout << event << "   bco : intt:0x" << std::hex << bcointt << ", diff 0x" << bcodiff
              << " : " << std::dec << " " << evtcnt << " " << diff_inttgl1 << " gl1:0x" << std::hex << bco_gl1 << std::dec << std::endl;
  }

  event++;

  h_bunch_gl1->Fill(bunch_gl1);

  std::set<uint> vUnique[8];
  std::map<uint, int> vchipbco[8];  // key, ihit

  // loop rawhits to remove copy hit
  uint nhit = rawhitmap->get_nhits();
  for (uint ihit = 0; ihit < nhit; ihit++)
  {
    InttRawHit* hit = rawhitmap->get_hit(ihit);

    int ifelix = hit->get_packetid() - 3001;
    uint bco = hit->get_FPHX_BCO();                         // 7bit
    uint64_t bcofull = (hit->get_bco() & 0xFFFFFFFFFFULL);  // 7bit

    uint ladder = hit->get_fee();               // 0-13 (4bit)
    uint chip = (hit->get_chip_id() - 1) % 26;  // 0-26 (5bit)
    uint chan = hit->get_channel_id();          // 0-127 (7bit)
    uint adc = hit->get_adc();                  // 3bit

    // check the difference between strobeBCO(bcofull) and gl1BCO
    int64_t bcogl1diff = bcofull - bco_gl1;
    // std::cout<<"bco-gl1diff " <<bcogl1diff<<std::endl;
    h_bcogl1diff_felix[ifelix]->Fill(bcogl1diff);

    // lad[25-22]+chip[21-17]+chan[16-10]+adc[9-7]+bco[6-0]
    uint key = ((ladder & 0xFU) << 22U) | ((chip & 0x1FU) << 17U) | ((chan & 0x7FU) << 10U) | ((adc & 0x7U) << 7U) | (bco & 0x7FU);

    auto ret = vUnique[ifelix].insert(key);
    if (ret.second)  // insertion successfull --> entry did not exist
    {
      uint chipbcokey = ((ladder & 0xFU) << 22U) | ((chip & 0x1FU) << 17U) | (bco & 0x7FU);
      vchipbco[ifelix].insert(std::make_pair(chipbcokey, ihit));  // no ADC info

      h_bco[ifelix]->Fill(ladder * 26 + chip + 0.5, bco + 0.5);
      h_hit[ifelix]->Fill(ladder * 26 + chip + 0.5, chan + 0.5);
    }

    // std::cout<<"    hit : "<<ihit<<" "<<ifelix<<" 0x"<<std::hex<<bco<<std::dec<<std::endl;
  }

  ////////////////////////////
  // felix by felix analysis
  std::map<int, int> vbcodiff_felix[8];
  for (int ifelix = 0; ifelix < 8; ifelix++)
  {
    for (auto val : vchipbco[ifelix])
    {
      uint bco = (val.first) & 0x7FU;
      h_bco_felix[ifelix]->Fill(bco);

      InttRawHit* hit = rawhitmap->get_hit(val.second);

      uint64_t bcofull = (hit->get_bco() & 0xFFFFFFFFFFULL);  // 7bit

      // stream mode
      // uint64_t bcofull_reco = bco + bcointt;
      int bcofull_reco = bco + bcofull;

      int bcointtgl1_diff = bcofull_reco - bco_gl1;
      int bcointthit_diff = bcofull_reco - bcointt;

      h_bcoreco_diff[ifelix]->Fill(bcointtgl1_diff);
      h_bcorecointt_diff[ifelix]->Fill(bcointthit_diff);

      h_bunch_bco[ifelix]->Fill(bco, bunch_gl1);
      h_bunch_evt_bcodiff[ifelix]->Fill(bcointtgl1_diff, bunch_gl1);

      auto bco_itr = vbcodiff_felix[ifelix].find(bcointtgl1_diff);
      if (bco_itr == vbcodiff_felix[ifelix].end())
      {
        vbcodiff_felix[ifelix].insert(std::make_pair(bcointtgl1_diff, 1));
        h_bcoreco_evt_diff[ifelix]->Fill(bcointtgl1_diff);  // event by event for each felix
      }
      else
      {
        bco_itr->second += 1;
      }
    }
  }

  ////////////////////////////
  // all felix combined analysis
  std::map<int, int> vbcodiff_all;

  for (auto& ifelix : vbcodiff_felix)
  {
    for (auto& val : ifelix)
    {
      int bcointtgl1_diff = val.first;
      // int count           = val.second;

      // std::cout<<"             recobco diff : "<<bcointtgl1_diff<<" "<<count<<" "<<ifelix<<std::endl;

      auto bco_all_itr = vbcodiff_all.find(bcointtgl1_diff);
      if (bco_all_itr == vbcodiff_all.end())
      {
        vbcodiff_all.insert(std::make_pair(bcointtgl1_diff, 1));
        h_bcoreco_evt_diff_all->Fill(bcointtgl1_diff);

        if (bcointtgl1_diff == 23)  // should use variable
        {
          h_bunch_all->Fill(bunch_gl1);
        }
      }
      else
      {
        bco_all_itr->second += 1;
      }
    }
  }

  prebcointt = bcointt;
}
  return Fun4AllReturnCodes::EVENT_OK;
}

/**
 * End the module and finish any data collection. Clean up any remaining
 * loose ends
 */
int InttStreamQA::End(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() > 1)
  {
    std::cout << "Ending InttStreamQA analysis package" << std::endl;
  }

  /*
    TFile *froot = TFile::Open("streamrecof4a.root","recreate");
    for(int i=0; i<8; i++){
      h_bco[i]->Write();
      h_hit[i]->Write();
      h_bco_felix[i]->Write();
      //h_bunch[i]->Write();
      h_bcoreco_diff[i]->Write();
      h_bcoreco_evt_diff[i]->Write();
      h_bcorecointt_diff[i]->Write();
      h_bcogl1diff_felix[i]->Write();

      //h_bunch_strb[i]->Write();
      h_bunch_evt_bcodiff[i]->Write();
      h_bunch_bco[i]->Write();
      //h_bcoprediff[i]->Write();
    }
    h_bcoreco_evt_diff_all->Write();
    h_bunch_all->Write();
    h_bunch_gl1->Write();
    h_bcointtgl1_diff->Write();
    froot->Close();
  */

  return 0;
}

std::string InttStreamQA::getHistoPrefix() const
{
  return std::string("h_") + Name() + std::string("_");
}

void InttStreamQA::createHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  std::string sname, stitle;
  for (int i = 0; i < 8; i++)
  {
    sname = (getHistoPrefix() + "bco_" + std::to_string(i));
    stitle = ("bco_" + std::to_string(i) + ";bco;ladder*chip");
    h_bco[i] = new TH2F(sname.c_str(), stitle.c_str(), 26 * 14, 0, 26 * 14, 140, -7, 133);
    hm->registerHisto(h_bco[i]);

    sname = (getHistoPrefix() + "hit_" + std::to_string(i));
    stitle = ("hit_" + std::to_string(i) + ";bco;ladder*chip");
    h_hit[i] = new TH2F(sname.c_str(), stitle.c_str(), 26 * 14, 0, 26 * 14, 128, 0, 128);
    hm->registerHisto(h_hit[i]);

    // RecoBCO - GL1BCO
    sname = (getHistoPrefix() + "bcoreco_diff_" + std::to_string(i));
    stitle = ("bcoreco diff_" + std::to_string(i));
    h_bcoreco_diff[i] = new TH1F(sname.c_str(), stitle.c_str(), 540, -270, 270);
    hm->registerHisto(h_bcoreco_diff[i]);

    sname = (getHistoPrefix() + "bcoreco_evt_diff_" + std::to_string(i));
    stitle = ("bcoreco evt diff_" + std::to_string(i));
    h_bcoreco_evt_diff[i] = new TH1F(sname.c_str(), stitle.c_str(), 540, -270, 270);
    hm->registerHisto(h_bcoreco_evt_diff[i]);

    // RecoBCO - StrobeBCO, same as FPHX BCO
    sname = (getHistoPrefix() + "bcorecointt_diff_" + std::to_string(i));
    stitle = ("bcoreco intt diff_" + std::to_string(i));
    h_bcorecointt_diff[i] = new TH1F(sname.c_str(), stitle.c_str(), 540, -270, 270);
    hm->registerHisto(h_bcorecointt_diff[i]);

    // ChipBCO
    sname = (getHistoPrefix() + "bco_felix_" + std::to_string(i));
    stitle = ("bco_felix_" + std::to_string(i));
    h_bco_felix[i] = new TH1F(sname.c_str(), stitle.c_str(), 128, 0, 128);
    hm->registerHisto(h_bco_felix[i]);

    // StrobeBCO - GL1BCO
    sname = (getHistoPrefix() + "bcogl1diff_felix_" + std::to_string(i));
    stitle = ("bcogl1diff_felix_" + std::to_string(i));
    h_bcogl1diff_felix[i] = new TH1F(sname.c_str(), stitle.c_str(), 1024, -512, 512);
    hm->registerHisto(h_bcogl1diff_felix[i]);

    // h_bunch[i] = new TH1F((getHistoPrefix()+Form("bunch_%d", i)).c_str(), Form("bunch @ trigger_%d", i), 150, -15, 135);
    // hm->registerHisto(h_bunch[i]);

    // h_bunch_strb[i] = new TH1F((getHistoPrefix()+Form("bunch_strb_%d", i)).c_str(), Form("bunch @ strobe_%d", i), 150, -15, 135);
    // hm->registerHisto(h_bunch_strb[i]);

    // RecoBCO - Gl1BCO, vs  bunch, to see the peak
    sname = (getHistoPrefix() + "bunch_evt_bcodiff_" + std::to_string(i));
    stitle = ("bunch @ strobe_" + std::to_string(i));
    h_bunch_evt_bcodiff[i] = new TH2F(sname.c_str(), stitle.c_str(), 750, -250, 500, 150, -15, 135);
    hm->registerHisto(h_bunch_evt_bcodiff[i]);

    // ChipBCO vs Bunch to check the linear correlation
    sname = (getHistoPrefix() + "bunch_bco_" + std::to_string(i));
    stitle = ("bunch vs BCO " + std::to_string(i));
    h_bunch_bco[i] = new TH2F(sname.c_str(), stitle.c_str(), 150, 0, 150, 150, -15, 135);
    hm->registerHisto(h_bunch_bco[i]);

    // StrobeBCO - Prev StrobeBCO
    // h_bcoprediff[i] = new TH1F((getHistoPrefix()+Form("bcoprediff_%d", i)).c_str(), Form("BCO - PreBCO %d", i), 1000, 0, 1000);
  }

  h_bunch_all = new TH1F((getHistoPrefix() + "bunch_all").c_str(), "bunch @ evt all felix", 150, -15, 135);
  h_bunch_gl1 = new TH1F((getHistoPrefix() + "bunch_gl1").c_str(), "bunch @ gl1", 150, -15, 135);
  hm->registerHisto(h_bunch_all);
  hm->registerHisto(h_bunch_gl1);

  h_bcoreco_evt_diff_all = new TH1F((getHistoPrefix() + "h_bcoreco_evt_diff_all").c_str(), "bcoreco evt diff_all", 540, -270, 270);
  h_bcointtgl1_diff = new TH1F((getHistoPrefix() + "h_bcointtgl1_diff").c_str(), "bco intt gl1 diff_", 540, -270, 270);
  hm->registerHisto(h_bcoreco_evt_diff_all);
  hm->registerHisto(h_bcointtgl1_diff);
}
