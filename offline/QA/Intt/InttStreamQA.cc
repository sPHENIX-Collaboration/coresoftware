#include "InttStreamQA.h"

/// Fun4All includes
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <ffarawobjects/Gl1Packet.h>
#include <ffarawobjects/InttRawHitContainer.h>
#include <ffarawobjects/InttRawHit.h>

#include <qautils/QAHistManagerDef.h>


#include <TH1.h>
#include <TH2.h>
#include <TFile.h>

#include <set>
#include <map>

using namespace std;


/**
 * Constructor of module
 */
InttStreamQA::InttStreamQA(const std::string &name)
  : SubsysReco(name)
{
}

/**
 * Destructor of module
 */
InttStreamQA::~InttStreamQA()
{
}

/**
 * Initialize the module and prepare looping over events
 */
int InttStreamQA::Init(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 5)
  {
    std::cout << "Beginning Init in InttStreamQA" << std::endl;
  }

  return 0;
}

int InttStreamQA::InitRun(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() > 5)
  {
    std::cout << "Beginning InitRun in InttStreamQA" << std::endl;
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
  static int event=0;
  static int nskip=0;

  Gl1Packet* gl1  = findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");

  uint64_t bco_gl1  = (gl1 !=nullptr) ? gl1->getBCO() : std::numeric_limits<uint64_t>::max();
  //int      evt_gl1  = (gl1 !=nullptr) ? gl1->getEvtSequence() : -1;
  int      bunch_gl1= (gl1!=nullptr) ? gl1->getBunchNumber() : -1;

  //uint64_t trig =  (gl1!=nullptr) ? gl1->getLiveVector() : 0;
  //uint64_t trig =  (gl1!=nullptr) ? gl1->getScaledVector() : 0;



  InttRawHitContainer* rawhitmap = findNode::getClass<InttRawHitContainer>(topNode, "INTTRAWHIT");
  if(rawhitmap==nullptr) {
    cout<<"rawhit  is null"<<endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  uint64_t bcointt  = (rawhitmap->get_nhits()>0)
                           ? rawhitmap->get_hit(0)->get_bco()
                           : std::numeric_limits<uint64_t>::max();


  static uint64_t prebcointt  = 0;

  if( bco_gl1 == std::numeric_limits<uint64_t>::max() ){
    cout<<"StreamQA bco is max. not valid"<<endl;
  }
  if( bcointt == std::numeric_limits<uint64_t>::max() ){
    cout<<"StreamQA inttbco is max. no intt data valid. skip nth: "<<nskip<<endl;
    nskip++;
  }
 
  // to 40bit
  bco_gl1  &= 0xFFFFFFFFFF;
  bcointt  &= 0xFFFFFFFFFF;

  uint64_t bcodiff = bcointt - prebcointt;

  int evtcnt= (rawhitmap->get_nhits()>0)
                           ? rawhitmap->get_hit(0)->get_event_counter()
                           : std::numeric_limits<uint64_t>::max();

  //////////////////////////////

  int64_t diff_inttgl1 = bcointt - bco_gl1;
  h_bcointtgl1_diff->Fill(diff_inttgl1);

  cout<<event<<"   bco : intt:0x"<<hex<<bcointt<<", diff 0x"<<bcodiff
             <<" : "<<dec<<" "<<evtcnt<<" "<<diff_inttgl1<<" gl1:0x"<<hex<<bco_gl1<<dec<<endl;

  event++;



  h_bunch_gl1->Fill(bunch_gl1);
  

  set<uint> vUnique[8];
  map<uint, int> vchipbco[8]; // key, ihit

  // loop rawhits to remove copy hit
  uint nhit = rawhitmap->get_nhits();
  for(uint ihit=0; ihit<nhit; ihit++){
    InttRawHit *hit = rawhitmap->get_hit(ihit);

    int      ifelix = hit->get_packetid() - 3001;
    int      bco    = hit->get_FPHX_BCO(); // 7bit
    uint64_t bcofull= (hit->get_bco()&0xFFFFFFFFFF); // 7bit

    int ladder = hit->get_fee();        // 0-13 (4bit)
    int chip   = (hit->get_chip_id()-1)%26;    // 0-26 (5bit)
    int chan   = hit->get_channel_id(); // 0-127 (7bit)
    int adc    = hit->get_adc();        // 3bit

    // check the difference between strobeBCO(bcofull) and gl1BCO
    int64_t bcogl1diff = bcofull - bco_gl1;
    //cout<<"bco-gl1diff " <<bcogl1diff<<endl;
    h_bcogl1diff_felix[ifelix]->Fill(bcogl1diff);

    // lad[25-22]+chip[21-17]+chan[16-10]+adc[9-7]+bco[6-0]
    uint key = ((ladder&0xF)<<22)|((chip&0x1F)<<17)|((chan&0x7F)<<10)|((adc&0x7)<<7)|(bco&0x7F) ;

    if(vUnique[ifelix].find(key)==vUnique[ifelix].end()) {
      vUnique[ifelix].insert(key);

      uint chipbcokey = ((ladder&0xF)<<22)|((chip&0x1F)<<17)|(bco&0x7F) ;
      vchipbco[ifelix].insert(std::make_pair(chipbcokey, ihit)); // no ADC info

      h_bco[ifelix]->Fill(ladder*26 + chip+0.5, bco+0.5);
      h_hit[ifelix]->Fill(ladder*26 + chip+0.5, chan+0.5);
    }

    //cout<<"    hit : "<<ihit<<" "<<ifelix<<" 0x"<<hex<<bco<<dec<<endl;
  }

  ////////////////////////////
  // felix by felix analysis
  map<int, int> vbcodiff_felix[8];
  for(int ifelix=0; ifelix<8; ifelix++){
   
    for(auto val : vchipbco[ifelix]){
      int bco = (val.first)&0x7F;
      h_bco_felix[ifelix]->Fill(bco);

      InttRawHit *hit = rawhitmap->get_hit(val.second);

      uint64_t bcofull= (hit->get_bco()&0xFFFFFFFFFF); // 7bit

      // stream mode
      //uint64_t bcofull_reco = bco + bcointt;
      int bcofull_reco = bco + bcofull;

      int  bcointtgl1_diff = bcofull_reco - bco_gl1;
      int  bcointthit_diff = bcofull_reco - bcointt;

      h_bcoreco_diff[ifelix]->Fill(bcointtgl1_diff);
      h_bcorecointt_diff[ifelix]->Fill(bcointthit_diff);

      h_bunch_bco[ifelix]->Fill(bco, bunch_gl1);
      h_bunch_evt_bcodiff[ifelix]->Fill(bcointtgl1_diff, bunch_gl1);

      auto bco_itr = vbcodiff_felix[ifelix].find(bcointtgl1_diff);
      if(bco_itr==vbcodiff_felix[ifelix].end()) {
        vbcodiff_felix[ifelix].insert(make_pair(bcointtgl1_diff,1));
        h_bcoreco_evt_diff[ifelix]->Fill(bcointtgl1_diff); // event by event for each felix
      } else {
        bco_itr->second += 1;
      }
    }
  }

  ////////////////////////////
  // all felix combined analysis
  map<int, int> vbcodiff_all;

  for(int ifelix=0; ifelix<8; ifelix++){
    for(auto& val : vbcodiff_felix[ifelix]){
      int bcointtgl1_diff = val.first;
      int count           = val.second;

      cout<<"             recobco diff : "<<bcointtgl1_diff<<" "<<count<<" "<<ifelix<<endl;

      auto bco_all_itr = vbcodiff_all.find(bcointtgl1_diff);
      if(bco_all_itr==vbcodiff_all.end()) {
        vbcodiff_all.insert(make_pair(bcointtgl1_diff,1));
        h_bcoreco_evt_diff_all->Fill(bcointtgl1_diff);
        
        if(bcointtgl1_diff==23) // should use variable
        {
          h_bunch_all->Fill(bunch_gl1);
        }

      } else {
        bco_all_itr->second += 1;
      }

    }
  }


  prebcointt  = bcointt;

  return Fun4AllReturnCodes::EVENT_OK;
}

/**
 * End the module and finish any data collection. Clean up any remaining
 * loose ends
 */
int InttStreamQA::End(PHCompositeNode * /*topNode*/)
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

std::string InttStreamQA::getHistoPrefix() const { 
  return std::string("h_") + Name() + std::string("_"); 
}

void InttStreamQA::createHistos(){

  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  for(int i=0; i<8; i++){
    h_bco[i] = new TH2F((getHistoPrefix()+Form("bco_%d", i)).c_str(), Form("bco_%d;bco;ladder*chip", i), 26*14, 0, 26*14, 140, -7, 133);
    h_hit[i] = new TH2F((getHistoPrefix()+Form("hit_%d", i)).c_str(), Form("hit_%d;bco;ladder*chip", i), 26*14, 0, 26*14, 128,  0, 128);
    hm->registerHisto(h_bco[i]);
    hm->registerHisto(h_hit[i]);

   
    // RecoBCO - GL1BCO
    h_bcoreco_diff[i]     = new TH1F((getHistoPrefix()+Form("bcoreco_diff_%d", i)).c_str(),     Form("bcoreco diff_%d", i),     540, -270, 270);
    h_bcoreco_evt_diff[i] = new TH1F((getHistoPrefix()+Form("bcoreco_evt_diff_%d", i)).c_str(), Form("bcoreco evt diff_%d", i), 540, -270, 270);
    hm->registerHisto(h_bcoreco_diff[i]);
    hm->registerHisto(h_bcoreco_evt_diff[i]);

    // RecoBCO - StrobeBCO, same as FPHX BCO
    h_bcorecointt_diff[i] = new TH1F((getHistoPrefix()+Form("bcorecointt_diff_%d", i)).c_str(), Form("bcoreco intt diff_%d", i), 540, -270, 270);
    hm->registerHisto(h_bcorecointt_diff[i]);

    // ChipBCO
    h_bco_felix[i]= new TH1F((getHistoPrefix()+Form("bco_felix_%d", i)).c_str(), Form("bcofelix_%d", i), 128,  0, 128);
    hm->registerHisto(h_bco_felix[i]);

    // StrobeBCO - GL1BCO
    h_bcogl1diff_felix[i]= new TH1F((getHistoPrefix()+Form("bcogl1diff_felix_%d", i)).c_str(), Form("bcogl1diff_felix_%d", i), 1024, -512, 512);
    hm->registerHisto(h_bcogl1diff_felix[i]);


    //h_bunch[i] = new TH1F((getHistoPrefix()+Form("bunch_%d", i)).c_str(), Form("bunch @ trigger_%d", i), 150, -15, 135);
    //hm->registerHisto(h_bunch[i]);

    //h_bunch_strb[i] = new TH1F((getHistoPrefix()+Form("bunch_strb_%d", i)).c_str(), Form("bunch @ strobe_%d", i), 150, -15, 135);
    //hm->registerHisto(h_bunch_strb[i]);

    // RecoBCO - Gl1BCO, vs  bunch, to see the peak
    h_bunch_evt_bcodiff[i] = new TH2F((getHistoPrefix()+Form("bunch_evt_bcodiff_%d", i)).c_str(), Form("bunch @ strobe_%d", i),
                                      750, -250, 500, 150, -15, 135);
    hm->registerHisto(h_bunch_evt_bcodiff[i]);

    //ChipBCO vs Bunch to check the linear correlation
    h_bunch_bco[i] = new TH2F((getHistoPrefix()+Form("bunch_bco_%d", i)).c_str(), Form("bunch vs BCO %d", i),
                                      150, 0, 150, 150, -15, 135);
    hm->registerHisto(h_bunch_bco[i]);

    // StrobeBCO - Prev StrobeBCO
    //h_bcoprediff[i] = new TH1F((getHistoPrefix()+Form("bcoprediff_%d", i)).c_str(), Form("BCO - PreBCO %d", i), 1000, 0, 1000);
  }

  h_bunch_all = new TH1F((getHistoPrefix()+"bunch_all").c_str(), "bunch @ evt all felix", 150, -15, 135);
  h_bunch_gl1 = new TH1F((getHistoPrefix()+"bunch_gl1").c_str(), "bunch @ gl1", 150, -15, 135);
  hm->registerHisto(h_bunch_all);
  hm->registerHisto(h_bunch_gl1);

  h_bcoreco_evt_diff_all = new TH1F((getHistoPrefix()+"h_bcoreco_evt_diff_all").c_str(), "bcoreco evt diff_all", 540, -270, 270);
  h_bcointtgl1_diff      = new TH1F((getHistoPrefix()+"h_bcointtgl1_diff").c_str(),      "bco intt gl1 diff_",   540, -270, 270);
  hm->registerHisto(h_bcoreco_evt_diff_all);
  hm->registerHisto(h_bcointtgl1_diff);
}

