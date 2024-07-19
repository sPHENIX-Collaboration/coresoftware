#include "TpcRawHitQA.h"

#include <qautils/QAHistManagerDef.h>
#include <qautils/QAUtil.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TMath.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <cassert>

//____________________________________________________________________________..
TpcRawHitQA::TpcRawHitQA(const std::string &name)
  : SubsysReco(name)
{
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
  }

  h_bco = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "bco").c_str()));

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcRawHitQA::process_event(PHCompositeNode * /*unused*/)
{

  //std::vector < TpcRawHit* > hits;
  std::vector < uint64_t > bcos;
  std::vector < int > sectors;
  std::vector < uint16_t > fees;

  //hits.clear();
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
      int32_t packet_id = hit->get_packetid();
      int ep = (packet_id - 4000) % 10;
      int sector = (packet_id - 4000 - ep) / 10;
      uint16_t fee = hit->get_fee();
      //hits.push_back( hit );
      bcos.push_back( bco );
      sectors.push_back( sector );
      fees.push_back( fee );
    }
  }

  // if no raw hit is found, skip this event
  if( raw_hit_num == 0 ) 
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  float nhit_sectors[24] = {0}; 
  float nhit_sectors_fees[24][26] = {{0}}; 
  
  for (int i=0; i < (int)raw_hit_num; i++)
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
      std::cout << nhit_sectors_fees[s][f] << std::endl;
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
    auto h = new TH1F(std::string(getHistoPrefix() + "nhits_sec" + std::to_string(s)).c_str(), 
                      std::string("Number of Hits in Sector " + std::to_string(s) + ";Number of Hits;Entries").c_str(),100,0,10000);
    hm->registerHisto(h);
    auto h2 = new TH2F(std::string(getHistoPrefix() + "nhits_sec" + std::to_string(s) + "_fees").c_str(),
                       std::string("Sector " + std::to_string(s) + " Fee Hit Distribution;FEE;Number of Hits").c_str(),26,-0.5,25.5,100,0,1000);
    hm->registerHisto(h2);
  }
  
  {
    auto h = new TH1F(std::string(getHistoPrefix() + "bco").c_str(), "BCO distribution;BCO;Entries",100,0,TMath::Power( 2, 40 ));
    hm->registerHisto(h);
  }
}
