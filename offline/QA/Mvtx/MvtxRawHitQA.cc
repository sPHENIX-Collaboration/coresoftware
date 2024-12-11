#include "MvtxRawHitQA.h"

#include <qautils/QAHistManagerDef.h>
#include <qautils/QAUtil.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHPointerListIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TMath.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <cassert>

//____________________________________________________________________________..
MvtxRawHitQA::MvtxRawHitQA(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int MvtxRawHitQA::InitRun(PHCompositeNode *topNode)
{
  createHistos();

 PHNodeIterator trkr_itr(topNode);
  PHCompositeNode *mvtx_node = dynamic_cast<PHCompositeNode *>(
      trkr_itr.findFirst("PHCompositeNode", "MVTX"));  
  if(!mvtx_node)
  {
    std::cout << PHWHERE << " No MVTX node found, exit" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHNodeIterator mvtx_itr(mvtx_node);
  PHPointerListIterator<PHNode> iter(mvtx_itr.ls());
  PHNode *thisnode;
  while((thisnode = iter()))
  {
    if(thisnode->getType() !="PHIODataNode")
    {
      continue;
    }
    // only want the raw hits, not the header nodes
    if((thisnode->getName()).find("HEADER") != std::string::npos)
    {
      continue;
    }
    PHIODataNode<MvtxRawHitContainer> *theNode = static_cast<PHIODataNode<MvtxRawHitContainer> *>(thisnode);
    if(theNode)
    {
      std::cout << PHWHERE << " Found Mvtx Raw hit container node " << theNode->getName() << std::endl;
      auto cont = (MvtxRawHitContainer*)theNode->getData();
      if(cont)
      {
        m_rawhit_containers.push_back(cont);
      }
    }
  }

  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  h_nhits_layer0 = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nhits_layer0").c_str()));
  h_nhits_layer1 = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nhits_layer1").c_str()));
  h_nhits_layer2 = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nhits_layer2").c_str()));

  h_bco = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "bco").c_str()));
  h_strobe_bc = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "strobe_bc").c_str()));
  h_chip_bc = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "chip_bc").c_str()));

  h_nhits_stave_chip_layer0 = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "nhits_stave_chip_layer0").c_str()));
  h_nhits_stave_chip_layer1 = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "nhits_stave_chip_layer1").c_str()));
  h_nhits_stave_chip_layer2 = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "nhits_stave_chip_layer2").c_str()));

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MvtxRawHitQA::process_event(PHCompositeNode * /*unused*/)
{

  std::vector < MvtxRawHit* > hits;
  std::vector < uint64_t > bcos;
  std::vector < uint32_t > strobe_bcs;
  std::vector < uint32_t > chip_bcs;
  std::vector < uint8_t > layers;
  std::vector < uint8_t > staves;
  std::vector < uint8_t > chips;
  std::vector < uint16_t > rows;
  std::vector < uint16_t > cols;

  hits.clear();
  bcos.clear();
  strobe_bcs.clear();
  chip_bcs.clear();
  layers.clear();
  staves.clear();
  chips.clear();
  rows.clear();
  cols.clear();

  unsigned int raw_hit_num = 0;
  for(auto& rawhitcont : m_rawhit_containers)
  {
  if (rawhitcont)
  {
    raw_hit_num = rawhitcont->get_nhits();
    for (unsigned int i = 0; i < raw_hit_num; i++)
    {
      auto hit = rawhitcont->get_hit(i);
      auto bco = hit->get_bco();
      auto strobe_bc = hit->get_strobe_bc();
      auto chip_bc = hit->get_chip_bc();
      auto layer = hit->get_layer_id();
      auto stave = hit->get_stave_id();
      auto chip = hit->get_chip_id();
      auto row = hit->get_row();
      auto col = hit->get_col();
      hits.push_back( hit );
      bcos.push_back( bco );
      strobe_bcs.push_back( strobe_bc );
      chip_bcs.push_back( chip_bc );
      layers.push_back( layer );
      staves.push_back( stave );
      chips.push_back( chip );
      rows.push_back( row );
      cols.push_back( col );
    }
  }

  // if no raw hit is found, skip this event
  if( raw_hit_num == 0 ) 
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  int nhit_layer0=0;
  int nhit_layer1=0;
  int nhit_layer2=0;
  for (int i=0; i<(int)raw_hit_num; i++)
  {
    if (layers[i]==0)
    {
      nhit_layer0++;
    }
    if (layers[i]==1)
    {
      nhit_layer1++;
    }
    if (layers[i]==2)
    {
      nhit_layer2++;
    }
  }
  h_nhits_layer0->Fill(nhit_layer0);
  h_nhits_layer1->Fill(nhit_layer1);
  h_nhits_layer2->Fill(nhit_layer2);

  h_bco->Fill(bcos[0]);
  h_strobe_bc->Fill(strobe_bcs[0]);
  h_chip_bc->Fill(chip_bcs[0]);

  for (int i=0; i<(int)raw_hit_num; i++)
  {
    if (layers[i]==0)
    {
      h_nhits_stave_chip_layer0->Fill(chips[i],staves[i]);
    }
    if (layers[i]==1)
    {
      h_nhits_stave_chip_layer1->Fill(chips[i],staves[i]);
    }
    if (layers[i]==2)
    {
      h_nhits_stave_chip_layer2->Fill(chips[i],staves[i]);
    }
  }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MvtxRawHitQA::EndRun(const int /*runnumber*/)
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MvtxRawHitQA::End(PHCompositeNode * /*unused*/) { return Fun4AllReturnCodes::EVENT_OK; }

std::string MvtxRawHitQA::getHistoPrefix() const { return std::string("h_") + Name() + std::string("_"); }

void MvtxRawHitQA::createHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "nhits_layer0").c_str(), "Number of hits in layer 0;Number of hits;Entries",100,0,10000);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "nhits_layer1").c_str(), "Number of hits in layer 1;Number of hits;Entries",100,0,10000);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "nhits_layer2").c_str(), "Number of hits in layer 2;Number of hits;Entries",100,0,10000);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "bco").c_str(), "BCO distribution;BCO;Entries",100,0,TMath::Power( 2, 40 ));
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "strobe_bc").c_str(), "Strobe BC distribution;Strobe BC;Entries",100,0,4000);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "chip_bc").c_str(), "Chip BC distribution;Chip BC;Entries",100,0,500);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "nhits_stave_chip_layer0").c_str(), "Hitmap in layer 0;ChipID;StaveID;Number of hits",9,0,9,12,0,12);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "nhits_stave_chip_layer1").c_str(), "Hitmap in layer 1;ChipID;StaveID;Number of hits",9,0,9,16,0,16);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "nhits_stave_chip_layer2").c_str(), "Hitmap in layer 2;ChipID;StaveID;Number of hits",9,0,9,20,0,20);
    hm->registerHisto(h);
  }
}
