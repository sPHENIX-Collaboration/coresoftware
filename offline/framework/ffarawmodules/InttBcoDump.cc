#include "InttBcoDump.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <fun4all/Fun4AllHistoManager.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>

#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

#include <iostream>  // for operator<<, endl, basic_ost...
#include <set>
#include <utility>  // for pair
#include <vector>   // for vector

//____________________________________________________________________________..
InttBcoDump::InttBcoDump(const std::string &name)
  : SubsysReco(name)
{
}
//____________________________________________________________________________..
int InttBcoDump::InitRun(PHCompositeNode * /*topNode*/)
{
  if (outfilename.empty())
  {
    std::cout << "no output filename given" << std::endl;
    gSystem->Exit(1);
  }

  outfile = new TFile(outfilename.c_str(), "RECREATE");
  outfile->SetCompressionSettings(505);  // ZSTD
  ttree = new TTree("bco", "bco");
  ttree->Branch("id", &m_id);
  ttree->Branch("evt", &m_evt);
  ttree->Branch("nfees", &m_nfees);
  ttree->Branch("bco", &m_bco);
  ttree->Branch("bcodiff", &m_bcodiff);
  ttree->SetAutoFlush(100000);
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InttBcoDump::process_event(PHCompositeNode *topNode)
{
  Event *evt = findNode::getClass<Event>(topNode, "PRDF");
  if (!evt)
  {
    std::cout << "No Event found" << std::endl;
    exit(1);
  }
  if (evt->getEvtType() == ENDRUNEVENT)
  {
    std::cout << "End run flag for INTT found, remaining INTT data is corrupted" << std::endl;
    delete evt;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  //  evt->identify();
  int EventSequence = evt->getEvtSequence();
  std::vector<Packet *> pktvec = evt->getPacketVector();
  std::map<int, std::set<uint64_t>> bcoset;
  for (auto *packet : pktvec)
  {
    int nbcos = packet->iValue(0, "NR_BCOS");
    for (int i = 0; i < nbcos; i++)
    {
      uint64_t bco = packet->lValue(i, "BCOLIST");
      int nfees = packet->iValue(i, "NR_FEES");
      bcoTaggedFees[bco] = nfees;
      for (int j = 0; j < nfees; j++)
      {
        int fee = packet->iValue(i, j, "FEELIST");
        lastbco.insert(std::make_pair(fee, 0));

        bcoset[fee].insert(bco);
      }
    }

    delete packet;
  }
  pktvec.clear();

  for (auto &mapiter : bcoset)
  {
    if (!mapiter.second.empty())
    {
      for (const auto &bco : mapiter.second)
      {
        uint64_t prevbco = lastbco[mapiter.first];
        if (prevbco > 0 && prevbco != bco)
        {
          int64_t diffbco = bco - prevbco;

          m_id = mapiter.first;
          m_evt = EventSequence;
          m_bco = bco;
          m_nfees = bcoTaggedFees[bco];
          m_bcodiff = diffbco;

          ttree->Fill();
        }
        lastbco[mapiter.first] = bco;
      }
      mapiter.second.clear();
    }
  }
  bcoset.clear();
  bcoTaggedFees.clear();
  lastbco.clear();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InttBcoDump::End(PHCompositeNode * /*topNode*/)
{
  outfile->cd();
  ttree->Write();
  outfile->Close();
  delete outfile;
  return Fun4AllReturnCodes::EVENT_OK;
}
