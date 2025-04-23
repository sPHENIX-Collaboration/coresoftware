#include "TpcBcoDump.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <fun4all/Fun4AllHistoManager.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>

#include <Event/Event.h>
#include <Event/packet.h>

#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

#include <iostream>  // for operator<<, endl, basic_ost...
#include <set>
#include <utility>  // for pair
#include <vector>   // for vector

//____________________________________________________________________________..
TpcBcoDump::TpcBcoDump(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int TpcBcoDump::InitRun(PHCompositeNode * /*topNode*/)
{
  if (outfilename.empty())
  {
    std::cout << "no output filename given" << std::endl;
    gSystem->Exit(1);
  }

  outfile = new TFile(outfilename.c_str(), "RECREATE");
  ttree = new TTree("bco", "bco");
  ttree->Branch("id", &m_id);
  ttree->Branch("evt", &m_evt);
  ttree->Branch("bco", &m_bco);
  ttree->Branch("bcodiff", &m_bcodiff);
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcBcoDump::process_event(PHCompositeNode *topNode)
{
  Event *evt = findNode::getClass<Event>(topNode, "PRDF");
  if (!evt)
  {
    std::cout << "No Event found" << std::endl;
    exit(1);
  }
  //  evt->identify();
  int EventSequence = evt->getEvtSequence();
  std::vector<Packet *> pktvec = evt->getPacketVector();
  std::map<int, std::set<uint64_t>> bcoset;
  for (auto *packet : pktvec)
  {
    int packetid = packet->getIdentifier();
    lastbco.insert(std::make_pair(packetid, 0));
    int numBCOs = packet->lValue(0, "N_TAGGER");
    for (int j = 0; j < numBCOs; j++)
    {
      const auto is_lvl1 = static_cast<uint8_t>(packet->lValue(j, "IS_LEVEL1_TRIGGER"));
      if (is_lvl1)
      {
        uint64_t bco = packet->lValue(j, "BCO");
        bcoset[packetid].insert(bco);
      }
    }
    delete packet;
  }

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
          m_bcodiff = diffbco;

          ttree->Fill();
        }
        lastbco[mapiter.first] = bco;
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcBcoDump::End(PHCompositeNode * /*topNode*/)
{
  outfile->cd();
  ttree->Write();
  outfile->Close();
  delete outfile;
  return Fun4AllReturnCodes::EVENT_OK;
}
