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
#include <TNtuple.h>
#include <TSystem.h>

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

  outTfile = new TFile(outfilename.c_str(), "RECREATE");
  ntup = new TNtuple("bco", "bco", "id:evt:bco:bcodiff");
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
  for (auto packet : pktvec)
  {
    int packetid = packet->getIdentifier();
    lastbco.insert(std::make_pair(packetid,0));
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
      for (auto &bco : mapiter.second)
      {
        uint64_t prevbco = lastbco[mapiter.first];
        if (prevbco > 0 && prevbco != bco)
        {
          int diffbco = bco - prevbco;
          ntup->Fill(mapiter.first, EventSequence, bco, diffbco);
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
  outTfile->cd();
  ntup->Write();
  outTfile->Close();
  delete outTfile;
  return Fun4AllReturnCodes::EVENT_OK;
}
