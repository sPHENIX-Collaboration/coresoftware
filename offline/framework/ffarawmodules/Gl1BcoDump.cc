#include "Gl1BcoDump.h"

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
Gl1BcoDump::Gl1BcoDump(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int Gl1BcoDump::InitRun(PHCompositeNode * /*topNode*/)
{
  if (outfilename.empty())
  {
    std::cout << "no output filename given" << std::endl;
    gSystem->Exit(1);
  }

  outTfile = new TFile(outfilename.c_str(), "RECREATE");
  ntup = new TTree("bco", "bco");
  ntup->Branch("id",&m_id);
  ntup->Branch("evt",&m_evt);
  ntup->Branch("bco",&m_bco);
  ntup->Branch("bcodiff",&m_bcodiff);
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int Gl1BcoDump::process_event(PHCompositeNode *topNode)
{
  Event *evt = findNode::getClass<Event>(topNode, "PRDF");
  if (!evt)
  {
    std::cout << "No Event found" << std::endl;
    exit(1);
  }
  //  evt->identify();
  int EventSequence = evt->getEvtSequence();
  Packet *packet = evt->getPacket(14001);
  if (!packet)
  {
    // evt->identify();
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  uint64_t bco = (uint64_t) packet->lValue(0, "BCO") & 0xFFFFFFFFFFU;
  if (lastbco > 0)
  {
    int64_t diffbco = bco - lastbco;
    m_id = 14001;
    m_evt = EventSequence;
    m_bco = bco;
    m_bcodiff = diffbco;
    ntup->Fill();
  }
  lastbco = bco;
  delete packet;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int Gl1BcoDump::End(PHCompositeNode * /*topNode*/)
{
  outTfile->cd();
  ntup->Write();
  outTfile->Close();
  delete outTfile;
  return Fun4AllReturnCodes::EVENT_OK;
}
