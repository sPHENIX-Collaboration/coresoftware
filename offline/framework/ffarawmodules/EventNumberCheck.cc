#include "EventNumberCheck.h"

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>

#include <Event/Event.h>
#include <Event/oncsEvent.h>

#include <TSystem.h>

#include <iostream>  // for operator<<, endl, basic_ost...
#include <utility>   // for pair
#include <vector>    // for vector

//____________________________________________________________________________..
EventNumberCheck::EventNumberCheck(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int EventNumberCheck::Init(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int EventNumberCheck::process_event(PHCompositeNode *topNode)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  Event *evt = findNode::getClass<Event>(topNode, m_MyPrdfNode);
  evt->identify();
  int eventno = evt->getEvtSequence();
  int nw = evt->getPacketList(plist, 10000);
  if (nw >= 10000)
  {
    std::cout << "Packet array too small, need " << nw << " entries" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }
  auto insert_chk = m_EventSeen.insert(eventno);
  if (!insert_chk.second)
  {
    std::cout << "event " << eventno << " exists already"
              << " counter: " << se->EventNumber() << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void EventNumberCheck::CheckFem(int nw)
{
  std::set<int> femclkcemc, clkcemc;
  std::set<int> femclkmbd, clkmbd;
  static int ifirst = 1;
  for (int i = 0; i < nw; i++)
  {
    int pktid = plist[i]->getIdentifier();
    if (pktid > 2000)
    {
      clkcemc.insert(plist[i]->iValue(0, "CLOCK"));
    }
    else
    {
      clkmbd.insert(plist[i]->iValue(0, "CLOCK"));
    }
    if (Verbosity() > 1)
    {
      std::cout << "packet " << plist[i]->getIdentifier() << ", evt nr "
                << plist[i]->iValue(0, "EVTNR") << ", bclk 0x" << std::hex
                << plist[i]->iValue(0, "CLOCK") << std::dec << std::endl;
    }
    for (int j = 0; j < plist[i]->iValue(0, "NRMODULES"); j++)
    {
      if (Verbosity() > 1)
      {
        std::cout << "FEM " << j << ", Clock 0x" << std::hex
                  << plist[i]->iValue(j, "FEMCLOCK") << std::dec << std::endl;
      }
      if (pktid > 2000)
      {
        femclkcemc.insert(plist[i]->iValue(j, "FEMCLOCK"));
      }
      else
      {
        femclkmbd.insert(plist[i]->iValue(j, "FEMCLOCK"));
      }
    }
    delete plist[i];
  }
  if (femclkcemc.size() > 1)
  {
    std::cout << "CEMC FEM clock mismatch, saw " << std::hex << std::endl;
    for (auto iter : femclkcemc)
    {
      std::cout << iter << std::endl;
    }
    std::cout << std::dec;
  }
  if (femclkmbd.size() > 1)
  {
    std::cout << "MBD FEM clock mismatch, saw " << std::hex << std::endl;
    for (auto iter : femclkmbd)
    {
      std::cout << iter << std::endl;
    }
    std::cout << std::dec;
  }
  if (clkcemc.size() > 1)
  {
    std::cout << "CEMC Packet clock mismatch, saw " << std::hex << std::endl;
    for (auto iter : clkcemc)
    {
      std::cout << iter << std::endl;
    }
    std::cout << std::dec;
  }
  if (clkmbd.size() > 1)
  {
    std::cout << "MBD Packet clock mismatch, saw " << std::hex << std::endl;
    for (auto iter : clkmbd)
    {
      std::cout << iter << std::endl;
    }
    std::cout << std::dec;
  }

  int femclockcemc = *(clkcemc.begin());
  int femclockmbd = *(clkmbd.begin());
  if (previous_event_clkdiff != femclockcemc - femclockmbd)
  {
    if (ifirst)
    {
      ifirst = 0;
    }
    else
    {
      std::cout << "clock diff changed at event "  //<< eventno
                << " from " << previous_event_clkdiff << " to "
                << femclockcemc - femclockmbd << std::endl;
    }
    previous_event_clkdiff = femclockcemc - femclockmbd;
  }
  return;
}
