#include "BcoRangeCheck.h"

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <ffarawobjects/Gl1Packet.h>
#include <ffarawobjects/InttRawHit.h>
#include <ffarawobjects/InttRawHitContainer.h>
#include <ffarawobjects/MicromegasRawHit.h>
#include <ffarawobjects/MicromegasRawHitContainer.h>
#include <ffarawobjects/MvtxRawHit.h>
#include <ffarawobjects/MvtxRawHitContainer.h>
#include <ffarawobjects/TpcRawHit.h>
#include <ffarawobjects/TpcRawHitContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>

#include <Event/packet.h>

#include <TSystem.h>

#include <bitset>
#include <iostream>  // for operator<<, endl, basic_ost...
#include <ranges>
#include <set>
#include <utility>  // for pair
#include <vector>   // for vector

//____________________________________________________________________________..
BcoRangeCheck::BcoRangeCheck(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int BcoRangeCheck::Init(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int BcoRangeCheck::process_event(PHCompositeNode *topNode)
{
  Gl1Packet *gl1cont = findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");
  if (!gl1cont)
  {
    std::cout << "could not find GL1RAWHIT node" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  InttRawHitContainer *inttcont = findNode::getClass<InttRawHitContainer>(topNode, "INTTRAWHIT");
  if (!inttcont)
  {
    std::cout << "could not find INTTRAWHIT node " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  MvtxRawHitContainer *mvtxcont = findNode::getClass<MvtxRawHitContainer>(topNode, "MVTXRAWHIT");
  if (!mvtxcont)
  {
    std::cout << "could not find MVTXRAWHIT node " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  TpcRawHitContainer *tpccont = findNode::getClass<TpcRawHitContainer>(topNode, "TPCRAWHIT");
  if (!tpccont)
  {
    std::cout << "could not find TPCRAWHIT node " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  MicromegasRawHitContainer *micromegascont = findNode::getClass<MicromegasRawHitContainer>(topNode, "MICROMEGASRAWHIT");
  if (!micromegascont)
  {
    std::cout << "could not find TPOTRAWHIT node " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  uint64_t gl1bco = (gl1cont->getBCO() & 0xFFFFFFFFFFU);  // 40 bit intt bco
                                                          //  inttcont->identify();
  std::set<uint64_t> inttbcoset;
  for (unsigned int i = 0; i < inttcont->get_nhits(); i++)
  {
    InttRawHit *inh = inttcont->get_hit(i);
    inttbcoset.insert(inh->get_bco());
  }
  std::set<uint64_t> mvtxbcoset;
  for (unsigned int i = 0; i < mvtxcont->get_nhits(); i++)
  {
    MvtxRawHit *inh = mvtxcont->get_hit(i);
    mvtxbcoset.insert(inh->get_bco());
  }
  std::set<uint64_t> tpcbcoset;
  for (unsigned int i = 0; i < tpccont->get_nhits(); i++)
  {
    TpcRawHit *inh = tpccont->get_hit(i);
    if (inh->get_bco() > 0)
    {
      tpcbcoset.insert(inh->get_bco());
    }
  }
  std::set<uint64_t> micromegasbcoset;
  for (unsigned int i = 0; i < micromegascont->get_nhits(); i++)
  {
    MicromegasRawHit *inh = micromegascont->get_hit(i);
    micromegasbcoset.insert(inh->get_bco());
  }
  if (!inttbcoset.empty())
  {
    std::cout << "inttrange: 0x" << std::hex << *inttbcoset.rbegin() << " to 0x" << *inttbcoset.begin()
              << std::dec << " diff: " << *inttbcoset.rbegin() - *inttbcoset.begin() << std::endl;
  }
  if (!mvtxbcoset.empty())
  {
    std::cout << "mvtxrange: 0x" << std::hex << *mvtxbcoset.rbegin() << " to 0x" << *mvtxbcoset.begin()
              << std::dec << " diff: " << *mvtxbcoset.rbegin() - *mvtxbcoset.begin() << std::endl;
  }
  if (!tpcbcoset.empty())
  {
    std::cout << "tpcrange: 0x" << std::hex << *tpcbcoset.rbegin() << " to 0x" << *tpcbcoset.begin()
              << std::dec << " diff: " << *tpcbcoset.rbegin() - *tpcbcoset.begin() << std::endl;
  }
  if (!micromegasbcoset.empty())
  {
    std::cout << "micromegasrange: 0x" << std::hex << *micromegasbcoset.rbegin() << " to 0x" << *micromegasbcoset.begin()
              << std::dec << " diff: " << *micromegasbcoset.rbegin() - *micromegasbcoset.begin() << std::endl;
  }
  std::cout << "Gl1: 0x" << std::hex << gl1bco << std::dec << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int BcoRangeCheck::End(PHCompositeNode * /*topNode*/)
{
  std::multimap<int, uint64_t> scoremap;
  for (auto &iter : diffcnt)
  {
    std::cout << "diff " << iter.first << " count " << iter.second << std::endl;
    scoremap.insert(std::make_pair(iter.second, iter.first));
  }
  int i = 0;
  for (auto &iter : std::ranges::reverse_view(scoremap))
  {
    std::cout << "high score " << iter.first << " for diff " << iter.second
              << std::endl;
    i++;
    if (i > 3)
    {
      break;
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
