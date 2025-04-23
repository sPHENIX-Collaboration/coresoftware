#include "InttGl1Check.h"

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <ffarawobjects/Gl1Packet.h>
#include <ffarawobjects/InttRawHit.h>
#include <ffarawobjects/InttRawHitContainer.h>

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
InttGl1Check::InttGl1Check(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int InttGl1Check::Init(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InttGl1Check::process_event(PHCompositeNode *topNode)
{
  InttRawHitContainer *inttcont = findNode::getClass<InttRawHitContainer>(topNode, m_EvtNodeName);
  Gl1Packet *gl1cont = findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");
  if (!gl1cont)
  {
    std::cout << "could not find GL1RAWHIT node" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if (!inttcont)
  {
    std::cout << "could not find node " << m_EvtNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  uint64_t gl1bco = (gl1cont->getBCO() & 0xFFFFFFFFFFU);  // 40 bit intt bco
                                                          //  inttcont->identify();
  std::set<uint64_t> inttbcoset;
  std::set<int64_t> diffgl1;
  //  std::map<uint64_t, std::tuple<uint64_t, uint64_t, uint64_t>> clocktuplemap;
  for (unsigned int i = 0; i < inttcont->get_nhits(); i++)
  {
    InttRawHit *inh = inttcont->get_hit(i);
    //      int64_t diffbck = inh->get_bco() - gl1bco;
    uint64_t fphxbco = inh->get_FPHX_BCO();
    uint64_t inttgl1 = inh->get_bco() + fphxbco;
    int64_t diffbck = inttgl1 - gl1bco;
    if (diffbck < 0)
    {
      continue;
    }
    diffgl1.insert(diffbck);
    //      int bco_diff = (inh->get_FPHX_BCO() - (inh->get_bco() & 0x7fU) + 128) % 128;
    // std::cout << "gl1 bco: 0x" << std::hex << gl1bco
    // 		<< ", intt: 0x" << inh->get_bco()
    // 		<< ", fphx: 0x" << fphxbco
    // 		<< ", intt corr: 0x" << inttgl1
    // 		<< ", diff: 0x" << diffbck
    // 		<< ", bcodiff: 0x" << bco_diff
    // 		<< std::dec << std::endl;
    diffcnt[diffbck]++;
    // std::cout << "b_intt: " << b_intt << std::endl;
    // std::cout << "b_cntt: " << b_cntt << std::endl;
    // std::cout << "b_fphx: " << b_fphx<< std::endl;
    // std::cout << "b_sums: " << b_inttgl1 << std::endl;
    //      clocktuplemap[gl1bco] = std::make_tuple(inh->get_bco(), chopbco, fphxbco,
    //      std::cout << "diffgl1: 0x" << std::hex << diffbck << std::dec << std::endl;
    inttbcoset.insert(inh->get_bco());
  }
  // for (auto &iter : diffgl1)
  // {
  //   std::cout << "diff with gl1: 0x" << std::hex << iter << std::dec << std::endl;
  // }
  std::set<uint64_t> diffs;
  for (const auto &iter : inttbcoset)
  {
    uint64_t refbco = iter;
    for (const auto &iter1 : inttbcoset)
    {
      if (iter1 <= refbco)
      {
        continue;
      }
      diffs.insert(iter1 - refbco);
    }
  }
  // for (auto &iter : diffs)
  // {
  //   std::cout << "diff: 0x" << std::hex <<  iter << std::dec << std::endl;
  // }
  // for (auto &iter : diffcnt)
  // {
  //   std::cout << "diff " << iter.first << " count " << iter.second << std::endl;
  // }
  return Fun4AllReturnCodes::EVENT_OK;
}

int InttGl1Check::End(PHCompositeNode * /*topNode*/)
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
