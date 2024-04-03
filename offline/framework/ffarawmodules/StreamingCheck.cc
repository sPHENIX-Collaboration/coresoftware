#include "StreamingCheck.h"

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <ffarawobjects/Gl1RawHit.h>
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

#include <TSystem.h>

#include <iostream>  // for operator<<, endl, basic_ost...
#include <set>
#include <utility>  // for pair
#include <vector>   // for vector

//____________________________________________________________________________..
StreamingCheck::StreamingCheck(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int StreamingCheck::Init(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int StreamingCheck::process_event(PHCompositeNode *topNode)
{
  Gl1RawHit *gl1rawhit = findNode::getClass<Gl1RawHit>(topNode, "GL1RAWHIT");
  if (!gl1rawhit)
  {
    std::cout << "could not find node GL1RAWHIT" << std::endl;
    exit(1);
  }
  uint64_t refBCO = gl1rawhit->get_bco();
  uint64_t refBCO_40Bits = refBCO & 0xFFFFFFFFFFU;
  TpcRawHitContainer *tpccont = findNode::getClass<TpcRawHitContainer>(topNode, "TPCRAWHIT");
  if (tpccont)
  {
    if (Verbosity() > 0)
    {
      tpccont->identify();
    }
    for (unsigned int i = 0; i < tpccont->get_nhits(); i++)
    {
      TpcRawHit *inh = tpccont->get_hit(i);
      if (inh->get_gtm_bco() < refBCO_40Bits)
      {
        std::cout << "bco mismatch (too small), gl1: 0x" << std::hex << refBCO_40Bits
                  << ", tpc: 0x" << inh->get_gtm_bco() << std::endl;
      }
      if (inh->get_gtm_bco() > refBCO_40Bits + tpc_bcorange)
      {
        std::cout << "bco mismatch (too large), gl1: 0x" << std::hex << refBCO_40Bits
                  << ", tpc: 0x" << inh->get_gtm_bco() << std::endl;
      }
    }
  }
  InttRawHitContainer *inttcont = findNode::getClass<InttRawHitContainer>(topNode, "INTTRAWHIT");
  if (inttcont)
  {
    if (Verbosity() > 0)
    {
      inttcont->identify();
    }
    for (unsigned int i = 0; i < inttcont->get_nhits(); i++)
    {
      InttRawHit *inh = inttcont->get_hit(i);
      if (inh->get_bco() != refBCO_40Bits)
      {
        std::cout << "bco mismatch (too small), gl1: 0x" << std::hex << refBCO_40Bits
                  << ", intt: 0x" << inh->get_bco() << std::endl;
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
