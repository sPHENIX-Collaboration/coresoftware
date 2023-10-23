#include "DumpTpcRawHitContainer.h"

#include <phool/PHIODataNode.h>

#include <ffarawobjects/TpcRawHit.h>
#include <ffarawobjects/TpcRawHitContainer.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<TpcRawHitContainer>;

DumpTpcRawHitContainer::DumpTpcRawHitContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpTpcRawHitContainer::process_Node(PHNode *myNode)
{
  TpcRawHitContainer *tpcrawhitcontainer = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    tpcrawhitcontainer = thisNode->getData();
  }
  if (tpcrawhitcontainer)
  {
    unsigned int nhits = tpcrawhitcontainer->get_nhits();
    *fout << "size: " << tpcrawhitcontainer->get_nhits() << std::endl;
    for (unsigned int ihit = 0; ihit < nhits; ihit++)
    {
      TpcRawHit *rawhit = tpcrawhitcontainer->get_hit(ihit);
      *fout << "bco: " << rawhit->get_bco() << std::endl;
      *fout << "gtm_bco: " << rawhit->get_gtm_bco() << std::endl;
      *fout << "packetid: " << rawhit->get_packetid() << std::endl;
      *fout << "fee: " << rawhit->get_fee() << std::endl;
      *fout << "sampaaddress: " << rawhit->get_sampaaddress() << std::endl;
      *fout << "sampachannel: " << rawhit->get_sampachannel() << std::endl;
      auto nsamples = rawhit->get_samples();
      *fout << "samples: " << rawhit->get_samples() << std::endl;
      for (auto isamp = 0; isamp < nsamples; isamp++)
      {
	*fout << "adc[" << isamp << "] =  " <<  rawhit->get_adc(isamp) << std::endl;
      }
    }
  }
  return 0;
}
