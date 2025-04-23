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
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);  // NOLINT(cppcoreguidelines-pro-type-static-cast-downcast)
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
      *fout << "channel:" << rawhit->get_channel() << std::endl;
      *fout << "sampaaddress: " << rawhit->get_sampaaddress() << std::endl;
      *fout << "sampachannel: " << rawhit->get_sampachannel() << std::endl;
      *fout << "type:" << rawhit->get_type() << std::endl;
      *fout << "userword:" << rawhit->get_userword() << std::endl;
      *fout << "checksum:" << rawhit->get_checksum() << std::endl;
      *fout << "checksumerror:" << rawhit->get_checksumerror() << std::endl;
      *fout << "parity:" << rawhit->get_parity() << std::endl;
      *fout << "parityerror:" << rawhit->get_parityerror() << std::endl;
      auto nsamples = rawhit->get_samples();
      *fout << "samples: " << rawhit->get_samples() << std::endl;
      for (auto isamp = 0; isamp < nsamples; isamp++)
      {
        *fout << "adc[" << isamp << "] =  " << rawhit->get_adc(isamp) << std::endl;
      }
    }
  }
  return 0;
}
