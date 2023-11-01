#include "DumpInttRawHitContainer.h"

#include <phool/PHIODataNode.h>

#include <ffarawobjects/InttRawHit.h>
#include <ffarawobjects/InttRawHitContainer.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<InttRawHitContainer>;

DumpInttRawHitContainer::DumpInttRawHitContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpInttRawHitContainer::process_Node(PHNode *myNode)
{
  InttRawHitContainer *inttrawhitcontainer = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    inttrawhitcontainer = thisNode->getData();
  }
  if (inttrawhitcontainer)
  {
    unsigned int nhits = inttrawhitcontainer->get_nhits();
    *fout << "size: " << inttrawhitcontainer->get_nhits() << std::endl;
    for (unsigned int ihit = 0; ihit < nhits; ihit++)
    {
      InttRawHit *rawhit = inttrawhitcontainer->get_hit(ihit);
      *fout << "bco: " << rawhit->get_bco() << std::endl;
      *fout << "packetid: " << rawhit->get_packetid() << std::endl;
      *fout << "word: " << rawhit->get_word() << std::endl;
      *fout << "fee: " << rawhit->get_fee() << std::endl;
      *fout << "channel_id: " << rawhit->get_channel_id() << std::endl;
      *fout << "chip_id: " << rawhit->get_chip_id() << std::endl;
      *fout << "adc: " << rawhit->get_adc() << std::endl;
      *fout << "FPHX_BCO: " << rawhit->get_FPHX_BCO() << std::endl;
      *fout << "full_FPHX: " << rawhit->get_full_FPHX() << std::endl;
      *fout << "full_ROC: " << rawhit->get_full_ROC() << std::endl;
      *fout << "amplitude: " << rawhit->get_amplitude() << std::endl;
    }
  }
  return 0;
}
