#include "DumpMvtxRawHitContainer.h"

#include <phool/PHIODataNode.h>

#include <ffarawobjects/MvtxRawHit.h>
#include <ffarawobjects/MvtxRawHitContainer.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<MvtxRawHitContainer>;

DumpMvtxRawHitContainer::DumpMvtxRawHitContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpMvtxRawHitContainer::process_Node(PHNode *myNode)
{
  MvtxRawHitContainer *mvtxrawhitcontainer = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    mvtxrawhitcontainer = thisNode->getData();
  }
  if (mvtxrawhitcontainer)
  {
    unsigned int nhits = mvtxrawhitcontainer->get_nhits();
    *fout << "size: " << mvtxrawhitcontainer->get_nhits() << std::endl;
    for (unsigned int ihit = 0; ihit < nhits; ihit++)
    {
      MvtxRawHit *rawhit = mvtxrawhitcontainer->get_hit(ihit);
      *fout << "bco: " << rawhit->get_bco() << std::endl;
      *fout << "strobe_bc: " << rawhit->get_strobe_bc() << std::endl;
      *fout << "chip_bc: " << rawhit->get_chip_bc() << std::endl;
      *fout << "layer_id: " << rawhit->get_layer_id() << std::endl;
      *fout << "stave_id: " << rawhit->get_stave_id() << std::endl;
      *fout << "chip_id: " << rawhit->get_chip_id() << std::endl;
      *fout << "row: " << rawhit->get_row() << std::endl;
      *fout << "col: " << rawhit->get_col() << std::endl;
    }
  }
  return 0;
}
