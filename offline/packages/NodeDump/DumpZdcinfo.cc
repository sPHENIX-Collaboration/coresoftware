#include "DumpZdcinfo.h"

#include <zdcinfo/Zdcinfo.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<Zdcinfo>;

DumpZdcinfo::DumpZdcinfo(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpZdcinfo::process_Node(PHNode *myNode)
{
  Zdcinfo *zdcinfo = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    zdcinfo = thisNode->getData();
  }
  if (zdcinfo && zdcinfo->isValid())
  {
    *fout << "Zdcinfo->get_zvertex: " << zdcinfo->get_zvertex() << std::endl;
    for (int i=0; i<2; i++)
    {
      *fout << "Zdcinfo->get_zdc_energy(" << i << "): " << zdcinfo->get_zdc_energy(i) << std::endl;
    *fout << "Zdcinfo->get_radius(" << i << "): " << zdcinfo->get_radius(i) << std::endl;
    }
  }
  return 0;
}
