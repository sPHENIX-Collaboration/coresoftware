#include "DumpEpInfo.h"

#include <eventplane/EpInfo.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<EpInfo>;

DumpEpInfo::DumpEpInfo(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpEpInfo::process_Node(PHNode *myNode)
{
  EpInfo *epinfo = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    epinfo = thisNode->getData();
  }
  if (epinfo)
  {
    *fout << "Max Order: " << epinfo->MaxOrder() << std::endl;
    // order starts at 1
    for (unsigned int i = 1; i < epinfo->MaxOrder(); ++i)
    {
      std::pair<double, double> tmp = epinfo->RawQ(i);
      *fout << "RawQ( " << i << "): " << tmp.first << ", " << tmp.second << std::endl;
      *fout << "SWRaw( " << i << "): " << epinfo->SWRaw(i) << std::endl;
      *fout << "RawPsi( " << i << "): " << epinfo->RawPsi(i) << std::endl;
    }
  }
  return 0;
}
