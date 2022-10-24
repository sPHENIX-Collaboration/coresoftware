#include "DumpCentralityInfo.h"

#include "DumpObject.h"

#include <centrality/CentralityInfo.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<CentralityInfo>;

DumpCentralityInfo::DumpCentralityInfo(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpCentralityInfo::process_Node(PHNode *myNode)
{
  CentralityInfo *centralityinfo = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    centralityinfo = thisNode->getData();
  }
  if (centralityinfo)
  {
    *fout << "CentralityInfo->isValid(): " << centralityinfo->isValid() << std::endl;
    if (centralityinfo->isValid())
    {
      for (int i = 0; i < 100; i++)
      {
        if (centralityinfo->has_quantity(static_cast<CentralityInfo::PROP>(i)))
        {
          *fout << "get_quantity(" << i << "): " << centralityinfo->get_quantity(static_cast<CentralityInfo::PROP>(i)) << std::endl;
        }
        if (centralityinfo->has_centile(static_cast<CentralityInfo::PROP>(i)))
        {
          *fout << "get_centile(" << i << "): " << centralityinfo->get_centile(static_cast<CentralityInfo::PROP>(i)) << std::endl;
        }
      }
    }
  }
  return 0;
}
