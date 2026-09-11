#include "DumpLaserEventInfo.h"

#include <phool/PHIODataNode.h>

#include <tpc/LaserEventInfo.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<LaserEventInfo>;

DumpLaserEventInfo::DumpLaserEventInfo(const std::string &NodeName)
  : DumpObject(NodeName)
{
}

int DumpLaserEventInfo::process_Node(PHNode *myNode)
{
  LaserEventInfo *info{nullptr};
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);  // NOLINT(cppcoreguidelines-pro-type-static-cast-downcast)
  if (thisNode)
  {
    info = thisNode->getData();
  }
  if (info)
  {
    *fout << "isLaserEvent: " << info->isLaserEvent() << std::endl;
    *fout << "isGl1LaserEvent: " << info->isGl1LaserEvent() << std::endl;
    *fout << "isGl1LaserPileupEvent: " << info->isGl1LaserPileupEvent() << std::endl;
    for (const bool side : {false, true})
    {
      *fout << "getPeakSample(" << side << "): " << info->getPeakSample(side) << std::endl;
      *fout << "getPeakWidth(" << side << "): " << info->getPeakWidth(side) << std::endl;
    }
  }
  return 0;
}
