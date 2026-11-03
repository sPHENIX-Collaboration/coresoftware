#include "DumpTriggerRunInfo.h"

#include "DumpObject.h"

#include <calotrigger/TriggerRunInfo.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<TriggerRunInfo>;

DumpTriggerRunInfo::DumpTriggerRunInfo(const std::string &NodeName)
  : DumpObject(NodeName)
{
  WriteRunEvent(0);  // do not write info for each event
  return;
}

int DumpTriggerRunInfo::process_Node(PHNode *myNode)
{
  TriggerRunInfo *trigrun = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);  // NOLINT(cppcoreguidelines-pro-type-static-cast-downcast)
  if (thisNode)
  {
    trigrun = thisNode->getData();
  }
  if (!node_written && trigrun)
  {
    for (int j = 0; j<64; j++)
    {
      *fout << "TriggerRunInfo->getPrescaleByBit(" << j<< "): " << trigrun->getPrescaleByBit(j) << std::endl;
      *fout << "TriggerRunInfo->getInitialPrescaleByBit(" << j<< "): " << trigrun->getInitialPrescaleByBit(j) << std::endl;
      *fout << "TriggerRunInfo->getScalersByBit(" << j<< "): " << trigrun->getScalersByBit(j) << std::endl;
      *fout << "TriggerRunInfo->getLiveScalersByBit(" << j<< "): " << trigrun->getLiveScalersByBit(j) << std::endl;
      *fout << "TriggerRunInfo->getRawScalersByBit(" << j<< "): " << trigrun->getRawScalersByBit(j) << std::endl;
      *fout << "TriggerRunInfo->getTriggerName(" << j<< "): " << trigrun->getTriggerName(j) << std::endl;
    }
    node_written = 1;
  }
  return 0;
}
