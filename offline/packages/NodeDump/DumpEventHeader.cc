#include "DumpEventHeader.h"

#include "DumpObject.h"

#include <ffaobjects/EventHeader.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<EventHeader>;

DumpEventHeader::DumpEventHeader(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpEventHeader::process_Node(PHNode *myNode)
{
  EventHeader *eventheader = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    eventheader = thisNode->getData();
  }
  if (eventheader)
  {
    *fout << "EventHeader->isValid(): " << eventheader->isValid() << std::endl;
    if (eventheader->isValid())
    {
      *fout << "get_EvtSequence(): " << eventheader->get_EvtSequence() << std::endl;
      *fout << "get_EvtType(): " << eventheader->get_EvtType() << std::endl;
      *fout << "get_BunchCrossing(): " << eventheader->get_BunchCrossing() << std::endl;
      *fout << "get_ncoll(): " << eventheader->get_ncoll() << std::endl;
      *fout << "get_npart(): " << eventheader->get_npart() << std::endl;
      *fout << "get_TimeStamp(): " << eventheader->get_TimeStamp() << std::endl;
      *fout << "get_ImpactParameter(): " << eventheader->get_ImpactParameter() << std::endl;
      *fout << "get_EventPlaneAngle(): " << eventheader->get_EventPlaneAngle() << std::endl;
      *fout << "get_eccentricity(): " << eventheader->get_eccentricity() << std::endl;
    }
  }
  return 0;
}
