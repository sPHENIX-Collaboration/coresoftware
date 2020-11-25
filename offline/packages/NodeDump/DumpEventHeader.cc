#include "DumpEventHeader.h"

#include "DumpObject.h"

#include <ffaobjects/EventHeader.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using namespace std;

typedef PHIODataNode<EventHeader> MyNode_t;

DumpEventHeader::DumpEventHeader(const string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpEventHeader::process_Node(PHNode *myNode)
{
  EventHeader *eventheader = 0;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    eventheader = thisNode->getData();
  }
  if (eventheader)
  {
    *fout << "EventHeader->isValid(): " << eventheader->isValid() << endl;
    if (eventheader->isValid())
    {
      *fout << "get_EvtSequence(): " << eventheader->get_EvtSequence() << endl;
      *fout << "get_EvtType(): " << eventheader->get_EvtType() << endl;
      *fout << "get_BunchCrossing(): " << eventheader->get_BunchCrossing() << endl;
      *fout << "get_ncoll(): " << eventheader->get_ncoll() << endl;
      *fout << "get_npart(): " << eventheader->get_npart() << endl;
      *fout << "get_ImpactParameter(): " << eventheader->get_ImpactParameter() << endl;
      *fout << "get_EventPlaneAngle(): " << eventheader->get_EventPlaneAngle() << endl;
      *fout << "get_eccentricity(): " << eventheader->get_eccentricity() << endl;
    }
  }
  return 0;
}
