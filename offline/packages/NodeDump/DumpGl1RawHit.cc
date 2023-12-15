#include "DumpGl1RawHit.h"

#include <ffarawobjects/Gl1RawHit.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<Gl1RawHit>;

DumpGl1RawHit::DumpGl1RawHit(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpGl1RawHit::process_Node(PHNode *myNode)
{
  Gl1RawHit *gl1rawhit = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    gl1rawhit = thisNode->getData();
  }
  if (gl1rawhit)
  {
    *fout << "get_bco: " << gl1rawhit->get_bco() << std::endl;
  }
  return 0;
}
