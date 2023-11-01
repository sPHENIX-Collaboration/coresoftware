#include "DumpMbdVertexMap.h"

#include <mbd/MbdVertex.h>
#include <mbd/MbdVertexMap.h>

#include <phool/PHIODataNode.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<MbdVertexMap>;

DumpMbdVertexMap::DumpMbdVertexMap(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpMbdVertexMap::process_Node(PHNode *myNode)
{
  MbdVertexMap *mbdvertexmap = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    mbdvertexmap = thisNode->getData();
  }
  if (mbdvertexmap)
  {
    MbdVertexMap::ConstIter biter_beg = mbdvertexmap->begin();
    MbdVertexMap::ConstIter biter_end = mbdvertexmap->end();
    *fout << "size: " << mbdvertexmap->size() << std::endl;
    for (MbdVertexMap::ConstIter biter = biter_beg; biter != biter_end; ++biter)
    {
      *fout << "id: " << biter->second->get_id() << std::endl;
      *fout << "t: " << biter->second->get_t() << std::endl;
      *fout << "t_err: " << biter->second->get_t_err() << std::endl;
      *fout << "z: " << biter->second->get_z() << std::endl;
      *fout << "z_err: " << biter->second->get_z_err() << std::endl;
    }
  }
  return 0;
}
