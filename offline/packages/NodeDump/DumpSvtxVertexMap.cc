#include "DumpSvtxVertexMap.h"

#include <phool/PHIODataNode.h>

#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<SvtxVertexMap>;

DumpSvtxVertexMap::DumpSvtxVertexMap(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpSvtxVertexMap::process_Node(PHNode *myNode)
{
  SvtxVertexMap *svtxvertexmap = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    svtxvertexmap = thisNode->getData();
  }
  if (svtxvertexmap)
  {
    SvtxVertexMap::ConstIter hiter;
    *fout << "size: " << svtxvertexmap->size() << std::endl;
    for (hiter = svtxvertexmap->begin(); hiter != svtxvertexmap->end(); hiter++)
    {
      *fout << "id: 0x" << std::hex << hiter->second->get_id() << std::dec << std::endl;
      *fout << "t0: " << hiter->second->get_t0() << std::endl;
      *fout << "x: " << hiter->second->get_x() << std::endl;
      *fout << "y: " << hiter->second->get_y() << std::endl;
      *fout << "z: " << hiter->second->get_z() << std::endl;
      *fout << "chisq: " << hiter->second->get_chisq() << std::endl;
      *fout << "ndof: " << hiter->second->get_ndof() << std::endl;
    }
  }
  return 0;
}
