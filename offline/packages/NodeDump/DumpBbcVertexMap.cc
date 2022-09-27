#include "DumpBbcVertexMap.h"

#include <g4bbc/BbcVertex.h>
#include <g4bbc/BbcVertexMap.h>

#include <phool/PHIODataNode.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<BbcVertexMap>;

DumpBbcVertexMap::DumpBbcVertexMap(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpBbcVertexMap::process_Node(PHNode *myNode)
{
  BbcVertexMap *bbcvertexmap = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    bbcvertexmap = thisNode->getData();
  }
  if (bbcvertexmap)
  {
    BbcVertexMap::ConstIter biter_beg = bbcvertexmap->begin();
    BbcVertexMap::ConstIter biter_end = bbcvertexmap->end();
    *fout << "size: " << bbcvertexmap->size() << std::endl;
    for (BbcVertexMap::ConstIter biter = biter_beg; biter != biter_end; ++biter)
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
