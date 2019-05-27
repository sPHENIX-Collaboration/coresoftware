#include "DumpBbcVertexMap.h"

#include <g4bbc/BbcVertex.h>
#include <g4bbc/BbcVertexMap.h>

#include <phool/PHIODataNode.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using namespace std;

typedef PHIODataNode<BbcVertexMap> MyNode_t;

DumpBbcVertexMap::DumpBbcVertexMap(const string &NodeName)
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
    *fout << "size: " << bbcvertexmap->size() << endl;
    for (BbcVertexMap::ConstIter biter = biter_beg; biter != biter_end; ++biter)
    {
      *fout << "id: " << biter->second->get_id() << endl;
      *fout << "t: " << biter->second->get_t() << endl;
      *fout << "t_err: " << biter->second->get_t_err() << endl;
      *fout << "z: " << biter->second->get_z() << endl;
      *fout << "z_err: " << biter->second->get_z_err() << endl;
    }
  }
  return 0;
}
