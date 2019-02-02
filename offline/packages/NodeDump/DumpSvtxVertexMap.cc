#include "DumpSvtxVertexMap.h"

#include <phool/PHIODataNode.h>

#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <string>

using namespace std;

typedef PHIODataNode<SvtxVertexMap> MyNode_t;

DumpSvtxVertexMap::DumpSvtxVertexMap(const string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpSvtxVertexMap::process_Node(PHNode *myNode)
{
  SvtxVertexMap *svtxvertexmap = NULL;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    svtxvertexmap = thisNode->getData();
  }
  if (svtxvertexmap)
  {
    SvtxVertexMap::ConstIter hiter;
    *fout << "size: " << svtxvertexmap->size() << endl;
    for (hiter = svtxvertexmap->begin(); hiter != svtxvertexmap->end(); hiter++)
    {
      *fout << "id: 0x" << hex << hiter->second->get_id() << dec << endl;
      *fout << "t0: " << hiter->second->get_t0() << endl;
      *fout << "x: " << hiter->second->get_x() << endl;
      *fout << "y: " << hiter->second->get_y() << endl;
      *fout << "z: " << hiter->second->get_z() << endl;
      *fout << "chisq: " << hiter->second->get_chisq() << endl;
      *fout << "ndof: " << hiter->second->get_ndof() << endl;
    }
  }
  return 0;
}
