#include "DumpGlobalVertexMap.h"

#include <g4vertex/GlobalVertex.h>
#include <g4vertex/GlobalVertexMap.h>

#include <phool/PHIODataNode.h>

#include <climits>
#include <string>

using namespace std;

typedef PHIODataNode<GlobalVertexMap> MyNode_t;

DumpGlobalVertexMap::DumpGlobalVertexMap(const string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpGlobalVertexMap::process_Node(PHNode *myNode)
{
  GlobalVertexMap *globalvertexmap = NULL;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    globalvertexmap = thisNode->getData();
  }
  if (globalvertexmap)
  {
    GlobalVertexMap::ConstIter viter_beg = globalvertexmap->begin();
    GlobalVertexMap::ConstIter viter_end = globalvertexmap->end();
    *fout << "size: " << globalvertexmap->size() << endl;
    for (GlobalVertexMap::ConstIter viter = viter_beg; viter != viter_end; ++viter)
    {
      *fout << "id: " << viter->second->get_id() << endl;
      *fout << "t: " << viter->second->get_t() << endl;
      *fout << "t_err: " << viter->second->get_t_err() << endl;
      *fout << "x: " << viter->second->get_x() << endl;
      *fout << "y: " << viter->second->get_y() << endl;
      *fout << "z: " << viter->second->get_z() << endl;
      *fout << "chisq: " << viter->second->get_chisq() << endl;
      *fout << "ndor: " << viter->second->get_ndof() << endl;
      for (int i = 0; i < 3; i++)
      {
        for (int j = 0; j < 3; j++)
        {
          *fout << "err[" << i << "][" << j << "]: " << viter->second->get_error(i, j) << endl;
        }
      }
      GlobalVertex::ConstVtxIter vtxbegin = viter->second->begin_vtxids();
      GlobalVertex::ConstVtxIter vtxend = viter->second->end_vtxids();
      for (GlobalVertex::ConstVtxIter vtxiter = vtxbegin; vtxiter != vtxend; ++vtxiter)
      {
        *fout << "type " << vtxiter->first << " id: " << vtxiter->second << endl;
      }
    }
  }
  return 0;
}
