#include "DumpGlobalVertexMap.h"

#include <g4vertex/GlobalVertex.h>
#include <g4vertex/GlobalVertexMap.h>

#include <phool/PHIODataNode.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<GlobalVertexMap>;

DumpGlobalVertexMap::DumpGlobalVertexMap(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpGlobalVertexMap::process_Node(PHNode *myNode)
{
  GlobalVertexMap *globalvertexmap = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    globalvertexmap = thisNode->getData();
  }
  if (globalvertexmap)
  {
    GlobalVertexMap::ConstIter viter_beg = globalvertexmap->begin();
    GlobalVertexMap::ConstIter viter_end = globalvertexmap->end();
    *fout << "size: " << globalvertexmap->size() << std::endl;
    for (GlobalVertexMap::ConstIter viter = viter_beg; viter != viter_end; ++viter)
    {
      *fout << "id: " << viter->second->get_id() << std::endl;
      *fout << "t: " << viter->second->get_t() << std::endl;
      *fout << "t_err: " << viter->second->get_t_err() << std::endl;
      *fout << "x: " << viter->second->get_x() << std::endl;
      *fout << "y: " << viter->second->get_y() << std::endl;
      *fout << "z: " << viter->second->get_z() << std::endl;
      *fout << "chisq: " << viter->second->get_chisq() << std::endl;
      *fout << "ndor: " << viter->second->get_ndof() << std::endl;
      for (int i = 0; i < 3; i++)
      {
        for (int j = 0; j < 3; j++)
        {
          *fout << "err[" << i << "][" << j << "]: " << viter->second->get_error(i, j) << std::endl;
        }
      }
      GlobalVertex::ConstVtxIter vtxbegin = viter->second->begin_vtxids();
      GlobalVertex::ConstVtxIter vtxend = viter->second->end_vtxids();
      for (GlobalVertex::ConstVtxIter vtxiter = vtxbegin; vtxiter != vtxend; ++vtxiter)
      {
        *fout << "type " << vtxiter->first << " id: " << vtxiter->second << std::endl;
      }
    }
  }
  return 0;
}
