#include "DumpGlobalVertexMap.h"

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

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
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);  // NOLINT(cppcoreguidelines-pro-type-static-cast-downcast)
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
      *fout << "first: " << viter->first << std::endl;
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
      *fout << "number of v1 vertices: " << viter->second->size_vtxids() << std::endl;
      GlobalVertex::ConstVtxIter vtxbegin = viter->second->begin_vtxids();
      GlobalVertex::ConstVtxIter vtxend = viter->second->end_vtxids();
      for (GlobalVertex::ConstVtxIter vtxiter = vtxbegin; vtxiter != vtxend; ++vtxiter)
      {
        *fout << "type " << vtxiter->first << " id: " << vtxiter->second << std::endl;
      }

      *fout << "number of v2 vertices: " << viter->second->size_vtxs() << std::endl;
      for (GlobalVertex::ConstVertexIter vtxiter = viter->second->begin_vertexes();
           vtxiter != viter->second->end_vertexes(); ++vtxiter)
      {
        *fout << "Vertex Type: " << vtxiter->first << std::endl;
        for (const auto *iter : vtxiter->second)
        {
          *fout << "id: " << iter->get_id() << std::endl;
          *fout << "beam_crossing: " << iter->get_beam_crossing() << std::endl;
          *fout << "t: " << iter->get_t() << std::endl;
          *fout << "t_err: " << iter->get_t_err() << std::endl;
          *fout << "x: " << iter->get_x() << std::endl;
          *fout << "y: " << iter->get_y() << std::endl;
          *fout << "z: " << iter->get_z() << std::endl;
          *fout << "chisq: " << iter->get_chisq() << std::endl;
          *fout << "ndof: " << iter->get_ndof() << std::endl;
          for (int i = 0; i < 3; i++)
          {
            for (int j = 0; j < 3; j++)
            {
              *fout << "err[" << i << "][" << j << "]: " << iter->get_error(i, j) << std::endl;
            }
          }
        }
      }
    }
  }
  return 0;
}
