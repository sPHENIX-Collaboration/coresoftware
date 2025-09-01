#include "DumpTruthVertexMap.h"

#include <globalvertex/TruthVertex.h>
#include <globalvertex/TruthVertexMap.h>

#include <phool/PHIODataNode.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<TruthVertexMap>;

DumpTruthVertexMap::DumpTruthVertexMap(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpTruthVertexMap::process_Node(PHNode *myNode)
{
  const auto original_precision = (*fout).precision();
  TruthVertexMap *truthvertexmap = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);  // NOLINT(cppcoreguidelines-pro-type-static-cast-downcast)
  if (thisNode)
  {
    truthvertexmap = thisNode->getData();
  }
  if (truthvertexmap)
  {
    (*fout).precision(std::numeric_limits<float>::max_digits10);
    TruthVertexMap::ConstIter biter_beg = truthvertexmap->begin();
    TruthVertexMap::ConstIter biter_end = truthvertexmap->end();
    *fout << "size: " << truthvertexmap->size() << std::endl;
    for (TruthVertexMap::ConstIter biter = biter_beg; biter != biter_end; ++biter)
    {
      *fout << "id: " << biter->second->get_id() << std::endl;
      *fout << "t: " << biter->second->get_t() << std::endl;
      *fout << "t_err: " << biter->second->get_t_err() << std::endl;
      *fout << "x: " << biter->second->get_x() << std::endl;
      *fout << "x_err: " << biter->second->get_x_err() << std::endl;
      *fout << "y: " << biter->second->get_y() << std::endl;
      *fout << "y_err: " << biter->second->get_y_err() << std::endl;
      *fout << "z: " << biter->second->get_z() << std::endl;
      *fout << "z_err: " << biter->second->get_z_err() << std::endl;
    }
    (*fout).precision(original_precision);
  }
  return 0;
}
