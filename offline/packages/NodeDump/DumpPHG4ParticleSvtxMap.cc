#include "DumpPHG4ParticleSvtxMap.h"

#include <phool/PHIODataNode.h>

#include <trackbase_historic/PHG4ParticleSvtxMap.h>

#include <map>
#include <ostream>
#include <set>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<PHG4ParticleSvtxMap>;

DumpPHG4ParticleSvtxMap::DumpPHG4ParticleSvtxMap(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpPHG4ParticleSvtxMap::process_Node(PHNode *myNode)
{
  PHG4ParticleSvtxMap *phg4particlesvtxmap = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    phg4particlesvtxmap = thisNode->getData();
  }
  if (phg4particlesvtxmap)
  {
    *fout << "size " << phg4particlesvtxmap->size() << std::endl;
    for (auto &iter : *phg4particlesvtxmap)
    {
      *fout << "Cluster: " << std::hex << iter.first << std::dec << std::endl;

      for (auto &iter2 : iter.second)
      {
        *fout << "weight: " << iter2.first << std::endl;
        for (unsigned int iter3 : iter2.second)
        {
          *fout << "track id " << iter3 << std::endl;
        }
      }
    }
  }
  return 0;
}
