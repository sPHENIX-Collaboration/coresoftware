#include "DumpSvtxClusterMap.h"

#include <phool/PHIODataNode.h>

#include <trackbase_historic/SvtxClusterMap.h>
#include <trackbase_historic/SvtxCluster.h>


#include <string>

using namespace std;

typedef PHIODataNode<SvtxClusterMap> MyNode_t;

DumpSvtxClusterMap::DumpSvtxClusterMap(const string &NodeName): DumpObject(NodeName)
{
  return ;
}

int DumpSvtxClusterMap::process_Node(PHNode *myNode)
{
  SvtxClusterMap *svtxclustermap = NULL;
  MyNode_t *thisNode = static_cast <MyNode_t *> (myNode);
  if (thisNode)
    {
      svtxclustermap = thisNode->getData();
    }
  if (svtxclustermap)
    {
      SvtxClusterMap::ConstIter hiter;
      *fout << "size: " << svtxclustermap->size() << endl;
      for (hiter = svtxclustermap->begin(); hiter != svtxclustermap->end(); hiter++)
        {
          *fout << "id: 0x" << hex << hiter->second->get_id() << dec << endl;
          *fout << "layer: " << hiter->second->get_layer() << endl;
          *fout << "x: " << hiter->second->get_x() << endl;
          *fout << "y: " << hiter->second->get_y() << endl;
          *fout << "z: " << hiter->second->get_z() << endl;
          *fout << "e: " << hiter->second->get_e() << endl;
          *fout << "adc: " << hiter->second->get_adc() << endl;
        }
    }
  return 0;
}

