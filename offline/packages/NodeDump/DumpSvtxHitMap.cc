#include "DumpSvtxHitMap.h"

#include <phool/PHIODataNode.h>

#include <trackbase_historic/SvtxHitMap.h>
#include <trackbase_historic/SvtxHit.h>


#include <string>

using namespace std;

typedef PHIODataNode<SvtxHitMap> MyNode_t;

DumpSvtxHitMap::DumpSvtxHitMap(const string &NodeName): DumpObject(NodeName)
{
  return ;
}

int DumpSvtxHitMap::process_Node(PHNode *myNode)
{
  SvtxHitMap *svtxhitmap = NULL;
  MyNode_t *thisNode = static_cast <MyNode_t *> (myNode);
  if (thisNode)
    {
      svtxhitmap = thisNode->getData();
    }
  if (svtxhitmap)
    {
      SvtxHitMap::ConstIter hiter;
      *fout << "size: " << svtxhitmap->size() << endl;
      for (hiter = svtxhitmap->begin(); hiter != svtxhitmap->end(); hiter++)
        {
          *fout << "id: 0x" << hex << hiter->second->get_id() << dec << endl;
          *fout << "layer: " << hiter->second->get_layer() << endl;
          *fout << "adc: " << hiter->second->get_adc() << endl;
          *fout << "e: " << hiter->second->get_e() << endl;
          *fout << "cellid: " << hiter->second->get_cellid() << endl;
        }
    }
  return 0;
}

