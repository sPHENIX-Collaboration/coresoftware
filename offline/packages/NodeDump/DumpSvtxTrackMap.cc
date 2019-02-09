#include "DumpSvtxTrackMap.h"

#include <phool/PHIODataNode.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <string>

using namespace std;

typedef PHIODataNode<SvtxTrackMap> MyNode_t;

DumpSvtxTrackMap::DumpSvtxTrackMap(const string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpSvtxTrackMap::process_Node(PHNode *myNode)
{
  SvtxTrackMap *svtxtrackmap = NULL;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    svtxtrackmap = thisNode->getData();
  }
  if (svtxtrackmap)
  {
    SvtxTrackMap::ConstIter hiter;
    *fout << "size: " << svtxtrackmap->size() << endl;
    for (hiter = svtxtrackmap->begin(); hiter != svtxtrackmap->end(); hiter++)
    {
      *fout << "id: 0x" << hex << hiter->second->get_id() << dec << endl;
      *fout << "positive_charge: " << hiter->second->get_positive_charge() << endl;
      *fout << "charge: " << hiter->second->get_charge() << endl;
      *fout << "chisq: " << hiter->second->get_chisq() << endl;
      *fout << "ndf: " << hiter->second->get_ndf() << endl;
      *fout << "dca: " << hiter->second->get_dca() << endl;
      *fout << "dca2d: " << hiter->second->get_dca2d() << endl;
      *fout << "dca2d_error: " << hiter->second->get_dca2d_error() << endl;
      *fout << "x: " << hiter->second->get_x() << endl;
      *fout << "y: " << hiter->second->get_y() << endl;
      *fout << "z: " << hiter->second->get_z() << endl;
      *fout << "px: " << hiter->second->get_px() << endl;
      *fout << "py: " << hiter->second->get_py() << endl;
      *fout << "pz: " << hiter->second->get_pz() << endl;
    }
  }
  return 0;
}
