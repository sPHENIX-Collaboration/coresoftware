#include "DumpSvtxTrackMap.h"

#include <phool/PHIODataNode.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<SvtxTrackMap>;

DumpSvtxTrackMap::DumpSvtxTrackMap(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpSvtxTrackMap::process_Node(PHNode *myNode)
{
  SvtxTrackMap *svtxtrackmap = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    svtxtrackmap = thisNode->getData();
  }
  if (svtxtrackmap)
  {
    SvtxTrackMap::ConstIter hiter;
    *fout << "size: " << svtxtrackmap->size() << std::endl;
    for (hiter = svtxtrackmap->begin(); hiter != svtxtrackmap->end(); hiter++)
    {
      *fout << "id: 0x" << std::hex << hiter->second->get_id() << std::dec << std::endl;
      *fout << "positive_charge: " << hiter->second->get_positive_charge() << std::endl;
      *fout << "charge: " << hiter->second->get_charge() << std::endl;
      *fout << "chisq: " << hiter->second->get_chisq() << std::endl;
      *fout << "ndf: " << hiter->second->get_ndf() << std::endl;
      *fout << "dca: " << hiter->second->get_dca() << std::endl;
      *fout << "dca2d: " << hiter->second->get_dca2d() << std::endl;
      *fout << "dca2d_error: " << hiter->second->get_dca2d_error() << std::endl;
      *fout << "x: " << hiter->second->get_x() << std::endl;
      *fout << "y: " << hiter->second->get_y() << std::endl;
      *fout << "z: " << hiter->second->get_z() << std::endl;
      *fout << "px: " << hiter->second->get_px() << std::endl;
      *fout << "py: " << hiter->second->get_py() << std::endl;
      *fout << "pz: " << hiter->second->get_pz() << std::endl;
      for (SvtxTrack::ConstStateIter trkstates = hiter->second->begin_states();
           trkstates != hiter->second->end_states();
           ++trkstates)
      {
        *fout << "trackstate path length: " << trkstates->first << std::endl;
        *fout << "trackstate name: " << trkstates->second->get_name() << std::endl;
        *fout << "trackstate x: " << trkstates->second->get_x() << std::endl;
        *fout << "trackstate y: " << trkstates->second->get_y() << std::endl;
        *fout << "trackstate z: " << trkstates->second->get_z() << std::endl;
        *fout << "trackstate px: " << trkstates->second->get_px() << std::endl;
        *fout << "trackstate py: " << trkstates->second->get_py() << std::endl;
        *fout << "trackstate pz: " << trkstates->second->get_pz() << std::endl;
        for (int i = 0; i < 6; i++)
        {
          for (int j = 0; j < 6; j++)
          {
            *fout << "trkstate covar[" << i << ", " << j << "]: "
                  << trkstates->second->get_error(i, j) << std::endl;
          }
        }
      }
    }
  }
  return 0;
}
