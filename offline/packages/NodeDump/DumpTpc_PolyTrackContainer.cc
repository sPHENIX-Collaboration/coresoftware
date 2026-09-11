#include "DumpTpc_PolyTrackContainer.h"

#include <phool/PHIODataNode.h>

#include <tpctrackreco/Tpc_PolyTrack.h>
#include <tpctrackreco/Tpc_PolyTrackContainer.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<Tpc_PolyTrackContainer>;

DumpTpc_PolyTrackContainer::DumpTpc_PolyTrackContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
}

int DumpTpc_PolyTrackContainer::process_Node(PHNode *myNode)
{
  Tpc_PolyTrackContainer *container{nullptr};
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);  // NOLINT(cppcoreguidelines-pro-type-static-cast-downcast)
  if (thisNode)
  {
    container = thisNode->getData();
  }
  if (container)
  {
    unsigned int is = container->size();
    *fout << "size: " << is << std::endl;
    for (unsigned int i = 0; i < is; ++i)
    {
      const Tpc_PolyTrack *track = container->get_track(i);
      *fout << "get_event(" << i << "): " << track->get_event() << std::endl;
      *fout << "get_track_id(" << i << "): " << track->get_track_id() << std::endl;
      *fout << "get_source_assembled_track_id(" << i << "): " << track->get_source_assembled_track_id() << std::endl;
      *fout << "get_fit_status(" << i << "): " << track->get_fit_status() << std::endl;
      *fout << "get_nclusters(" << i << "): " << track->get_nclusters() << std::endl;
      *fout << "get_x(" << i << "): " << track->get_x() << std::endl;
      *fout << "get_y(" << i << "): " << track->get_y() << std::endl;
      *fout << "get_z(" << i << "): " << track->get_z() << std::endl;
      *fout << "get_px(" << i << "): " << track->get_px() << std::endl;
      *fout << "get_py(" << i << "): " << track->get_py() << std::endl;
      *fout << "get_pz(" << i << "): " << track->get_pz() << std::endl;
      *fout << "get_charge(" << i << "): " << track->get_charge() << std::endl;
      *fout << "get_chi2(" << i << "): " << track->get_chi2() << std::endl;
      *fout << "get_ndf(" << i << "): " << track->get_ndf() << std::endl;
      *fout << "get_dedx(" << i << "): " << track->get_dedx() << std::endl;
      *fout << "get_seed_x0(" << i << "): " << track->get_seed_x0() << std::endl;
      *fout << "get_seed_y0(" << i << "): " << track->get_seed_y0() << std::endl;
      *fout << "get_helix_x0(" << i << "): " << track->get_helix_x0() << std::endl;
      *fout << "get_helix_y0(" << i << "): " << track->get_helix_y0() << std::endl;
      *fout << "get_seed_z0(" << i << "): " << track->get_seed_z0() << std::endl;
      *fout << "get_seed_phi(" << i << "): " << track->get_seed_phi() << std::endl;
      *fout << "get_seed_slope(" << i << "): " << track->get_seed_slope() << std::endl;
      *fout << "get_seed_q_over_r(" << i << "): " << track->get_seed_q_over_r() << std::endl;
      for (unsigned int row = 0; row < 6; ++row)
      {
        for (unsigned int column = 0; column < 6; ++column)
        {
          *fout << "get_cov(" << i << ", " << row << ", " << column << "): " << track->get_cov(row, column) << std::endl;
        }
      }
      *fout << "size_cluster_keys(" << i << "): " << track->size_cluster_keys() << std::endl;
      for (unsigned int j = 0; j < track->size_cluster_keys(); ++j)
      {
        *fout << std::hex << "get_cluster_key(" << i << ", " << j << "): 0x" << track->get_cluster_key(j) << std::dec << std::endl;
      }
    }
  }
  return 0;
}
