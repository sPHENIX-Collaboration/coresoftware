#include "DumpTpc_PolyClusterContainer.h"

#include <phool/PHIODataNode.h>

#include <tpctrackreco/Tpc_PolyCluster.h>
#include <tpctrackreco/Tpc_PolyClusterContainer.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<Tpc_PolyClusterContainer>;

DumpTpc_PolyClusterContainer::DumpTpc_PolyClusterContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
}

int DumpTpc_PolyClusterContainer::process_Node(PHNode *myNode)
{
  Tpc_PolyClusterContainer *container{nullptr};
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
      const Tpc_PolyCluster *cluster = container->get_cluster(i);
      *fout << "get_event(" << i << "): " << cluster->get_event() << std::endl;
      *fout << "get_cluster_id(" << i << "): " << cluster->get_cluster_id() << std::endl;
      *fout << "get_source_assembled_track_id(" << i << "): " << cluster->get_source_assembled_track_id() << std::endl;
      *fout << std::hex << "get_trkr_cluster_key(" << i << "): 0x" << cluster->get_trkr_cluster_key() << std::dec << std::endl;
      *fout << "get_side(" << i << "): " << cluster->get_side() << std::endl;
      *fout << "get_centroid_x(" << i << "): " << cluster->get_centroid_x() << std::endl;
      *fout << "get_centroid_y(" << i << "): " << cluster->get_centroid_y() << std::endl;
      *fout << "get_centroid_z(" << i << "): " << cluster->get_centroid_z() << std::endl;
      *fout << "get_rms_x(" << i << "): " << cluster->get_rms_x() << std::endl;
      *fout << "get_rms_y(" << i << "): " << cluster->get_rms_y() << std::endl;
      *fout << "get_rms_z(" << i << "): " << cluster->get_rms_z() << std::endl;
      *fout << "get_adc(" << i << "): " << cluster->get_adc() << std::endl;
      *fout << "get_phi_width(" << i << "): " << cluster->get_phi_width() << std::endl;
      *fout << "get_time_width(" << i << "): " << cluster->get_time_width() << std::endl;
      *fout << "get_phase(" << i << "): " << cluster->get_phase() << std::endl;
      *fout << "size_hits(" << i << "): " << cluster->size_hits() << std::endl;
      for (unsigned int j = 0; j < cluster->size_hits(); ++j)
      {
        const auto hit = cluster->get_hit_index(j);
        *fout << std::hex << "get_hit_index(" << i << ", " << j << "): hitsetkey 0x" << hit.first
              << ", hitkey 0x" << hit.second << std::dec << std::endl;
        *fout << "get_hit_x(" << i << ", " << j << "): " << cluster->get_hit_x(j) << std::endl;
        *fout << "get_hit_y(" << i << ", " << j << "): " << cluster->get_hit_y(j) << std::endl;
        *fout << "get_hit_z(" << i << ", " << j << "): " << cluster->get_hit_z(j) << std::endl;
      }
    }
  }
  return 0;
}
