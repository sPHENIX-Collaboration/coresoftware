#include "DumpTpc_PolyTrackVertexContainer.h"

#include <phool/PHIODataNode.h>

#include <tpctrackreco/Tpc_PolyTrackVertex.h>
#include <tpctrackreco/Tpc_PolyTrackVertexContainer.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<Tpc_PolyTrackVertexContainer>;

DumpTpc_PolyTrackVertexContainer::DumpTpc_PolyTrackVertexContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
}

int DumpTpc_PolyTrackVertexContainer::process_Node(PHNode *myNode)
{
  Tpc_PolyTrackVertexContainer *container{nullptr};
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
      const Tpc_PolyTrackVertex *vertex = container->get_vertex(i);
      *fout << "get_track_id(" << i << "): " << vertex->get_track_id() << std::endl;
      *fout << "get_source_assembled_track_id(" << i << "): " << vertex->get_source_assembled_track_id() << std::endl;
      *fout << "get_dca2d(" << i << "): " << vertex->get_dca2d() << std::endl;
      *fout << "get_z0(" << i << "): " << vertex->get_z0() << std::endl;
      *fout << "get_pca_valid(" << i << "): " << vertex->get_pca_valid() << std::endl;
      *fout << "get_pca_x(" << i << "): " << vertex->get_pca_x() << std::endl;
      *fout << "get_pca_y(" << i << "): " << vertex->get_pca_y() << std::endl;
      *fout << "get_pca_z(" << i << "): " << vertex->get_pca_z() << std::endl;
      *fout << "get_pca_radius(" << i << "): " << vertex->get_pca_radius() << std::endl;
      *fout << "get_pca_phi(" << i << "): " << vertex->get_pca_phi() << std::endl;
    }

    *fout << "get_collision_vertex_valid: " << container->get_collision_vertex_valid() << std::endl;
    *fout << "get_collision_min_clusters: " << container->get_collision_min_clusters() << std::endl;
    *fout << "get_collision_vertex_count: " << container->get_collision_vertex_count() << std::endl;
    for (unsigned int i = 0; i < container->get_collision_vertex_count(); ++i)
    {
      *fout << "get_collision_x(" << i << "): " << container->get_collision_x(i) << std::endl;
      *fout << "get_collision_y(" << i << "): " << container->get_collision_y(i) << std::endl;
      *fout << "get_collision_z(" << i << "): " << container->get_collision_z(i) << std::endl;
      *fout << "get_collision_z_rms(" << i << "): " << container->get_collision_z_rms(i) << std::endl;
      *fout << "get_collision_ntracks(" << i << "): " << container->get_collision_ntracks(i) << std::endl;
    }
  }
  return 0;
}
