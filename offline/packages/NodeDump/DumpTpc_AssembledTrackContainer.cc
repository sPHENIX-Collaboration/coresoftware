#include "DumpTpc_AssembledTrackContainer.h"

#include <phool/PHIODataNode.h>

#include <tpctrackreco/Tpc_AssembledTrack.h>
#include <tpctrackreco/Tpc_AssembledTrackContainer.h>

#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<Tpc_AssembledTrackContainer>;

DumpTpc_AssembledTrackContainer::DumpTpc_AssembledTrackContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpTpc_AssembledTrackContainer::process_Node(PHNode *myNode)
{
  Tpc_AssembledTrackContainer *assembledtrack_container = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);  // NOLINT(cppcoreguidelines-pro-type-static-cast-downcast)
  if (thisNode)
  {
    assembledtrack_container = thisNode->getData();
  }
  if (assembledtrack_container)
  {
    *fout << "size: " << assembledtrack_container->size() << std::endl;
    for (unsigned int i = 0; i < assembledtrack_container->size(); i++)
    {
      Tpc_AssembledTrack *track = assembledtrack_container->get_track(i);
      *fout << "get_event(" << i << "): " << track->get_event() << std::endl;
      *fout << "get_track_id(" << i << "): " << track->get_track_id() << std::endl;
      *fout << "get_side(" << i << "): " << track->get_side() << std::endl;

      // --- Topology ---
      *fout << "get_nsegments(" << i << "): " << track->get_nsegments() << std::endl;
      *fout << "get_nblobs(" << i << "): " << track->get_nblobs() << std::endl;
      *fout << "get_nrawhits(" << i << "): " << track->get_nrawhits() << std::endl;
      *fout << "get_first_layer(" << i << "): " << track->get_first_layer() << std::endl;
      *fout << "get_last_layer(" << i << "): " << track->get_last_layer() << std::endl;
      *fout << "get_first_sector(" << i << "): " << track->get_first_sector() << std::endl;
      *fout << "get_last_sector(" << i << "): " << track->get_last_sector() << std::endl;
      *fout << "get_first_region(" << i << "): " << track->get_first_region() << std::endl;
      *fout << "get_last_region(" << i << "): " << track->get_last_region() << std::endl;

      // --- Full-track fit in global sector-unwrapped coordinates ---
      // global_phi is in radians.  It is built from sector + local pad using the
      // same local pad/phi calibration constants as Tpc_ModuleTrackReco.
      *fout << "get_phi_slope(" << i << "): " << track->get_phi_slope() << std::endl;
      *fout << "get_phi_intercept(" << i << "): " << track->get_phi_intercept() << std::endl;
      *fout << "get_tbin_slope(" << i << "): " << track->get_tbin_slope() << std::endl;
      *fout << "get_tbin_intercept(" << i << "): " << track->get_tbin_intercept() << std::endl;
      *fout << "get_chi2_phi(" << i << "): " << track->get_chi2_phi() << std::endl;
      *fout << "get_chi2_tbin(" << i << "): " << track->get_chi2_tbin() << std::endl;
      *fout << "get_ndof_phi(" << i << "): " << track->get_ndof_phi() << std::endl;
      *fout << "get_ndof_tbin(" << i << "): " << track->get_ndof_tbin() << std::endl;
      *fout << "get_vertex_valid(" << i << "): " << track->get_vertex_valid() << std::endl;
      *fout << "get_vertex_x(" << i << "): " << track->get_vertex_x() << std::endl;
      *fout << "get_vertex_y(" << i << "): " << track->get_vertex_y() << std::endl;
      *fout << "get_vertex_r(" << i << "): " << track->get_vertex_r() << std::endl;
      *fout << "get_vertex_phi(" << i << "): " << track->get_vertex_phi() << std::endl;
      *fout << "get_vertex_tbin(" << i << "): " << track->get_vertex_tbin() << std::endl;
      *fout << "get_vertex_npairs(" << i << "): " << track->get_vertex_npairs() << std::endl;
      *fout << "get_vertex_quality(" << i << "): " << track->get_vertex_quality() << std::endl;
      *fout << "get_seed_valid(" << i << "): " << track->get_seed_valid() << std::endl;
      *fout << "get_seed_x(" << i << "): " << track->get_seed_x() << std::endl;
      *fout << "get_seed_y(" << i << "): " << track->get_seed_y() << std::endl;
      *fout << "get_seed_z(" << i << "): " << track->get_seed_z() << std::endl;
      *fout << "get_seed_px(" << i << "): " << track->get_seed_px() << std::endl;
      *fout << "get_seed_py(" << i << "): " << track->get_seed_py() << std::endl;
      *fout << "get_seed_pz(" << i << "): " << track->get_seed_pz() << std::endl;
      for (unsigned int j = 0; j < 6; j++)
      {
        for (unsigned int k = 0; k < 6; k++)
        {
          *fout << "trk " << i << " get_seed_cov(" << j << "," << k << "): " << track->get_seed_cov(j, k) << std::endl;
        }
      }
      *fout << "size_source_tracks(" << i << "): " << track->size_source_tracks() << std::endl;
      for (unsigned int j = 0; j < track->size_source_tracks(); j++)
      {
        *fout << "get_source_track_id(" << i << "," << j << "): " << track->get_source_track_id(j) << std::endl;
        *fout << "get_source_region(" << i << "," << j << "): " << track->get_source_region(j) << std::endl;
        *fout << "get_source_sector(" << i << "," << j << "): " << track->get_source_sector(j) << std::endl;
      }
      *fout << "size_hit_indices(" << i << "): " << track->size_hit_indices() << std::endl;
      auto hitindices = track->get_hit_indices();
      for (auto iter : hitindices)
      {
        *fout << std::hex << "hitsetkey 0x" << iter.first << ", hitkey 0x" << iter.second << std::dec << std::endl;
      }
    }
  }
  return 0;
}
