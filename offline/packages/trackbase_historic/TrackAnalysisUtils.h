#ifndef TRACKBASEHISTORIC_TRACKANALYSISUTILS_H
#define TRACKBASEHISTORIC_TRACKANALYSISUTILS_H

#include <trackbase/TrkrDefs.h>
#include <Acts/Definitions/Algebra.hpp>

#include <utility>

class SvtxTrack;
class TrkrClusterContainer;
class PHG4TpcCylinderGeomContainer;

namespace TrackAnalysisUtils
{
  /// Returns DCA as .first and uncertainty on DCA as .second
  using DCA = std::pair<float, float>;
  using DCAPair = std::pair<DCA, DCA>;

  DCAPair get_dca(SvtxTrack* track, Acts::Vector3& vertex);

  std::vector<TrkrDefs::cluskey> get_cluster_keys(SvtxTrack* track);

  float get_dEdx(SvtxTrack* track, TrkrClusterContainer* cluster_map, PHG4TpcCylinderGeomContainer* geom_container);

};  // namespace TrackAnalysisUtils

#endif
