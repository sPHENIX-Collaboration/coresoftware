#ifndef TRACKBASEHISTORIC_TRACKANALYSISUTILS_H
#define TRACKBASEHISTORIC_TRACKANALYSISUTILS_H

#include <trackbase/TrkrDefs.h>
#include <Acts/Definitions/Algebra.hpp>

#include <utility>

class SvtxTrack;
class TrackSeed;
class ActsGeometry;
class TrkrClusterContainer;
namespace TrackAnalysisUtils
{
  /// Returns DCA as .first and uncertainty on DCA as .second
  using DCA = std::pair<float, float>;
  using DCAPair = std::pair<DCA, DCA>;

  DCAPair get_dca(SvtxTrack* track, Acts::Vector3& vertex);

  std::vector<TrkrDefs::cluskey> get_cluster_keys(SvtxTrack* track);
  float calc_dedx_calib(TrackSeed* tpcseed, TrkrClusterContainer* cluster_map,
                        ActsGeometry* tgeometry,
                        // this is the thickness from R1[0],R1[1],R2[2], and R4[3]. You need
                        // to pass these from the geometry object, which keeps the dependencies
                        // of this helper class minimal. This will also help us catch any changes
                        // when/if the tpc geometry changes in the future. This is to get us going
                        float thickness_per_region[4]);
  float calc_dedx(TrackSeed* tpcseed, TrkrClusterContainer* clustermap, ActsGeometry* tgeometry,
                  float thickness_per_region[4]);

};  // namespace TrackAnalysisUtils

#endif
