#ifndef _TRACKBASEHISTORIC_TRACKANALYSISUTILS_H
#define _TRACKBASEHISTORIC_TRACKANALYSISUTILS_H

#include <utility>
#include <Acts/Definitions/Algebra.hpp>
#include <trackbase/TrkrDefs.h>

class SvtxTrack;

namespace TrackAnalysisUtils
{
  /// Returns DCA as .first and uncertainty on DCA as .second
  using DCA = std::pair<float, float>;
  using DCAPair = std::pair<DCA, DCA>;

  DCAPair get_dca(SvtxTrack* track, Acts::Vector3& vertex);

  std::vector<TrkrDefs::cluskey> get_cluster_keys(SvtxTrack* track);

};

#endif
