#ifndef _TRACKBASEHISTORIC_TRACKANALYSISUTILS_H
#define _TRACKBASEHISTORIC_TRACKANALYSISUTILS_H

#include <utility>
#include <Acts/Definitions/Algebra.hpp>

class SvtxTrack;

namespace TrackAnalysisUtils
{
  /// Returns DCA as .first and uncertainty on DCA as .second
  using DCA = std::pair<float, float>;
  using DCAPair = std::pair<DCA, DCA>;

  DCAPair get_dca(SvtxTrack* track, Acts::Vector3& vertex);

};

#endif
