#ifndef _TRACKBASEHISTORIC_TRACKANALYSISUTILS_H
#define _TRACKBASEHISTORIC_TRACKANALYSISUTILS_H

#include <utility>

class SvtxTrack;
class SvtxVertex;

class TrackAnalysisUtils
{
 public:
  using DCA = std::pair<float, float>;
  using DCAPair = std::pair<DCA, DCA>;

  TrackAnalysisUtils() = default;
  ~TrackAnalysisUtils() {}

  DCAPair get_dca(SvtxTrack* track, SvtxVertex* svtxVertex);

 private:
};

#endif
