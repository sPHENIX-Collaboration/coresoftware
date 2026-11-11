#ifndef QA_TPCCALIBRATIONQA_DISTORTION_H
#define QA_TPCCALIBRATIONQA_DISTORTION_H

#include "TpcQAInternal.h"

namespace TpcCalibrationQADetail
{
  // Histogram postprocessing; no event nodes or conditions services are required.
  class TpcDistortionQA
  {
   public:
    static Fit fit(TH1 *histogram, int minEntries, bool background = false);
    static Curves analyze(Hists &hs, const Sample &sample, const Config &cfg,
                          TDirectory *dir, std::ostream &csv, std::ostream &trends);
  };
}  // namespace TpcCalibrationQADetail
#endif
