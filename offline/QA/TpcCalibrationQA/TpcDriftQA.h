#ifndef QA_TPCCALIBRATIONQA_DRIFT_H
#define QA_TPCCALIBRATIONQA_DRIFT_H

#include "TpcQAInternal.h"

namespace TpcCalibrationQADetail
{
  class TpcDriftQA
  {
   public:
    static DriftInput load(const Sample &sample, std::ostream &availability,
                           double defaultRawShift = 40.);
    static void apply_reference(DriftInput &d, const Sample &sample, const char *path);
    static Curves analyze(DriftInput &d, const Sample &sample, TDirectory *dir,
                          std::ostream &csv, ObjectIndex &outputIndex);
  };
}  // namespace TpcCalibrationQADetail
#endif
