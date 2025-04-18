
#ifndef _TGEO_DETECTOR_WITH_OPTIONS_H_
#define _TGEO_DETECTOR_WITH_OPTIONS_H_

#include "IBaseDetector.h"
#include <ActsExamples/TGeoDetector/TGeoDetector.hpp>
#include <ActsExamples/Utilities/OptionsFwd.hpp>

namespace ActsExamples {

class TGeoDetectorWithOptions : public IBaseDetector {
 public:
  TGeoDetector m_detector;

  void addOptions(
      boost::program_options::options_description& opt) const override;

  auto finalize(const boost::program_options::variables_map& vm,
                std::shared_ptr<const Acts::IMaterialDecorator> mdecorator)
      -> std::pair<TrackingGeometryPtr, ContextDecorators> override;
};
}  // namespace ActsExamples

#endif  // _TGEO_DETECTOR_WITH_OPTIONS_H_