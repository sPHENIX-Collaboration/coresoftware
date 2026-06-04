
#ifndef _TGEO_DETECTOR_WITH_OPTIONS_H_
#define _TGEO_DETECTOR_WITH_OPTIONS_H_

#include "IBaseDetector.h"
#include <ActsExamples/TGeoDetector/TGeoDetector.hpp>

namespace ActsExamples {

class TGeoDetectorWithOptions : public IBaseDetector {
 public:
  TGeoDetectorWithOptions(TGeoDetector::Config config) : m_detector(config) {}
  TGeoDetector m_detector;

  void addOptions(
      boost::program_options::options_description& opt) const override;
       
};
}  // namespace ActsExamples

#endif  // _TGEO_DETECTOR_WITH_OPTIONS_H_