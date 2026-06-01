
#ifndef TRACKBASE_CALIBRATOR_H
#define TRACKBASE_CALIBRATOR_H

#include "ActsSourceLink.h"
#include "TrkrDefs.h"
#include "alignmentTransformationContainer.h"

#include <ActsExamples/EventData/Measurement.hpp>
#include <Acts/EventData/MultiTrajectory.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>


class Calibrator
{
 public:
  using MeasurementContainer = ActsExamples::MeasurementContainer;

  void calibrate(const MeasurementContainer& measurements,
                 const Acts::GeometryContext& gctx,
                 const Acts::CalibrationContext&,
                 const Acts::SourceLink& sourceLink,
                 Acts::VectorMultiTrajectory::TrackStateProxy& trackState) const;
  virtual ~Calibrator() = default;
};

class CalibratorAdapter
{
 public:
  /// Construct using a user-provided container to chose measurements from.
  CalibratorAdapter(const Calibrator& calibrator,
                    const Calibrator::MeasurementContainer& measurements)
    : m_calibrator{calibrator}
    , m_measurements{measurements}
  {
  }

  CalibratorAdapter() = delete;

  void calibrate(const Acts::GeometryContext& gctx,
                 const Acts::CalibrationContext& cctx,
                 const Acts::SourceLink& sourceLink,
                 Acts::VectorMultiTrajectory::TrackStateProxy trackState) const;

 private:
  const Calibrator& m_calibrator;
  const Calibrator::MeasurementContainer& m_measurements;
};

#endif
