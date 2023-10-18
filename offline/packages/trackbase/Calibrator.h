
#ifndef TRACKBASE_CALIBRATOR_H
#define TRACKBASE_CALIBRATOR_H

#include <Acts/EventData/Measurement.hpp>
#include <Acts/EventData/MultiTrajectory.hpp>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <Acts/EventData/VectorMultiTrajectory.hpp>
#pragma GCC diagnostic pop

#include "ActsSourceLink.h"
#include "TrkrDefs.h"
#include "alignmentTransformationContainer.h"

class Calibrator
{
 public:
  using Measurement = ::Acts::BoundVariantMeasurement;
  /// Container of measurements.
  ///
  /// In contrast to the source links, the measurements themself must not be
  /// orderable. The source links stored in the measurements are treated
  /// as opaque here and no ordering is enforced on the stored measurements.
  using MeasurementContainer = std::vector<Measurement>;

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
