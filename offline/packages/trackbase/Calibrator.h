
#ifndef TRACKBASE_CALIBRATOR_H
#define TRACKBASE_CALIBRATOR_H

#include <Acts/EventData/Measurement.hpp>
#include <Acts/EventData/MultiTrajectory.hpp>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <Acts/EventData/VectorMultiTrajectory.hpp>
#pragma GCC diagnostic pop

#include <Acts/EventData/SourceLink.hpp>

#include "TrkrDefs.h"
#include "alignmentTransformationContainer.h"
#include "ActsSourceLink.h"

class Calibrator
{
  using Measurement = ::Acts::BoundVariantMeasurement;
  /// Container of measurements.
  ///
  /// In contrast to the source links, the measurements themself must not be
  /// orderable. The source links stored in the measurements are treated
  /// as opaque here and no ordering is enforced on the stored measurements.
  using MeasurementContainer = std::vector<Measurement>;

 public:
  /// Construct an invalid calibrator. Required to allow copying.
  Calibrator() = default;
  /// Construct using a user-provided container to chose measurements from.
  Calibrator(const MeasurementContainer& measurements)
    : m_measurements(&measurements)
  {
  }

  /// Find the measurement corresponding to the source link.
  ///
  /// @tparam parameters_t Track parameters type
  /// @param gctx The geometry context (unused)
  /// @param trackState The track state to calibrate

  void calibrate(const Acts::GeometryContext& gctx,
                 Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::TrackStateProxy trackState) const
  {
    const auto& sourceLink =
        static_cast<const ActsSourceLink&>(trackState.uncalibrated());
    assert(m_measurements and
           "Undefined measurement container in DigitizedCalibrator");
    assert((sourceLink.index() < m_measurements->size()) and
           "Source link index is outside the container bounds");
    std::visit(
        [&](const auto& uncalibmeas) {
          std::array<Acts::BoundIndices, 2> indices;
          indices[0] = Acts::BoundIndices::eBoundLoc0;
          indices[1] = Acts::BoundIndices::eBoundLoc1;

          Acts::ActsVector<2> loc;
          loc(0) = uncalibmeas.parameters()[Acts::eBoundLoc0];
          loc(1) = uncalibmeas.parameters()[Acts::eBoundLoc1];

          auto cov = uncalibmeas.covariance();
          const auto& cluskey = sourceLink.cluskey();
          const auto trkrid = TrkrDefs::getTrkrId(cluskey);
          const double misalignmentFactor = gctx.get<alignmentTransformationContainer*>()->getMisalignmentFactor(trkrid);

          Acts::ActsSymMatrix<2> expandedCov = Acts::ActsSymMatrix<2>::Zero();

          for (int i = 0; i < cov.rows(); i++)
          {
            for (int j = 0; j < cov.cols(); j++)
            {
              expandedCov(i, j) = cov(i, j) * misalignmentFactor;
            }
          }

          Acts::Measurement<Acts::BoundIndices, 2> meas(
              uncalibmeas.sourceLink(), indices,
              loc, expandedCov);

	  trackState.allocateCalibrated(meas.size());
          trackState.setCalibrated(meas);

    },
        (*m_measurements)[sourceLink.index()]);
  }
    
 private:
  // use pointer so the calibrator is copyable and default constructible.
  const MeasurementContainer* m_measurements = nullptr;
};

#endif
