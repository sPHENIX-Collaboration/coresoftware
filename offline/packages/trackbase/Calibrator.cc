#include "Calibrator.h"

void Calibrator::calibrate(const Calibrator::MeasurementContainer& measurements,
                           const Acts::GeometryContext& gctx,
                           const Acts::CalibrationContext& /*unused*/,
                           const Acts::SourceLink& sourceLink,
                           Acts::VectorMultiTrajectory::TrackStateProxy& trackState) const
{
  trackState.setUncalibratedSourceLink(Acts::SourceLink{sourceLink});
  const ActsSourceLink sl = sourceLink.get<ActsSourceLink>();
  const ActsExamples::ConstVariableBoundMeasurementProxy measurement =
      measurements.getMeasurement(sl.index());

  Acts::visit_measurement(measurement.size(), [&](auto N) -> void
                          {
    constexpr std::size_t kMeasurementSize = decltype(N)::value;
    const ActsExamples::ConstFixedBoundMeasurementProxy<kMeasurementSize> fixedMeasurement =
        static_cast<ActsExamples::ConstFixedBoundMeasurementProxy<kMeasurementSize>>(
            measurement);
    const auto cov = fixedMeasurement.covariance();
    const TrkrDefs::cluskey cluskey = sl.cluskey();
    const uint8_t layer = TrkrDefs::getLayer(cluskey);
    const double misalignmentFactor = gctx.get<alignmentTransformationContainer*>()->getMisalignmentFactor(layer);

    Acts::ActsSquareMatrix<kMeasurementSize> expandedCov = Acts::ActsSquareMatrix<kMeasurementSize>::Zero();

    for (int i = 0; i < cov.rows(); i++)
    {
      for (int j = 0; j < cov.cols(); j++)
      {
        expandedCov(i, j) = cov(i, j) * misalignmentFactor;
      }
    }
    trackState.allocateCalibrated(fixedMeasurement.parameters().eval(),
                                  expandedCov.eval());
    trackState.setProjectorSubspaceIndices(fixedMeasurement.subspaceIndices()); });
}

void CalibratorAdapter::calibrate(
    const Acts::GeometryContext& gctx,
    const Acts::CalibrationContext& cctx,
    const Acts::SourceLink& sourceLink,
    Acts::VectorMultiTrajectory::TrackStateProxy trackState) const
{
  return m_calibrator.calibrate(m_measurements,
                                gctx,
                                cctx,
                                sourceLink,
                                trackState);
}
