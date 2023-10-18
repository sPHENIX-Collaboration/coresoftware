#include "Calibrator.h"

void Calibrator::calibrate(const Calibrator::MeasurementContainer& measurements,
                           const Acts::GeometryContext& gctx,
                           const Acts::CalibrationContext&,
                           const Acts::SourceLink& sourceLink,
                           Acts::VectorMultiTrajectory::TrackStateProxy& trackState) const
{
  trackState.setUncalibratedSourceLink(sourceLink);
  const ActsSourceLink sl = sourceLink.get<ActsSourceLink>();
  const ActsSourceLink::Index index = sl.index();
  std::visit(
      [&](const auto& uncalibmeas)
      {
        std::array<Acts::BoundIndices, 2> indices;
        indices[0] = Acts::BoundIndices::eBoundLoc0;
        indices[1] = Acts::BoundIndices::eBoundLoc1;

        Acts::ActsVector<2> loc;
        loc(0) = uncalibmeas.parameters()[Acts::eBoundLoc0];
        loc(1) = uncalibmeas.parameters()[Acts::eBoundLoc1];

        auto cov = uncalibmeas.covariance();
        const TrkrDefs::cluskey cluskey = sl.cluskey();
        const uint8_t layer = TrkrDefs::getLayer(cluskey);
        const double misalignmentFactor = gctx.get<alignmentTransformationContainer*>()->getMisalignmentFactor(layer);

        Acts::ActsSquareMatrix<2> expandedCov = Acts::ActsSquareMatrix<2>::Zero();

        for (int i = 0; i < cov.rows(); i++)
        {
          for (int j = 0; j < cov.cols(); j++)
          {
            expandedCov(i, j) = cov(i, j) * misalignmentFactor;
          }
        }

        Acts::Measurement<Acts::BoundIndices, 2> meas(sourceLink,
                                                      indices,
                                                      loc, expandedCov);

        trackState.allocateCalibrated(meas.size());
        trackState.setCalibrated(meas);
      },
      (measurements)[index]);
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
