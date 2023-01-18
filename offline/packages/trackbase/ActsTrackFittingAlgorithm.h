
#ifndef TRACKBASE_ACTSTRACKFITTINGALGORITHM_H
#define TRACKBASE_ACTSTRACKFITTINGALGORITHM_H

#include "ActsSourceLink.h"
#include "Calibrator.h"
#include "ResidualOutlierFinder.h"

#include <Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wunused-value"
#include <Acts/TrackFitting/KalmanFitter.hpp>
#pragma GCC diagnostic pop

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <Acts/EventData/VectorMultiTrajectory.hpp>
#pragma GCC diagnostic pop
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include <Acts/Propagator/MultiEigenStepperLoop.hpp>
#pragma GCC diagnostic pop

#include <Acts/EventData/TrackParameters.hpp>

#include <functional>
#include <memory>
#include <vector>

namespace Acts
{
  class TrackingGeometry;
}

class ActsTrackFittingAlgorithm final
{
 public:
  using TrackParameters = ::Acts::BoundTrackParameters;
  using Measurement = ::Acts::BoundVariantMeasurement;
  using MeasurementContainer = std::vector<Measurement>;

  /// Track fitter function that takes input measurements, initial trackstate
  /// and fitter options and returns some track-fitter-specific result.
  using TrackFitterOptions = Acts::KalmanFitterOptions<Acts::VectorMultiTrajectory>;

  using TrackFitterResult = Acts::Result<Acts::KalmanFitterResult<Acts::VectorMultiTrajectory>>;

  /// General options that do not depend on the fitter type, but need to be
  /// handed over by the algorithm
  struct GeneralFitterOptions
  {
    std::reference_wrapper<const Acts::GeometryContext> geoContext;
    std::reference_wrapper<const Acts::MagneticFieldContext> magFieldContext;
    std::reference_wrapper<const Acts::CalibrationContext> calibrationContext;
    std::reference_wrapper<const Calibrator> calibrator;
    const Acts::Surface* referenceSurface = nullptr;
    Acts::LoggerWrapper logger;
    Acts::PropagatorPlainOptions propOptions;
  };

  /// Fit function that takes the above parameters and runs a fit
  /// @note This is separated into a virtual interface to keep compilation units
  /// small
  class TrackFitterFunction
  {
   public:
    virtual ~TrackFitterFunction() = default;
    virtual TrackFitterResult operator()(
        const std::vector<std::reference_wrapper<const ActsSourceLink>>&,
        const TrackParameters&, const GeneralFitterOptions&,
        std::shared_ptr<Acts::VectorMultiTrajectory>& trajectory) const = 0;

    virtual void outlierFinder(const ResidualOutlierFinder&) {}
  };

  /// Fit function that takes the above parameters plus a sorted surface
  /// sequence for the DirectNavigator to follow
  /// @note This is separated into a virtual interface to keep compilation units
  /// small
  class DirectedTrackFitterFunction
  {
   public:
    virtual ~DirectedTrackFitterFunction() = default;

    virtual TrackFitterResult operator()(
        const std::vector<std::reference_wrapper<const ActsSourceLink>>&,
        const TrackParameters&, const GeneralFitterOptions&,
        const std::vector<const Acts::Surface*>&,
        std::shared_ptr<Acts::VectorMultiTrajectory>& trajectory) const = 0;
  };

  struct Config
  {
    bool directNavigation;
    /// Type erased fitter function.
    std::shared_ptr<TrackFitterFunction> fit;
    /// Type erased direct navigation fitter function
    std::shared_ptr<DirectedTrackFitterFunction> dFit;
  };

  /// Constructor of the fitting algorithm
  ///
  /// @param config is the config struct to configure the algorihtm
  /// @param level is the logging level
  ActsTrackFittingAlgorithm(Config config, Acts::Logging::Level level);

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

  /// Create the track fitter function implementation.
  ///
  /// The magnetic field is intentionally given by-value since the variant
  /// contains shared_ptr anyways.
  /// @param trackingGeometry
  /// @param multipleScattering correct for MCS (mainly for debugging)
  /// @param energyLoss correct for e-loss
  /// @param reverseFilteringMomThreshold at which threshold
  /// @param freeToBoundCorrection Correction for non-linearity effect during transform from free to bound
  static std::shared_ptr<TrackFitterFunction> makeKalmanFitterFunction(
      std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
      std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
      bool multipleScattering = true, bool energyLoss = true,
      double reverseFilteringMomThreshold = 0.0,
      Acts::FreeToBoundCorrection freeToBoundCorrection =
          Acts::FreeToBoundCorrection());

  static std::shared_ptr<DirectedTrackFitterFunction> makeKalmanFitterFunction(
      std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
      bool multipleScattering = true, bool energyLoss = true,
      double reverseFilteringMomThreshold = 0.0,
      Acts::FreeToBoundCorrection freeToBoundCorrection =
          Acts::FreeToBoundCorrection());

 private:
  /// Helper function to call correct FitterFunction
  TrackFitterResult fitTrack(
      const std::vector<std::reference_wrapper<
          const ActsSourceLink>>& sourceLinks,
      const TrackParameters& initialParameters,
      const GeneralFitterOptions& options,
      const std::vector<const Acts::Surface*>& surfSequence,
      std::shared_ptr<Acts::VectorMultiTrajectory>& trajectory) const;

  Config m_cfg;
};

inline ActsTrackFittingAlgorithm::TrackFitterResult
ActsTrackFittingAlgorithm::fitTrack(
    const std::vector<std::reference_wrapper<
        const ActsSourceLink>>& sourceLinks,
    const TrackParameters& initialParameters,
    const ActsTrackFittingAlgorithm::GeneralFitterOptions& options,
    const std::vector<const Acts::Surface*>& surfSequence,
    std::shared_ptr<Acts::VectorMultiTrajectory>& trajectory) const
{
  if (m_cfg.directNavigation)
  {
    return (*m_cfg.dFit)(sourceLinks, initialParameters, options, surfSequence, trajectory);
  }

  return (*m_cfg.fit)(sourceLinks, initialParameters, options, trajectory);
}

#endif
