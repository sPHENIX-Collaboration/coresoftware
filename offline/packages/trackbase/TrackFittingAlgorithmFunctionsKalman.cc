#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#include <Acts/Definitions/TrackParametrization.hpp>
#pragma GCC diagnostic pop
#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wunused-value"
#include <Acts/Propagator/EigenStepper.hpp>
#pragma GCC diagnostic pop

#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>
#include <Acts/Utilities/Helpers.hpp>

#include "ActsTrackFittingAlgorithm.h"

namespace
{
  using Updater = Acts::GainMatrixUpdater;
  using Smoother = Acts::GainMatrixSmoother;
  using Stepper = Acts::EigenStepper<>;
  using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
  using Fitter = Acts::KalmanFitter<Propagator, Acts::VectorMultiTrajectory>;
  using DirectPropagator = Acts::Propagator<Stepper, Acts::DirectNavigator>;
  using DirectFitter = Acts::KalmanFitter<DirectPropagator, Acts::VectorMultiTrajectory>;

  struct SimpleReverseFilteringLogic
  {
    SimpleReverseFilteringLogic() = default;
    double momentumThreshold = 0.0;

    bool doBackwardFiltering(
        Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::ConstTrackStateProxy trackState) const
    {
      auto momentum = fabs(1 / trackState.filtered()[Acts::eBoundQOverP]);
      return (momentum <= momentumThreshold);
    }
  };

  template <typename TrackFitterFunction>
  auto makeKfOptions(const TrackFitterFunction& f,
                     ActsTrackFittingAlgorithm::GeneralFitterOptions options)
  {
    Acts::KalmanFitterExtensions<Acts::VectorMultiTrajectory> extensions;
    // cppcheck-suppress constStatement
    extensions.updater.connect<&Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(&f.kfUpdater);
    // cppcheck-suppress constStatement
    extensions.smoother.connect<&Acts::GainMatrixSmoother::operator()<Acts::VectorMultiTrajectory>>(&f.kfSmoother);
    extensions.reverseFilteringLogic
        .connect<&SimpleReverseFilteringLogic::doBackwardFiltering>(
            &f.reverseFilteringLogic);

    Acts::KalmanFitterOptions<Acts::VectorMultiTrajectory> kfOptions(
        options.geoContext, options.magFieldContext, options.calibrationContext,
        extensions, options.propOptions,
        &(*options.referenceSurface));

    kfOptions.multipleScattering = f.multipleScattering;
    kfOptions.energyLoss = f.energyLoss;
    kfOptions.freeToBoundCorrection = f.freeToBoundCorrection;

    return kfOptions;
  }

  struct TrackFitterFunctionImpl
    : public ActsTrackFittingAlgorithm::TrackFitterFunction
  {
    Fitter trackFitter;
    ActsSourceLink::SurfaceAccessor m_slSurfaceAccessor;

    Acts::GainMatrixUpdater kfUpdater;
    Acts::GainMatrixSmoother kfSmoother;
    SimpleReverseFilteringLogic reverseFilteringLogic;

    bool multipleScattering;
    bool energyLoss;
    Acts::FreeToBoundCorrection freeToBoundCorrection;

    TrackFitterFunctionImpl(Fitter&& f,
                            const Acts::TrackingGeometry& trkGeo)
      : trackFitter(std::move(f))
      , m_slSurfaceAccessor{trkGeo}
      , multipleScattering(true)
      , energyLoss(true)
      , freeToBoundCorrection(Acts::FreeToBoundCorrection())
    {
    }

    ActsTrackFittingAlgorithm::TrackFitterResult operator()(
        const std::vector<Acts::SourceLink>& sourceLinks,
        const ActsTrackFittingAlgorithm::TrackParameters& initialParameters,
        const ActsTrackFittingAlgorithm::GeneralFitterOptions& options,
        const CalibratorAdapter& calibrator,
        ActsTrackFittingAlgorithm::TrackContainer& tracks) const override
    {
      auto kfOptions = makeKfOptions(*this, options);
      kfOptions.extensions.surfaceAccessor
          .connect<&ActsSourceLink::SurfaceAccessor::operator()>(
              &m_slSurfaceAccessor);
      kfOptions.extensions.calibrator
          .connect<&CalibratorAdapter::calibrate>(
              &calibrator);
      return trackFitter.fit(sourceLinks.begin(), sourceLinks.end(),
                             initialParameters, kfOptions, tracks);
    }
  };

  struct DirectedFitterFunctionImpl
    : public ActsTrackFittingAlgorithm::DirectedTrackFitterFunction
  {
    DirectFitter fitter;
    ActsSourceLink::SurfaceAccessor m_slSurfaceAccessor;

    Acts::GainMatrixUpdater kfUpdater;
    Acts::GainMatrixSmoother kfSmoother;
    SimpleReverseFilteringLogic reverseFilteringLogic;

    bool multipleScattering;
    bool energyLoss;
    Acts::FreeToBoundCorrection freeToBoundCorrection;

    DirectedFitterFunctionImpl(DirectFitter&& f,
                               const Acts::TrackingGeometry& trkGeo)
      : fitter(std::move(f))
      , m_slSurfaceAccessor{trkGeo}
      , multipleScattering(true)
      , energyLoss(true)
      , freeToBoundCorrection(Acts::FreeToBoundCorrection())
    {
    }

    ActsTrackFittingAlgorithm::TrackFitterResult operator()(
        const std::vector<Acts::SourceLink>& sourceLinks,
        const ActsTrackFittingAlgorithm::TrackParameters& initialParameters,
        const ActsTrackFittingAlgorithm::GeneralFitterOptions& options,
        const std::vector<const Acts::Surface*>& sSequence,
        const CalibratorAdapter& calibrator,
        ActsTrackFittingAlgorithm::TrackContainer& tracks) const override
    {
      auto kfOptions = makeKfOptions(*this, options);
      kfOptions.extensions.calibrator
          .connect<&CalibratorAdapter::calibrate>(
              &calibrator);
      kfOptions.extensions.surfaceAccessor
          .connect<&ActsSourceLink::SurfaceAccessor::operator()>(
              &m_slSurfaceAccessor);
      return fitter.fit(sourceLinks.begin(), sourceLinks.end(), initialParameters,
                        kfOptions, sSequence, tracks);
    };
  };

}  // namespace

struct sPHENIXTrackFitterFunctionImpl : public TrackFitterFunctionImpl
{
  ResidualOutlierFinder m_oFinder;

  bool use_OF = false;

  void outlierFinder(const ResidualOutlierFinder& finder) override
  {
    m_oFinder = finder;
    use_OF = true;
  }

  sPHENIXTrackFitterFunctionImpl(Fitter&& f,
                                 const Acts::TrackingGeometry& trkGeo)
    : TrackFitterFunctionImpl(std::move(f), trkGeo)
  {
  }

  ActsTrackFittingAlgorithm::TrackFitterResult operator()(
      const std::vector<Acts::SourceLink>& sourceLinks,
      const ActsTrackFittingAlgorithm::TrackParameters& initialParameters,
      const ActsTrackFittingAlgorithm::GeneralFitterOptions& options,
      const CalibratorAdapter& calibrator,
      ActsTrackFittingAlgorithm::TrackContainer& tracks)
      const override
  {
    auto kfOptions = makeKfOptions(*this, options);
    kfOptions.extensions.calibrator.connect<&CalibratorAdapter::calibrate>(
        &calibrator);
    kfOptions.extensions.surfaceAccessor
        .connect<&ActsSourceLink::SurfaceAccessor::operator()>(
            &m_slSurfaceAccessor);
    if (use_OF)
    {
      kfOptions.extensions.outlierFinder.connect<&ResidualOutlierFinder::operator()>(&m_oFinder);
    }

    return trackFitter.fit(sourceLinks.begin(), sourceLinks.end(),
                           initialParameters, kfOptions, tracks);
  };
};

std::shared_ptr<ActsTrackFittingAlgorithm::TrackFitterFunction>
ActsTrackFittingAlgorithm::makeKalmanFitterFunction(
    const std::shared_ptr<const Acts::TrackingGeometry>& trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
    bool multipleScattering, bool energyLoss,
    double reverseFilteringMomThreshold,
    Acts::FreeToBoundCorrection freeToBoundCorrection,
    const Acts::Logger& logger)
{
  Stepper stepper(std::move(magneticField));
  const auto& geo = *trackingGeometry;

  Acts::Navigator::Config cfg{trackingGeometry};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg, logger.cloneWithSuffix("Navigator"));
  Propagator propagator(std::move(stepper), std::move(navigator),
                        logger.cloneWithSuffix("Propagator"));
  Fitter trackFitter(std::move(propagator), logger.cloneWithSuffix("Fitter"));

  // build the fitter functions. owns the fitter object.
  auto fitterFunction =
      std::make_shared<sPHENIXTrackFitterFunctionImpl>(std::move(trackFitter), geo);
  fitterFunction->multipleScattering = multipleScattering;
  fitterFunction->energyLoss = energyLoss;
  fitterFunction->reverseFilteringLogic.momentumThreshold =
      reverseFilteringMomThreshold;
  fitterFunction->freeToBoundCorrection = freeToBoundCorrection;

  return fitterFunction;
}

std::shared_ptr<
    ActsTrackFittingAlgorithm::DirectedTrackFitterFunction>
ActsTrackFittingAlgorithm::makeDirectedKalmanFitterFunction(
    const std::shared_ptr<const Acts::TrackingGeometry>& trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
    bool multipleScattering, bool energyLoss,
    double reverseFilteringMomThreshold,
    Acts::FreeToBoundCorrection freeToBoundCorrection,
    const Acts::Logger& logger)
{
  // construct all components for the fitter
  const Stepper stepper(std::move(magneticField));
  const auto& geo = *trackingGeometry;

  Acts::DirectNavigator navigator{
      logger.cloneWithSuffix("DirectNavigator")};
  DirectPropagator propagator(stepper, std::move(navigator),
                              logger.cloneWithSuffix("DirectPropagator"));
  DirectFitter fitter(std::move(propagator),
                      logger.cloneWithSuffix("DirectFitter"));

  // build the fitter functions. owns the fitter object.
  auto fitterFunction =
      std::make_shared<DirectedFitterFunctionImpl>(std::move(fitter), geo);
  fitterFunction->multipleScattering = multipleScattering;
  fitterFunction->energyLoss = energyLoss;
  fitterFunction->reverseFilteringLogic.momentumThreshold =
      reverseFilteringMomThreshold;
  fitterFunction->freeToBoundCorrection = freeToBoundCorrection;

  return fitterFunction;
}
