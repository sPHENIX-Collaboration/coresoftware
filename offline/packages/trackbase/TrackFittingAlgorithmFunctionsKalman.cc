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

  template <typename TrackFitterFunktion>
  auto makeKfOptions(
      const TrackFitterFunktion& f,
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
        extensions, options.logger, options.propOptions,
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

    Acts::GainMatrixUpdater kfUpdater;
    Acts::GainMatrixSmoother kfSmoother;
    SimpleReverseFilteringLogic reverseFilteringLogic;

    bool multipleScattering;
    bool energyLoss;
    Acts::FreeToBoundCorrection freeToBoundCorrection;

    TrackFitterFunctionImpl(Fitter&& f)
      : trackFitter(std::move(f))
      , multipleScattering(true)
      , energyLoss(true)
      , freeToBoundCorrection(Acts::FreeToBoundCorrection())
    {
    }

    ActsTrackFittingAlgorithm::TrackFitterResult operator()(
        const std::vector<std::reference_wrapper<
            const ActsSourceLink>>& sourceLinks,
        const ActsTrackFittingAlgorithm::TrackParameters& initialParameters,
        const ActsTrackFittingAlgorithm::GeneralFitterOptions& options,
        std::shared_ptr<Acts::VectorMultiTrajectory>& trajectory)
        const override
    {
      auto kfOptions = makeKfOptions(*this, options);
      kfOptions.extensions.calibrator
          .connect<&Calibrator::calibrate>(
              &options.calibrator.get());
      return trackFitter.fit(sourceLinks.begin(), sourceLinks.end(),
                             initialParameters, kfOptions, trajectory);
    };
  };

  struct DirectedFitterFunctionImpl
    : public ActsTrackFittingAlgorithm::DirectedTrackFitterFunction
  {
    DirectFitter fitter;

    Acts::GainMatrixUpdater kfUpdater;
    Acts::GainMatrixSmoother kfSmoother;
    SimpleReverseFilteringLogic reverseFilteringLogic;

    bool multipleScattering;
    bool energyLoss;
    Acts::FreeToBoundCorrection freeToBoundCorrection;

    DirectedFitterFunctionImpl(DirectFitter&& f)
      : fitter(std::move(f))
      , multipleScattering(true)
      , energyLoss(true)
      , freeToBoundCorrection(Acts::FreeToBoundCorrection())
    {
    }

    ActsTrackFittingAlgorithm::TrackFitterResult operator()(
        const std::vector<std::reference_wrapper<
            const ActsSourceLink>>& sourceLinks,
        const ActsTrackFittingAlgorithm::TrackParameters& initialParameters,
        const ActsTrackFittingAlgorithm::GeneralFitterOptions& options,
        const std::vector<const Acts::Surface*>& sSequence,
        std::shared_ptr<Acts::VectorMultiTrajectory>& trajectory) const override
    {
      auto kfOptions = makeKfOptions(*this, options);
      kfOptions.extensions.calibrator
          .connect<&Calibrator::calibrate>(
              &options.calibrator.get());
      return fitter.fit(sourceLinks.begin(), sourceLinks.end(), initialParameters,
                        kfOptions, sSequence, trajectory);
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

  sPHENIXTrackFitterFunctionImpl(Fitter&& f)
    : TrackFitterFunctionImpl(std::move(f))
  {
  }

  ActsTrackFittingAlgorithm::TrackFitterResult operator()(
      const std::vector<std::reference_wrapper<
          const ActsSourceLink>>& sourceLinks,
      const ActsTrackFittingAlgorithm::TrackParameters& initialParameters,
      const ActsTrackFittingAlgorithm::GeneralFitterOptions& options,
      std::shared_ptr<Acts::VectorMultiTrajectory>& trajectory)
      const override
  {
    auto kfOptions = makeKfOptions(*this, options);
    kfOptions.extensions.calibrator.connect<&Calibrator::calibrate>(
        &options.calibrator.get());

    if (use_OF)
    {
      kfOptions.extensions.outlierFinder.connect<&ResidualOutlierFinder::operator()>(&m_oFinder);
    }

    return trackFitter.fit(sourceLinks.begin(), sourceLinks.end(),
                           initialParameters, kfOptions, trajectory);
  };
};

std::shared_ptr<ActsTrackFittingAlgorithm::TrackFitterFunction>
ActsTrackFittingAlgorithm::makeKalmanFitterFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
    bool multipleScattering, bool energyLoss,
    double reverseFilteringMomThreshold,
    Acts::FreeToBoundCorrection freeToBoundCorrection)
{
  Stepper stepper(std::move(magneticField));
  Acts::Navigator::Config cfg{trackingGeometry};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg);
  Propagator propagator(std::move(stepper), std::move(navigator));
  Fitter trackFitter(std::move(propagator));

  // build the fitter functions. owns the fitter object.
  auto fitterFunction =
      std::make_shared<sPHENIXTrackFitterFunctionImpl>(std::move(trackFitter));
  fitterFunction->multipleScattering = multipleScattering;
  fitterFunction->energyLoss = energyLoss;
  fitterFunction->reverseFilteringLogic.momentumThreshold =
      reverseFilteringMomThreshold;
  fitterFunction->freeToBoundCorrection = freeToBoundCorrection;

  return fitterFunction;
}

std::shared_ptr<
    ActsTrackFittingAlgorithm::DirectedTrackFitterFunction>
ActsTrackFittingAlgorithm::makeKalmanFitterFunction(
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
    bool multipleScattering, bool energyLoss,
    double reverseFilteringMomThreshold,
    Acts::FreeToBoundCorrection freeToBoundCorrection)
{
  // construct all components for the fitter
  Stepper stepper(std::move(magneticField));
  Acts::DirectNavigator navigator;
  DirectPropagator propagator(std::move(stepper), navigator);
  DirectFitter fitter(std::move(propagator));

  // build the fitter functions. owns the fitter object.
  auto fitterFunction =
      std::make_shared<DirectedFitterFunctionImpl>(std::move(fitter));
  fitterFunction->multipleScattering = multipleScattering;
  fitterFunction->energyLoss = energyLoss;
  fitterFunction->reverseFilteringLogic.momentumThreshold =
      reverseFilteringMomThreshold;
  fitterFunction->freeToBoundCorrection = freeToBoundCorrection;

  return fitterFunction;
}
