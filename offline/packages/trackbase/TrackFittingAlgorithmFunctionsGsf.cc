
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/SharedBField.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wunused-value"
#include <Acts/Propagator/MultiEigenStepperLoop.hpp>
#pragma GCC diagnostic pop

#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>
#include <Acts/TrackFitting/GaussianSumFitter.hpp>
#include <Acts/Utilities/Helpers.hpp>
#include <ActsExamples/MagneticField/MagneticField.hpp>

#include "ActsTrackFittingAlgorithm.h"

using namespace ActsExamples;

namespace {

template <typename FitterFunction>
auto makeGsfOptions(
    const FitterFunction& f,
    const ActsTrackFittingAlgorithm::GeneralFitterOptions& options) {
  Acts::GsfExtensions extensions;
  extensions.updater.connect<&Acts::GainMatrixUpdater::operator()>(&f.updater);

  Acts::GsfOptions gsfOptions{options.geoContext,
                              options.magFieldContext,
                              options.calibrationContext,
                              extensions,
                              options.logger,
                              options.propOptions,
                              &(*options.referenceSurface),
                              f.maxComponents,
                              f.abortOnError,
                              f.disableAllMaterialHandling};

  return gsfOptions;
}

template <typename Fitter>
struct GsfFitterFunctionImpl
    : public ActsTrackFittingAlgorithm::TrackFitterFunction {
  Fitter trackFitter;
  Acts::GainMatrixUpdater updater;

  /// These default to the values set in the header file defining
  /// makeGsfFitterFunction
  std::size_t maxComponents = 4;
  bool abortOnError = true;
  bool disableAllMaterialHandling = false;

  GsfFitterFunctionImpl(Fitter&& f) : trackFitter(std::move(f)) {}

  ActsTrackFittingAlgorithm::TrackFitterResult operator()(
      const std::vector<std::reference_wrapper<
          const ActsSourceLink>>& sourceLinks,
      const ActsExamples::TrackParameters& initialParameters,
      const ActsTrackFittingAlgorithm::GeneralFitterOptions& options)
      const override {
    auto gsfOptions = makeGsfOptions(*this, options);
    gsfOptions.extensions.calibrator
        .template connect<&Calibrator::calibrate>(
            &options.calibrator.get());

    return trackFitter.fit(sourceLinks.begin(), sourceLinks.end(),
                           initialParameters, gsfOptions);
  }
};

template <typename Fitter>
struct DirectedFitterFunctionImpl
    : public ActsTrackFittingAlgorithm::DirectedTrackFitterFunction {
  Fitter trackFitter;
  Acts::GainMatrixUpdater updater;

  std::size_t maxComponents = 4;
  bool abortOnError = true;
  bool disableAllMaterialHandling = false;

  DirectedFitterFunctionImpl(Fitter&& f) : trackFitter(std::move(f)) {}

  ActsTrackFittingAlgorithm::TrackFitterResult operator()(
      const std::vector<std::reference_wrapper<
          const ActsSourceLink>>& sourceLinks,
      const ActsExamples::TrackParameters& initialParameters,
      const ActsTrackFittingAlgorithm::GeneralFitterOptions& options,
      const std::vector<const Acts::Surface*>& sSequence) const override {
    auto gsfOptions = makeGsfOptions(*this, options);
    gsfOptions.extensions.calibrator
        .template connect<&Calibrator::calibrate>(
            &options.calibrator.get());

    return trackFitter.fit(sourceLinks.begin(), sourceLinks.end(),
                           initialParameters, gsfOptions, sSequence);
  }
};
}  // namespace


std::shared_ptr<ActsTrackFittingAlgorithm::TrackFitterFunction>
ActsTrackFittingAlgorithm::makeGsfFitterFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
    std::size_t maxComponents, bool abortOnError,
    bool disableAllMaterialHandling) {
  Acts::MultiEigenStepperLoop stepper(std::move(magneticField));
  Acts::Navigator::Config cfg{trackingGeometry};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg);
  Acts::Propagator propagator(std::move(stepper), std::move(navigator));
  Acts::GaussianSumFitter<decltype(propagator)> trackFitter(
      std::move(propagator));

  // build the fitter functions. owns the fitter object.
  auto fitterFunction =
      std::make_shared<GsfFitterFunctionImpl<decltype(trackFitter)>>(
          std::move(trackFitter));
  fitterFunction->maxComponents = maxComponents;
  fitterFunction->abortOnError = abortOnError;
  fitterFunction->disableAllMaterialHandling = disableAllMaterialHandling;

  return fitterFunction;
}

std::shared_ptr<ActsTrackFittingAlgorithm::DirectedTrackFitterFunction>
ActsTrackFittingAlgorithm::makeGsfFitterFunction(
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
    std::size_t maxComponents, bool abortOnError,
    bool disableAllMaterialHandling) {
  Acts::MultiEigenStepperLoop stepper(std::move(magneticField));
  Acts::DirectNavigator navigator;
  Acts::Propagator propagator(std::move(stepper), navigator);
  Acts::GaussianSumFitter<decltype(propagator)> trackFitter(
      std::move(propagator));

  // build the fitter functions. owns the fitter object.
  auto fitterFunction =
      std::make_shared<DirectedFitterFunctionImpl<decltype(trackFitter)>>(
          std::move(trackFitter));
  fitterFunction->maxComponents = maxComponents;
  fitterFunction->abortOnError = abortOnError;
  fitterFunction->disableAllMaterialHandling = disableAllMaterialHandling;

  return fitterFunction;
}
