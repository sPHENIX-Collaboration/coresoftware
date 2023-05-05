#include "ActsGsfTrackFittingAlgorithm.h"

std::shared_ptr<ActsTrackFittingAlgorithm::TrackFitterFunction>
ActsGsfTrackFittingAlgorithm::makeGsfFitterFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
    BetheHeitlerApprox betheHeitlerApprox, std::size_t maxComponents,
    Acts::FinalReductionMethod finalReductionMethod, bool abortOnError,
    bool disableAllMaterialHandling) {
  MultiStepper stepper(std::move(magneticField), finalReductionMethod);

  // Standard fitter
  Acts::Navigator::Config cfg{trackingGeometry};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg);
  Propagator propagator(std::move(stepper), std::move(navigator));
  Fitter trackFitter(std::move(propagator),
                     BetheHeitlerApprox(betheHeitlerApprox));

  // build the fitter functions. owns the fitter object.
  auto fitterFunction = std::make_shared<GsfFitterFunctionImpl>(
      std::move(trackFitter));
  fitterFunction->maxComponents = maxComponents;
  fitterFunction->abortOnError = abortOnError;
  fitterFunction->disableAllMaterialHandling = disableAllMaterialHandling;

  return fitterFunction;
}
