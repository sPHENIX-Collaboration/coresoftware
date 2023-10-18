#include "ActsGsfTrackFittingAlgorithm.h"

std::shared_ptr<ActsTrackFittingAlgorithm::TrackFitterFunction>
ActsGsfTrackFittingAlgorithm::makeGsfFitterFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
    BetheHeitlerApprox betheHeitlerApprox, std::size_t maxComponents,
    double weightCutoff,
    Acts::MixtureReductionMethod finalReductionMethod, bool abortOnError,
    bool disableAllMaterialHandling, const Acts::Logger& logger)
{
  MultiStepper stepper(magneticField, finalReductionMethod,
                       logger.cloneWithSuffix("GSFStep"));
  const auto& geo = *trackingGeometry;

  // Standard fitter
  Acts::Navigator::Config cfg{trackingGeometry};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg, logger.cloneWithSuffix("GSFNavigator"));
  Propagator propagator(std::move(stepper), std::move(navigator),
                        logger.cloneWithSuffix("GSFPropagator"));
  Fitter trackFitter(std::move(propagator),
                     BetheHeitlerApprox(betheHeitlerApprox),
                     logger.cloneWithSuffix("GSFFitter"));

  // build the fitter functions. owns the fitter object.
  auto fitterFunction = std::make_shared<GsfFitterFunctionImpl>(
      std::move(trackFitter), geo);
  fitterFunction->maxComponents = maxComponents;
  fitterFunction->weightCutoff = weightCutoff;
  fitterFunction->abortOnError = abortOnError;
  fitterFunction->disableAllMaterialHandling = disableAllMaterialHandling;
  fitterFunction->reductionMethod = finalReductionMethod;

  return fitterFunction;
}
