#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wunused-value"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#include <Acts/Definitions/TrackParametrization.hpp>
#pragma GCC diagnostic ignored "-Wshadow"
#include <Acts/TrackFitting/GaussianSumFitter.hpp>
#pragma GCC diagnostic pop

#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/TrackFitting/BetheHeitlerApprox.hpp>
#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>
#include <Acts/Utilities/Helpers.hpp>

#include "ActsTrackFittingAlgorithm.h"

namespace
{
  /// This type is used in the Examples framework for the Bethe-Heitler
  /// approximation
  using BetheHeitlerApprox = Acts::Experimental::AtlasBetheHeitlerApprox<6, 5>;

  using MultiStepper = Acts::MultiEigenStepperLoop<>;
  using Propagator = Acts::Propagator<MultiStepper, Acts::Navigator>;
  using DirectPropagator = Acts::Propagator<MultiStepper, Acts::DirectNavigator>;

  using Fitter =
      Acts::Experimental::GaussianSumFitter<Propagator,
                                            BetheHeitlerApprox,
                                            Acts::VectorMultiTrajectory>;
  using DirectFitter =
      Acts::Experimental::GaussianSumFitter<DirectPropagator,
                                            BetheHeitlerApprox,
                                            Acts::VectorMultiTrajectory>;

  struct GsfFitterFunctionImpl
    : public ActsTrackFittingAlgorithm::TrackFitterFunction
  {
    Fitter fitter;

    Acts::GainMatrixUpdater updater;

    std::size_t maxComponents = 0;
    double weightCutoff = 0;
    bool abortOnError = false;
    bool disableAllMaterialHandling = false;

    GsfFitterFunctionImpl(Fitter&& f)
      : fitter(std::move(f))
    {
    }

    auto makeGsfOptions(
        const ActsTrackFittingAlgorithm::GeneralFitterOptions& options)
        const
    {
      Acts::Experimental::GsfExtensions<Acts::VectorMultiTrajectory> extensions;
      // cppcheck-suppress constStatement
      extensions.updater.connect<&Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(&updater);

      Acts::Experimental::GsfOptions<Acts::VectorMultiTrajectory> gsfOptions{
          options.geoContext,
          options.magFieldContext,
          options.calibrationContext,
          extensions,
          options.logger,
          options.propOptions,
          &(*options.referenceSurface),
          maxComponents,
          abortOnError,
          disableAllMaterialHandling};

      gsfOptions.extensions.calibrator
          .template connect<&Calibrator::calibrate>(
              &options.calibrator.get());

      return gsfOptions;
    }

    ActsTrackFittingAlgorithm::TrackFitterResult operator()(
        const std::vector<std::reference_wrapper<
            const ActsSourceLink>>& sourceLinks,
        const ActsTrackFittingAlgorithm::TrackParameters& initialParameters,
        const ActsTrackFittingAlgorithm::GeneralFitterOptions& options,
        std::shared_ptr<Acts::VectorMultiTrajectory>& trajectory) const override
    {
      const auto gsfOptions = makeGsfOptions(options);
      return fitter.fit(sourceLinks.begin(), sourceLinks.end(), initialParameters,
                        gsfOptions, trajectory);
    }
  };

}  // namespace

// Have a separate class befriend the main class to ensure that GSF specific
// track fitting headers only stay here to avoid library clashes
class ActsGsfTrackFittingAlgorithm
{
 public:
  friend class ActsTrackFittingAlgorithm;

  std::shared_ptr<ActsTrackFittingAlgorithm::TrackFitterFunction>
  makeGsfFitterFunction(
      std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
      std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
      BetheHeitlerApprox betheHeitlerApprox, std::size_t maxComponents,
      Acts::FinalReductionMethod finalReductionMethod, bool abortOnError,
      bool disableAllMaterialHandling);
};
