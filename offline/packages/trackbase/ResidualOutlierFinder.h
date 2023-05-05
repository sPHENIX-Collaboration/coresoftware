#ifndef TRACKBASE_RESIDUALOUTLIERFINDER_H
#define TRACKBASE_RESIDUALOUTLIERFINDER_H

#include <Acts/Definitions/Units.hpp>
#include <Acts/EventData/Measurement.hpp>
#include <Acts/EventData/MeasurementHelpers.hpp>
#include <Acts/EventData/MultiTrajectory.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>

struct ResidualOutlierFinder
{
  int verbosity = 0;
  std::map<long unsigned int, float> chi2Cuts;

  bool operator()(Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::ConstTrackStateProxy state) const
  {
    // can't determine an outlier w/o a measurement or predicted parameters
    if (!state.hasCalibrated() || !state.hasPredicted())
    {
      return false;
    }
  
    const auto predicted = state.predicted();
    const auto predictedCovariance = state.predictedCovariance();
    double chi2 = std::numeric_limits<float>::max();
    
    auto fullCalibrated = state
      .template calibrated<Acts::MultiTrajectoryTraits::MeasurementSizeMax>().data();
    auto fullCalibratedCovariance = state
      .template calibratedCovariance<Acts::MultiTrajectoryTraits::MeasurementSizeMax>().data();

    chi2 = Acts::visit_measurement(state.calibratedSize(), [&](auto N) -> double {
	constexpr size_t kMeasurementSize = decltype(N)::value;
	typename Acts::TrackStateTraits<kMeasurementSize, true>::Measurement calibrated{
	  fullCalibrated};

	typename Acts::TrackStateTraits<kMeasurementSize, true>::MeasurementCovariance
	  calibratedCovariance{fullCalibratedCovariance};

	using ParametersVector = Acts::ActsVector<kMeasurementSize>;
	const auto H = state.projector().template topLeftCorner<kMeasurementSize, Acts::eBoundSize>().eval();
	ParametersVector res;
	res = calibrated - H * predicted;
	chi2 = (res.transpose() * ((calibratedCovariance + H * predictedCovariance * H.transpose())).inverse() * res).eval()(0, 0);
	
	return chi2;
      });

    if (verbosity > 2)
    {
      auto distance = Acts::visit_measurement(state.calibratedSize(), [&](auto N) {
      constexpr size_t kMeasurementSize = decltype(N)::value;
      auto residuals =
          state.template calibrated<kMeasurementSize>() -
          state.projector()
      .template topLeftCorner<kMeasurementSize, Acts::eBoundSize>() *
              state.predicted();
      auto cdistance = residuals.norm();
      return cdistance;
    });
      std::cout << "Measurement has distance, chi2 "
                << distance << ", " << chi2
                << std::endl;
    }

    auto volume = state.referenceSurface().geometryId().volume();
    auto layer = state.referenceSurface().geometryId().layer();

    bool outlier = false;
    float chi2cut = chi2Cuts.find(volume)->second;
    if (chi2 > chi2cut)
    {
      outlier = true;
    }

    if (verbosity > 2)
    {
      std::cout << "Meas vol id and layer " << volume << ", " << layer
                << " and chi2cut "
                << chi2cut << " so returning outlier : " << outlier
                << std::endl;
    }

    return outlier;
  }
};

#endif
