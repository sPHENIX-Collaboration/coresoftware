#ifndef TRACKBASE_RESIDUALOUTLIERFINDER_H
#define TRACKBASE_RESIDUALOUTLIERFINDER_H

#include <TFile.h>
#include <TH2.h>
#include <TNtuple.h>
#include <Acts/Definitions/Units.hpp>
#include <Acts/EventData/Measurement.hpp>
#include <Acts/EventData/MeasurementHelpers.hpp>
#include <Acts/EventData/MultiTrajectory.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>
struct ResidualOutlierFinder
{
  int verbosity = 0;
  std::map<long unsigned int, float> chi2Cuts;
  std::string outfilename = "OutlierFinder.root";
  TH2 *hChi2 = new TH2F("h_layerChi2", ";layer;#chi^{2}", 60, 0, 60, 1000, 0, 500);

  TH2 *hDistance = new TH2F("h_layerDistance", ";layer;distance", 60, 0, 60, 1000, 0, 500);

  TNtuple *tree = new TNtuple("ntp_outlierfinder", "tree with predicted states and cluster info",
                              "sphenixlayer:layer:volume:distance:chi2:predgx:predgy:predgz:predlx:predlz:clusgx:clusgy:clusgz:cluslx:cluslz");
  void outfileName(const std::string& name)
  {
    outfilename = name;
  }
  void Write()
  {
    TFile *file = new TFile(outfilename.c_str(), "recreate");

    file->cd();
    hChi2->Write();
    hDistance->Write();
    tree->Write();
    file->ls();
  }
  bool operator()(Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::ConstTrackStateProxy state) const
  {
    // can't determine an outlier w/o a measurement or predicted parameters
    if (!state.hasCalibrated() || !state.hasPredicted())
    {
      return false;
    }
    
    auto sourceLink = state.getUncalibratedSourceLink().template get<ActsSourceLink>();
    const auto& cluskey = sourceLink.cluskey();

    const auto predicted = state.predicted();
    const auto predictedCovariance = state.predictedCovariance();
    float chi2 = std::numeric_limits<float>::max();

    auto fullCalibrated = state
                              .template calibrated<Acts::MultiTrajectoryTraits::MeasurementSizeMax>()
                              .data();
    auto fullCalibratedCovariance = state
                                        .template calibratedCovariance<Acts::MultiTrajectoryTraits::MeasurementSizeMax>()
                                        .data();

    chi2 = Acts::visit_measurement(state.calibratedSize(), [&](auto N) -> double
                                   {
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
	
	return chi2; });

    float distance = Acts::visit_measurement(state.calibratedSize(), [&](auto N)
                                             {
      constexpr size_t kMeasurementSize = decltype(N)::value;
      auto residuals =
          state.template calibrated<kMeasurementSize>() -
          state.projector()
      .template topLeftCorner<kMeasurementSize, Acts::eBoundSize>() *
              state.predicted();
      auto cdistance = residuals.norm();
      return cdistance; });

    if (verbosity > 2)
    {
      std::cout << "Measurement has distance, chi2 "
                << distance << ", " << chi2
                << std::endl;
    }

    auto volume = state.referenceSurface().geometryId().volume();
    auto layer = state.referenceSurface().geometryId().layer();

    bool outlier = false;
    float chi2cut = chi2Cuts.find(volume)->second;
    int sphenixlayer = TrkrDefs::getLayer(cluskey);
    hChi2->Fill(sphenixlayer, chi2);
    hDistance->Fill(sphenixlayer, distance);

    float data[] = {
        (float) sphenixlayer, (float) layer, (float) volume, distance, chi2,
        (float) predicted[Acts::eFreePos0], (float) predicted[Acts::eFreePos1], (float) predicted[Acts::eFreePos2],
        (float) predicted[Acts::eBoundLoc0], (float) predicted[Acts::eBoundLoc1],
        (float) fullCalibrated[Acts::eFreePos0], (float) fullCalibrated[Acts::eFreePos1], (float) fullCalibrated[Acts::eFreePos2],
        (float) fullCalibrated[Acts::eBoundLoc0], (float) fullCalibrated[Acts::eBoundLoc1]};
    tree->Fill(data);

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
