#ifndef TRACKBASE_RESIDUALOUTLIERFINDER_H
#define TRACKBASE_RESIDUALOUTLIERFINDER_H

#include <phool/phool.h>

#include <TFile.h>
#include <TH2.h>
#include <TNtuple.h>

#include <Acts/Definitions/Units.hpp>

#include <ActsExamples/EventData/Measurement.hpp>

#include <Acts/EventData/MeasurementHelpers.hpp>
#include <Acts/EventData/MultiTrajectory.hpp>
#include <Acts/EventData/MultiTrajectoryHelpers.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>

#include <Acts/Utilities/TrackHelpers.hpp>

struct ResidualOutlierFinder
{
  ActsGeometry* m_tGeometry = nullptr;
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

    double chi2 = Acts::calculatePredictedChi2(state);

    float distance = Acts::visit_measurement(state.calibratedSize(), [&](auto N)
                                             {
      constexpr size_t kMeasurementSize = decltype(N)::value;
      auto [residual, residualCovariance] =
          calculatePredictedResidual<kMeasurementSize>(state);
      return residual.norm(); });

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
    if(!m_tGeometry)
    {
      std::cout << PHWHERE << "no geometry set in residual outlier finder" << std::endl;
      exit(1);
    }
    const auto predicted = state.predicted();
    auto fullCalibrated = state
                              .template calibrated<Acts::kMeasurementSizeMax>()
                              .data();
    Acts::FreeVector freeParams =
        Acts::transformBoundToFreeParameters(state.referenceSurface(),
                                                     m_tGeometry->geometry().getGeoContext(),
                                                     predicted);
    Acts::Vector2 local(fullCalibrated[Acts::eBoundLoc0], fullCalibrated[Acts::eBoundLoc1]);
    Acts::Vector3 global = state.referenceSurface().localToGlobal(
        m_tGeometry->geometry().getGeoContext(), local,
        Acts::Vector3(1, 1, 1));
    float data[] = {
        (float) sphenixlayer, (float) layer, (float) volume, distance, (float)chi2,
        (float) freeParams[Acts::eFreePos0], (float) freeParams[Acts::eFreePos1], (float) freeParams[Acts::eFreePos2],
        (float) predicted[Acts::eBoundLoc0], (float) predicted[Acts::eBoundLoc1],
        (float) global[Acts::eFreePos0], (float) global[Acts::eFreePos1], (float) global[Acts::eFreePos2],
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
