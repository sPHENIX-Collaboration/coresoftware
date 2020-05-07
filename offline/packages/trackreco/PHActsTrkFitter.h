/*!
 *  \file		PHActsTrkFitter.h
 *  \brief		Refit SvtxTracks with Acts.
 *  \details	Refit SvtxTracks with Acts
 *  \author		Tony Frawley <afrawley@fsu.edu>
 */

#ifndef TRACKRECO_ACTSTRKFITTER_H
#define TRACKRECO_ACTSTRKFITTER_H

#include "PHTrackFitting.h"
#include "PHActsSourceLinks.h"

#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <Acts/EventData/MeasurementHelpers.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>

#include <ACTFW/Fitting/TrkrClusterFittingAlgorithm.hpp>

#include <memory>
#include <string>

namespace FW
{
  namespace Data
  {
    class TrkrClusterSourceLink;
  }
}  // namespace FW

class ActsTrack;
class MakeActsGeometry;
class SvtxTrack;
class SvtxTrackMap;
using SourceLink = FW::Data::TrkrClusterSourceLink;

class PHActsTrkFitter : public PHTrackFitting
{
 public:
  /// Default constructor
  PHActsTrkFitter(const std::string& name = "PHActsTrkFitter");

  /// Destructor
  ~PHActsTrkFitter();

  /// End, write and close files
  int End(PHCompositeNode*);

  /// Get and create nodes
  int Setup(PHCompositeNode* topNode);

  /// Process each event by calling the fitter
  int Process();

 private:
  /// Event counter
  int m_event;

  /// Get all the nodes
  int getNodes(PHCompositeNode*);

  /// Create new nodes
  int createNodes(PHCompositeNode*);

  Acts::BoundSymMatrix rotateCovarianceLocalToGlobal(const Acts::KalmanFitterResult<SourceLink>& fitOutput);

  /// Convert the acts track fit result to an svtx track
  void updateSvtxTrack(const Acts::KalmanFitterResult<SourceLink>& fitOutput, const unsigned int trackKey);

  /// Map of acts tracks and track key created by PHActsTracks
  std::map<unsigned int, ActsTrack>* m_actsProtoTracks;

  /// Options that Acts::Fitter needs to run from MakeActsGeometry
  ActsTrackingGeometry *m_tGeometry;

  /// Configuration containing the fitting function instance
  FW::TrkrClusterFittingAlgorithm::Config fitCfg;

  /// TrackMap containing SvtxTracks
  SvtxTrackMap *m_trackMap;

};

#endif
