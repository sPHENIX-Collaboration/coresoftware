#ifndef TRACKRECO_PHACTSSILICONMICROMEGASFITTER_H
#define TRACKRECO_PHACTSSILICONMICROMEGASFITTER_H

#include "PHTrackFitting.h"
#include "PHActsSourceLinks.h"
#include "ActsTrack.h"

#include <trackbase/TrkrDefs.h>

#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <Acts/EventData/MeasurementHelpers.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>

#include <ActsExamples/Fitting/TrkrClusterFittingAlgorithm.hpp>
#include <ActsExamples/EventData/TrkrClusterMultiTrajectory.hpp>

#include <memory>
#include <string>

class TFile;
class TH1;
class TH2;

using SourceLink = ActsExamples::TrkrClusterSourceLink;
using SourceLinkVec = std::vector<SourceLink>;
using FitResult = Acts::KalmanFitterResult<SourceLink>;
using Trajectory = ActsExamples::TrkrClusterMultiTrajectory;
using Measurement = Acts::Measurement<ActsExamples::TrkrClusterSourceLink,
                                      Acts::BoundIndices,
                                      Acts::eBoundLoc0,
                                      Acts::eBoundLoc1>;
using SurfaceVec = std::vector<const Acts::Surface*>;

class PHActsSiliconMicromegasFitter : public PHTrackFitting
{

 public:
  PHActsSiliconMicromegasFitter(
      const std::string &name = "PHActsSiliconMicromegasFitter");

  ~PHActsSiliconMicromegasFitter(){}
  
  int End(PHCompositeNode *topNode);
  int Setup(PHCompositeNode *topNode);
  int Process();
  int ResetEvent(PHCompositeNode *topNode);
  
 protected:

 private:

 int getNodes(PHCompositeNode *topNode);
 int createNodes(PHCompositeNode *topNode);
 SourceLinkVec getSiliconMMsSls(SourceLinkVec trackSls, 
			        SurfaceVec &surfaces);
 void checkSurfaceVec(SurfaceVec &surfaces);
 void updateSvtxTrack(Trajectory traj, const unsigned int trackKey,
		      Acts::Vector3D vertex);

 int m_event;
 SvtxTrackMap *m_trackMap;
 
 std::map<TrkrDefs::cluskey, unsigned int> *m_hitIdClusKey;
std::map<const unsigned int, Trajectory> *m_actsFitResults;
 std::map<unsigned int, ActsTrack> *m_actsProtoTracks;
 ActsTrackingGeometry *m_tGeometry;
 ActsExamples::TrkrClusterFittingAlgorithm::Config m_fitCfg;

};

#endif 
