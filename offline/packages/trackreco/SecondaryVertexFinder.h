// Tell emacs that this is a C++ source
//  -*- C++ -*-.

/*!
 *  \file		  SecondaryVertexFinder
 *  \brief		Class for determining primary vertices usin g fitted tracks
 *  \author	 Tony Frawley <afrawley@fsu.edu>
 */

#ifndef SECONDARYVERTEXFINDER_H
#define SECONDARYVERTEXFINDER_H

#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/TpcDefs.h>

#include <tpc/TpcClusterZCrossingCorrection.h>
#include <tpc/TpcDistortionCorrection.h>

#include <Acts/Definitions/Algebra.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Acts/Propagator/Propagator.hpp>
#pragma GCC diagnostic pop

#include <Acts/Utilities/Result.hpp>
#include <Acts/Surfaces/CylinderSurface.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <ActsExamples/EventData/Trajectories.hpp>

#include <string>
#include <vector>
#include <map>
#include <set>

#include <Eigen/Dense>

class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class SvtxVertexMap;
class TrkrCluster;
class TrackSeed;
class ActsGeometry;
class TrkrClusterContainer;

using BoundTrackParam = const Acts::BoundTrackParameters;
using BoundTrackParamResult = Acts::Result<BoundTrackParam>;
using SurfacePtr = std::shared_ptr<const Acts::Surface>;
using Trajectory = ActsExamples::Trajectories;

class SecondaryVertexFinder : public SubsysReco
{
 public:

  SecondaryVertexFinder(const std::string &name = "SecondaryVertexFinder");

  ~SecondaryVertexFinder() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 void setTrackDcaCut(const double cut) {_track_dcacut = cut;}
 void setTwoTrackDcaCut(const double cut) {_two_track_dcacut = cut;}
 void setTrackQualityCut(double cut) {_qual_cut = cut;}
 void setRequireMVTX(bool set) {_require_mvtx = set;}
 void setNmvtxRequired(unsigned int n) {_nmvtx_required = n;}

 private:

  int GetNodes(PHCompositeNode* topNode);
//  int CreateNodes(PHCompositeNode* topNode);

void findPcaTwoTracks(SvtxTrack *tr1, SvtxTrack *tr2,
					      double &dca, Eigen::Vector3d &PCA1, Eigen::Vector3d &PCA2);
  double findTwoTrackPCA(SvtxTrack *track1, SvtxTrack *track2, Eigen::Vector3d &PCA1, Eigen::Vector3d &PCA2);  
  double dcaTwoLines(const Eigen::Vector3d &p1, const Eigen::Vector3d &v1, 
		     const Eigen::Vector3d &p2, const Eigen::Vector3d &v2, 
		     Eigen::Vector3d &PCA1, Eigen::Vector3d &PCA2);
  std::vector<float> fitClusters(TrackSeed *tracklet);
  void getTrackletClusters(TrackSeed *tracklet, std::vector<Eigen::Vector3d>& global_vec, std::vector<TrkrDefs::cluskey>& cluskey_vec);
  std::vector<double> circle_circle_intersection(double r0, double x0, double y0, double r1, double x1, double y1 );
  void makeTpcGlobalCorrections(TrkrDefs::cluskey cluster_key, short int crossing, Eigen::Vector3d &global);
  Acts::BoundTrackParameters makeTrackParams(SvtxTrack* track);
  BoundTrackParamResult propagateTrack(const Acts::BoundTrackParameters& params, const Eigen::Vector3d PCA);
  void updateSvtxTrack(SvtxTrack* track, const Acts::BoundTrackParameters& params);
  Acts::Vector3 getVertex(SvtxTrack* track);

  SvtxTrackMap *_track_map{nullptr};
  SvtxTrack *_track{nullptr};  
  SvtxVertexMap *_svtx_vertex_map{nullptr};
  TrkrClusterContainer *_cluster_map;
  ActsGeometry *_tGeometry;

  TpcClusterZCrossingCorrection _clusterCrossingCorrection;
  TpcDistortionCorrectionContainer* _dcc_static{nullptr};
  TpcDistortionCorrectionContainer* _dcc_average{nullptr};
  TpcDistortionCorrectionContainer* _dcc_fluctuation{nullptr};

 TpcDistortionCorrection _distortionCorrection;

  double _track_dcacut = 0.10;  
  double _two_track_dcacut = 0.050;  // 500 microns 
  double _qual_cut = 5.0;
  bool _require_mvtx = false;
  unsigned int _nmvtx_required = 3; 
  double _track_pt_cut = 0.0;

  std::multimap<unsigned int, unsigned int> _vertex_track_map;
  using matrix_t = Eigen::Matrix<double,3,3>;
  std::multimap<unsigned int, std::pair<unsigned int, double>> _track_pair_map;
  // Eigen::Vector3d is an Eigen::Matrix<double,3,1>  
  std::multimap<unsigned int, std::pair<unsigned int, std::pair<Eigen::Vector3d,
  Eigen::Vector3d>>>  _track_pair_pca_map;
  std::map<unsigned int, Eigen::Vector3d> _vertex_position_map;
  std::map<unsigned int, matrix_t> _vertex_covariance_map;
  std::set<unsigned int> _vertex_set;
  
};

#endif //SECONDARYVERTEXFINDER_H
