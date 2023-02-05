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
class TNtuple;
class TH2D;
class TFile;

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

  void setTrackDcaCut(const double cutxy, const double cutz) {_track_dcaxy_cut = cutxy; _track_dcaz_cut = cutz;}
 void setTwoTrackDcaCut(const double cut) {_two_track_dcacut = cut;}
 void setMinPathCut(const double cut) {_min_path_cut = cut;}
 void setTrackQualityCut(double cut) {_qual_cut = cut;}
 void setRequireMVTX(bool set) {_require_mvtx = set;}
 void setOutfileName(std::string filename) {outfile = filename;}
 void setDecayParticleMass(double mass) {_decaymass = mass;}

 private:

  int GetNodes(PHCompositeNode* topNode);

  bool  hasSiliconSeed(SvtxTrack* tr); 
  void outputTrackDetails(SvtxTrack *tr);
  void get_dca(SvtxTrack* track, float& dca3dxy, float& dca3dz, float& dca3dxysigma, float& dca3dzsigma);
  std::vector<float> fitClusters(TrackSeed *tracklet);
  void getTrackletClusters(TrackSeed *tracklet, std::vector<Eigen::Vector3d>& global_vec,
			   std::vector<TrkrDefs::cluskey>& cluskey_vec);
  bool circle_circle_intersection(double r0, double x0, double y0, double r1, double x1, double y1, std::vector<double>& intersectionXY);
  void makeTpcGlobalCorrections(TrkrDefs::cluskey cluster_key, short int crossing, Eigen::Vector3d &global);
  Acts::BoundTrackParameters makeTrackParams(SvtxTrack* track);
  BoundTrackParamResult propagateTrack(const Acts::BoundTrackParameters& params, const Eigen::Vector3d PCA);
  void updateSvtxTrack(SvtxTrack* track, const Acts::BoundTrackParameters& params);
  Acts::Vector3 getVertex(SvtxTrack* track);
  bool projectTrackToCylinder(SvtxTrack* track, double Radius, Eigen::Vector3d& pos, Eigen::Vector3d& mom);
  void findPcaTwoLines(Eigen::Vector3d pos1, Eigen::Vector3d mom1, Eigen::Vector3d pos2, Eigen::Vector3d mom2,
		       double &dca, Eigen::Vector3d &PCA1, Eigen::Vector3d &PCA2);
  bool findTwoTrackIntersection(SvtxTrack *track1, SvtxTrack *track2, Eigen::Vector2d& intersect1,  Eigen::Vector2d& intersect2);
  BoundTrackParamResult propagateTrack(const Acts::BoundTrackParameters& params, const SurfacePtr& targetSurf);
  void fillNtp(SvtxTrack *track1, SvtxTrack *track2, double dca3dxy1, double dca3dz1, double dca3dxy2, double dca3dz2,  Eigen::Vector3d vpos1,  Eigen::Vector3d vpos2, Acts::Vector3 pca_rel1, Acts::Vector3 pca_rel2, double pair_dca, double invariantMass, double invariantPt, double path, int has_silicon_1, int has_siilicon_2);

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

 // these are minimal cuts used to make the ntuple
 // They can be tightened later when analyzing the ntuple
 double _track_dcaxy_cut = 0.020;  
 double _track_dcaz_cut = 0.020;  
 double _two_track_dcacut = 0.5;  // 5000 microns 
 double _qual_cut = 10.0;
 bool _require_mvtx = false;
 double _track_pt_cut = 0.0;
 double _min_path_cut = 0.2;
 double _max_intersection_radius = 40.0;  // discard intersections at greater than 40 cm radius
 double _projected_track_z_cut = 0.5;
 double _decaymass = 0.13957;  // pion, default
	      
  TH2D *recomass{nullptr};
  TNtuple *ntp{nullptr};
  std::string outfile;

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
