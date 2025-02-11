// Tell emacs that this is a C++ source
//  -*- C++ -*-.

/*!
 *  \file		  PHSimpleVertexFinder
 *  \brief		Class for determining primary vertices usin g fitted tracks
 *  \author	 Tony Frawley <afrawley@fsu.edu>
 */

#ifndef PHSIMPLEVERTEXFINDER_H
#define PHSIMPLEVERTEXFINDER_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/ActsGeometry.h>

#include <map>
#include <set>
#include <string>
#include <vector>

#include <Eigen/Dense>

class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class SvtxVertexMap;
class TrkrCluster;
class TrackVertexCrossingAssoc;
class TrackSeed;
class TrkrClusterContainer;
class ActsGeometry;

class PHSimpleVertexFinder : public SubsysReco
{
 public:
  PHSimpleVertexFinder(const std::string &name = "PHSimpleVertexFinder");

  ~PHSimpleVertexFinder() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  void setBeamLineCut(const double cut) { _beamline_xy_cut = cut; }
  void setDcaCut(const double cut) { _base_dcacut = cut; }
  void setTrackQualityCut(double cut) { _qual_cut = cut; }
  void setRequireMVTX(bool set) { _require_mvtx = set; }
  void setNmvtxRequired(unsigned int n) { _nmvtx_required = n; }
  void setTrackPtCut(const double cut) { _track_pt_cut = cut; }
  // void setUseTrackCovariance(bool set) {_use_track_covariance = set;}
  void setOutlierPairCut(const double cut) { _outlier_cut = cut; }
  void setTrackMapName(const std::string &name) { _track_map_name = name; }
  void setVertexMapName(const std::string &name) { _vertex_map_name = name; }
  void zeroField(const bool flag) { _zero_field = flag; }
  void setTrkrClusterContainerName(std::string &name){ m_clusterContainerName = name; }

 private:
  int GetNodes(PHCompositeNode *topNode);
  int CreateNodes(PHCompositeNode *topNode);

  void checkDCAs(SvtxTrackMap *track_map);
  void checkDCAsZF(SvtxTrackMap *track_map);
  void checkDCAs();

  void getTrackletClusterList(TrackSeed* tracklet, std::vector<TrkrDefs::cluskey>& cluskey_vec);
  
  void findDcaTwoTracks(SvtxTrack *tr1, SvtxTrack *tr2);
  double dcaTwoLines(const Eigen::Vector3d &p1, const Eigen::Vector3d &v1,
                     const Eigen::Vector3d &p2, const Eigen::Vector3d &v2,
                     Eigen::Vector3d &PCA1, Eigen::Vector3d &PCA2);
  std::vector<std::set<unsigned int>> findConnectedTracks();
  void removeOutlierTrackPairs();
  double getMedian(std::vector<double> &v);
  double getAverage(std::vector<double> &v);

  SvtxTrackMap *_track_map{nullptr};
  TrkrClusterContainer* _cluster_map{nullptr};
  SvtxVertexMap *_svtx_vertex_map{nullptr};
  ActsGeometry* _tGeometry{nullptr};

  double _base_dcacut = 0.0080;  // 80 microns
  double _active_dcacut = 0.080;
  double _beamline_xy_cut = 0.2;  // must be within 2 mm of beam line
  double _qual_cut = 10.0;
  bool _require_mvtx = true;
  unsigned int _nmvtx_required = 3;
  double _track_pt_cut = 0.0;
  double _outlier_cut = 0.015;
  //name of TRKR_CLUSTER Container
  std::string m_clusterContainerName = "TRKR_CLUSTER";

  bool _zero_field = false;     // fit straight lines if true

  std::string _track_map_name = "SvtxTrackMap";
  std::string _vertex_map_name = "SvtxVertexMap";
  std::multimap<unsigned int, unsigned int> _vertex_track_map;
  using matrix_t = Eigen::Matrix<double, 3, 3>;
  std::multimap<unsigned int, std::pair<unsigned int, double>> _track_pair_map;
  // Eigen::Vector3d is an Eigen::Matrix<double,3,1>
  std::multimap<unsigned int, std::pair<unsigned int, std::pair<Eigen::Vector3d,
                                                                Eigen::Vector3d>>>
      _track_pair_pca_map;
  std::map<unsigned int, Eigen::Vector3d> _vertex_position_map;
  std::map<unsigned int, matrix_t> _vertex_covariance_map;
  std::set<unsigned int> _vertex_set;

  TrackVertexCrossingAssoc *_track_vertex_crossing_map{nullptr};
};

#endif  // PHSIMPLEVERTEXFINDER_H
