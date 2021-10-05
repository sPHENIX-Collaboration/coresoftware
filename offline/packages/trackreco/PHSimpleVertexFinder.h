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
class TpcSeedTrackMap;

class PHSimpleVertexFinder : public SubsysReco
{
 public:

  PHSimpleVertexFinder(const std::string &name = "PHSimpleVertexFinder");

  ~PHSimpleVertexFinder() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
 void setBeamLineCut(const double cut) {beamline_xy_cut = cut;}
 void setDcaCut(const double cut) {dcacut = cut;}
 void setTrackQualityCut(double cut) {qual_cut = cut;}
 void setRequireMVTX(bool set) {require_mvtx = set;}

 private:

  int GetNodes(PHCompositeNode* topNode);
  int CreateNodes(PHCompositeNode* topNode);
  
  void checkDCAs();

  void findDcaTwoTracks(SvtxTrack *tr1, SvtxTrack *tr2);  
  double dcaTwoLines(const Eigen::Vector3d &p1, const Eigen::Vector3d &v1, 
		     const Eigen::Vector3d &p2, const Eigen::Vector3d &v2, 
		     Eigen::Vector3d &PCA1, Eigen::Vector3d &PCA2);

  SvtxTrackMap *_track_map{nullptr};
  SvtxTrack *_track{nullptr};
  
  SvtxVertexMap *_svtx_vertex_map{nullptr};
  
  double dcacut = 0.0080;  // 80 microns 
  double beamline_xy_cut = 0.2;  // must be within 2 mm of beam line
  double qual_cut = 6.0;
  bool require_mvtx = true;
  
  std::multimap<unsigned int, unsigned int> _vertex_track_map;
  using matrix_t = Eigen::Matrix<float,3,3>;
  std::multimap<unsigned int, matrix_t> _vertex_cov_map;
  std::multimap<unsigned int, std::pair<unsigned int, double>> _track_pair_map;
  //std::multimap<unsigned int, std::pair<unsigned int, Eigen::Vector3d>> _track_pca_map;
  std::multimap<unsigned int, std::pair<unsigned int, std::pair<Eigen::Vector3d,
  Eigen::Vector3d>>>  _track_pair_pca_map;
  std::map<unsigned int, Eigen::Vector3d> _vertex_position_map;
  std::map<unsigned int, matrix_t> _vertex_covariance_map;
  std::set<unsigned int> _vertex_set;
  
};

#endif // PHSIMPLEVERTEXFINDER_H
