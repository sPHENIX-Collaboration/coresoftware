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

 private:

  int GetNodes(PHCompositeNode* topNode);

  double findDcaTwoTracks(SvtxTrack *tr1, SvtxTrack *tr2);
double dcaTwoLines(const Eigen::Vector3d &p1, const Eigen::Vector3d &v1, const Eigen::Vector3d &p2, const Eigen::Vector3d &v2, 
		     Eigen::Vector3d &PCA1, Eigen::Vector3d &PCA2);

SvtxTrackMap *_track_map{nullptr};
SvtxTrack *_track{nullptr};

SvtxVertexMap *_svtx_vertex_map{nullptr};

double dcacut = 0.0080;  // 80 microns 

std::multimap<unsigned int, unsigned int> _vertex_track_map;
std::multimap<unsigned int, std::pair<unsigned int, double>> _track_pair_map;
std::multimap<unsigned int, std::pair<unsigned int, Eigen::Vector3d>> _track_pca_map;
std::map<unsigned int, Eigen::Vector3d> _vertex_position_map;
std::set<unsigned int> _connected;
std::set<unsigned int> _vertex_set;

};

#endif // PHSIMPLEVERTEXFINDER_H
