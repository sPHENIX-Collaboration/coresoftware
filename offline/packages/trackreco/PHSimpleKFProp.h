/*!
 *  \file		  PHTrackPropagating.h
 *  \brief		Base class for track seeding
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef TRACKRECO_PHSIMPLEKFPROP_H
#define TRACKRECO_PHSIMPLEKFPROP_H

// PHENIX includes
#include "PHTrackPropagating.h"
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/SvtxTrack_v1.h>
#include <phfield/PHField.h>
#include "nanoflann.hpp"
#include "ALICEKF.h"

// STL includes
#include <string>
#include <vector>
#include <memory>

#include <Eigen/Core>

// forward declarations
class PHCompositeNode;

class TrkrClusterContainer;
class SvtxVertexMap;
class SvtxTrackMap;
class AssocInfoContainer;

/// \class PHTrackPropagating
///
/// \brief Base class for track seeding
///
class PHSimpleKFProp : public PHTrackPropagating
{
 public:
  PHSimpleKFProp(const std::string &name = "PHSimpleKFProp");
  virtual ~PHSimpleKFProp() {}

  //int InitRun(PHCompositeNode *topNode) override;
  //int process_event(PHCompositeNode *topNode) override;
  //int End(PHCompositeNode *topNode) override;
  //void set_track_map_name(const std::string &map_name) { _track_map_name = map_name; }
  //void SetUseTruthClusters(bool setit){_use_truth_clusters = setit;}
  void set_field_dir(const double rescale)
  {
    std::cout << "rescale: " << rescale << std::endl;
    _fieldDir = 1;
    if(rescale > 0)
      _fieldDir = -1;
  }
  void set_max_window(double s){_max_dist = s;}
 protected:
  /// setup interface for trackers, called in InitRun, setup things like pointers to nodes.
  /// overrided in derived classes
  int Setup(PHCompositeNode *topNode) override;

  /// process event interface for trackers, called in process_event.
  /// implemented in derived classes
  int Process() override;

  ///
  int End() override;


  //SvtxClusterMap *_cluster_map;
  //TrkrClusterContainer *_cluster_map;
  //SvtxVertexMap *_vertex_map;
  //SvtxTrackMap *_track_map;
  //AssocInfoContainer *_assoc_container;
  PHField* _field_map;

  //std::string _track_map_name;

  //bool _use_truth_clusters = false;

 private:
  /// fetch node pointers
  int get_nodes(PHCompositeNode *topNode);
  std::vector<double> radii;
  std::vector<double> _vertex_x;
  std::vector<double> _vertex_y;
  std::vector<double> _vertex_z;
  std::vector<double> _vertex_xerr;
  std::vector<double> _vertex_yerr;
  std::vector<double> _vertex_zerr;
  std::vector<double> _vertex_ids;
  double _Bzconst = 10*0.000299792458f;
  //double _Bz = 1.4*_Bzconst;
  double _max_dist = .05;
  size_t _min_clusters_per_track = 20;
  double _fieldDir = -1;
  double _max_sin_phi = 1.;
  void PrepareKDTrees();
  std::vector<TrkrDefs::cluskey> PropagateTrack(SvtxTrack* track);
  template <typename T>
  struct KDPointCloud
  {
    KDPointCloud<T>(){}
    std::vector<std::vector<T>> pts;
    inline size_t kdtree_get_point_count() const
    {
      return pts.size();
    }
    inline T kdtree_distance(const T* p1, const size_t idx_p2, size_t /*size*/) const
    {
      const T d0 = p1[0] - pts[idx_p2][0];
      const T d1 = p1[1] - pts[idx_p2][1];
      const T d2 = p1[2] - pts[idx_p2][2];
      return d0 * d0 + d1 * d1 + d2 * d2;
    }
    inline T kdtree_get_pt(const size_t idx, int dim) const
    {
      if (dim == 0)
        return pts[idx][0];
      else if (dim == 1)
        return pts[idx][1];
      else
        return pts[idx][2];
    }
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /*bb*/) const
    {
      return false;
    }
  };
  std::vector<std::shared_ptr<KDPointCloud<double>>> _ptclouds;
  std::vector<std::shared_ptr<nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, KDPointCloud<double>>,
                                                KDPointCloud<double>,3>>> _kdtrees;
  std::shared_ptr<ALICEKF> fitter;
  double get_Bz(double x, double y, double z);
  void publishSeeds(std::vector<SvtxTrack_v1>);
  void MoveToVertex();
  void MoveToFirstTPCCluster();
  void line_fit(std::vector<std::pair<double,double>> points, double &A, double &B);
  void line_fit_clusters(std::vector<TrkrCluster*> clusters, double &A, double &B);
  void CircleFitByTaubin(std::vector<std::pair<double,double>> points, double &R, double &X0, double &Y0);
};

#endif
