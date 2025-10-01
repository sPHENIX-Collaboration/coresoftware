#pragma once

#include <fun4all/SubsysReco.h>
#include <tpc/TpcClusterZCrossingCorrection.h>
#include <tpc/TpcGlobalPositionWrapper.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/TrkrDefs.h>

#include <Eigen/Dense>

#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>

#include <array>

class PHCompositeNode;
class ActsGeometry;
class TrackSeed;
class TrackSeedContainer;
class TrkrCluster;
class TrkrClusterContainer;
class SvtxTrackMap;
class SvtxAlignmentStateMap;

class TPad;

class WeightedFitter : public SubsysReco
{
 public:
  class ClusterFitPoint
  {
   public:
    TrkrDefs::cluskey cluster_key;

    // All are expected in global coordinates
    Eigen::Vector3d pos;  // cluster position
    Eigen::Vector3d o;    // sensor center
    Eigen::Vector3d x;    // sensor local x axis
    Eigen::Vector3d y;    // sensor local y axis
    Eigen::Vector3d z;    // sensor local z axis

    // Cluster uncertainties
    double sigma_x{1.0};  // Local x ~ Global rphi
    double sigma_y{1.0};  // Local y ~ Global z
  };

  class FitErrorCalculator
  {
   public:
    FitErrorCalculator(std::vector<ClusterFitPoint> const&);
    ~FitErrorCalculator() = default;

    double operator()(double const*) const;

   private:
    std::vector<ClusterFitPoint> m_points;
  };

  WeightedFitter(std::string const& = "WeightedFitter");
  virtual ~WeightedFitter() override;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  void set_track_seed_container_node_name(std::string const& name) { m_track_seed_container_node_name = name; }
  void set_trkr_cluster_container_node_name(std::string const& name) { m_trkr_cluster_container_node_name = name; }

  // minimum required MVTX clusters to use the seed
  void set_min_mvtx(int const& num_mvtx) { m_num_mvtx = num_mvtx; }
  // minimum required INTT clusters to use the seed
  void set_min_intt(int const& num_intt) { m_num_intt = num_intt; }
  // minimum required TPC clusters to use the seed
  void set_min_tpc(int const& num_tpc) { m_num_tpc = num_tpc; }

 private:
  void get_nodes(PHCompositeNode*);

  bool get_points(TrackSeed const*);
  bool do_fit();
  bool add_track(TrackSeed*);

  void modify_point_tpc(TrkrDefs::cluskey, TrkrCluster*, ClusterFitPoint&);

  void draw();
  void draw_cluster_xy(TPad*, ClusterFitPoint const&);
  void draw_fit_xy(TPad*);

  std::string m_geometry_node_name{"ActsGeometry"};
  ActsGeometry* m_geometry{nullptr};

  std::string m_track_seed_container_node_name{"SvtxTrackSeedContainer"};
  TrackSeedContainer* m_track_seed_container{nullptr};

  std::string m_trkr_cluster_container_node_name{"TRKR_CLUSTER"};
  TrkrClusterContainer* m_trkr_cluster_container{nullptr};

  std::string m_track_map_node_name{"WeightedFitterTrackMap"};
  SvtxTrackMap* m_track_map{nullptr};

  std::string m_alignment_map_node_name{"WeightedFitterAlignmentStateMap"};
  SvtxAlignmentStateMap* m_alignment_map{nullptr};

  int m_track_id{0};
  std::vector<ClusterFitPoint> m_points;
  ROOT::Math::Minimizer* m_minimizer{nullptr};

  int m_num_mvtx{3};  // minimum required MVTX clusters to use the seed
  int m_num_intt{2};  // minimum required INTT clusters to use the seed
  int m_num_tpc{30};  // minimum required TPC clusters to use the seed

  TpcGlobalPositionWrapper m_globalPositionWrapper;
  ClusterErrorPara m_cluster_error_para;
};
