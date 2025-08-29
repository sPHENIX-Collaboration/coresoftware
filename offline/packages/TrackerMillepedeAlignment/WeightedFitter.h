#pragma once

#include <fun4all/SubsysReco.h>
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
class TrkrClusterContainer;
class SvtxTrackMap;
class SvtxAlignmentStateMap;

class TPad;

class ClusterFitPoint {
public:
	TrkrDefs::cluskey cluster_key;

	// All are expected in global coordinates
	Eigen::Vector3d pos; // cluster position
	Eigen::Vector3d o; // sensor center
	Eigen::Vector3d x; // sensor local x axis
	Eigen::Vector3d y; // sensor local x axis
	Eigen::Vector3d z; // sensor local x axis

	// Cluster uncertainties
	double sigma_x{1.0};
	double sigma_y{1.0};
};

class FitErrorCalculator {
public:
	FitErrorCalculator(std::vector<ClusterFitPoint> const&);
	~FitErrorCalculator() = default;

	double operator()(double const*) const;

private:
	std::vector<ClusterFitPoint> m_points;
};

class WeightedFitter : public SubsysReco {
public:
	WeightedFitter(std::string const& = "WeightedFitter");
	virtual ~WeightedFitter() override;

	int InitRun(PHCompositeNode*) override;
	int process_event(PHCompositeNode*) override;

	void set_track_seed_container_node_name (std::string const& name) { m_track_seed_container_node_name = name; }
	void set_trkr_cluster_container_node_name (std::string const& name) { m_trkr_cluster_container_node_name = name; }

private:
	void get_nodes(PHCompositeNode*);

	bool get_points(TrackSeed const*);
	bool do_fit();
	bool add_track(TrackSeed*);

	void draw();
	void draw_cluster_xy(TPad*, ClusterFitPoint const&);
	void draw_fit_xy(TPad*);

	std::string m_geometry_node_name{"ActsGeometry"};
	ActsGeometry* m_geometry{};

	std::string m_track_seed_container_node_name{"SvtxTrackSeedContainer"};
	TrackSeedContainer* m_track_seed_container{};

	std::string m_trkr_cluster_container_node_name{"TRKR_CLUSTER"};
	TrkrClusterContainer* m_trkr_cluster_container{};

	std::string m_track_map_node_name{"WeightedFitterTrackMap"};
	SvtxTrackMap* m_track_map{};

	std::string m_alignment_map_node_name{"WeightedFitterAlignmentStateMap"};
	SvtxAlignmentStateMap* m_alignment_map{};

	int m_track_id{0};
	std::vector<ClusterFitPoint> m_points;
	ROOT::Math::Minimizer* m_minimizer{nullptr};
};
