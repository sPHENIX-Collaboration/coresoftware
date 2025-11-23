#pragma once

#include <fun4all/SubsysReco.h>

// #include <trackbase/ClusterErrorPara.h>
#include <trackbase/TrkrDefs.h>

// #include <tpc/TpcClusterZCrossingCorrection.h>
#include <tpc/TpcGlobalPositionWrapper.h>

#include <Eigen/Dense>

#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>

#include <array>
#include <functional>
#include <type_traits>
#include <set>

class PHCompositeNode;
class ActsGeometry;
class SvtxTrackMap;
class SvtxAlignmentStateMap;
class SvtxVertexMap;
class TrackSeed;
class TrackSeedContainer;
class TrkrCluster;
class TrkrClusterContainer;
class WeightedTrackMap;
class WeightedTrack;

class TFile;
class TNtuple;

class WeightedFitter : public SubsysReco {
public:
	enum which_track_e { k_silicon_tracks, k_tpc_tracks, k_svtx_tracks, };

	WeightedFitter(std::string const& = "WeightedFitter");
	virtual ~WeightedFitter() override;

	int InitRun(PHCompositeNode*) override;
	int process_event(PHCompositeNode*) override;
	int End(PHCompositeNode*) override;

	/// Input setting methods
	void set_which_tracks (which_track_e which_tracks) { m_which_tracks = which_tracks; }
	void set_trkr_cluster_container_node_name (std::string const& name) { m_trkr_cluster_container_node_name = name; }
	void set_svtx_track_seed_container_node_name (std::string const& name) { m_svtx_track_seed_container_node_name = name; }
	void set_silicon_track_seed_container_node_name (std::string const& name) { m_silicon_track_seed_container_node_name = name; }
	void set_tpc_track_seed_container_node_name (std::string const& name) { m_tpc_track_seed_container_node_name = name; }

	/// Set the minimum number of clusters for the track to be used
	void set_min_num_mvtx (int num_mvtx) { m_min_num_mvtx = num_mvtx; }
	void set_min_num_intt (int num_intt) { m_min_num_intt = num_intt; }
	void set_min_num_tpc (int num_tpc) { m_min_num_tpc = num_tpc; }
	void set_min_num_tpot (int num_tpot) { m_min_num_tpot = num_tpot; }

	void use_vertex (bool const& use_vertex = true) { m_use_vertex = use_vertex; }
	void set_vertex_node_name (std::string const& name) { m_vertex_map_node_name = name; }

	/// If this method is called such that the internal m_included_layers member is non-empty,
	/// only those layers which have been explicitly added will be used in the fit
	void include_layer (int layer) { m_included_layers.insert(layer); m_excluded_layers.erase(layer); }
	/// Explicitly exclude layers
	void exclude_layer (int layer) { m_excluded_layers.insert(layer); m_included_layers.erase(layer); }

	/// Set the type of tracks to use for the fit via a template argument
	/// By default, WeightedFitterZeroField tracks are used
	template <typename T>
	void set_track_type() {
		// Helps clarify compilation failures to users
		static_assert (std::is_convertible_v<T*, WeightedTrack*>, "Template argument must be derived from trackbase_historic/WeightedTrack");
		m_track_factory = []{ return new T; };
	}

	/// Set the method to produce "new" (owned by an instance of this class) tracks
	/// Useful if some configuration beyond default construction, otherwise, see set_track_type
	void set_track_factory (std::function<WeightedTrack*(void)> track_factory) { m_track_factory = track_factory; }

	/// Output setting methods
	void set_track_map_node_name (std::string const& name) { m_track_map_node_name = name; }
	void set_alignment_map_node_name (std::string const& name) { m_alignment_map_node_name = name; }

	void set_ntuple_file_name (std::string const& name) { m_ntuple_file_name = name; }
	void reassign_cluster_sides (bool reassign_sides = true) { m_reassign_sides = reassign_sides; }

private:
	void get_nodes(PHCompositeNode*);
	void make_nodes(PHCompositeNode*);
	void make_ntuple();

	bool get_cluster_keys (TrackSeed*);
	bool get_points();
	bool do_fit();
	bool refit_with_vertex();
	bool add_track();

	std::string m_ntuple_file_name{}; // empty--by default do not make it
	TFile* m_file{};
	TNtuple* m_ntuple{};

	std::string m_geometry_node_name{"ActsGeometry"};
	ActsGeometry* m_geometry{};

	std::string m_trkr_cluster_container_node_name{"TRKR_CLUSTER"};
	TrkrClusterContainer* m_trkr_cluster_container{};

	which_track_e m_which_tracks = k_svtx_tracks;

	std::string m_silicon_track_seed_container_node_name{"SiliconTrackSeedContainer"};
	TrackSeedContainer* m_silicon_track_seed_container{};
	TrackSeed* m_silicon_seed{};

	std::string m_tpc_track_seed_container_node_name{"TpcTrackSeedContainer"};
	TrackSeedContainer* m_tpc_track_seed_container{};
	TrackSeed* m_tpc_seed{};

	std::string m_svtx_track_seed_container_node_name{"SvtxTrackSeedContainer"};
	TrackSeedContainer* m_svtx_track_seed_container{};

	bool m_use_vertex{false};
	std::string m_vertex_map_node_name{"SvtxVertexMap"};
	SvtxVertexMap* m_vertex_map{};

	// Nodes that this class will populate
	std::string m_track_map_node_name{"WeightedFitterTrackMap"};
	SvtxTrackMap* m_track_map{};

	std::string m_alignment_map_node_name{"WeightedFitterAlignmentStateMap"};
	SvtxAlignmentStateMap* m_alignment_map{};

	// Layers to be included--if populated, only those layers will be included
	std::set<int> m_included_layers;
	// Layers to be excluded
	std::set<int> m_excluded_layers;

	int m_track_id{0};
	int m_crossing{SHRT_MAX};
	std::vector<TrkrDefs::cluskey> m_cluster_keys;

	ROOT::Math::Minimizer* m_minimizer{};
	WeightedTrack* m_weighted_track{};
	std::function<WeightedTrack*(void)> m_track_factory;

	bool m_reassign_sides{false};
	int m_num_mvtx{0};
	int m_num_intt{0};
	int m_num_tpc{0};
	int m_num_tpot{0};
	int m_side{-1};

	int m_min_num_mvtx{0}; // 3
	int m_min_num_intt{0}; // 2
	int m_min_num_tpc{0}; // 20
	int m_min_num_tpot{0}; // 0

	TpcGlobalPositionWrapper m_global_position_wrapper;
	// ClusterErrorPara m_cluster_error_para;
	// TpcClusterZCrossingCorrection m_clusterCrossingCorrection; // unused for now, maybe to get uncorrected positions later
};

