#include "WeightedFitter.h"

#include <trackbase/ActsGeometry.h>
#include <trackbase/ActsSurfaceMaps.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainerv4.h>

#include <trackbase_historic/SvtxAlignmentState_v1.h>
#include <trackbase_historic/SvtxAlignmentStateMap_v1.h>
#include <trackbase_historic/SvtxTrack_v4.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/SvtxTrackSeed_v1.h>
#include <trackbase_historic/SvtxTrackSeed_v2.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeed_v1.h>
#include <trackbase_historic/TrackSeed_v2.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/WeightedTrack.h>
#include <trackbase_historic/WeightedTrackZeroField.h>
#include <trackbase_historic/WeightedTrackMap.h>

#include <globalvertex/SvtxVertexMap_v1.h>
#include <globalvertex/SvtxVertex_v2.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/phool.h> // PHWHERE

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>

#include <TFile.h>
#include <TNtuple.h>

#include <iostream>
#include <functional>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

namespace { // anonymous
	template <typename T>
	T sqr (T const& t) { return t*t; }
}

WeightedFitter::WeightedFitter (
	std::string const& name
) : SubsysReco (name),
	m_minimizer(ROOT::Math::Factory::CreateMinimizer("Minuit2"))
{
	// Setup Minimizer
	m_minimizer->SetMaxFunctionCalls(1E8); // For Minuit/Minuit2
	m_minimizer->SetMaxIterations(1E6); // For GSL Minimizers--probably moot to call
	m_minimizer->SetTolerance(1E-8);
	m_minimizer->SetPrintLevel(0);

	set_track_type<WeightedTrackZeroField>();
}

WeightedFitter::~WeightedFitter (
) {
	delete m_weighted_track;
	delete m_minimizer;
}

void
WeightedFitter::get_nodes (
	PHCompositeNode* top_node
) {
	std::vector<std::string> missing_node_names;

	m_geometry = findNode::getClass<ActsGeometry>(top_node, m_geometry_node_name);
	if (!m_geometry) { missing_node_names.push_back(m_geometry_node_name); }

	m_trkr_cluster_container = findNode::getClass<TrkrClusterContainer>(top_node, m_trkr_cluster_container_node_name);
	if (!m_trkr_cluster_container) { missing_node_names.push_back(m_trkr_cluster_container_node_name); }

	if (m_which_tracks == k_svtx_tracks) {
		m_svtx_track_seed_container = findNode::getClass<TrackSeedContainer>(top_node, m_svtx_track_seed_container_node_name);
		if (!m_svtx_track_seed_container) { missing_node_names.push_back(m_svtx_track_seed_container_node_name); }
	}

	if (m_which_tracks == k_svtx_tracks || m_which_tracks == k_silicon_tracks) {
		m_silicon_track_seed_container = findNode::getClass<TrackSeedContainer>(top_node, m_silicon_track_seed_container_node_name);
		if (!m_silicon_track_seed_container) { missing_node_names.push_back(m_silicon_track_seed_container_node_name); }
	}

	if (m_which_tracks == k_svtx_tracks || m_which_tracks == k_tpc_tracks) {
		m_tpc_track_seed_container = findNode::getClass<TrackSeedContainer>(top_node, m_tpc_track_seed_container_node_name);
		if (!m_tpc_track_seed_container) { missing_node_names.push_back(m_tpc_track_seed_container_node_name); }
	}

	if (m_use_vertex) {
		m_vertex_map = findNode::getClass<SvtxVertexMap>(top_node, m_vertex_map_node_name);
		if (!m_vertex_map) { missing_node_names.push_back(m_vertex_map_node_name); }
	}

	if (!missing_node_names.empty()) {
		std::stringstream what;
		what
			<< PHWHERE
			<< " Couldn't get node(s): ";
		for (auto const& name : missing_node_names) { what << "\n " << name; }
		throw std::runtime_error(what.str());
	}
}

void
WeightedFitter::make_nodes (
	PHCompositeNode* top_node
) {
	PHNodeIterator dst_itr(top_node);
	PHCompositeNode* dst_node = dynamic_cast<PHCompositeNode*>(dst_itr.findFirst("PHCompositeNode", "DST"));
	if (!dst_node) {
		std::stringstream what;
		what
			<< PHWHERE
			<< " Couldn't get node: DST";
		throw std::runtime_error(what.str());
	}

	PHCompositeNode* svtx_node = dynamic_cast<PHCompositeNode*>(dst_itr.findFirst("PHCompositeNode", "SVTX"));
	if (!svtx_node) {
		svtx_node = new PHCompositeNode("SVTX");
		dst_node->addNode(svtx_node);
	}

	m_track_map = findNode::getClass<SvtxTrackMap>(top_node, m_track_map_node_name);
	if (!m_track_map) {
		m_track_map = new SvtxTrackMap_v2;
		PHIODataNode<PHObject>* track_map_node = new PHIODataNode<PHObject>(m_track_map, m_track_map_node_name, "PHObject");
		svtx_node->addNode(track_map_node);
	}

	m_alignment_map = findNode::getClass<SvtxAlignmentStateMap>(top_node, m_alignment_map_node_name);
	if (!m_alignment_map) {
		m_alignment_map = new SvtxAlignmentStateMap_v1;
		PHIODataNode<PHObject>* alignment_map_node = new PHIODataNode<PHObject>(m_alignment_map, m_alignment_map_node_name, "PHObject");
		svtx_node->addNode(alignment_map_node);
	}
}

void
WeightedFitter::make_ntuple (
) {
	if (m_ntuple_file_name.empty()) {
		if (m_file) { m_file->Close(); }
		m_file = nullptr;

		delete m_ntuple;
		m_ntuple = nullptr;

		return;
	}

	m_file = TFile::Open(m_ntuple_file_name.c_str(), "RECREATE");
	if (!m_file) {
		std::stringstream what;
		what
			<< PHWHERE
			<< " Couldn't create file"
			<< m_ntuple_file_name;
		throw std::runtime_error(what.str());
	}

	delete m_ntuple;
	m_ntuple = new TNtuple (
		"ntp", "ntp",
		"event:track:"
		"fitstatus:"
		"nmaps:nintt:ntpc:cluslayer:"
		"clusstave:cluschip:clusstrobe:"
		"clusladderz:clusladderphi:clustimebucket:"
		"clusside:clussector:"
		"cluslx:cluslz:cluselx:cluselz:"
		"clusgx:clusgy:clusgz:"
		"statelx:statelz:stateelx:stateelz:"
		"stategx:stategy:stategz:"
		"dXdx0:dXdy0:dXdtheta:dXdphi:"
		"dYdx0:dYdy0:dYdtheta:dYdphi:"
		"dca_xy:vtxgx:vtxgy:vtxgz:"
	);
	m_ntuple->SetDirectory(m_file);
}

int
WeightedFitter::InitRun (
	PHCompositeNode* top_node
) {
	m_global_position_wrapper.loadNodes(top_node);
	m_global_position_wrapper.set_suppressCrossing(true); // Suppresses print states in the case the crossing is SHRT_MAX

	try {
		make_nodes(top_node);
		make_ntuple();
		return Fun4AllReturnCodes::EVENT_OK;
	} catch (std::exception const& e) {
		std::cout
			<< e.what()
			<< std::endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}
}

int
WeightedFitter::End (
	PHCompositeNode* /*top_node*/
) {
	if (m_ntuple && m_file) {
		m_ntuple->Write();
		m_file->Write();
		m_file->Close();
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

int
WeightedFitter::process_event (
	PHCompositeNode* top_node
) {
	try {
		get_nodes(top_node);
	} catch (std::exception const& e) {
		std::cout
			<< e.what()
			<< std::endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	m_track_id = 0;
	std::string track_seed_container_node_name{};
	TrackSeedContainer* track_seed_container{};
	switch (m_which_tracks) {
	case k_silicon_tracks:
		track_seed_container_node_name = m_silicon_track_seed_container_node_name;
		track_seed_container = m_silicon_track_seed_container;
		break;
	case k_tpc_tracks:
		track_seed_container_node_name = m_tpc_track_seed_container_node_name;
		track_seed_container = m_tpc_track_seed_container;
		break;
	case k_svtx_tracks:
		track_seed_container_node_name = m_svtx_track_seed_container_node_name;
		track_seed_container = m_svtx_track_seed_container;
		break;
	}

	if (!track_seed_container) {
		if (Verbosity()) {
			std::cout
				<< PHWHERE
				<< " TrackSeedContainer " << track_seed_container_node_name << " is null"
				<< std::endl;
		}
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	if (Verbosity()) {
		std::cout
			<< PHWHERE
			<< " TrackSeedContainer " << track_seed_container_node_name
			<< " identify(): ";
		track_seed_container->identify();
		std::cout << std::endl;
	}

	for (auto* track_seed_ptr : *track_seed_container) {
		if (!track_seed_ptr) { continue; }

		if (get_cluster_keys(track_seed_ptr)) {
			if (1 < Verbosity()) { std::cout << PHWHERE << "continuing" << std::endl; }
			continue;
		}
		if (get_points()) {
			if (1 < Verbosity()) { std::cout << PHWHERE << "continuing" << std::endl; }
			continue;
		}
		if (do_fit()) {
			if (1 < Verbosity()) { std::cout << PHWHERE << "continuing" << std::endl; }
			continue;
		}
		if (m_use_vertex && refit_with_vertex()) {
			if (1 < Verbosity()) { std::cout << PHWHERE << "continuing" << std::endl; }
			continue;
		}
		if (add_track()) {
			if (1 < Verbosity()) { std::cout << PHWHERE << "continuing" << std::endl; }
			continue;
		}
		++m_track_id;
	}

	if (Verbosity()) {
		std::cout
			<< PHWHERE
			<< " processed " << m_track_id << " tracks"
			<< std::endl;
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

bool
WeightedFitter::get_cluster_keys (
	TrackSeed* track_seed
) {
	if (1 < Verbosity()) {
		std::cout
			<< PHWHERE
			<< " track seed identify: ";
		track_seed->identify();
		std::cout << std::endl;
	}

	switch (m_which_tracks) {
	case k_silicon_tracks:
		m_silicon_seed = track_seed;
		m_tpc_seed = nullptr;
		break;
	case k_tpc_tracks:
		m_silicon_seed = nullptr;
		m_tpc_seed = track_seed;
		break;
	case k_svtx_tracks:
		m_silicon_seed = m_silicon_track_seed_container->get(track_seed->get_silicon_seed_index());
		m_tpc_seed = m_tpc_track_seed_container->get(track_seed->get_tpc_seed_index());
	}
	m_crossing = m_silicon_seed ? m_silicon_seed->get_crossing() : SHRT_MAX;

	m_cluster_keys.clear();
	for (auto const* seed : {m_silicon_seed, m_tpc_seed}) {
		if (!seed) { continue; }
        std::copy(seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter(m_cluster_keys));
	}

	return false;
}

bool
WeightedFitter::get_points (
) {
	m_num_mvtx = 0;
	m_num_intt = 0;
	m_num_tpc = 0;

	delete m_weighted_track;
	m_weighted_track = m_track_factory();

	// Assign all clusters (e.g., Si clusters) the same side
	// Choose mode as the value to set all points to
	// Print variations with high enough verbosity
	std::map<int, int> sides;

	for (auto const& cluster_key : m_cluster_keys) {
		TrkrCluster* cluster = m_trkr_cluster_container->findCluster(cluster_key);
		if (!cluster) { continue; }

		Surface const surf = m_geometry->maps().getSurface(cluster_key, cluster);
		if (!surf) { continue; }
	
		auto local_to_global_transform = surf->transform(m_geometry->geometry().getGeoContext()); // in mm
		local_to_global_transform.translation() /= Acts::UnitConstants::cm; // converted to cm
		Eigen::Vector3d local_pos = Eigen::Vector3d { cluster->getLocalX(), cluster->getLocalY(), 0.0 }; // in cm

		ClusterFitPoint point;
		point.cluster_key = cluster_key;
		point.cluster_position = (local_to_global_transform * local_pos);
		point.cluster_errors = { cluster->getRPhiError(), cluster->getZError() };
		point.sensor_local_to_global_transform = local_to_global_transform;

		// Modify position, error using helper classes
		// Note that these calls are effectively NOOPs for non-TPC points
		point.cluster_position = m_global_position_wrapper.getGlobalPositionDistortionCorrected(cluster_key, cluster, m_crossing);

		// note the "clusterRadius" argument (supplied as 0.0) is unused
		auto rphi_z_error_pair = ClusterErrorPara::get_clusterv5_modified_error(cluster, 0.0, cluster_key);
		point.cluster_errors = Eigen::Vector2d {std::sqrt(rphi_z_error_pair.first), std::sqrt(rphi_z_error_pair.second) };

		m_weighted_track->push_back(point);

		TrkrDefs::TrkrId trkr_id = static_cast<TrkrDefs::TrkrId>(TrkrDefs::getTrkrId(cluster_key));
		switch (trkr_id) {
			case TrkrDefs::mvtxId:
				++m_num_mvtx;
				break;
			case TrkrDefs::inttId:
				++m_num_intt;
				break;
			case TrkrDefs::tpcId:
				++m_num_tpc;
				++sides[TpcDefs::getSide(cluster_key)];
				break;
			default:
				continue;
		}
	}

	if (m_reassign_sides && !sides.empty()) {
		m_side = sides[0] < sides[1] ? 1 : 0;
	}

	return false;
}

bool
WeightedFitter::do_fit (
) {
	try {
		m_weighted_track->configure_minimizer(*m_minimizer);
	} catch (std::exception const& e) {
		if (Verbosity()) {
			std::cout
				<< PHWHERE << "\n"
				<< e.what()
				<< std::endl;
		}
		return true;
	}

	std::function<double(double const*)> lambda = [&](double const* params) { return m_weighted_track->get_squared_error(params); };
	ROOT::Math::Functor functor(lambda, m_weighted_track->get_n_parameters());
	m_minimizer->SetFunction(functor);

	bool fit_succeeded = m_minimizer->Minimize();
	if (1 < Verbosity()) {
		std::cout
			<< PHWHERE
			<< " fit " << (fit_succeeded ? "succeeded" : "failed")
			<< " status: " << m_minimizer->Status()
			<< std::endl;
	}

	double const* params = m_minimizer->X();
	m_weighted_track->set_parameters(params);

	return !fit_succeeded;
}

bool
WeightedFitter::refit_with_vertex (
) {
	if (!m_vertex_map) {
		std::stringstream what;
		what
			<< PHWHERE
			<< " Couldn't get node: "
			<< m_vertex_map_node_name;
		throw std::runtime_error(what.str());
	}

	if (m_vertex_map->empty()) {
		if (1 < Verbosity()) {
			std::cout
				<< PHWHERE
				<< " vertex map is empty"
				<< std::endl;
		}
		return true;
	}

	SvtxVertex* closest_vertex{};
	double min_dca = std::numeric_limits<double>::max();
	for (auto const& [vertex_id, svtx_vertex] : *m_vertex_map) {
		if (!svtx_vertex) { continue; }

		Eigen::Vector3d vertex_pos {
			svtx_vertex->get_x(),
			svtx_vertex->get_y(),
			svtx_vertex->get_z(),
		};

		double dca = (m_weighted_track->get_pca(vertex_pos) - vertex_pos).norm();
		if (dca < min_dca) {
			min_dca = dca;
			closest_vertex = svtx_vertex;
		}
	}

	if (!closest_vertex) {
		std::stringstream what;
		what
			<< PHWHERE
			<< "no vertex after loop";
		throw std::runtime_error(what.str());
	}

	m_weighted_track->use_vertex();
	for (int i = 0; i < 3; ++i) {
		m_weighted_track->vertex_position(i) = closest_vertex->get_position(i);
		for (int j = 0; j < 3; ++j) {
			m_weighted_track->vertex_covariance(i, j) = closest_vertex->get_error(i, j);
		}
	}

	if (1 < Verbosity()) {
		std::cout
			<< PHWHERE
			<< " vertex: " << m_weighted_track->get_vertex_position().transpose()
			<< " track pca: " << m_weighted_track->get_pca(m_weighted_track->get_vertex_position()).transpose()
			<< " dca: " << min_dca
			<< std::endl;
	}

	return do_fit();
}

bool
WeightedFitter::add_track (
) {
	double const* params = m_minimizer->X();
	m_weighted_track->set_parameters(params);
	Eigen::Vector3d slope = m_weighted_track->get_slope_at_path_length(0); // Straight line tracks have uniform slope, we can pass any value for path_length here

	SvtxTrack_v4 fitted_track;
	fitted_track.set_id(m_track_id);
	fitted_track.set_silicon_seed(m_silicon_seed);
	fitted_track.set_tpc_seed(m_tpc_seed);
	fitted_track.set_crossing(m_crossing);
	fitted_track.set_charge(m_tpc_seed ? m_tpc_seed->get_charge() : 0);
	fitted_track.set_x(params[0]);
	fitted_track.set_y(params[1]);
	fitted_track.set_z(0.0);
	fitted_track.set_px(slope(0));
	fitted_track.set_py(slope(1));
	fitted_track.set_pz(slope(2));

	SvtxAlignmentStateMap::StateVec alignment_states;
	for (auto const& point : *m_weighted_track) {
		Acts::Vector3 intersection = m_weighted_track->get_intersection(point.sensor_local_to_global_transform);

		SvtxTrackState_v1 svtx_track_state(intersection.norm());
		svtx_track_state.set_x(intersection(0));
		svtx_track_state.set_y(intersection(1));
		svtx_track_state.set_z(intersection(2));
		svtx_track_state.set_px(slope(0));
		svtx_track_state.set_py(slope(1));
		svtx_track_state.set_pz(slope(2));
		fitted_track.insert_state(&svtx_track_state);

		auto alignment_state = std::make_unique<SvtxAlignmentState_v1>();
		SvtxAlignmentState::GlobalMatrix global_derivative_matrix = SvtxAlignmentState::GlobalMatrix::Zero();
		SvtxAlignmentState::LocalMatrix local_derivative_matrix = SvtxAlignmentState::LocalMatrix::Zero();

		Eigen::Vector2d residual = point.get_residuals(intersection);
		Eigen::Matrix<double, 2, 3> projection = m_weighted_track->get_projection(point.sensor_local_to_global_transform);

		// PARTIAL derivatives of global alignment parameters
		// (the projection vectors take care of the implicit intersection dependency)
		std::array<Eigen::Vector3d, 6> global_derivatives {
			Eigen::Vector3d {1.0, 0.0, 0.0}.cross(intersection - point.sensor_local_to_global_transform.translation()),
			Eigen::Vector3d {0.0, 1.0, 0.0}.cross(intersection - point.sensor_local_to_global_transform.translation()),
			Eigen::Vector3d {0.0, 0.0, 1.0}.cross(intersection - point.sensor_local_to_global_transform.translation()),
			Eigen::Vector3d {1.0, 0.0, 0.0},
			Eigen::Vector3d {0.0, 1.0, 0.0},
			Eigen::Vector3d {0.0, 0.0, 1.0},
		};

		for (int i = 0; i < 4; ++i) {
			local_derivative_matrix.col(i) = m_weighted_track->get_residual_derivative(i, point.sensor_local_to_global_transform);
		}

		for (int i = 0; i < 6; ++i) {
			global_derivative_matrix.col(i) = projection * global_derivatives[i];
		}

		alignment_state->set_cluster_key(point.cluster_key);
		alignment_state->set_residual(residual);
		alignment_state->set_local_derivative_matrix(local_derivative_matrix);
		alignment_state->set_global_derivative_matrix(global_derivative_matrix);
		alignment_states.push_back(alignment_state.release());

		if (m_ntuple) {
			if (1 < Verbosity()) {
				std::cout
					<< " Filling ntuple"
					<< std::endl;
			}

			float layer = TrkrDefs::getLayer(point.cluster_key);
			float stave{-1};
			float chip{-1};
			float strobe{-1};
			float ladder_z{-1};
			float ladder_phi{-1};
			float time_bucket{-1};
			float side{-1};
			float sector{-1};
			TrkrDefs::TrkrId trkr_id = static_cast<TrkrDefs::TrkrId>(TrkrDefs::getTrkrId(point.cluster_key));
			switch (trkr_id) {
				case TrkrDefs::mvtxId:
					stave = MvtxDefs::getStaveId(point.cluster_key);
					chip = MvtxDefs::getChipId(point.cluster_key);
					strobe = MvtxDefs::getStrobeId(point.cluster_key);
					break;
				case TrkrDefs::inttId:
					ladder_z = InttDefs::getLadderZId(point.cluster_key);
					ladder_phi = InttDefs::getLadderPhiId(point.cluster_key);
					time_bucket = InttDefs::getTimeBucketId(point.cluster_key);
					break;
				case TrkrDefs::tpcId:
					side = TpcDefs::getSide(point.cluster_key);
					sector = TpcDefs::getSectorId(point.cluster_key);
					break;
				default:
					continue;
			}
			if (m_reassign_sides) { side = m_side; }

			Eigen::Affine3d global_to_local_transform = point.sensor_local_to_global_transform.inverse();
			Eigen::Vector3d cluster_local_position =  global_to_local_transform * point.cluster_position;
			Eigen::Vector3d intersection_local_position = global_to_local_transform * intersection;

			float dca_xy = std::numeric_limits<float>::quiet_NaN();
			Eigen::Vector3f vertex = {
				std::numeric_limits<float>::quiet_NaN(),
				std::numeric_limits<float>::quiet_NaN(),
				std::numeric_limits<float>::quiet_NaN(),
			};

			if (m_use_vertex) {
				vertex(0) = (float)m_weighted_track->vertex_position(0);
				vertex(1) = (float)m_weighted_track->vertex_position(1);
				vertex(2) = (float)m_weighted_track->vertex_position(2);

				// displacement between the associated vertex and the point of closest approach on the track
				Eigen::Vector3d delta = m_weighted_track->get_pca(m_weighted_track->get_vertex_position()) - m_weighted_track->get_vertex_position();
				dca_xy = std::sqrt(sqr(delta(0)) + sqr(delta(1)));
			}

			float ntp_data[] = {
				(float)Fun4AllServer::instance()->EventNumber(), (float)m_track_id,
				(float)m_minimizer->Status(),
				(float)m_num_mvtx, (float)m_num_intt, (float)m_num_tpc, layer,
				stave, chip, strobe,
				ladder_z, ladder_phi, time_bucket,
				side, sector,
				(float)cluster_local_position(0), (float)cluster_local_position(1), (float)point.cluster_errors(0), (float)point.cluster_errors(1),
				(float)point.cluster_position(0), (float)point.cluster_position(1), (float)point.cluster_position(2),
				(float)intersection_local_position(0), (float)intersection_local_position(1), (float)point.cluster_errors(0), (float)point.cluster_errors(1),
				(float)intersection(0), (float)intersection(1), (float)intersection(2),
				(float)local_derivative_matrix(0, 0), (float)local_derivative_matrix(0, 1), (float)local_derivative_matrix(0, 2), (float)local_derivative_matrix(0, 3),
				(float)local_derivative_matrix(1, 0), (float)local_derivative_matrix(1, 1), (float)local_derivative_matrix(1, 2), (float)local_derivative_matrix(1, 3),
				dca_xy, vertex(0), vertex(1), vertex(2),
			};
			m_ntuple->Fill(ntp_data);
		}

		if (2 < Verbosity()) {
			std::cout
				<< "\tcluster_key:  " << point.cluster_key << "\n"
				<< "\tintersection: " << intersection.transpose() << "\n"
				<< "\tpos:          " << point.cluster_position.transpose() << "\n"
				<< "\to:            " << point.sensor_local_to_global_transform.translation().transpose() << "\n"
				<< "\tx:            " << point.sensor_local_to_global_transform.rotation().col(0).transpose() << "\n"
				<< "\ty:            " << point.sensor_local_to_global_transform.rotation().col(1).transpose() << "\n"
				<< "\tz:            " << point.sensor_local_to_global_transform.rotation().col(2).transpose() << "\n"
				<< "\tproj_x:       " << projection.row(0) << "\n"
				<< "\tproj_y:       " << projection.row(1) << "\n"
				<< std::endl;
		}
	}

	// cuts...
	// return true;

	m_track_map->insertWithKey(&fitted_track, m_track_id);
	m_alignment_map->insertWithKey(m_track_id, alignment_states);

	return false;
}

