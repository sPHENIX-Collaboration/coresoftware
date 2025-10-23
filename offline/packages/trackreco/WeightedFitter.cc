#include "WeightedFitter.h"

#include <trackbase/ActsGeometry.h>
#include <trackbase/ActsSurfaceMaps.h>
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
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeed_v1.h>
#include <trackbase_historic/TrackSeed_v2.h>
#include <trackbase_historic/SvtxTrackSeed_v1.h>
#include <trackbase_historic/SvtxTrackSeed_v2.h>
#include <trackbase_historic/TrackSeedContainer.h>

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
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

namespace {
	template <typename T>
	T sqr (T const& t) { return t * t; }
}

double
WeightedFitter::FitErrorCalculator::operator() (
	double const* params
) const {
	set_parameters(params);

	double error{0};

	if (m_vertex.use) {
		Eigen::Vector3d displacement = get_pca(m_vertex.pos) - m_vertex.pos;
		error += displacement.transpose() * m_vertex.cov * displacement;
	}

	for (auto const& point : m_points) {
		// in-plane displacement (in global coordinates)
		Eigen::Vector3d displacement = get_intersection(point) - point.pos;

		// Square error
		error += sqr(displacement.dot(point.x) / point.sigma_x);
		error += sqr(displacement.dot(point.y) / point.sigma_y);
	}
	return error;
}

void
WeightedFitter::FitErrorCalculator::set_parameters (
	double const* params
) const {
	m_cos_theta = std::cos(params[2]);
	m_sin_theta = std::sin(params[2]);
	m_cos_phi = std::cos(params[3]);
	m_sin_phi = std::sin(params[3]);

	m_intercept = Eigen::Vector3d {
		params[0], // x-component of xy plane intercept
		params[1], // y-component of xy plane intercept
		0,
	};

	m_slope = Eigen::Vector3d {
		m_sin_theta * m_cos_phi,
		m_sin_theta * m_sin_phi,
		m_cos_theta
	};
}

Eigen::Vector3d
WeightedFitter::FitErrorCalculator::get_intersection (
	WeightedFitter::ClusterFitPoint const& point
) const {
	// solve s
	// (m_slope * s + m_intercept - point.pos).dot(point.z) = 0
	double s = (point.pos - m_intercept).dot(point.z) / m_slope.dot(point.z);
	return m_slope * s + m_intercept;
}

Eigen::Vector3d
WeightedFitter::FitErrorCalculator::get_pca (
	Eigen::Vector3d const& point
) const {
	// solve s
	// (m_slope * s + m_intercept - point).dot(m_slope) = 0
	double s = (point - m_intercept).dot(m_slope); // m_slope is a unit vector
	return m_slope * s + m_intercept;
}

WeightedFitter::WeightedFitter (
	std::string const& name
) : SubsysReco (
	name
) {
	// Setup Minimizer
	m_minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2");
	m_minimizer->SetMaxFunctionCalls(1E8); // For Minuit/Minuit2
	m_minimizer->SetMaxIterations(1E6); // For GSL Minimizers--probably moot to call
	m_minimizer->SetTolerance(1E-8);
	m_minimizer->SetPrintLevel(0);

	m_fit_error_calculator = new FitErrorCalculator;
}

WeightedFitter::~WeightedFitter (
) {
	delete m_fit_error_calculator;
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

	if (missing_node_names.size()) {
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
		if (m_file) m_file->Close();
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
		"Ndx0:Ndy0:Ndtheta:Ndphi:"
		"dca:vtxgx:vtxgy:vtxgz:"
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
	PHCompositeNode*
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

	for (auto track_seed_ptr : *track_seed_container) {
		if (!track_seed_ptr) continue;

		if (get_cluster_keys(track_seed_ptr)) continue;
		if (get_points()) continue;
		if (do_fit()) continue;
		if (m_use_vertex && refit_with_vertex()) continue;
		if (add_track()) continue;
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
		if (!seed) continue;
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

	m_fit_error_calculator->m_vertex.reset();
	m_fit_error_calculator->m_points.clear();

	// Assign all clusters (e.g., Si clusters) the same side
	// Choose mode as the value to set all points to
	// Print variations with high enough verbosity
	std::map<int, int> sides;

	for (auto const& cluster_key : m_cluster_keys) {
		TrkrCluster* cluster = m_trkr_cluster_container->findCluster(cluster_key);
		if (!cluster) continue;

		Surface const surf = m_geometry->maps().getSurface(cluster_key, cluster);
		if (!surf) continue;
	
		auto local_to_global_transform = surf->transform(m_geometry->geometry().getGeoContext()); // in mm

		Eigen::Vector3d local_pos = Eigen::Vector3d {
			cluster->getLocalX(),
			cluster->getLocalY(),
			0.0
		} * Acts::UnitConstants::cm; // converted to mm

		ClusterFitPoint point {
			.cluster_key = cluster_key,
			.layer = TrkrDefs::getLayer(point.cluster_key),
			.pos = (local_to_global_transform * local_pos) / Acts::UnitConstants::cm,
			.o = local_to_global_transform.translation() / Acts::UnitConstants::cm,
			.x = local_to_global_transform.rotation().col(0), 
			.y = local_to_global_transform.rotation().col(1),
			.z = local_to_global_transform.rotation().col(2),
			.sigma_x = cluster->getRPhiError(),
			.sigma_y = cluster->getZError(),
		};

		// Modify position, error using helper classes
		// Note that these calls are effectively NOOPs for non-TPC points
		point.pos = m_global_position_wrapper.getGlobalPositionDistortionCorrected(cluster_key, cluster, m_crossing);

		// note the "clusterRadius" argument (supplied as 0.0) is unused
		auto rphi_z_error_pair = m_cluster_error_para.get_clusterv5_modified_error(cluster, 0.0, cluster_key);
		point.sigma_x = sqrt(rphi_z_error_pair.first);
		point.sigma_y = sqrt(rphi_z_error_pair.second);

		TrkrDefs::TrkrId trkr_id = static_cast<TrkrDefs::TrkrId>(TrkrDefs::getTrkrId(cluster_key));
		switch (trkr_id) {
			case TrkrDefs::mvtxId:
				point.stave = MvtxDefs::getStaveId(cluster_key);
				point.chip = MvtxDefs::getChipId(cluster_key);
				point.strobe = MvtxDefs::getStrobeId(cluster_key);
				++m_num_mvtx;
				break;
			case TrkrDefs::inttId:
				point.ladder_z = InttDefs::getLadderZId(cluster_key);
				point.ladder_phi = InttDefs::getLadderPhiId(cluster_key);
				point.time_bucket = InttDefs::getTimeBucketId(cluster_key);
				++m_num_intt;
				break;
			case TrkrDefs::tpcId: {
				point.side = TpcDefs::getSide(cluster_key);
				point.sector = TpcDefs::getSectorId(cluster_key);
				++sides[point.side];
				++m_num_tpc;
				break;
			}
			default:
				continue;
		}

		m_fit_error_calculator->m_points.push_back (point);
	}

	if (m_reassign_sides && !sides.empty()) {
		int side = sides[0] < sides[1] ? 1 : 0;
		for (auto& point : m_fit_error_calculator->m_points) {
			if (Verbosity() && (point.side != -1) && (point.side != side)) {
				std::cout
					<< PHWHERE
					<< " track with mixed TPC sides"
					<< " (at layer " << point.layer << ", mode is " << side << ", got " << point.side << ")"
					<< std::endl;
			}
			point.side = side;
		}
	}

	return false;
}

bool
WeightedFitter::do_fit (
) {
	if (m_fit_error_calculator->m_points.size() < 2) {
		if (1 < Verbosity()) {
			std::cout
				<< PHWHERE
				<< " too few points to perform fit"
				<< std::endl;
		}
		return true;
	}

	// initial slope
	Eigen::Vector3d slope = (m_fit_error_calculator->m_points[1].pos - m_fit_error_calculator->m_points[0].pos).normalized();
	double phi = std::atan2(slope(1), slope(0));
	double theta = std::acos(slope(2));

	// initlal intercept
	// Solve for s (slope * s + m_points[0].pos).dot({0, 0, 1}) = 0
	double s = -m_fit_error_calculator->m_points[0].pos(2) / slope(2);
	double x = (slope * s + m_fit_error_calculator->m_points[0].pos)(0);
	double y = (slope * s + m_fit_error_calculator->m_points[0].pos)(1);

	m_minimizer->SetVariable(0, "x",     x,     1E-6);
	m_minimizer->SetVariable(1, "y",     y,     1E-6);
	m_minimizer->SetVariable(2, "theta", theta, 1E-6);
	m_minimizer->SetVariable(3, "phi",   phi,   1E-6);

	// Bounds for angular parameters
	m_minimizer->SetVariableLimits(2,          0,     3.1416);
	m_minimizer->SetVariableLimits(3, phi - 1.58, phi + 1.58);

	ROOT::Math::Functor get_fit_error(*m_fit_error_calculator, &FitErrorCalculator::operator(), 4);
	m_minimizer->SetFunction(get_fit_error);

	bool fit_succeeded = m_minimizer->Minimize();
	if (1 < Verbosity()) {
		std::cout
			<< PHWHERE
			<< " fit " << (fit_succeeded ? "succeeded" : "failed")
			<< " status: " << m_minimizer->Status()
			<< std::endl;
	}

	double const* params = m_minimizer->X();
	m_fit_error_calculator->set_parameters(params);

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
		if (!svtx_vertex) continue;

		Eigen::Vector3d vertex_pos {
			svtx_vertex->get_x(),
			svtx_vertex->get_y(),
			svtx_vertex->get_z(),
		};

		double dca = (m_fit_error_calculator->get_pca(vertex_pos) - vertex_pos).norm();
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

	m_fit_error_calculator->m_vertex.use = true;
	m_fit_error_calculator->m_vertex.pos = Eigen::Vector3d {
		closest_vertex->get_x(),
		closest_vertex->get_y(),
		closest_vertex->get_z(),
	};
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			m_fit_error_calculator->m_vertex.cov(i, j) = closest_vertex->get_error(i, j);
		}
	}

	if (1 < Verbosity()) {
		std::cout
			<< PHWHERE
			<< " vertex: " << m_fit_error_calculator->m_vertex.pos.transpose()
			<< " track pca: " << m_fit_error_calculator->get_pca(m_fit_error_calculator->m_vertex.pos).transpose()
			<< " dca: " << min_dca
			<< std::endl;
	}

	return do_fit();
}

bool
WeightedFitter::add_track (
) {
	double const* params = m_minimizer->X();
	m_fit_error_calculator->set_parameters(params);

	// xy-intercept
	Acts::Vector3 b {
		params[0], // x-component of xy plane intercept
		params[1], // y-component of xy plane intercept
		0
	};

	// slope (spherical-polar parameterized)
	double cos_theta = std::cos(params[2]), sin_theta = std::sin(params[2]);
	double cos_phi = std::cos(params[3]), sin_phi = std::sin(params[3]);
	Acts::Vector3 m {
		sin_theta * cos_phi,
		sin_theta * sin_phi,
		cos_theta
	};

	SvtxTrack_v4 fitted_track;
	fitted_track.set_id(m_track_id);
	fitted_track.set_silicon_seed(m_silicon_seed);
	fitted_track.set_tpc_seed(m_tpc_seed);
	fitted_track.set_crossing(m_crossing);
	fitted_track.set_charge(m_tpc_seed ? m_tpc_seed->get_charge() : 0);
	fitted_track.set_x(params[0]);
	fitted_track.set_y(params[1]);
	fitted_track.set_z(0.0);
	fitted_track.set_px(m(0));
	fitted_track.set_py(m(1));
	fitted_track.set_pz(m(2));

	SvtxAlignmentStateMap::StateVec alignment_states;
	for (auto const& point : m_fit_error_calculator->m_points) {
		Acts::Vector3 intersection = m_fit_error_calculator->get_intersection(point);
		double s = (point.pos - b).dot(point.z) / m.dot(point.z);

		SvtxTrackState_v1 svtx_track_state(intersection.norm());
		svtx_track_state.set_x(intersection(0));
		svtx_track_state.set_y(intersection(1));
		svtx_track_state.set_z(intersection(2));
		svtx_track_state.set_px(m(0));
		svtx_track_state.set_py(m(1));
		svtx_track_state.set_pz(m(2));
		fitted_track.insert_state(&svtx_track_state);

		// Move this implementation into the fitter helper class eventually
		auto alignment_state = std::make_unique<SvtxAlignmentState_v1>();
		SvtxAlignmentState::GlobalMatrix global_derivative_matrix = SvtxAlignmentState::GlobalMatrix::Zero();
		SvtxAlignmentState::LocalMatrix local_derivative_matrix = SvtxAlignmentState::LocalMatrix::Zero();

		// TODO: move these to the helper class implementation

		Acts::Vector2 residual {
			point.x.dot(intersection - point.pos),
			point.y.dot(intersection - point.pos),
		};

		// Follows from the Atlas paper (Global \Chi^{2} approach to the Alignment of the ATLAS Silicon Tracking Detectors)
		std::array<Acts::Vector3, 2> proj = {
			point.x - (m.dot(point.x) / m.dot(point.z)) * point.z,
			point.y - (m.dot(point.y) / m.dot(point.z)) * point.z,
		};

		// PARTIAL derivatives of track parameters
		// (the projection vectors take care of the implicit intersection dependency)
		std::array<Eigen::Vector3d, 4> local_derivatives {
			Eigen::Vector3d {1.0, 0.0, 0.0},
			Eigen::Vector3d {0.0, 1.0, 0.0},
			Eigen::Vector3d { cos_theta * cos_phi, cos_theta * sin_phi, -sin_theta} * s,
			Eigen::Vector3d {-sin_theta * sin_phi, sin_theta * cos_phi, 0.0} * s,
		};

		// PARTIAL derivatives of alignment parameters
		// (the projection vectors take care of the implicit intersection dependency)
		std::array<Eigen::Vector3d, 6> global_derivatives {
			Eigen::Vector3d {1.0, 0.0, 0.0}.cross(intersection - point.o),
			Eigen::Vector3d {0.0, 1.0, 0.0}.cross(intersection - point.o),
			Eigen::Vector3d {0.0, 0.0, 1.0}.cross(intersection - point.o),
			Eigen::Vector3d {1.0, 0.0, 0.0},
			Eigen::Vector3d {0.0, 1.0, 0.0},
			Eigen::Vector3d {0.0, 0.0, 1.0},
		};

		for (int i = 0; i < 4; ++i) {
			local_derivative_matrix(0, i) = proj[0].dot(local_derivatives[i]); // X residual partial derivative
			local_derivative_matrix(1, i) = proj[1].dot(local_derivatives[i]); // Y residual partial derivative
		}

		for (int i = 0; i < 6; ++i) {
			global_derivative_matrix(0, i) = proj[0].dot(global_derivatives[i]); // X residual partial derivative
			global_derivative_matrix(1, i) = proj[1].dot(global_derivatives[i]); // Y residual partial derivative
		}

		std::array<double, 4> numeric_partial_derivatives{};
		for (int i = 0; i < 4; ++i) {
			static double const epsilon = 1.0E-8;

			double temp_params[4]{};
			for (int j = 0; j < 4; ++j) temp_params[j] = params[j];

			temp_params[i] = params[i] + epsilon;
			m_fit_error_calculator->set_parameters(temp_params);
			Eigen::Vector3d temp = m_fit_error_calculator->get_intersection(point) - point.pos;
			numeric_partial_derivatives[i] += sqr(point.x.dot(temp) / point.sigma_x);
			numeric_partial_derivatives[i] += sqr(point.y.dot(temp) / point.sigma_y);

			temp_params[i] = params[i] - epsilon;
			m_fit_error_calculator->set_parameters(temp_params);
			temp = m_fit_error_calculator->get_intersection(point) - point.pos;
			numeric_partial_derivatives[i] -= sqr(point.x.dot(temp) / point.sigma_x);
			numeric_partial_derivatives[i] -= sqr(point.y.dot(temp) / point.sigma_y);

			numeric_partial_derivatives[i] /= 2.0 * epsilon;
		}
		m_fit_error_calculator->set_parameters(params);

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

			float dca = std::numeric_limits<float>::quiet_NaN();
			Eigen::Vector3f vertex = {
				std::numeric_limits<float>::quiet_NaN(),
				std::numeric_limits<float>::quiet_NaN(),
				std::numeric_limits<float>::quiet_NaN(),
			};

			if (m_fit_error_calculator->m_vertex.use) {
				vertex(0) = (float)m_fit_error_calculator->m_vertex.pos(0);
				vertex(1) = (float)m_fit_error_calculator->m_vertex.pos(1);
				vertex(2) = (float)m_fit_error_calculator->m_vertex.pos(2);
				dca = (float)(m_fit_error_calculator->get_pca(m_fit_error_calculator->m_vertex.pos) - m_fit_error_calculator->m_vertex.pos).norm();
			}

			float ntp_data[] = {
				(float)Fun4AllServer::instance()->EventNumber(), (float)m_track_id,
				(float)m_minimizer->Status(),
				(float)m_num_mvtx, (float)m_num_intt, (float)m_num_tpc, (float)point.layer,
				(float)point.stave, (float)point.chip, (float)point.strobe,
				(float)point.ladder_z, (float)point.ladder_phi, (float)point.time_bucket,
				(float)point.side, (float)point.sector,
				(float)point.x.dot(point.pos - point.o), (float)point.y.dot(point.pos - point.o), (float)point.sigma_x, (float)point.sigma_y,
				(float)point.pos(0), (float)point.pos(1), (float)point.pos(2),
				(float)point.x.dot(intersection - point.o), (float)point.y.dot(intersection - point.o), (float)point.sigma_x, (float)point.sigma_y,
				(float)intersection(0), (float)intersection(1), (float)intersection(2),
				(float)local_derivative_matrix(0, 0), (float)local_derivative_matrix(0, 1), (float)local_derivative_matrix(0, 2), (float)local_derivative_matrix(0, 3),
				(float)local_derivative_matrix(1, 0), (float)local_derivative_matrix(1, 1), (float)local_derivative_matrix(1, 2), (float)local_derivative_matrix(1, 3),
				(float)numeric_partial_derivatives[0],  (float)numeric_partial_derivatives[1], (float)numeric_partial_derivatives[2], (float)numeric_partial_derivatives[3], 
				dca, vertex(0), vertex(1), vertex(2),
			};
			m_ntuple->Fill(ntp_data);
		}

		if (2 < Verbosity()) {
			std::cout
				<< "\tcluster_key:  " << point.cluster_key << "\n"
				<< "\tintersection: " << intersection.transpose() << "\n"
				<< "\tpos:          " << point.pos.transpose() << "\n"
				<< "\to:            " << point.o.transpose() << "\n"
				<< "\tx:            " << point.x.transpose() << "\n"
				<< "\ty:            " << point.y.transpose() << "\n"
				<< "\tz:            " << point.z.transpose() << "\n"
				<< "\tproj_x:       " << proj[0].transpose() << "\n"
				<< "\tproj_y:       " << proj[1].transpose() << "\n"
				<< std::endl;
		}
	}

	// cuts...
	// return true;

	m_track_map->insertWithKey(&fitted_track, m_track_id);
	m_alignment_map->insertWithKey(m_track_id, alignment_states);

	return false;
}

