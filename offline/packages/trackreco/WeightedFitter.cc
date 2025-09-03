#include "WeightedFitter.h"

#include <trackbase/ActsGeometry.h>
#include <trackbase/ActsSurfaceMaps.h>
// #include <trackbase/InttDefs.h>
// #include <trackbase/MvtxDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainerv4.h>

#include <trackbase_historic/SvtxAlignmentState_v1.h>
#include <trackbase_historic/SvtxAlignmentStateMap_v1.h>
#include <trackbase_historic/SvtxTrack_v4.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/TrackSeed_v2.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/phool.h> // PHWHERE

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>

#include <TAxis.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TLine.h>
#include <TPad.h>

#include <array>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace {
	template <typename T>
	T sqr(T const& t) { return t * t; }
}

WeightedFitter::FitErrorCalculator::FitErrorCalculator (
	std::vector<WeightedFitter::ClusterFitPoint> const& points
) : m_points{points} {}

double
WeightedFitter::FitErrorCalculator::operator() (
	double const* params
) const {
	// xy-intercept
	Eigen::Vector3d b {
		params[0], // x-component of xy plane intercept
		params[1], // y-component of xy plane intercept
		0,
	};

	// compute trig functions (reused often)
	double cos_theta = std::cos(params[2]), sin_theta = std::sin(params[2]);
	double cos_phi = std::cos(params[3]), sin_phi = std::sin(params[3]);

	// slope (spherical-polar parameterized)
	Eigen::Vector3d m {
		sin_theta * cos_phi,
		sin_theta * sin_phi,
		cos_theta
	};

	double error{0};
	for (auto const& point : m_points) {
		// Get intersection, solve s
		// (m * s + b - point.pos).dot(point.z) = 0
		double s = (point.pos - b).dot(point.z) / m.dot(point.z);

		// in-plane displacement (in global coordinates)
		Eigen::Vector3d displacement = m * s + b - point.pos;

		// Square error
		error += sqr(displacement.dot(point.x) / point.sigma_x);
		error += sqr(displacement.dot(point.y) / point.sigma_y);
	}

	return error;
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
	m_minimizer->SetTolerance(1E-4);
	m_minimizer->SetPrintLevel(0);
}

WeightedFitter::~WeightedFitter (
) {
	delete m_minimizer;
}

void
WeightedFitter::get_nodes (
	PHCompositeNode* top_node
) {
	m_geometry = findNode::getClass<ActsGeometry>(top_node, m_geometry_node_name);
	if (!m_geometry) {
		std::stringstream what;
		what
			<< PHWHERE
			<< " Couldn't get node: "
			<< m_geometry_node_name;
		throw std::runtime_error(what.str());
	}

	m_trkr_cluster_container = findNode::getClass<TrkrClusterContainer>(top_node, m_trkr_cluster_container_node_name);
	if (!m_trkr_cluster_container) {
		std::stringstream what;
		what
			<< PHWHERE
			<< " Couldn't get node: "
			<< m_trkr_cluster_container_node_name;
		throw std::runtime_error(what.str());
	}

	m_track_seed_container = findNode::getClass<TrackSeedContainer>(top_node, m_track_seed_container_node_name);
	if (!m_track_seed_container) {
		std::stringstream what;
		what
			<< PHWHERE
			<< " Couldn't get node: "
			<< m_track_seed_container_node_name;
		throw std::runtime_error(what.str());
	}

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

int
WeightedFitter::InitRun (
	PHCompositeNode* top_node
) {
	try {
		get_nodes(top_node);
		return Fun4AllReturnCodes::EVENT_OK;
	} catch (std::exception const& e) {
		std::cout
			<< e.what()
			<< std::endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}
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
	for (
		auto track_seed_itr = m_track_seed_container->begin();
		track_seed_itr != m_track_seed_container->end();
		++track_seed_itr
	) {
		if (!*track_seed_itr) continue;

		if (get_points(*track_seed_itr)) continue;
		if (do_fit()) continue;
		if (add_track(*track_seed_itr)) continue;
		++m_track_id;

		// draw(); // leaks memory in current implementation
		// break;
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
WeightedFitter::get_points (
	TrackSeed const* track_seed
) {
	int n_mvtx{0}, n_intt{0}, n_tpc{0};
	m_points.clear();	

	if (Verbosity()) {
		std::cout
			<< PHWHERE
			<< " tracklet clusters: " << track_seed->size_cluster_keys()
			<< std::endl;
	}

	for (
		auto cluster_key_itr = track_seed->begin_cluster_keys();
		cluster_key_itr != track_seed->end_cluster_keys();
		++cluster_key_itr
	) {
		TrkrCluster* cluster = m_trkr_cluster_container->findCluster(*cluster_key_itr);
		if (!cluster) continue;

		// Ignore crazy clusters
		if (1.0 < cluster->getRPhiError()) continue;
		if (1.0 < cluster->getZError()) continue;

		Surface const surf = m_geometry->maps().getSurface(*cluster_key_itr, cluster);
		if (!surf) continue;
	
		auto local_to_global_transform = surf->transform(m_geometry->geometry().getGeoContext()); // in mm

		Eigen::Vector3d local_pos = Eigen::Vector3d {
			cluster->getLocalX(),
			cluster->getLocalY(),
			0.0
		} * Acts::UnitConstants::cm; // convernted to mm

		ClusterFitPoint point {
			.cluster_key = *cluster_key_itr,
			.pos = (local_to_global_transform * local_pos) / Acts::UnitConstants::cm,
			.o = local_to_global_transform.translation() / Acts::UnitConstants::cm,
			.x = local_to_global_transform.rotation().col(0), 
			.y = local_to_global_transform.rotation().col(1),
			.z = -local_to_global_transform.rotation().col(2), // deliberate but idk why it's a left-handed coordinate system
			.sigma_x = cluster->getRPhiError(),
			.sigma_y = cluster->getZError(),
		};

		switch (static_cast<TrkrDefs::TrkrId>(TrkrDefs::getTrkrId(*cluster_key_itr))) {
			case TrkrDefs::mvtxId:
				++n_mvtx;
				break;
			case TrkrDefs::inttId:
				++n_intt;
				break;
			case TrkrDefs::tpcId: {
				++n_tpc;
				modify_point_tpc(*cluster_key_itr, cluster, point);
				break;
			}
			default:
				continue;
		}

		m_points.push_back (point);
	}

	if (n_mvtx < m_num_mvtx) { return true; }
	if (n_intt < m_num_intt) { return true; }
	if (n_tpc < m_num_tpc) { return true; }

	return false;
}

bool
WeightedFitter::do_fit (
) {
	// Get initial values for parameters from the first 2 m_points
	if (m_points.size() < 2) { return true; }

	// slope
	Eigen::Vector3d m = (m_points[1].pos - m_points[0].pos).normalized();
	double phi = std::atan2(m(1), m(0));
	double theta = std::acos(m(2));

	// intercept
	// Solve for s (m * s + m_points[0].pos).dot({0, 0, 1}) = 0
	double s = -m_points[0].pos(2) / m(2);
	double x = (m * s + m_points[0].pos)(0);
	double y = (m * s + m_points[0].pos)(1);

	FitErrorCalculator fit_error(m_points);
	ROOT::Math::Functor get_fit_error (fit_error, &FitErrorCalculator::operator(), 4);
	m_minimizer->SetFunction(get_fit_error);

	m_minimizer->SetVariable(0, "x",     x,     1E-6);
	m_minimizer->SetVariable(1, "y",     y,     1E-6);
	m_minimizer->SetVariable(2, "theta", theta, 1E-6);
	m_minimizer->SetVariable(3, "phi",   phi,   1E-6);

	// Bounds for angular parameters
	m_minimizer->SetVariableLimits(2,       0, 3.1416);
	m_minimizer->SetVariableLimits(3, -3.1416, 3.1416);

	return !m_minimizer->Minimize();
}

bool
WeightedFitter::add_track (
	TrackSeed* track_seed
) {
	double const* params = m_minimizer->X();

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

	TrackSeed_v2 fitted_seed;
	for (
		auto cluster_key_itr = track_seed->begin_cluster_keys();
		cluster_key_itr != track_seed->end_cluster_keys();
		++cluster_key_itr
	) {
		fitted_seed.insert_cluster_key(*cluster_key_itr);
	}

	fitted_seed.set_qOverR(1.0);
	fitted_seed.set_phi(track_seed->get_phi());
	fitted_seed.set_X0(b(0)); // TODO: check this
	fitted_seed.set_Y0(b(1)); // TODO: check this
	fitted_seed.set_Z0(b(2)); // TODO: check this
	fitted_seed.set_slope(sin_theta / cos_theta); // TODO: check this

	SvtxTrack_v4 fitted_track;
	fitted_track.set_id(m_track_id);
	fitted_track.set_silicon_seed(track_seed);
	// fitted_track.set_tpc_seed();
	fitted_track.set_crossing(track_seed->get_crossing());
	fitted_track.set_charge(track_seed->get_charge());
	fitted_track.set_x(params[0]);
	fitted_track.set_y(params[1]);
	fitted_track.set_z(0.0);
	fitted_track.set_px(fitted_seed.get_p() * m(0));
	fitted_track.set_py(fitted_seed.get_p() * m(1));
	fitted_track.set_pz(fitted_seed.get_p() * m(2));

	SvtxAlignmentStateMap::StateVec alignment_states;
	for (auto const& point : m_points) {
		// Get intersection, solve s
		// (m * s + b - point.pos).dot(point.z) = 0
		double s = (point.pos - b).dot(point.z) / m.dot(point.z);
		Acts::Vector3 intersection = m * s + b; // global coordinates

    	SvtxTrackState_v1 svtx_track_state(intersection.norm());
		svtx_track_state.set_x(intersection(0));
		svtx_track_state.set_y(intersection(1));
		svtx_track_state.set_z(intersection(2));
		svtx_track_state.set_px(fitted_seed.get_p() * m(0));
		svtx_track_state.set_py(fitted_seed.get_p() * m(1));
		svtx_track_state.set_pz(fitted_seed.get_p() * m(2));
		fitted_track.insert_state(&svtx_track_state);

		// Move this implementation into the fitter helper class eventually
		auto alignment_state = std::make_unique<SvtxAlignmentState_v1>();
		SvtxAlignmentState::GlobalMatrix global_derivative_matrix = SvtxAlignmentState::GlobalMatrix::Zero();
		SvtxAlignmentState::LocalMatrix local_derivative_matrix = SvtxAlignmentState::LocalMatrix::Zero();

		// TODO: move these to the helper class implementation

		Acts::Vector2 residual {
			point.x.dot(intersection - point.o),
			point.y.dot(intersection - point.o),
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
			Eigen::Vector3d { cos_theta * cos_phi, cos_theta * sin_phi, -sin_theta},
			Eigen::Vector3d {-sin_theta * sin_phi, sin_theta * cos_phi, 0.0},
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

		alignment_state->set_cluster_key(point.cluster_key);
		alignment_state->set_residual(residual);
		alignment_state->set_local_derivative_matrix(local_derivative_matrix);
		alignment_state->set_global_derivative_matrix(global_derivative_matrix);
		alignment_states.push_back(alignment_state.release());

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

void
WeightedFitter::modify_point_tpc (
	TrkrDefs::cluskey cluster_key,
	TrkrCluster* cluster,
	WeightedFitter::ClusterFitPoint& point
) {
	// See HelicalFitter::convertTimeToZ

	// must convert local Y from cluster average time of arival to local cluster z position
	double const drift_velocity = m_geometry->get_drift_velocity();
	double const zdriftlength = cluster->getLocalY() * drift_velocity;
	double const surfCenterZ = 52.89;          // 52.89 is where G4 thinks the surface center is
	double zloc = surfCenterZ - zdriftlength;  // converts z drift length to local z position in the TPC in north
	unsigned int const side = TpcDefs::getSide(cluster_key);
	if (side == 0)
	{
	  zloc = -zloc;
	}
	point.pos.z() = zloc;  // in cm

	// See HelicalFitter::makeTpcGlobalCorrections

	// make all corrections to global position of TPC cluster
	// crossing argument is set to 0 (in HelicalFitter implementation)
	// TODO: maybe retrieve from cluster later?
	point.pos.z() = m_clusterCrossingCorrection.correctZ(point.pos.z(), side, 0);

	// apply distortion corrections
	point.pos = m_globalPositionWrapper.applyDistortionCorrections(point.pos);
}

void
WeightedFitter::draw (
) {
	TCanvas* canvas = new TCanvas (
		"canvas", "canvas",
		1600, 800
	);
	canvas->Divide(2, 1);

	// xy event display
	{
		canvas->cd(1);
		gPad->SetFillStyle(4000);
		double* x = new double[m_points.size()];
		double* y = new double[m_points.size()];
		for (std::size_t s = 0; s < m_points.size(); ++s) {
			x[s] = m_points[s].pos.x();
			y[s] = m_points[s].pos.y();
		}
		TGraph* graph_xy = new TGraph (
			m_points.size(), x, y
		);
		// 15cm is good for Si Only
		// 80cm is good for Gas detectors also
		double range = 80;
		graph_xy->GetXaxis()->SetLimits(-range, range);
		graph_xy->GetXaxis()->SetRangeUser(-range, range);
		graph_xy->GetYaxis()->SetLimits(-range, range);
		graph_xy->GetYaxis()->SetRangeUser(-range, range);
		graph_xy->SetMarkerStyle(kFullCircle);
		graph_xy->SetMarkerSize(1);
		graph_xy->SetMarkerColor(kBlack);
		graph_xy->Draw("AP");

		draw_fit_xy(canvas);
		for (auto const& point : m_points) {
			draw_cluster_xy(canvas, point);
		}
	}

	canvas->Update();
	canvas->Show();
	canvas->SaveAs("EventDisplay.png");
}

void
WeightedFitter::draw_fit_xy (
	TPad* pad
) {
	double x_min = -120;
	double x_max =  120;

	double const* params = m_minimizer->X();

	if (std::cos(params[3]) < 0) {
		x_max = 0;
	} else {
		x_min = 0;
	}

	pad->cd();
	TF1* func = new TF1 (
		"model",
		"[0] + [1] * x",
		x_min, x_max
	);

	double m = std::tan(params[3]);
	func->SetParameter(0, params[1] - m * params[0]);
	func->SetParameter(1, m); 

	func->SetLineStyle(kSolid);
	func->SetLineWidth(3);
	func->SetLineColor(kBlue);
	func->DrawCopy("LSAME");
}

void
WeightedFitter::draw_cluster_xy (
	TPad* pad,
	ClusterFitPoint const& cluster
) {
	pad->cd();

	std::array<Eigen::Vector3d, 2> endpoint = {
		cluster.pos - 5.0E+3 * cluster.sigma_x * cluster.x,
		cluster.pos + 5.0E+3 * cluster.sigma_x * cluster.x,
	};

	TLine* line = new TLine (
		endpoint[0].x(), endpoint[0].y(),
		endpoint[1].x(), endpoint[1].y()
	);
	line->SetLineStyle(kSolid);
	line->SetLineWidth(3);
	line->SetLineColor(kBlack);
	line->Draw("SAME");
}

