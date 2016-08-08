/*!
 *  \file		PHG4TrackFastSim.cc
 *  \brief		Kalman Filter based on smeared truth PHG4Hit
 *  \details	Kalman Filter based on smeared truth PHG4Hit
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#include <cmath>
#include <map>
#include <utility>
#include <fun4all/Fun4AllReturnCodes.h>
#include <GenFit/AbsMeasurement.h>
#include <GenFit/EventDisplay.h>
#include <GenFit/MeasuredStateOnPlane.h>
#include <GenFit/RKTrackRep.h>
#include <GenFit/StateOnPlane.h>
#include <GenFit/Track.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>

#include <TMath.h>
#include <TMatrixF.h>
#include <TRandom.h>
#include <TString.h>
#include <g4hough/SvtxTrackMap.h>
#include <g4hough/SvtxTrackMap_v1.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <phgenfit/Fitter.h>
#include <phgenfit/PlanarMeasurement.h>
#include <phgenfit/Track.h>
#include <phgeom/PHGeomUtility.h>
#include <phgenfit/SpacepointMeasurement.h>
#include "SvtxTrack.h"
#include "SvtxTrack_FastSim.h"

#include "PHG4TrackFastSim.h"


#define LogDebug(exp)		std::cout<<"DEBUG: "	<<__FILE__<<": "<<__LINE__<<": "<< exp <<"\n"
#define LogError(exp)		std::cout<<"ERROR: "	<<__FILE__<<": "<<__LINE__<<": "<< exp <<"\n"
#define LogWarning(exp)	std::cout<<"WARNING: "	<<__FILE__<<": "<<__LINE__<<": "<< exp <<"\n"

#define _DEBUG_MODE_ 0

//#define _N_DETECTOR_LAYER 5

using namespace std;

PHG4TrackFastSim::PHG4TrackFastSim(const std::string &name) :
		SubsysReco(name),
		_detector_type(Vertical_Plane),
		_truth_container(NULL),
		_sub_top_node_name("SVTX"), /*_clustermap_out_name("SvtxClusterMap"),*/ _trackmap_out_name("SvtxTrackMap"),
		/*_clustermap_out(NULL),*/_trackmap_out(NULL),
		_fitter (NULL),_mag_field_file_name("/phenix/upgrades/decadal/fieldmaps/fsPHENIX.2d.root"),
		 _mag_field_re_scaling_factor(1.), _reverse_mag_field(false), _fit_alg_name("KalmanFitterRefTrack"), _do_evt_display(false),
		 _use_vertex_in_fitting(true), _vertex_xy_resolution(50E-4), _vertex_z_resolution(50E-4),
		 _phi_resolution(50E-4), _r_resolution(1.), _z_resolution(50E-4),
		 _pat_rec_hit_finding_eff(1.), _pat_rec_nosise_prob(0.),
		 _N_DETECTOR_LAYER(5)
		 {

	_event = -1;

	for(int i=0;i<_N_DETECTOR_LAYER;i++) {
		_phg4hits_names.push_back(Form("G4HIT_FGEM_%1d",i));
	}
}

PHG4TrackFastSim::~PHG4TrackFastSim() {
	delete _fitter;
}

/*
 * Init
 */
int PHG4TrackFastSim::Init(PHCompositeNode *topNode) {

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TrackFastSim::InitRun(PHCompositeNode *topNode) {

	_event = -1;


	CreateNodes(topNode);

	TGeoManager* tgeo_manager = PHGeomUtility::GetTGeoManager(topNode);

	//_fitter = new PHGenFit::Fitter("sPHENIX_Geo.root","sPHENIX.2d.root", 1.4 / 1.5);
	_fitter = PHGenFit::Fitter::getInstance(tgeo_manager,
			_mag_field_file_name.data(), (_reverse_mag_field) ? -1.*_mag_field_re_scaling_factor : _mag_field_re_scaling_factor, "KalmanFitterRefTrack", "RKTrackRep",
			_do_evt_display);

	if (!_fitter) {
		cerr << PHWHERE << endl;
		return Fun4AllReturnCodes::ABORTRUN;
	}

	_fitter->set_verbosity(verbosity);

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TrackFastSim::End(PHCompositeNode *topNode) {

	if(_do_evt_display) {
		_fitter->displayEvent();
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TrackFastSim::process_event(PHCompositeNode *topNode) {

	_event++;

	if(verbosity >= 2) std::cout<< "PHG4TrackFastSim::process_event: "<< _event << ".\n";

	GetNodes(topNode);

//	if(_clustermap_out)
//		_clustermap_out->empty();
//	else {
//		LogError("_clustermap_out not found!");
//		return Fun4AllReturnCodes::ABORTRUN;
//	}

	if(_trackmap_out)
		_trackmap_out->empty();
	else {
		LogError("_trackmap_out not found!");
		return Fun4AllReturnCodes::ABORTRUN;
	}

	vector<genfit::Track*> rf_gf_tracks;

	for (PHG4TruthInfoContainer::ConstIterator itr =
			_truth_container->GetPrimaryParticleRange().first;
			itr != _truth_container->GetPrimaryParticleRange().second; ++itr) {
		PHG4Particle* particle = itr->second;

		TVector3 seed_pos(0, 0, 0);
		TVector3 seed_mom(0, 0, 0);
		TMatrixDSym seed_cov(6);

		//! Create measurements
		std::vector<PHGenFit::Measurement*> measurements;

//		_use_vertex_in_fitting = true;

		PHGenFit::Measurement* vtx_meas = NULL;

		if (_use_vertex_in_fitting) {
			vtx_meas = VertexMeasurement(
					TVector3(0, 0, 0), _vertex_xy_resolution, _vertex_z_resolution);
			measurements.push_back(vtx_meas);
		}

		PseudoPatternRecognition(particle, measurements, seed_pos, seed_mom, seed_cov);

		if (measurements.size() < 3) {
			if(verbosity >= 2) {
				//LogWarning("measurements.size() < 3");
				std::cout<<"event: "<< _event << " : measurements.size() < 3" <<"\n";
			}
			continue;
		}

		//! Build TrackRep from particle assumption
		/*!
		 * mu+:	-13
		 * mu-:	13
		 * pi+:	211
		 * pi-:	-211
		 * e-:	11
		 * e+:	-11
		 */
		int pid = 13; //
		//SMART(genfit::AbsTrackRep) rep = NEW(genfit::RKTrackRep)(pid);
		genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pid);

		//rep->setDebugLvl(1); //DEBUG

		//! Initiallize track with seed from pattern recognition
		PHGenFit::Track* track = new PHGenFit::Track(rep, seed_pos, seed_mom,
				seed_cov);

		rf_gf_tracks.push_back(track->getGenFitTrack());

		//LogDEBUG;
		//! Add measurements to track
		track->addMeasurements(measurements);

		//LogDEBUG;
		//! Fit the track
		int fitting_err = _fitter->processTrack(track, false);

		if(fitting_err != 0)
		{
			if(verbosity >= 2) {
				//LogWarning("measurements.size() < 3");
				std::cout<<"event: "<< _event << " : fitting_err != 0, next track." <<"\n";
			}
			continue;
		}

		SvtxTrack* svtx_track_out = MakeSvtxTrack(track,particle->get_track_id());

		_trackmap_out->insert(svtx_track_out);
	} // Loop all primary particles


	//! add tracks to event display
	if(_do_evt_display)
		_fitter->getEventDisplay()->addEvent(rf_gf_tracks);
	else
		rf_gf_tracks.clear();

//	if(_trackmap_out->get(0)) {
//		_trackmap_out->get(0)->identify();
//		std::cout<<"DEBUG : "<< _trackmap_out->get(0)->get_px() <<"\n";
//		std::cout<<"DEBUG : "<< _trackmap_out->get(0)->get_truth_track_id() <<"\n";
//	}

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TrackFastSim::CreateNodes(PHCompositeNode *topNode) {
	// create nodes...
	PHNodeIterator iter(topNode);

	PHCompositeNode *dstNode = static_cast<PHCompositeNode*>(iter.findFirst(
			"PHCompositeNode", "DST"));
	if (!dstNode) {
		cerr << PHWHERE << "DST Node missing, doing nothing." << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	// Create the FGEM node
	PHCompositeNode* tb_node = dynamic_cast<PHCompositeNode*>(iter.findFirst(
			"PHCompositeNode", _sub_top_node_name.c_str()));
	if (!tb_node) {
		tb_node = new PHCompositeNode(_sub_top_node_name.c_str());
		dstNode->addNode(tb_node);
		if (verbosity > 0)
			cout << _sub_top_node_name.c_str() <<" node added" << endl;
	}


//	_clustermap_out = new SvtxClusterMap_v1;
//
//	PHIODataNode<PHObject>* clusters_node = new PHIODataNode<PHObject>(
//			_clustermap_out, _clustermap_out_name.c_str(), "PHObject");
//	tb_node->addNode(clusters_node);
//	if (verbosity > 0)
//		cout << _clustermap_out_name.c_str() <<" node added" << endl;

	_trackmap_out = new SvtxTrackMap_v1;

	PHIODataNode<PHObject>* tracks_node = new PHIODataNode<PHObject>(
			_trackmap_out, _trackmap_out_name.c_str(), "PHObject");
	tb_node->addNode(tracks_node);
	if (verbosity > 0)
		cout << _trackmap_out_name.c_str() << " node added" << endl;

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TrackFastSim::GetNodes(PHCompositeNode * topNode) {
	//DST objects
	//Truth container
	_truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode,
			"G4TruthInfo");
	if (!_truth_container && _event < 2) {
		cout << PHWHERE << " PHG4TruthInfoContainer node not found on node tree"
				<< endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	for(int i=0;i<_N_DETECTOR_LAYER;i++) {
		PHG4HitContainer* phg4hit = findNode::getClass<PHG4HitContainer>(topNode,
				_phg4hits_names[i].c_str());
		if (!phg4hit && _event < 2) {
			cout << PHWHERE << _phg4hits_names[i].c_str() <<" node not found on node tree"
					<< endl;
			return Fun4AllReturnCodes::ABORTEVENT;
		}

		_phg4hits.push_back(phg4hit);
	}


//	_clustermap_out = findNode::getClass<SvtxClusterMap>(topNode,
//			_clustermap_out_name.c_str());
//	if (!_clustermap_out && _event < 2) {
//		cout << PHWHERE << _clustermap_out_name.c_str() << " node not found on node tree"
//				<< endl;
//		return Fun4AllReturnCodes::ABORTEVENT;
//	}

	_trackmap_out = findNode::getClass<SvtxTrackMap>(topNode,
			_trackmap_out_name.c_str());
	if (!_trackmap_out && _event < 2) {
		cout << PHWHERE << _trackmap_out_name.c_str() << " node not found on node tree"
				<< endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TrackFastSim::PseudoPatternRecognition(
		const PHG4Particle* particle,
		std::vector<PHGenFit::Measurement*>& meas_out, TVector3& seed_pos,
		TVector3& seed_mom, TMatrixDSym& seed_cov, const bool do_smearing) {

	seed_pos.SetXYZ(0, 0, 0);
	seed_mom.SetXYZ(0, 0, 10);
	seed_cov.ResizeTo(6, 6);

	for (int i = 0; i < 3; i++) {
		seed_cov[i][i] = _phi_resolution * _phi_resolution;
	}

	for (int i = 3; i < 6; i++) {
		seed_cov[i][i] = 10;
	}

	if (particle) {
		TVector3 True_mom(particle->get_px(), particle->get_py(),
				particle->get_pz());

		seed_mom.SetXYZ(particle->get_px(), particle->get_py(),
				particle->get_pz());
		if(do_smearing)
		{
			const double momSmear = 3. / 180. * TMath::Pi();     // rad
			const double momMagSmear = 0.1;   // relative

			seed_mom.SetPhi(gRandom->Gaus(True_mom.Phi(), momSmear));
			seed_mom.SetTheta(gRandom->Gaus(True_mom.Theta(), momSmear));
			seed_mom.SetMag(
					gRandom->Gaus(True_mom.Mag(),
							momMagSmear * True_mom.Mag()));
		}
	}

	for (int ilayer = 0; ilayer < _N_DETECTOR_LAYER; ilayer++) {

		if (!_phg4hits[ilayer]) {
			LogError("No _phg4hits[i] found!");
			continue;
		}

#if _DEBUG_MODE_ == 1
		std::cout<<"DEBUG: ilayer: " << ilayer <<"; nsublayers: " <<_phg4hits[ilayer]->num_layers() <<" \n";
#endif

		for (PHG4HitContainer::LayerIter layerit =
				_phg4hits[ilayer]->getLayers().first;
				layerit != _phg4hits[ilayer]->getLayers().second; layerit++) {
			for (PHG4HitContainer::ConstIterator itr =
					_phg4hits[ilayer]->getHits(*layerit).first;
					itr != _phg4hits[ilayer]->getHits(*layerit).second; ++itr) {

				PHG4Hit * hit = itr->second;
				if (!hit) {
					LogDebug("No PHG4Hit Found!");
					continue;
				}

#if _DEBUG_MODE_ == 1
				std::cout<<"DEBUG: ilayer: " << ilayer <<"; sublayer: " <<*layerit << "; itr->first : "<< itr->first <<" \n";
				//hit->identify();
#endif

				if (hit->get_trkid() == particle->get_track_id()
						|| gRandom->Uniform(0, 1) < _pat_rec_nosise_prob) {

					PHGenFit::Measurement* meas = NULL;
					if (_detector_type == Vertical_Plane)
						meas = PHG4HitToMeasurementVerticalPlane(hit,
								_phi_resolution, _r_resolution);
					else if (_detector_type == Cylinder)
						meas = PHG4HitToMeasurementCylinder(hit,
								_phi_resolution, _z_resolution);
					else {
						LogError("Type not implemented!");
						return Fun4AllReturnCodes::ABORTEVENT;
					}
					if (gRandom->Uniform(0, 1) <= _pat_rec_hit_finding_eff) {
						meas_out.push_back(meas);
						//meas->getMeasurement()->Print(); //DEBUG
					}
				}
			}
		} /*Loop layers within one detector layer*/
	} /*Loop detector layers*/

	return Fun4AllReturnCodes::EVENT_OK;
}

SvtxTrack* PHG4TrackFastSim::MakeSvtxTrack(
		const PHGenFit::Track* phgf_track, const unsigned int truth_track_id) {

	double chi2 = phgf_track->get_chi2();
	double ndf = phgf_track->get_ndf();

	genfit::MeasuredStateOnPlane* gf_state = phgf_track->extrapolateToPlane(TVector3(0.,0.,0.), TVector3(0.,0.,1.));
	TVector3 mom = gf_state->getMom();
	TVector3 pos = gf_state->getPos();
	TMatrixDSym cov = gf_state->get6DCov();

//	SvtxTrack_v1* out_track = new SvtxTrack_v1(*static_cast<const SvtxTrack_v1*> (svtx_track));
//	SvtxTrack_v1* out_track = new SvtxTrack_v1();

	SvtxTrack* out_track = new SvtxTrack_FastSim();

	out_track->set_truth_track_id(truth_track_id);

	/*!
	 * FIXME: check the definition
	 *  1/p, u'/z', v'/z', u, v
	 *  u is defined as mom X beam line at POCA
	 *  so u is the dca2d direction
	 */
	double dca2d = gf_state->getState()[3];
	out_track->set_dca2d(dca2d);
	out_track->set_dca2d_error(gf_state->getCov()[3][3]);
	double dca3d = sqrt(
			dca2d*dca2d +
			gf_state->getState()[4]*gf_state->getState()[4]);
	out_track->set_dca(dca3d);


	out_track->set_chisq(chi2);
	out_track->set_ndf(ndf);
	out_track->set_charge((_reverse_mag_field) ? -1.*phgf_track->get_charge() : phgf_track->get_charge());

	out_track->set_px(mom.Px());
	out_track->set_py(mom.Py());
	out_track->set_pz(mom.Pz());

	out_track->set_x(pos.X());
	out_track->set_y(pos.Y());
	out_track->set_z(pos.Z());

	for(int i=0;i<6;i++)
	{
		for(int j=i;j<6;j++)
		{
			out_track->set_error(i,j,cov[i][j]);
		}
	}

	return out_track;
}



PHGenFit::PlanarMeasurement* PHG4TrackFastSim::PHG4HitToMeasurementVerticalPlane(
		const PHG4Hit* g4hit, const double phi_resolution,
		const double r_resolution) {

	PHGenFit::PlanarMeasurement* meas = NULL;

	TVector3 pos(g4hit->get_avg_x(), g4hit->get_avg_y(), g4hit->get_avg_z());

	TVector3 v(pos.X(), pos.Y(), 0);
	v = 1 / v.Mag() * v;

	TVector3 u = v.Cross(TVector3(0,0,1));
	u = 1/u.Mag() * u;

	double u_smear = gRandom->Gaus(0, phi_resolution);
	double v_smear = gRandom->Gaus(0, r_resolution);
	pos.SetX(g4hit->get_avg_x() + u_smear*u.X() + v_smear*v.X());
	pos.SetY(g4hit->get_avg_y() + u_smear*u.Y() + v_smear*v.Y());


	meas = new PHGenFit::PlanarMeasurement(pos, u, v, phi_resolution, r_resolution);

//	std::cout<<"------------\n";
//	pos.Print();
//	std::cout<<"at "<<istation<<" station, "<<ioctant << " octant \n";
//	u.Print();
//	v.Print();

	//dynamic_cast<PHGenFit::PlanarMeasurement*> (meas)->getMeasurement()->Print();

	return meas;
}


PHGenFit::PlanarMeasurement* PHG4TrackFastSim::PHG4HitToMeasurementCylinder(const PHG4Hit* g4hit, const double phi_resolution, const double z_resolution) {

	PHGenFit::PlanarMeasurement* meas = NULL;

	TVector3 pos(
			g4hit->get_avg_x(),
			g4hit->get_avg_y(),
			g4hit->get_avg_z());


	TVector3 v(0, 0, pos.Z());
	v = 1/v.Mag() * v;

	TVector3 u = v.Cross(TVector3(pos.X(), pos.Y(), 0));
	u = 1/u.Mag() * u;

	double u_smear = gRandom->Gaus(0, phi_resolution);
	double v_smear = gRandom->Gaus(0, z_resolution);
	pos.SetX(g4hit->get_avg_x() + u_smear*u.X() + v_smear*v.X());
	pos.SetY(g4hit->get_avg_y() + u_smear*u.Y() + v_smear*v.Y());


	meas = new PHGenFit::PlanarMeasurement(pos, u, v, phi_resolution, z_resolution);

//	std::cout<<"------------\n";
//	pos.Print();
//	std::cout<<"at "<<istation<<" station, "<<ioctant << " octant \n";
//	u.Print();
//	v.Print();

	//dynamic_cast<PHGenFit::PlanarMeasurement*> (meas)->getMeasurement()->Print();

	return meas;
}

PHGenFit::Measurement* PHG4TrackFastSim::VertexMeasurement(const TVector3 &vtx, const double dxy,
		const double dz) {
	PHGenFit::Measurement* meas = NULL;

	TMatrixDSym cov(3);
	cov.Zero();
	cov(0, 0) = dxy*dxy;
	cov(1, 1) = dxy*dxy;
	cov(2, 2) = dz*dz;

	TVector3 pos = vtx;

	double xy_smear = gRandom->Gaus(0, dxy);
	double z_smear = gRandom->Gaus(0, dz);
	pos.SetX(vtx.X() + xy_smear);
	pos.SetY(vtx.Y() + xy_smear);
	pos.SetZ(vtx.Z() + z_smear);


	meas = new PHGenFit::SpacepointMeasurement(pos, cov);

	return meas;
}

















