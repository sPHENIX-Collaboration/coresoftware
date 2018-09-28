/*!
 *  \file		PHG4TrackKalmanFitter.C
 *  \brief		Refit SvtxTracks with PHGenFit.
 *  \details	Refit SvtxTracks with PHGenFit.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#include "PHG4TrackKalmanFitter.h"
#include "SvtxCluster.h"
#include "SvtxClusterMap.h"
#include "SvtxTrackState_v1.h"
#include "SvtxHit_v1.h"
#include "SvtxHitMap.h"
#include "SvtxTrack.h"
#include "SvtxTrack_v1.h"
#include "SvtxVertex_v1.h"
#include "SvtxTrackMap.h"
#include "SvtxTrackMap_v1.h"
#include "SvtxVertexMap_v1.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>

#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CylinderGeom_MAPS.h>
#include <g4detectors/PHG4CylinderGeom_Siladders.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Particlev2.h>
#include <g4main/PHG4VtxPointv1.h>

#include <phgenfit/Fitter.h>
#include <phgenfit/PlanarMeasurement.h>
#include <phgenfit/Track.h>
#include <phgenfit/SpacepointMeasurement.h>

#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>

#include <phgeom/PHGeomUtility.h>
#include <phfield/PHFieldUtility.h>

#include <GenFit/FieldManager.h>
#include <GenFit/GFRaveVertex.h>
#include <GenFit/GFRaveVertexFactory.h>
#include <GenFit/MeasuredStateOnPlane.h>
#include <GenFit/RKTrackRep.h>
#include <GenFit/StateOnPlane.h>
#include <GenFit/Track.h>
#include <GenFit/KalmanFitterInfo.h>

//Rave
#include <rave/Version.h>
#include <rave/Track.h>
#include <rave/VertexFactory.h>
#include <rave/ConstantMagneticField.h>

//GenFit
#include <GenFit/GFRaveConverters.h>

#include <TClonesArray.h>
#include <TMatrixDSym.h>
#include <TTree.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <TRotation.h>



#include <iostream>
#include <map>
#include <utility>
#include <vector>
#include <memory>


#define LogDebug(exp)		std::cout<<"DEBUG: "  <<__FILE__<<": "<<__LINE__<<": "<< exp <<std::endl
#define LogError(exp)		std::cout<<"ERROR: "  <<__FILE__<<": "<<__LINE__<<": "<< exp <<std::endl
#define LogWarning(exp)	std::cout<<"WARNING: "<<__FILE__<<": "<<__LINE__<<": "<< exp <<std::endl

#define WILD_FLOAT -9999.

#define _DEBUG_MODE_ 0

//#define _DEBUG_

using namespace std;


class PHRaveVertexFactory {

public:
	//! ctor
	PHRaveVertexFactory(const int verbosity) {
		rave::ConstantMagneticField mfield(0., 0., 0.); // RAVE use Tesla
		_factory = new rave::VertexFactory(mfield, rave::VacuumPropagator(),
				"default", verbosity);

		IdGFTrackStateMap_.clear();
	}

	//! dotr
	~PHRaveVertexFactory() {
		clearMap();

		delete _factory;
	}

	void findVertices(std::vector<genfit::GFRaveVertex*>* vertices,
			const std::vector<genfit::Track*>& tracks, const bool use_beamspot =
					false) {

		clearMap();

		try {
			genfit::RaveToGFVertices(vertices,
					_factory->create(
							genfit::GFTracksToTracks(tracks, NULL,
									IdGFTrackStateMap_, 0), use_beamspot),
					IdGFTrackStateMap_);
		} catch (genfit::Exception & e) {
			std::cerr << e.what();
		}
	}

	void findVertices(std::vector<genfit::GFRaveVertex*>* vertices,
			const std::vector<genfit::Track*>& tracks,
			std::vector<genfit::MeasuredStateOnPlane*> & GFStates,
			const bool use_beamspot = false) {

		clearMap();

		try {
			genfit::RaveToGFVertices(vertices,
					_factory->create(
							genfit::GFTracksToTracks(tracks, &GFStates,
									IdGFTrackStateMap_, 0), use_beamspot),
					IdGFTrackStateMap_);
		} catch (genfit::Exception & e) {
			std::cerr << e.what();
		}
	}

private:
	void clearMap() {

		for (unsigned int i = 0; i < IdGFTrackStateMap_.size(); ++i)
			delete IdGFTrackStateMap_[i].state_;

		IdGFTrackStateMap_.clear();
	}

	std::map<int, genfit::trackAndState> IdGFTrackStateMap_;

	rave::VertexFactory* _factory;

};

/*
 * Constructor
 */
PHG4TrackKalmanFitter::PHG4TrackKalmanFitter(const string &name) :
		SubsysReco(name),
		_flags(NONE),
		_output_mode(PHG4TrackKalmanFitter::MakeNewNode),
		_over_write_svtxtrackmap(true),
		_over_write_svtxvertexmap(true),
		_fit_primary_tracks(false),
		_use_truth_vertex(false),
		_fitter( NULL),
		_track_fitting_alg_name("DafRef"),
		_primary_pid_guess(211),
		_fit_min_pT(0.1),
		_vertex_min_ndf(20),
		_vertex_finder(NULL),
		_vertexing_method("avf-smoothing:1"),
		_truth_container(NULL),
		_clustermap(NULL),
		_trackmap(NULL),
		_vertexmap(NULL),
		_trackmap_refit(NULL),
		_primary_trackmap(NULL),
		_vertexmap_refit(NULL),
		_do_eval(false),
		_eval_outname("PHG4TrackKalmanFitter_eval.root"),
		_eval_tree(NULL),
		_tca_particlemap(NULL),
		_tca_vtxmap(NULL),
		_tca_trackmap(NULL),
		_tca_vertexmap(NULL),
		_tca_trackmap_refit(NULL),
		_tca_primtrackmap(NULL),
		_tca_vertexmap_refit(NULL),
		_do_evt_display(false) {

	Verbosity(0);

	_event = 0;

	_cluster_eval_tree = NULL;
	_cluster_eval_tree_x = WILD_FLOAT;
	_cluster_eval_tree_y = WILD_FLOAT;
	_cluster_eval_tree_z = WILD_FLOAT;
	_cluster_eval_tree_gx = WILD_FLOAT;
	_cluster_eval_tree_gy = WILD_FLOAT;
	_cluster_eval_tree_gz = WILD_FLOAT;
}

/*
 * Init
 */
int PHG4TrackKalmanFitter::Init(PHCompositeNode *topNode) {

  //RCC add extrapolation ntuple definition here.
//	CreateNodes(topNode);

	return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * Init run
 */
int PHG4TrackKalmanFitter::InitRun(PHCompositeNode *topNode) {

	CreateNodes(topNode);

	TGeoManager* tgeo_manager = PHGeomUtility::GetTGeoManager(topNode);
	PHField * field = PHFieldUtility::GetFieldMapNode(nullptr, topNode);

	//_fitter = new PHGenFit::Fitter("sPHENIX_Geo.root","sPHENIX.2d.root", 1.4 / 1.5);
	_fitter = PHGenFit::Fitter::getInstance(tgeo_manager,
	    field, _track_fitting_alg_name,
			"RKTrackRep", _do_evt_display);
	_fitter->set_verbosity(verbosity);

	if (!_fitter) {
		cerr << PHWHERE << endl;
		return Fun4AllReturnCodes::ABORTRUN;
	}

	//LogDebug(genfit::FieldManager::getInstance()->getFieldVal(TVector3(0, 0, 0)).Z());

	_vertex_finder = new genfit::GFRaveVertexFactory(verbosity);
	//_vertex_finder->setMethod("kalman-smoothing:1"); //! kalman-smoothing:1 is the defaul method
	_vertex_finder->setMethod(_vertexing_method.data());
	//_vertex_finder->setBeamspot();

	//_vertex_finder = new PHRaveVertexFactory(verbosity);

	if (!_vertex_finder) {
		cerr << PHWHERE << endl;
		return Fun4AllReturnCodes::ABORTRUN;
	}

	if (_do_eval) {
		if(verbosity >= 1)
			cout << PHWHERE << " Opening file: " << _eval_outname << endl;
		PHTFileServer::get().open(_eval_outname, "RECREATE");
		init_eval_tree();
	}

	return Fun4AllReturnCodes::EVENT_OK;
}
/*
 * process_event():
 *  Call user instructions for every event.
 *  This function contains the analysis structure.
 *
 */
int PHG4TrackKalmanFitter::process_event(PHCompositeNode *topNode) {
	_event++;

	if(verbosity > 1)
		std::cout << PHWHERE << "(Heads up! You're running a non-stock PHG4TrackKalmanFitter) Events processed: " << _event << std::endl;
//	if (_event % 1000 == 0)
//		cout << PHWHERE << "Events processed: " << _event << endl;

	GetNodes(topNode);

	//! stands for Refit_GenFit_Tracks
	vector<genfit::Track*> rf_gf_tracks;
	rf_gf_tracks.clear();

	vector< std::shared_ptr<PHGenFit::Track> > rf_phgf_tracks;
	rf_phgf_tracks.clear();

	map<unsigned int, unsigned int> svtxtrack_genfittrack_map;

	if (_trackmap_refit)
		_trackmap_refit->empty();

	for (SvtxTrackMap::Iter iter = _trackmap->begin(); iter != _trackmap->end();
			++iter) {
		SvtxTrack* svtx_track = iter->second;
		if (!svtx_track)
			continue;
		if (!(svtx_track->get_pt() > _fit_min_pT))
			continue;

		//! stands for Refit_PHGenFit_Track
		std::shared_ptr<PHGenFit::Track> rf_phgf_track = ReFitTrack(topNode, svtx_track);

		if (rf_phgf_track) {
			svtxtrack_genfittrack_map[svtx_track->get_id()] =
					rf_phgf_tracks.size();
			rf_phgf_tracks.push_back(rf_phgf_track);
			if(rf_phgf_track->get_ndf() > _vertex_min_ndf)
				rf_gf_tracks.push_back(rf_phgf_track->getGenFitTrack());
		}

		//RCC do some studies of track extrapolation to the tpc:
	     
	       if (_do_eval) {
		 if (verbosity >= 2) std::cout << PHWHERE << "Starting extrapolation Eval"<<endl;

		if (rf_phgf_track){
		  //line_point and line_direction are the axis of the cylinder.
		  const TVector3 line_point=TVector3(0.,0.,0.);
		  const TVector3 line_direction=TVector3(0.,0.,1.);
		  const int inner_wall_r=20.0;//cm.  Need tocheck the native units RCC

		  int npoints=rf_phgf_track->getGenFitTrack()->getNumPointsWithMeasurement();
		   if (verbosity >= 2) std::cout << PHWHERE << npoints << " points in measured track RCC"<<endl;

		  //this ought to give the same result regardless of where we extrapolate from, since it should pick the best track ref, but... not clear if that's true.  pick the outermost point.  Note that currently this doesn't include extra material effects, so it won't work for points beyond the inner wall of the tpc -- there's an additional scattering growth to the covariance.
		  std::shared_ptr<genfit::MeasuredStateOnPlane> gf_cylinder_state = NULL;

		  bool covariance_okay=true;

		  try {
		    gf_cylinder_state = std::shared_ptr < genfit::MeasuredStateOnPlane
							  > (rf_phgf_track->extrapolateToCylinder(inner_wall_r,line_point,line_direction,npoints-1,1));
		  } catch (...) {
		    if (verbosity >= 2)
		      LogWarning("extrapolateToCylinder failed!");
		  }
		  if (!gf_cylinder_state) {
		    covariance_okay=false;
		    LogWarning("No Cylinder state RCC");
		    //std::cout << PHWHERE << "Cylinder State not filled RCC"<<endl;
		  }

		  TVector3 pos(NAN,NAN,NAN);
		  TVector3 mom(NAN,NAN,NAN);
		  //float pos_rphi_error = NAN;
		  //float pos_z_error  = NAN;

		  TMatrixF pos_in(3,1);
		  TMatrixF cov_in(3,3);
		  TMatrixF pos_out(3,1);
		  TMatrixF cov_out(3,3);
		  if(covariance_okay){
		    try{
		      TVectorD state6(6); // pos(3), mom(3)
		      TMatrixDSym cov6(6,6); //
		      
		      gf_cylinder_state->get6DStateCov(state6, cov6);
		    
		      //for extrapolation to cylinder, the normal direction is just the xy position of the hit (I think...)
		      //look at the other instance for confirmation
		      TVector3 vn(state6[0], state6[1], 0);
		      
		      pos_in[0][0] = state6[0];
		      pos_in[1][0] = state6[1];
		      pos_in[2][0] = state6[2];
		      pos.SetXYZ(state6[0],state6[1],state6[2]);
		      mom.SetXYZ(state6[3],state6[4],state6[5]);

		      for(int i=0;i<3;++i){
			for(int j=0;j<3;++j){
			  cov_in[i][j] = cov6[i][j];
			    }
		      }

		      if (verbosity > 4){
			cout << "cov_kalman:"<<endl;
			for (int j=0;j<3;j++){
			  for (int k=0;k<3;k++){
			    cout << cov_in[j][k] << '\t';
			  }
			  cout << endl;
			}
		      }
		      
		      pos_cov_XYZ_to_RZ(vn, pos_in, cov_in, pos_out, cov_out);
		      
		      //pos_rphi = pos_out[1][0];
		      //pos_z  = pos_out[2][0];
		      //pos_rphi_error = sqrt(cov_out[0][0]);
		      //pos_z_error  = sqrt(cov_out[2][2]);	
		    } catch (...) {
		      if (verbosity > 0)
			LogWarning("State Covariance unavailable!");
		      // std::cout << PHWHERE << "Extraction of State Cov failed RCC"<<endl;

		      covariance_okay=false;
		    }
		  }

		  TVector3 pos_linear[3];//the linear extrapolation from the last two hits.
		  TMatrixF cov_linout[3];//three of these, for the three covariance matrices in the linear approach(3,3);
		    float t=0; //path length from [1] to [2] in units of the distance between [0] and [1]

		  //also try a straight-line extrapolation: (rcchere)
		  try{
		    //std::shared_ptr <genfit::TrackPoint> trackpoint[2];
		    genfit::TrackPoint* trackpoint[2];
		    genfit::MeasuredStateOnPlane* trackstate[2];

		    TVectorD hit[3];//don't actually need the last one
			      
		    TMatrixF cov_linear[3];//(3,3);
		    TVectorD state_lin6(6); // pos(3), mom(3)
		    TMatrixDSym cov_lin6(6,6);
		    for (int i=0;i<3;i++){
		      hit[i].ResizeTo(3);
		      pos_linear[i].SetXYZ(NAN,NAN,NAN);
		      cov_linear[i].ResizeTo(3,3);
		      cov_linout[i].ResizeTo(3,3);

		    }
		    
      
		    //get the info for the last two points in the track, in order:
		    for (int i=0;i<2;i++){
		      trackpoint[i]= rf_phgf_track->getGenFitTrack()->getPointWithMeasurementAndFitterInfo(npoints-2+i);
		      trackstate[i]=trackpoint[i]->getFitterInfo()->getFittedState().clone();
		      //not sure this works right: pos_linear[i]=trackstate[i]->getPos();
		      //cov_lin6=trackstate[i]->get6DCov();  For reasons that aren't clear, this call does not work correctly.  Get as below.
		      trackstate[i]->get6DStateCov(state_lin6, cov_lin6);
		      (pos_linear[i]).SetXYZ(state_lin6[0],state_lin6[1],state_lin6[2]);

		      
		      for(int j=0;j<3;++j){ 
			  for(int k=0;k<3;++k){
			    //std::cout << PHWHERE << " filling cov_linear[ "<<j<<"]["<<k<<"]" << std::endl;
			    //std::cout << PHWHERE << " con_lin6[ "<<j<<"]["<<k<<"]="<<cov_lin6[j][k]<< std::endl;
			    (cov_linear[i])[j][k] = cov_lin6[j][k];
			  }
			}
		    }



		    //extrapolate from those last two points to the cylinder:
		    TVector3 delta=pos_linear[1]-pos_linear[0];
		    TMatrixF cov_delta=cov_linear[1]-cov_linear[0];
		    
		    //to extrapolate to a cylinder of radius r, we want to solve:
		    //(pos')^2=R^2, where pos'=delta*t+pos, a vector quantity.
		    //This leads to a quadratic equation: 
		    double quad_a=delta*delta;
		    double quad_b=2*delta*(pos_linear[1]);
		    double quad_c=(pos_linear[1]*pos_linear[1])-inner_wall_r*inner_wall_r;
		    //pick the positive root to keep going the direction we defined:
		    t=(-quad_b+TMath::Sqrt(quad_b*quad_b-4*quad_a*quad_c))/(2*quad_a);
		    //extrapolate out all the dimensions to the intersection.
		    pos_linear[2]=pos_linear[1]+delta*t;
		    cov_linear[2]=cov_linear[1]+cov_delta*t;

		    if (verbosity > 4){
		      
		      for(int i=0;i<3;i++){
			cout << "pos_linear " << i << ":"<<endl;
			for (int j=0;j<3;j++){
			  cout << (pos_linear[i])[j] << '\t';
			}
			cout  << "R="<< pos_linear[i].Perp() << " Phi="<<pos_linear[i].Phi() << " Z=" <<pos_linear[i].Z();
			cout << endl;
		      }
		      
		      cout << "pathlength=" << t << endl;
		      for(int i=0;i<3;i++){
			cout << "cov_linear " << i << ":"<<endl;
			for (int j=0;j<3;j++){
			  for (int k=0;k<3;k++){
			    cout << (cov_linear[i])[j][k] << '\t';
			  }
			  cout << endl;
			}
		      }
		      cout << "cov_delta:"<<endl;
		      for (int j=0;j<3;j++){
			for (int k=0;k<3;k++){
			  cout << (cov_delta)[j][k] << '\t';
			}
			cout << endl;
		      }
		    }
		    //rotate cov_linear to Rzphi for each of our three matrices:
		    for (int i=0;i<3;i++){
		      pos_cov_XYZ_to_RZ(pos_linear[2], pos_in, cov_linear[i], pos_out, cov_linout[i]);//don't need the pos rotations, but had to fill those in again anyway.
		      //std::cout << PHWHERE<< " cov_linout["<<i<<"] has " << cov_linout[i].GetNrows() << " rows and " << cov_linout[i].GetNcols() << " cols."<<endl;
		    }

		    
		  } catch (...) {
		    if (verbosity >= 2)
		      LogWarning("gf_cylinder_state_prev failed!");
		    //std::cout <<PHWHERE << " Problem in lienar extrapolation."<< endl;
		  }
		  
		  if (verbosity > 0){// and !covariance_okay){
		    std::cout << PHWHERE << "Covariance_okay="<< (covariance_okay?"true":"false") <<" RCC!" <<endl;
		  }
		  _kalman_extrapolation_eval_tree_phi=pos.Phi();
		  _kalman_extrapolation_eval_tree_z=pos.Z();
		  _kalman_extrapolation_eval_tree_r=pos.Perp();
		  _kalman_extrapolation_eval_tree_okay=covariance_okay;
		  if (covariance_okay){
		    _kalman_extrapolation_eval_tree_sigma_r=cov_out[0][0];
		    _kalman_extrapolation_eval_tree_sigma_rphi=cov_out[1][1];
		    _kalman_extrapolation_eval_tree_sigma_z=cov_out[2][2];
		    _kalman_extrapolation_eval_tree_sigma_r_rphi=cov_out[0][1];
		    _kalman_extrapolation_eval_tree_sigma_rphi_z=cov_out[1][2];
		    _kalman_extrapolation_eval_tree_sigma_z_r=cov_out[2][0];
		    //before rotation:
		    _kalman_extrapolation_eval_tree_covin_x=cov_in[0][0];
		    _kalman_extrapolation_eval_tree_covin_y=cov_in[1][1];
		    _kalman_extrapolation_eval_tree_covin_z=cov_in[2][2];
		    
		    //linear data:
		    _kalman_extrapolation_eval_tree_lin_pathlength=t;

		    //std::cout << PHWHERE <<endl;
		    _kalman_extrapolation_eval_tree_lin_pos0_x=(pos_linear[0]).X();
		    _kalman_extrapolation_eval_tree_lin_pos0_y=(pos_linear[0]).Y();
		    _kalman_extrapolation_eval_tree_lin_pos0_z=(pos_linear[0]).Z();
		    //std::cout << PHWHERE<<endl;
		    _kalman_extrapolation_eval_tree_lin_pos1_x=(pos_linear[1]).X();
		    _kalman_extrapolation_eval_tree_lin_pos1_y=(pos_linear[1]).Y();
		    _kalman_extrapolation_eval_tree_lin_pos1_z=(pos_linear[1]).Z();
		    //std::cout << PHWHERE<<endl;
		    _kalman_extrapolation_eval_tree_lin_pos2_x=(pos_linear[2]).X();
		    _kalman_extrapolation_eval_tree_lin_pos2_y=(pos_linear[2]).Y();
		    _kalman_extrapolation_eval_tree_lin_pos2_z=(pos_linear[2]).Z();
		    // std::cout << PHWHERE<<endl;
		    _kalman_extrapolation_eval_tree_lin_phi=(pos_linear[2]).Phi();
		    _kalman_extrapolation_eval_tree_lin_r=(pos_linear[2]).Perp();
		    _kalman_extrapolation_eval_tree_lin_z=(pos_linear[2]).Z();
		    //std::cout << PHWHERE<< " cov_linout[0] has " << cov_linout[0].GetNrows() << " rows and " << cov_linout[0].GetNcols() << " cols."<<endl;
		    _kalman_extrapolation_eval_tree_lin_sigma0_r=(cov_linout[0])[0][0];
		    _kalman_extrapolation_eval_tree_lin_sigma0_rphi=(cov_linout[0])[1][1];
		    _kalman_extrapolation_eval_tree_lin_sigma0_z=(cov_linout[0])[2][2];
		    //std::cout << PHWHERE<< " cov_linout[1] has " << cov_linout[1].GetNrows() << " rows and " << cov_linout[1].GetNcols() << " cols."<<endl;
		    _kalman_extrapolation_eval_tree_lin_sigma1_r=(cov_linout[1])[0][0];
		    _kalman_extrapolation_eval_tree_lin_sigma1_rphi=(cov_linout[1])[1][1];
		    _kalman_extrapolation_eval_tree_lin_sigma1_z=(cov_linout[1])[2][2];
		    //std::cout << PHWHERE<< " cov_linout[2] has " << cov_linout[2].GetNrows() << " rows and " << cov_linout[2].GetNcols() << " cols."<<endl;
		    _kalman_extrapolation_eval_tree_lin_sigma2_r=(cov_linout[2])[0][0];
		    _kalman_extrapolation_eval_tree_lin_sigma2_rphi=(cov_linout[2])[1][1];
		    _kalman_extrapolation_eval_tree_lin_sigma2_z=(cov_linout[2])[2][2];
		    //std::cout << PHWHERE<<endl;
		  }
		    
		    _kalman_extrapolation_eval_tree->Fill();
		    if (verbosity >= 2) std::cout << PHWHERE << "end of extrapolation eval."<<endl;
		  
		}
	       }
	}


	/*
	 * add tracks to event display
	 * needs to make copied for smart ptrs will be destroied even
	 * there are still references in TGeo::EventView
	 */
	if (_do_evt_display) {
		vector<genfit::Track*> copy;
		for(genfit::Track* t : rf_gf_tracks){
			copy.push_back(new genfit::Track(*t));
		}
		_fitter->getEventDisplay()->addEvent(copy);
	}

	//! find vertex using tracks
	std::vector<genfit::GFRaveVertex*> rave_vertices;
	rave_vertices.clear();
	if (rf_gf_tracks.size() >= 2) {
		//_vertex_finder->findVertices(&rave_vertices,rf_gf_tracks,rf_gf_states);
		try {
			_vertex_finder->findVertices(&rave_vertices, rf_gf_tracks);
		} catch (...) {
			if(verbosity > 1)
				std::cout << PHWHERE << "GFRaveVertexFactory::findVertices failed!";
		}
	}

	FillSvtxVertexMap(rave_vertices, rf_gf_tracks);

	for (SvtxTrackMap::Iter iter = _trackmap->begin(); iter != _trackmap->end();
			++iter) {
		std::shared_ptr<PHGenFit::Track> rf_phgf_track = NULL;

		if (svtxtrack_genfittrack_map.find(iter->second->get_id())
				!= svtxtrack_genfittrack_map.end()) {
			unsigned int itrack =
					svtxtrack_genfittrack_map[iter->second->get_id()];
			rf_phgf_track = rf_phgf_tracks[itrack];
		}

		if (rf_phgf_track) {

			//FIXME figure out which vertex to use.
			SvtxVertex* vertex = NULL;
			if (_over_write_svtxvertexmap) {
				if (_vertexmap->size() > 0)
					vertex = _vertexmap->get(0);
			} else {
				if (_vertexmap_refit->size() > 0)
					vertex = _vertexmap_refit->get(0);
			}

			//BEGIN DEBUG
//			vertex = NULL;
//
//			PHG4VtxPoint *truth_vtx = _truth_container->GetVtx(
//					_truth_container->GetPrimaryVertexIndex());
//			if(!truth_vtx) {
//				LogDebug("!truth_vtx");
//				return Fun4AllReturnCodes::ABORTEVENT;
//			}
//
//			LogDebug("");
//			truth_vtx->identify();
//
//			vertex = new SvtxVertex_v1();
//			vertex->set_x(truth_vtx->get_x());
//			vertex->set_y(truth_vtx->get_y());
//			vertex->set_z(truth_vtx->get_z());
//
//			for(int i=0;i<3;i++)
//				for(int j=0;j<3;j++)
//					vertex->set_error(i,j,0);
			//END DEBUG

//			SvtxTrack* rf_track = MakeSvtxTrack(iter->second, rf_phgf_track,
//					vertex);
#ifdef _DEBUG_
			cout<<__LINE__<<endl;
#endif
			std::shared_ptr<SvtxTrack> rf_track = MakeSvtxTrack(iter->second, rf_phgf_track,
					vertex);
#ifdef _DEBUG_
		cout<<__LINE__<<endl;
#endif
			if(!rf_track) {
				//if (_output_mode == OverwriteOriginalNode)
#ifdef _DEBUG_
						LogDebug("!rf_track, continue.");
#endif
				if (_over_write_svtxtrackmap)
					_trackmap->erase(iter->first);

				continue;
			}

//			delete vertex;//DEBUG

//			rf_phgf_tracks.push_back(rf_phgf_track);
//			rf_gf_tracks.push_back(rf_phgf_track->getGenFitTrack());

			if (!(_over_write_svtxtrackmap) || _output_mode == DebugMode)
				if (_trackmap_refit) {
					_trackmap_refit->insert(rf_track.get());
//					delete rf_track;
				}

			if (_over_write_svtxtrackmap
					|| _output_mode == DebugMode) {
				*(dynamic_cast<SvtxTrack_v1*>(iter->second)) =
						*(dynamic_cast<SvtxTrack_v1*>(rf_track.get()));
//				delete rf_track;
#ifdef _DEBUG_
		cout<<__LINE__<<endl;
#endif
			}
		} else {
			if (_over_write_svtxtrackmap)
				_trackmap->erase(iter->first);
		}
	}

#ifdef _DEBUG_
		cout<<__LINE__<<endl;
		cout << "about to clear tracks, if we're told to."<<endl;
		cout << "size: " << rf_phgf_tracks.size() << endl;
		cout << "_do_evt_display="<< (_do_evt_display?"true":"false") << endl;
#endif

	// Need to keep tracks if _do_evt_display
	if(!_do_evt_display) {
		rf_phgf_tracks.clear();
	}

#ifdef _DEBUG_
		cout<<__LINE__<<endl;
#endif
	/*!
	 * Fit track as primary track, This part need to be called after FillSvtxVertexMap
	 */
	if (_fit_primary_tracks && rave_vertices.size() > 0) {
		_primary_trackmap->empty();

		//FIXME figure out which vertex to use.
		SvtxVertex* vertex = NULL;
		if (_over_write_svtxvertexmap) {
			if (_vertexmap->size() > 0)
				vertex = _vertexmap->get(0);
		} else {
			if (_vertexmap_refit->size() > 0)
				vertex = _vertexmap_refit->get(0);
		}

		if (vertex) {
			for (SvtxTrackMap::ConstIter iter = _trackmap->begin();
					iter != _trackmap->end(); ++iter) {
				SvtxTrack* svtx_track = iter->second;
				if (!svtx_track)
					continue;
				if (!(svtx_track->get_pt() > _fit_min_pT))
					continue;
				/*!
				 * rf_phgf_track stands for Refit_PHGenFit_Track
				 */
				std::shared_ptr<PHGenFit::Track> rf_phgf_track = ReFitTrack(topNode, svtx_track,
						vertex);
				if (rf_phgf_track) {
//					//FIXME figure out which vertex to use.
//					SvtxVertex* vertex = NULL;
//					if (_vertexmap_refit->size() > 0)
//						vertex = _vertexmap_refit->get(0);
#ifdef _DEBUG_
			cout<<__LINE__<<endl;
#endif
					std::shared_ptr<SvtxTrack> rf_track = MakeSvtxTrack(svtx_track,
							rf_phgf_track, vertex);
#ifdef _DEBUG_
			cout<<__LINE__<<endl;
#endif
					//delete rf_phgf_track;
					if(!rf_track) {
#ifdef _DEBUG_
						LogDebug("!rf_track, continue.");
#endif
						continue;
					}
					_primary_trackmap->insert(rf_track.get());
				}
			}
		} else {
			LogError("No vertex in SvtxVertexMapRefit!");
		}
	}
#ifdef _DEBUG_
		cout<<__LINE__<<endl;
#endif
	for(genfit::GFRaveVertex *vertex: rave_vertices) {
		delete vertex;
	}
	rave_vertices.clear();

	if (_do_eval) {
		fill_eval_tree(topNode);
	}
#ifdef _DEBUG_
		cout<<__LINE__<<endl;
#endif
	return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * End
 */
int PHG4TrackKalmanFitter::End(PHCompositeNode *topNode) {

	if (_do_eval) {
		if(verbosity >= 1)
			cout << PHWHERE << " Writing to file: " << _eval_outname << endl;
		PHTFileServer::get().cd(_eval_outname);
		_eval_tree->Write();
		_cluster_eval_tree->Write();
		_kalman_extrapolation_eval_tree->Write();
		//put it on the same place as those trees, or make our own file per KalmanPatRec.
		//RCC add ntuple write.
	}

	if (_do_evt_display)
		_fitter->displayEvent();

	return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * dtor
 */
PHG4TrackKalmanFitter::~PHG4TrackKalmanFitter() {
	delete _fitter;
	delete _vertex_finder;
}

/*
 * fill_eval_tree():
 */
void PHG4TrackKalmanFitter::fill_eval_tree(PHCompositeNode *topNode) {
	//! Make sure to reset all the TTree variables before trying to set them.
	reset_eval_variables();

	int i = 0;
	for (PHG4TruthInfoContainer::ConstIterator itr =
			_truth_container->GetPrimaryParticleRange().first;
			itr != _truth_container->GetPrimaryParticleRange().second; ++itr)
		new ((*_tca_particlemap)[i++]) (PHG4Particlev2)(
				*dynamic_cast<PHG4Particlev2*>(itr->second));

	i = 0;
	for (PHG4TruthInfoContainer::ConstVtxIterator itr =
			_truth_container->GetPrimaryVtxRange().first;
			itr != _truth_container->GetPrimaryVtxRange().second; ++itr)
		new ((*_tca_vtxmap)[i++]) (PHG4VtxPointv1)(
				*dynamic_cast<PHG4VtxPointv1*>(itr->second));

	i = 0;
	for (SvtxTrackMap::ConstIter itr = _trackmap->begin();
			itr != _trackmap->end(); ++itr)
		new ((*_tca_trackmap)[i++]) (SvtxTrack_v1)(
				*dynamic_cast<SvtxTrack_v1*>(itr->second));

	i = 0;
	if (_vertexmap)
		for (SvtxVertexMap::ConstIter itr = _vertexmap->begin();
				itr != _vertexmap->end(); ++itr)
			new ((*_tca_vertexmap)[i++]) (SvtxVertex_v1)(
					*dynamic_cast<SvtxVertex_v1*>(itr->second));

	if (_trackmap_refit) {
		i = 0;
		for (SvtxTrackMap::ConstIter itr = _trackmap_refit->begin();
				itr != _trackmap_refit->end(); ++itr)
			new ((*_tca_trackmap_refit)[i++]) (SvtxTrack_v1)(
					*dynamic_cast<SvtxTrack_v1*>(itr->second));
	}

	if (_fit_primary_tracks) {
		i = 0;
		for (SvtxTrackMap::ConstIter itr = _primary_trackmap->begin();
				itr != _primary_trackmap->end(); ++itr)
			new ((*_tca_primtrackmap)[i++]) (SvtxTrack_v1)(
					*dynamic_cast<SvtxTrack_v1*>(itr->second));
	}

	if (_vertexmap_refit) {
		i = 0;
		for (SvtxVertexMap::ConstIter itr = _vertexmap_refit->begin();
				itr != _vertexmap_refit->end(); ++itr)
			new ((*_tca_vertexmap_refit)[i++]) (SvtxVertex_v1)(
					*dynamic_cast<SvtxVertex_v1*>(itr->second));
	}

	_eval_tree->Fill();

	return;
}

/*
 * init_eval_tree
 */
void PHG4TrackKalmanFitter::init_eval_tree() {
	if (!_tca_particlemap)
		_tca_particlemap = new TClonesArray("PHG4Particlev2");
	if (!_tca_vtxmap)
		_tca_vtxmap = new TClonesArray("PHG4VtxPointv1");

	if (!_tca_trackmap)
		_tca_trackmap = new TClonesArray("SvtxTrack_v1");
	if (!_tca_vertexmap)
		_tca_vertexmap = new TClonesArray("SvtxVertex_v1");
	if (!_tca_trackmap_refit)
		_tca_trackmap_refit = new TClonesArray("SvtxTrack_v1");
	if (_fit_primary_tracks)
		if (!_tca_primtrackmap)
			_tca_primtrackmap = new TClonesArray("SvtxTrack_v1");
	if (!_tca_vertexmap_refit)
		_tca_vertexmap_refit = new TClonesArray("SvtxVertex_v1");

	//! create TTree
	_eval_tree = new TTree("T", "PHG4TrackKalmanFitter Evaluation");

	_eval_tree->Branch("PrimaryParticle", _tca_particlemap);
	_eval_tree->Branch("TruthVtx", _tca_vtxmap);

	_eval_tree->Branch("SvtxTrack", _tca_trackmap);
	_eval_tree->Branch("SvtxVertex", _tca_vertexmap);
	_eval_tree->Branch("SvtxTrackRefit", _tca_trackmap_refit);
	if (_fit_primary_tracks)
		_eval_tree->Branch("PrimSvtxTrack", _tca_primtrackmap);
	_eval_tree->Branch("SvtxVertexRefit", _tca_vertexmap_refit);

	_cluster_eval_tree = new TTree("cluster_eval", "cluster eval tree");
	_cluster_eval_tree->Branch("x", &_cluster_eval_tree_x, "x/F");
	_cluster_eval_tree->Branch("y", &_cluster_eval_tree_y, "y/F");
	_cluster_eval_tree->Branch("z", &_cluster_eval_tree_z, "z/F");
	_cluster_eval_tree->Branch("gx", &_cluster_eval_tree_gx, "gx/F");
	_cluster_eval_tree->Branch("gy", &_cluster_eval_tree_gy, "gy/F");
	_cluster_eval_tree->Branch("gz", &_cluster_eval_tree_gz, "gz/F");

	//RCC add extrapolation branch to the output so we can check resolutions:
	_kalman_extrapolation_eval_tree = new TTree("kalman_extrapolation_eval","kalman extrapolation eval tree");
	_kalman_extrapolation_eval_tree->Branch("phi", &_kalman_extrapolation_eval_tree_phi, "phi/F");
	_kalman_extrapolation_eval_tree->Branch("z", &_kalman_extrapolation_eval_tree_z, "z/F");
	_kalman_extrapolation_eval_tree->Branch("r", &_kalman_extrapolation_eval_tree_r, "r/F");
	_kalman_extrapolation_eval_tree->Branch("okay", &_kalman_extrapolation_eval_tree_okay, "ok/O");
	_kalman_extrapolation_eval_tree->Branch("sigma_phi", &_kalman_extrapolation_eval_tree_sigma_rphi, "sigma_rphi/F");
	_kalman_extrapolation_eval_tree->Branch("sigma_z", &_kalman_extrapolation_eval_tree_sigma_z, "sigma_z/F");
	_kalman_extrapolation_eval_tree->Branch("sigma_r", &_kalman_extrapolation_eval_tree_sigma_r, "sigma_r/F");
	_kalman_extrapolation_eval_tree->Branch("sigma_phi", &_kalman_extrapolation_eval_tree_sigma_rphi_z, "sigma_phi_z/F");
	_kalman_extrapolation_eval_tree->Branch("sigma_z", &_kalman_extrapolation_eval_tree_sigma_z_r, "sigma_z_r/F");
	_kalman_extrapolation_eval_tree->Branch("sigma_r", &_kalman_extrapolation_eval_tree_sigma_r_rphi, "sigma_r_rphi/F");

	//before rotation:
	_kalman_extrapolation_eval_tree->Branch("sigma_x",&_kalman_extrapolation_eval_tree_covin_x);
	_kalman_extrapolation_eval_tree->Branch("sigma_y",&_kalman_extrapolation_eval_tree_covin_y);
	_kalman_extrapolation_eval_tree->Branch("sigma_z",&_kalman_extrapolation_eval_tree_covin_z);

	//linear data:	
	_kalman_extrapolation_eval_tree->Branch("lin_path",&_kalman_extrapolation_eval_tree_lin_pathlength);
	_kalman_extrapolation_eval_tree->Branch("lin_pos0_x",&_kalman_extrapolation_eval_tree_lin_pos0_x);
	_kalman_extrapolation_eval_tree->Branch("lin_pos0_y",&_kalman_extrapolation_eval_tree_lin_pos0_y);
	_kalman_extrapolation_eval_tree->Branch("lin_pos0_z",&_kalman_extrapolation_eval_tree_lin_pos0_z);
	_kalman_extrapolation_eval_tree->Branch("lin_pos1_x",&_kalman_extrapolation_eval_tree_lin_pos1_x);
	_kalman_extrapolation_eval_tree->Branch("lin_pos1_y",&_kalman_extrapolation_eval_tree_lin_pos1_y);
	_kalman_extrapolation_eval_tree->Branch("lin_pos1_z",&_kalman_extrapolation_eval_tree_lin_pos1_z);
	_kalman_extrapolation_eval_tree->Branch("lin_pos2_x",&_kalman_extrapolation_eval_tree_lin_pos2_x);
	_kalman_extrapolation_eval_tree->Branch("lin_pos2_y",&_kalman_extrapolation_eval_tree_lin_pos2_y);
	_kalman_extrapolation_eval_tree->Branch("lin_pos2_z",&_kalman_extrapolation_eval_tree_lin_pos2_z);
	_kalman_extrapolation_eval_tree->Branch("lin_phi",&_kalman_extrapolation_eval_tree_lin_phi);
	_kalman_extrapolation_eval_tree->Branch("lin_r",&_kalman_extrapolation_eval_tree_lin_r);
	_kalman_extrapolation_eval_tree->Branch("lin_z",&_kalman_extrapolation_eval_tree_lin_z);
	_kalman_extrapolation_eval_tree->Branch("lin_sig0_r",&_kalman_extrapolation_eval_tree_lin_sigma0_r);
	_kalman_extrapolation_eval_tree->Branch("lin_sig0_rphi",&_kalman_extrapolation_eval_tree_lin_sigma0_rphi);
	_kalman_extrapolation_eval_tree->Branch("lin_sig0_z",&_kalman_extrapolation_eval_tree_lin_sigma0_z);
	_kalman_extrapolation_eval_tree->Branch("lin_sig1_r",&_kalman_extrapolation_eval_tree_lin_sigma1_r);
	_kalman_extrapolation_eval_tree->Branch("lin_sig1_rphi",&_kalman_extrapolation_eval_tree_lin_sigma1_rphi);
	_kalman_extrapolation_eval_tree->Branch("lin_sig1_z",&_kalman_extrapolation_eval_tree_lin_sigma1_z);
	_kalman_extrapolation_eval_tree->Branch("lin_sig2_r",&_kalman_extrapolation_eval_tree_lin_sigma2_r);
	_kalman_extrapolation_eval_tree->Branch("lin_sig2_rphi",&_kalman_extrapolation_eval_tree_lin_sigma2_rphi);
	_kalman_extrapolation_eval_tree->Branch("lin_sig2_z",&_kalman_extrapolation_eval_tree_lin_sigma2_z);
}

/*
 * reset_eval_variables():
 *  Reset all the tree variables to their default values.
 *  Needs to be called at the start of every event
 */
void PHG4TrackKalmanFitter::reset_eval_variables() {
	_tca_particlemap->Clear();
	_tca_vtxmap->Clear();

	_tca_trackmap->Clear();
	_tca_vertexmap->Clear();
	_tca_trackmap_refit->Clear();
	if (_fit_primary_tracks)
		_tca_primtrackmap->Clear();
	_tca_vertexmap_refit->Clear();

	_cluster_eval_tree_x = WILD_FLOAT;
	_cluster_eval_tree_y = WILD_FLOAT;
	_cluster_eval_tree_z = WILD_FLOAT;
	_cluster_eval_tree_gx = WILD_FLOAT;
	_cluster_eval_tree_gy = WILD_FLOAT;
	_cluster_eval_tree_gz = WILD_FLOAT;

	_kalman_extrapolation_eval_tree_phi = WILD_FLOAT;
	 _kalman_extrapolation_eval_tree_z = WILD_FLOAT;
	 _kalman_extrapolation_eval_tree_r = WILD_FLOAT;
	 _kalman_extrapolation_eval_tree_okay = false;
	 _kalman_extrapolation_eval_tree_sigma_rphi = WILD_FLOAT;
	 _kalman_extrapolation_eval_tree_sigma_z = WILD_FLOAT;
	 _kalman_extrapolation_eval_tree_sigma_r = WILD_FLOAT;
	 _kalman_extrapolation_eval_tree_sigma_rphi_z = WILD_FLOAT;
	 _kalman_extrapolation_eval_tree_sigma_z_r = WILD_FLOAT;
	 _kalman_extrapolation_eval_tree_sigma_r_rphi = WILD_FLOAT;

	 //before rotation:
	_kalman_extrapolation_eval_tree_covin_x = WILD_FLOAT;
	_kalman_extrapolation_eval_tree_covin_y = WILD_FLOAT;
	_kalman_extrapolation_eval_tree_covin_z = WILD_FLOAT;

	//linear data:
	_kalman_extrapolation_eval_tree_lin_pos0_x = WILD_FLOAT;
	_kalman_extrapolation_eval_tree_lin_pos0_y = WILD_FLOAT;
	_kalman_extrapolation_eval_tree_lin_pos0_z = WILD_FLOAT;
	_kalman_extrapolation_eval_tree_lin_pos1_x = WILD_FLOAT;
	_kalman_extrapolation_eval_tree_lin_pos1_y = WILD_FLOAT;
	_kalman_extrapolation_eval_tree_lin_pos1_z = WILD_FLOAT;
	_kalman_extrapolation_eval_tree_lin_pos2_x = WILD_FLOAT;
	_kalman_extrapolation_eval_tree_lin_pos2_y = WILD_FLOAT;
	_kalman_extrapolation_eval_tree_lin_pos2_z = WILD_FLOAT;
	_kalman_extrapolation_eval_tree_lin_phi = WILD_FLOAT;
	_kalman_extrapolation_eval_tree_lin_r = WILD_FLOAT;
	_kalman_extrapolation_eval_tree_lin_z = WILD_FLOAT;
	_kalman_extrapolation_eval_tree_lin_sigma0_r = WILD_FLOAT;
	_kalman_extrapolation_eval_tree_lin_sigma0_rphi = WILD_FLOAT;
	_kalman_extrapolation_eval_tree_lin_sigma0_z = WILD_FLOAT;
	_kalman_extrapolation_eval_tree_lin_sigma1_r = WILD_FLOAT;
	_kalman_extrapolation_eval_tree_lin_sigma1_rphi = WILD_FLOAT;
	_kalman_extrapolation_eval_tree_lin_sigma1_z = WILD_FLOAT;
	_kalman_extrapolation_eval_tree_lin_sigma2_r = WILD_FLOAT;
	_kalman_extrapolation_eval_tree_lin_sigma2_rphi = WILD_FLOAT;
	_kalman_extrapolation_eval_tree_lin_sigma2_z = WILD_FLOAT;
}

int PHG4TrackKalmanFitter::CreateNodes(PHCompositeNode *topNode) {
	// create nodes...
	PHNodeIterator iter(topNode);

	PHCompositeNode *dstNode = static_cast<PHCompositeNode*>(iter.findFirst(
			"PHCompositeNode", "DST"));
	if (!dstNode) {
		cerr << PHWHERE << "DST Node missing, doing nothing." << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}
  PHNodeIterator iter_dst(dstNode);

	// Create the SVTX node
	PHCompositeNode* tb_node = dynamic_cast<PHCompositeNode*>(iter_dst.findFirst(
			"PHCompositeNode", "SVTX"));
	if (!tb_node) {
		tb_node = new PHCompositeNode("SVTX");
		dstNode->addNode(tb_node);
		if (verbosity > 0)
			cout << "SVTX node added" << endl;
	}

	if (!(_over_write_svtxtrackmap) || _output_mode == DebugMode) {
		_trackmap_refit = new SvtxTrackMap_v1;
		PHIODataNode<PHObject>* tracks_node = new PHIODataNode<PHObject>(
				_trackmap_refit, "SvtxTrackMapRefit", "PHObject");
		tb_node->addNode(tracks_node);
		if (verbosity > 0)
			cout << "Svtx/SvtxTrackMapRefit node added" << endl;
	}

	if (_fit_primary_tracks) {
		_primary_trackmap = new SvtxTrackMap_v1;
		PHIODataNode<PHObject>* primary_tracks_node =
				new PHIODataNode<PHObject>(_primary_trackmap, "PrimaryTrackMap",
						"PHObject");
		tb_node->addNode(primary_tracks_node);
		if (verbosity > 0)
			cout << "Svtx/PrimaryTrackMap node added" << endl;
	}

	if (!(_over_write_svtxvertexmap)) {
		_vertexmap_refit = new SvtxVertexMap_v1;
		PHIODataNode<PHObject>* vertexes_node = new PHIODataNode<PHObject>(
				_vertexmap_refit, "SvtxVertexMapRefit", "PHObject");
		tb_node->addNode(vertexes_node);
		if (verbosity > 0)
			cout << "Svtx/SvtxVertexMapRefit node added" << endl;
	} else if (!findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap")) {
		_vertexmap = new SvtxVertexMap_v1;
		PHIODataNode<PHObject>* vertexes_node = new PHIODataNode<PHObject>(
				_vertexmap, "SvtxVertexMap", "PHObject");
		tb_node->addNode(vertexes_node);
		if (verbosity > 0)
			cout << "Svtx/SvtxVertexMap node added" << endl;
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * GetNodes():
 *  Get all the all the required nodes off the node tree
 */
int PHG4TrackKalmanFitter::GetNodes(PHCompositeNode * topNode) {
	//DST objects
	//Truth container
	_truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode,
			"G4TruthInfo");
	if (!_truth_container && _event < 2) {
		cout << PHWHERE << " PHG4TruthInfoContainer node not found on node tree"
				<< endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	// Input Svtx Clusters
	_clustermap = findNode::getClass<SvtxClusterMap>(topNode, "SvtxClusterMap");
	if (!_clustermap && _event < 2) {
		cout << PHWHERE << " SvtxClusterMap node not found on node tree"
				<< endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	// Input Svtx Tracks
	_trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
	if (!_trackmap && _event < 2) {
		cout << PHWHERE << " SvtxClusterMap node not found on node tree"
				<< endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	// Input Svtx Vertices
	_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
	if (!_vertexmap && _event < 2) {
		cout << PHWHERE << " SvtxVertexrMap node not found on node tree"
				<< endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	// Output Svtx Tracks
	if (!(_over_write_svtxtrackmap) || _output_mode == DebugMode) {
		_trackmap_refit = findNode::getClass<SvtxTrackMap>(topNode,
				"SvtxTrackMapRefit");
		if (!_trackmap_refit && _event < 2) {
			cout << PHWHERE << " SvtxTrackMapRefit node not found on node tree"
					<< endl;
			return Fun4AllReturnCodes::ABORTEVENT;
		}
	}

	// Output Primary Svtx Tracks
	if (_fit_primary_tracks) {
		_primary_trackmap = findNode::getClass<SvtxTrackMap>(topNode,
				"PrimaryTrackMap");
		if (!_primary_trackmap && _event < 2) {
			cout << PHWHERE << " PrimaryTrackMap node not found on node tree"
					<< endl;
			return Fun4AllReturnCodes::ABORTEVENT;
		}
	}

	// Output Svtx Vertices
	if (!(_over_write_svtxvertexmap)) {
		_vertexmap_refit = findNode::getClass<SvtxVertexMap>(topNode,
				"SvtxVertexMapRefit");
		if (!_vertexmap_refit && _event < 2) {
			cout << PHWHERE << " SvtxVertexMapRefit node not found on node tree"
					<< endl;
			return Fun4AllReturnCodes::ABORTEVENT;
		}
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

//struct CompMeasurementByR {
//  bool operator() (PHGenFit::Measurement *m1,PHGenFit::Measurement *m2) {
//	  float x1 = m1->getMeasurement()
//
//	  return (i<j);}
//} myobject;

/*
 * fit track with SvtxTrack as input seed.
 * \param intrack Input SvtxTrack
 * \param invertex Input Vertex, if fit track as a primary vertex
 */
//PHGenFit::Track* PHG4TrackKalmanFitter::ReFitTrack(PHCompositeNode *topNode, const SvtxTrack* intrack,
std::shared_ptr<PHGenFit::Track> PHG4TrackKalmanFitter::ReFitTrack(PHCompositeNode *topNode, const SvtxTrack* intrack,
								   const SvtxVertex* invertex){//rcc temporary hack:, const bool use_svtx, const bool use_intt, const bool use_mvtx) {

  //rcc temporarily turning off the tpc:
  bool use_svtx=false;
  bool use_intt=true;
  bool use_maps=true;
	//std::shared_ptr<PHGenFit::Track> empty_track(NULL);

	if (!intrack) {
		cerr << PHWHERE << " Input SvtxTrack is NULL!" << endl;
		return NULL;
	}

	// get node containing the digitized hits
	SvtxHitMap* hitsmap =  findNode::getClass<SvtxHitMap>(topNode, "SvtxHitMap");
	if (!hitsmap) {
		cout << PHWHERE << "ERROR: Can't find node SvtxHitMap" << endl;
		return NULL;
	}

	PHG4CellContainer* cells_svtx = findNode::getClass<PHG4CellContainer>(topNode,
			"G4CELL_SVTX");

	PHG4CellContainer* cells_intt = findNode::getClass<PHG4CellContainer>(
			topNode, "G4CELL_SILICON_TRACKER");

	PHG4CellContainer* cells_maps = findNode::getClass<PHG4CellContainer>(
			topNode, "G4CELL_MAPS");

	if (!cells_svtx and !cells_intt and !cells_maps) {
		if (verbosity >= 0) {
			LogError("No PHG4CellContainer found!");
		}
		return nullptr;
	}

	PHG4CylinderGeomContainer* geom_container_intt = findNode::getClass<
			PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_SILICON_TRACKER");

	PHG4CylinderGeomContainer* geom_container_maps = findNode::getClass<
			PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MAPS");

	if (!cells_svtx && !cells_maps && !cells_intt) {
		cout << PHWHERE << "ERROR: Can't find any cell node!" << endl;
		return NULL;
	}

	// prepare seed
	TVector3 seed_mom(100, 0, 0);
	TVector3 seed_pos(0, 0, 0);
	TMatrixDSym seed_cov(6);
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			seed_cov[i][j] = 100.;
		}
	}

	// Create measurements
	std::vector<PHGenFit::Measurement*> measurements;

	/*!
	 * if fit track as a primary track
	 */

//	if(invertex and verbosity >= 2)
//	{
//		LogDebug(invertex->size_tracks());
//		LogDebug(invertex->get_chisq());
//		LogDebug(invertex->get_ndof());
//		for (unsigned int i = 0; i < 3; i++)
//			for (unsigned int j = 0; j < 3; j++)
//			{
//				LogDebug(invertex->get_error(i,j));
//			}
//
//	}
	/*!
	 *
	 */
#if _DEBUG_MODE_ == 1
	if (invertex
//			and invertex->size_tracks() == 1
	) {
		TRandom3 rand(0);
		double dxy = 0.0007;	//7 um
		double dz = 0.003;//30 um

		TVector3 pos(invertex->get_x(), invertex->get_y(), invertex->get_z());
		TMatrixDSym cov(3);

		// Use smeared position instead of reco'd one.
		double x = rand.Gaus(0, dxy);
		double y = rand.Gaus(0, dxy);
		double z = rand.Gaus(0, dz);
		pos.SetXYZ(x, y, z);

		for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
		cov[i][j] = 0;

		cov[0][0] = dxy * dxy;
		cov[1][1] = dxy * dxy;
		cov[2][2] =dz * dz;

		PHGenFit::Measurement* meas = new PHGenFit::SpacepointMeasurement(
				pos, cov);
		measurements.push_back(meas);
	}
#else

	//! 1000 is a arbitrary number for now
	const double vertex_chi2_over_dnf_cut = 1000;
	const double vertex_cov_element_cut = 10000; //arbitrary cut cm*cm

	if (invertex and invertex->size_tracks() > 1
			and invertex->get_chisq() / invertex->get_ndof()
					< vertex_chi2_over_dnf_cut) {
		TVector3 pos(invertex->get_x(), invertex->get_y(), invertex->get_z());
		TMatrixDSym cov(3);
		cov.Zero();
		bool is_vertex_cov_sane = true;
		for (unsigned int i = 0; i < 3; i++)
			for (unsigned int j = 0; j < 3; j++) {

				cov(i, j) = invertex->get_error(i, j);

				if (i == j) {
					if (!(invertex->get_error(i, j) > 0
							and invertex->get_error(i, j)
									< vertex_cov_element_cut))
						is_vertex_cov_sane = false;
				}
			}

		if (is_vertex_cov_sane) {
			PHGenFit::Measurement* meas = new PHGenFit::SpacepointMeasurement(
					pos, cov);
			measurements.push_back(meas);
//			if(verbosity >= 2)
//			{
//				meas->getMeasurement()->Print();
//			}
		}
	}
#endif

	// sort clusters with radius before fitting

	std::map<float, unsigned int> m_r_cluster_id;
	for (auto iter = intrack->begin_clusters();
			iter != intrack->end_clusters(); ++iter) {
		unsigned int cluster_id = *iter;
		SvtxCluster* cluster = _clustermap->get(cluster_id);
		float x = cluster->get_x();
		float y = cluster->get_y();
		float r = sqrt(x*x+y*y);
		m_r_cluster_id.insert(std::pair<float, unsigned int>(r, cluster_id));
	}

	for (auto iter = m_r_cluster_id.begin();
			iter != m_r_cluster_id.end();
			++iter) {

		unsigned int cluster_id = iter->second;
		SvtxCluster* cluster = _clustermap->get(cluster_id);
		if (!cluster) {
			LogError("No cluster Found!");
			continue;
		}

#ifdef _DEBUG_
		cout
		<< __LINE__
		<<": ID: " << cluster_id
		<<": layer: " << cluster->get_layer()
		<<endl;
#endif

		TVector3 pos(cluster->get_x(), cluster->get_y(), cluster->get_z());

		// DEBUG: BEGIN
		if (_do_eval) {
			PHG4HitContainer* phg4hits_svtx = findNode::getClass<
					PHG4HitContainer>(topNode, "G4HIT_SVTX");

			PHG4HitContainer* phg4hits_intt = findNode::getClass<
					PHG4HitContainer>(topNode, "G4HIT_SILICON_TRACKER");

			PHG4HitContainer* phg4hits_maps = findNode::getClass<
					PHG4HitContainer>(topNode, "G4HIT_MAPS");

			if (!phg4hits_svtx and !phg4hits_intt and !phg4hits_maps) {
				if (verbosity >= 0) {
					LogError("No PHG4HitContainer found!");
				}
				continue;
			}

			SvtxHit* svtxhit = hitsmap->find(*cluster->begin_hits())->second;

			PHG4Cell* cell = nullptr;
			if(cells_svtx) cell = cells_svtx->findCell(svtxhit->get_cellid());
			if(!cell && cells_intt) cell = cells_intt->findCell(svtxhit->get_cellid());
			if(!cell && cells_maps) cell = cells_maps->findCell(svtxhit->get_cellid());
			if(!cell){
				if(verbosity>=0)
					LogError("!cell");
				continue;
			}

			PHG4Hit *phg4hit = nullptr;
			if(phg4hits_svtx) phg4hit = phg4hits_svtx->findHit(cell->get_g4hits().first->first);
			if(!phg4hit and phg4hits_intt) phg4hit = phg4hits_intt->findHit(cell->get_g4hits().first->first);
			if(!phg4hit and phg4hits_maps) phg4hit = phg4hits_maps->findHit(cell->get_g4hits().first->first);

			if (!phg4hit) {
				if (verbosity >= 0)
					LogError("!phg4hit");
				continue;
			}

			TVector3 phg4hit_position(phg4hit->get_avg_x(),
					phg4hit->get_avg_y(), phg4hit->get_avg_z());
			TVector3 cluster_position(cluster->get_x(), cluster->get_y(),
					cluster->get_z());

			_cluster_eval_tree_x = cluster_position.X();
			_cluster_eval_tree_y = cluster_position.Y();
			_cluster_eval_tree_z = cluster_position.Z();
			_cluster_eval_tree_gx = phg4hit_position.X();
			_cluster_eval_tree_gy = phg4hit_position.Y();
			_cluster_eval_tree_gz = phg4hit_position.Z();

			_cluster_eval_tree->Fill();
		}

//		if (phg4hit_position.Perp() > 30.) {
//			pos.SetXYZ(phg4hit_position.X(), phg4hit_position.Y(),phg4hit_position.Z()); //DEBUG
//			//pos.SetPerp(phg4hit_position.Perp());
//			//pos.SetPhi(TMath::ATan2(phg4hit_position.Y(),phg4hit_position.X()));
//		}
//
//		if(phg4hit->get_trkid()!=1) {
//			continue;
//		}
		// DEBUG: END



		//TODO use u, v explicitly?
		TVector3 n(cluster->get_x(), cluster->get_y(), 0);

		unsigned int begin_hit_id = *(cluster->begin_hits());
		//LogDebug(begin_hit_id);
		SvtxHit* svtxhit = hitsmap->find(begin_hit_id)->second;
		//LogDebug(svtxhit->get_cellid());

		PHG4Cell* cell_svtx = nullptr;
		PHG4Cell* cell_intt = nullptr;
		PHG4Cell* cell_maps = nullptr;

		if(cells_svtx) cell_svtx = cells_svtx->findCell(svtxhit->get_cellid());
		if(cells_intt) cell_intt = cells_intt->findCell(svtxhit->get_cellid());
		if(cells_maps) cell_maps = cells_maps->findCell(svtxhit->get_cellid());
		if(!(cell_svtx or cell_intt or cell_maps)){
			if(verbosity>=0)
				LogError("!(cell_svtx or cell_intt or cell_maps)");
			continue;
		}


		//if the detector is turned off by the bools, avoid adding it to the measurement vector by short circuiting here (check variable name conventions! RCC)
		if ( (cell_maps and !use_maps) or (cell_intt and !use_intt) or (cell_svtx and !use_svtx) ) continue;

		//RCC moved the seed updating to not update unless we're including the hits in the measurement set.
		seed_mom.SetPhi(pos.Phi());
		seed_mom.SetTheta(pos.Theta());
		
		//17.4, 17.4, 17.4, 14.0, 14.0, 12.0, 11.5
		//float phi_tilt[7] = {0.304, 0.304, 0.304, 0.244, 0.244, 0.209, 0.201};

		unsigned int layer = cluster->get_layer();
		//std::cout << "cluster layer: " << layer << std::endl;
		if (cell_maps) {
			PHG4Cell* cell = cell_maps;

			int stave_index = cell->get_stave_index();
			int half_stave_index = cell->get_half_stave_index();
			int module_index = cell->get_module_index();
			int chip_index = cell->get_chip_index();

			double ladder_location[3] = { 0.0, 0.0, 0.0 };
			PHG4CylinderGeom_MAPS *geom =
					(PHG4CylinderGeom_MAPS*) geom_container_maps->GetLayerGeom(
							layer);
			// returns the center of the sensor in world coordinates - used to get the ladder phi location
			geom->find_sensor_center(stave_index, half_stave_index,
					module_index, chip_index, ladder_location);
			//n.Print();
			n.SetXYZ(ladder_location[0], ladder_location[1], 0);
			n.RotateZ(geom->get_stave_phi_tilt());
			//n.Print();
		} else if (cell_intt) {
			PHG4Cell* cell = cell_intt;
			PHG4CylinderGeom_Siladders* geom =
					(PHG4CylinderGeom_Siladders*) geom_container_intt->GetLayerGeom(
							layer);
			double hit_location[3] = { 0.0, 0.0, 0.0 };
			geom->find_segment_center(cell->get_ladder_z_index(),
					cell->get_ladder_phi_index(), hit_location);

			n.SetXYZ(hit_location[0], hit_location[1], 0);
			n.RotateZ(geom->get_strip_phi_tilt());
		}

		PHGenFit::Measurement* meas = new PHGenFit::PlanarMeasurement(pos, n,
				cluster->get_rphi_error(), cluster->get_z_error());

//		TMatrixF cov_uvn(3,3);
//		TMatrixF cov_xyz(3,3);
//
//		double er2 = 2500E-8/12.;
//		if(pos.Perp() >= 30.) er2 = 25./36./12.;
//
//		// u: r; v: r-phi; n: z
//		cov_uvn[0][0] = er2;
//		cov_uvn[0][1] = 0.;
//		cov_uvn[0][2] = 0.;
//
//		cov_uvn[1][0] = 0.;
//		cov_uvn[1][1] = cluster->get_rphi_error()*cluster->get_rphi_error();
//		cov_uvn[1][2] = 0.;
//
//		cov_uvn[2][0] = 0.;
//		cov_uvn[2][1] = 0.;
//		cov_uvn[2][2] = cluster->get_z_error()*cluster->get_z_error();
//		cov_uvn.Print();
//
//		TVector3 u(cluster->get_x(), cluster->get_y(), 0);
//		u = u.Unit();
//		TVector3 n(0, 0, 1);
//		TVector3 v = n.Cross(u);
//
//		TMatrixF R = get_rotation_matrix(TVector3(1, 0, 0), TVector3(0, 1, 0),
//				TVector3(0, 0, 1), u, v, n);
//
//		R.Print();
//		TMatrixF R_T(3,3);
//		R_T.Transpose(R);
//		cov_xyz = R * cov_uvn * R_T;
//		cov_xyz.Print();
//		TMatrixDSym cov_xyz_dsym(3);
//		for(int i=0;i<3;i++)
//			for(int j=0;j<3;j++)
//				cov_xyz_dsym[i][j] = cov_xyz[i][j];
//		cov_xyz_dsym.Print();
//
//		PHGenFit::Measurement* meas = new PHGenFit::SpacepointMeasurement(pos, cov_xyz_dsym);

		measurements.push_back(meas);
	}

	/*!
	 * mu+:	-13
	 * mu-:	13
	 * pi+:	211
	 * pi-:	-211
	 * e-:	11
	 * e+:	-11
	 */
	//TODO Add multiple TrackRep choices.
	//int pid = 211;
	genfit::AbsTrackRep* rep = new genfit::RKTrackRep(_primary_pid_guess);
	std::shared_ptr<PHGenFit::Track> track(new PHGenFit::Track(rep, seed_pos, seed_mom,
			seed_cov));

	track->addMeasurements(measurements);

	/*!
	 *  Fit the track
	 *  ret code 0 means 0 error or good status
	 */
	if (_fitter->processTrack(track.get(), false) != 0) {
		if (verbosity >= 1)
			LogWarning("Track fitting failed");
		//delete track;
		return NULL;
	}

	return track;
}

/*
 * Make SvtxTrack from PHGenFit::Track and SvtxTrack
 */
//SvtxTrack* PHG4TrackKalmanFitter::MakeSvtxTrack(const SvtxTrack* svtx_track,
std::shared_ptr<SvtxTrack> PHG4TrackKalmanFitter::MakeSvtxTrack(const SvtxTrack* svtx_track,
		const std::shared_ptr<PHGenFit::Track>& phgf_track, const SvtxVertex* vertex) {


	double chi2 = phgf_track->get_chi2();
	double ndf = phgf_track->get_ndf();

	TVector3 vertex_position(0, 0, 0);
	TMatrixF vertex_cov(3,3);
	double dvr2 = 0;
	double dvz2 = 0;

	if(_use_truth_vertex) {
		PHG4VtxPoint* first_point = _truth_container->GetPrimaryVtx(_truth_container->GetPrimaryVertexIndex());
		vertex_position.SetXYZ(first_point->get_x(), first_point->get_y(), first_point->get_z());
		if(verbosity > 1) {
			cout<<"Using: truth vertex: {" << vertex_position.X() << ", " << vertex_position.Y() << ", " << vertex_position.Z() << "} " <<endl;
		}
	} else if (vertex) {
		vertex_position.SetXYZ(vertex->get_x(), vertex->get_y(),
				vertex->get_z());
		dvr2 = vertex->get_error(0, 0) + vertex->get_error(1, 1);
		dvz2 = vertex->get_error(2, 2);

		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
				vertex_cov[i][j] = vertex->get_error(i,j);
	}

	//genfit::MeasuredStateOnPlane* gf_state_beam_line_ca = NULL;
	std::shared_ptr<genfit::MeasuredStateOnPlane> gf_state_beam_line_ca = NULL;
	try {
		gf_state_beam_line_ca = std::shared_ptr<genfit::MeasuredStateOnPlane>(phgf_track->extrapolateToLine(vertex_position,
				TVector3(0., 0., 1.)));
	} catch (...) {
		if (verbosity >= 2)
			LogWarning("extrapolateToLine failed!");
	}
	if(!gf_state_beam_line_ca) return NULL;

	/*!
	 *  1/p, u'/z', v'/z', u, v
	 *  u is defined as momentum X beam line at POCA of the beam line
	 *  v is alone the beam line
	 *  so u is the dca2d direction
	 */

	double u = gf_state_beam_line_ca->getState()[3];
	double v = gf_state_beam_line_ca->getState()[4];

	double du2 = gf_state_beam_line_ca->getCov()[3][3];
	double dv2 = gf_state_beam_line_ca->getCov()[4][4];

	//delete gf_state_beam_line_ca;

	//const SvtxTrack_v1* temp_track = static_cast<const SvtxTrack_v1*> (svtx_track);
//	SvtxTrack_v1* out_track = new SvtxTrack_v1(
//			*static_cast<const SvtxTrack_v1*>(svtx_track));
	std::shared_ptr < SvtxTrack_v1 > out_track = std::shared_ptr < SvtxTrack_v1
			> (new SvtxTrack_v1(*static_cast<const SvtxTrack_v1*>(svtx_track)));

	out_track->set_dca2d(u);
	out_track->set_dca2d_error(sqrt(du2 + dvr2));

	std::shared_ptr<genfit::MeasuredStateOnPlane> gf_state_vertex_ca = NULL;
	try {
		gf_state_vertex_ca = std::shared_ptr < genfit::MeasuredStateOnPlane
				> (phgf_track->extrapolateToPoint(vertex_position));
	} catch (...) {
		if (verbosity >= 2)
			LogWarning("extrapolateToPoint failed!");
	}
	if (!gf_state_vertex_ca) {
		//delete out_track;
		return NULL;
	}

	TVector3 mom = gf_state_vertex_ca->getMom();
	TVector3 pos = gf_state_vertex_ca->getPos();
	TMatrixDSym cov = gf_state_vertex_ca->get6DCov();

//	genfit::MeasuredStateOnPlane* gf_state_vertex_ca =
//			phgf_track->extrapolateToLine(vertex_position,
//					TVector3(0., 0., 1.));

	u = gf_state_vertex_ca->getState()[3];
	v = gf_state_vertex_ca->getState()[4];

	du2 = gf_state_vertex_ca->getCov()[3][3];
	dv2 = gf_state_vertex_ca->getCov()[4][4];


	double dca3d = sqrt(u * u + v * v);
	double dca3d_error = sqrt(du2 + dv2 + dvr2 + dvz2);

	out_track->set_dca(dca3d);
	out_track->set_dca_error(dca3d_error);

	/*!
	 * dca3d_xy, dca3d_z
	 */


	/*
	// Rotate from u,v,n to r: n X Z, Z': n X r, n using 5D state/cov
	// commented on 2017-10-09

	TMatrixF pos_in(3,1);
	TMatrixF cov_in(3,3);
	pos_in[0][0] = gf_state_vertex_ca->getState()[3];
	pos_in[1][0] = gf_state_vertex_ca->getState()[4];
	pos_in[2][0] = 0.;

	cov_in[0][0] = gf_state_vertex_ca->getCov()[3][3];
	cov_in[0][1] = gf_state_vertex_ca->getCov()[3][4];
	cov_in[0][2] = 0.;
	cov_in[1][0] = gf_state_vertex_ca->getCov()[4][3];
	cov_in[1][1] = gf_state_vertex_ca->getCov()[4][4];
	cov_in[1][2] = 0.;
	cov_in[2][0] = 0.;
	cov_in[2][1] = 0.;
	cov_in[2][2] = 0.;

	TMatrixF pos_out(3,1);
	TMatrixF cov_out(3,3);

	TVector3 vu = gf_state_vertex_ca->getPlane().get()->getU();
	TVector3 vv = gf_state_vertex_ca->getPlane().get()->getV();
	TVector3 vn = vu.Cross(vv);

	pos_cov_uvn_to_rz(vu, vv, vn, pos_in, cov_in, pos_out, cov_out);

	//! vertex cov in (u',v',n')
	TMatrixF vertex_cov_out(3,3);

	get_vertex_error_uvn(vu,vv,vn, vertex_cov, vertex_cov_out);

	float dca3d_xy = pos_out[0][0];
	float dca3d_z  = pos_out[1][0];

	float dca3d_xy_error = sqrt(cov_out[0][0] + vertex_cov_out[0][0]);
	float dca3d_z_error  = sqrt(cov_out[1][1] + vertex_cov_out[1][1]);

		//Begin DEBUG
//	LogDebug("rotation debug---------- ");
//	gf_state_vertex_ca->Print();
//	LogDebug("dca rotation---------- ");
//	pos_out = pos_in;
//	cov_out = cov_in;
//	pos_in.Print();
//	cov_in.Print();
//	pos_out.Print();
//	cov_out.Print();
//	cout
//		<<"dca3d_xy: "<<dca3d_xy <<" +- "<<dca3d_xy_error*dca3d_xy_error
//		<<"; dca3d_z: "<<dca3d_z<<" +- "<< dca3d_z_error*dca3d_z_error
//		<<"\n";
//	gf_state_vertex_ca->get6DCov().Print();
//	LogDebug("vertex rotation---------- ");
//	vertex_position.Print();
//	vertex_cov.Print();
//	vertex_cov_out.Print();
	//End DEBUG
	*/

	//
	// in: X, Y, Z; out; r: n X Z, Z X r, Z

	float dca3d_xy = NAN;
	float dca3d_z  = NAN;
	float dca3d_xy_error = NAN;
	float dca3d_z_error  = NAN;

	try{
		TMatrixF pos_in(3,1);
		TMatrixF cov_in(3,3);
		TMatrixF pos_out(3,1);
		TMatrixF cov_out(3,3);

		TVectorD state6(6); // pos(3), mom(3)
		TMatrixDSym cov6(6,6); //

		gf_state_vertex_ca->get6DStateCov(state6, cov6);

		TVector3 vn(state6[3], state6[4], state6[5]);

		// mean of two multivariate gaussians Pos - Vertex
		pos_in[0][0] = state6[0] - vertex_position.X();
		pos_in[1][0] = state6[1] - vertex_position.Y();
		pos_in[2][0] = state6[2] - vertex_position.Z();


		for(int i=0;i<3;++i){
			for(int j=0;j<3;++j){
				cov_in[i][j] = cov6[i][j] + vertex_cov[i][j];
			}
		}

		pos_cov_XYZ_to_RZ(vn, pos_in, cov_in, pos_out, cov_out);

		dca3d_xy = pos_out[0][0];
		dca3d_z  = pos_out[2][0];
		dca3d_xy_error = sqrt(cov_out[0][0]);
		dca3d_z_error  = sqrt(cov_out[2][2]);

#ifdef _DEBUG_
		cout<<__LINE__<<": Vertex: ----------------"<<endl;
		vertex_position.Print();
		vertex_cov.Print();

		cout<<__LINE__<<": State: ----------------"<<endl;
		state6.Print();
		cov6.Print();

		cout<<__LINE__<<": Mean: ----------------"<<endl;
		pos_in.Print();
		cout<<"===>"<<endl;
		pos_out.Print();

		cout<<__LINE__<<": Cov: ----------------"<<endl;
		cov_in.Print();
		cout<<"===>"<<endl;
		cov_out.Print();

		cout<<endl;
#endif

	} catch (...) {
		if (verbosity > 0)
			LogWarning("DCA calculationfailed!");
	}

	out_track->set_dca3d_xy(dca3d_xy);
	out_track->set_dca3d_z(dca3d_z);
	out_track->set_dca3d_xy_error(dca3d_xy_error);
	out_track->set_dca3d_z_error(dca3d_z_error);

	//if(gf_state_vertex_ca) delete gf_state_vertex_ca;

	out_track->set_chisq(chi2);
	out_track->set_ndf(ndf);
	out_track->set_charge(phgf_track->get_charge());

	out_track->set_px(mom.Px());
	out_track->set_py(mom.Py());
	out_track->set_pz(mom.Pz());

	out_track->set_x(pos.X());
	out_track->set_y(pos.Y());
	out_track->set_z(pos.Z());

	for (int i = 0; i < 6; i++) {
		for (int j = i; j < 6; j++) {
			out_track->set_error(i, j, cov[i][j]);
		}
	}

//	for (SvtxTrack::ConstClusterIter iter = svtx_track->begin_clusters();
//			iter != svtx_track->end_clusters(); ++iter) {
//		unsigned int cluster_id = *iter;
//		SvtxCluster* cluster = _clustermap->get(cluster_id);
//		if (!cluster) {
//			LogError("No cluster Found!");
//			continue;
//		}
//		//cluster->identify(); //DEBUG
//
//		//unsigned int l = cluster->get_layer();
//
//		TVector3 pos(cluster->get_x(), cluster->get_y(), cluster->get_z());
//
//		double radius = pos.Pt();
//
//		std::shared_ptr<genfit::MeasuredStateOnPlane> gf_state = NULL;
//		try {
//			gf_state = std::shared_ptr < genfit::MeasuredStateOnPlane
//					> (phgf_track->extrapolateToCylinder(radius,
//							TVector3(0, 0, 0), TVector3(0, 0, 1), 0));
//		} catch (...) {
//			if (verbosity >= 2)
//				LogWarning("Exrapolation failed!");
//		}
//		if (!gf_state) {
//			if (verbosity > 1)
//				LogWarning("Exrapolation failed!");
//			continue;
//		}
//
//		//SvtxTrackState* state = new SvtxTrackState_v1(radius);
//		std::shared_ptr<SvtxTrackState> state = std::shared_ptr<SvtxTrackState> (new SvtxTrackState_v1(radius));
//		state->set_x(gf_state->getPos().x());
//		state->set_y(gf_state->getPos().y());
//		state->set_z(gf_state->getPos().z());
//
//		state->set_px(gf_state->getMom().x());
//		state->set_py(gf_state->getMom().y());
//		state->set_pz(gf_state->getMom().z());
//
//		//gf_state->getCov().Print();
//
//		for (int i = 0; i < 6; i++) {
//			for (int j = i; j < 6; j++) {
//				state->set_error(i, j, gf_state->get6DCov()[i][j]);
//			}
//		}
//
//		out_track->insert_state(state.get());
//
//#ifdef _DEBUG_
//		cout
//		<<__LINE__
//		<<": " << radius <<" => "
//		<<sqrt(state->get_x()*state->get_x() + state->get_y()*state->get_y())
//		<<endl;
//#endif
//	}

#ifdef _DEBUG_
	cout << __LINE__ << endl;
#endif

	const genfit::Track *gftrack = phgf_track->getGenFitTrack();
	const genfit::AbsTrackRep *rep = gftrack->getCardinalRep();
	for(unsigned int id = 0; id< gftrack->getNumPointsWithMeasurement();++id) {
#ifdef _DEBUG_
	  cout << __LINE__ << " Trying to get trackpoint with id="<<id<<" of " << gftrack->getNumPointsWithMeasurement() << " rep ptr="<< (void*)rep << endl;
#endif
		genfit::TrackPoint *trpoint = gftrack->getPointWithMeasurementAndFitterInfo(id, gftrack->getCardinalRep());

		if(!trpoint) {
			if (verbosity > 1)
				LogWarning("!trpoint");
			continue;
		}
#ifdef _DEBUG_
	cout << __LINE__ << endl;
#endif
		genfit::KalmanFitterInfo* kfi = static_cast<genfit::KalmanFitterInfo*>( trpoint->getFitterInfo(rep) );
		if(!kfi) {
			if (verbosity > 1)
				LogWarning("!kfi");
			continue;
		}
#ifdef _DEBUG_
	cout << __LINE__ << endl;
#endif
		std::shared_ptr<const genfit::MeasuredStateOnPlane> gf_state = NULL;
		try {
			//gf_state = std::shared_ptr <genfit::MeasuredStateOnPlane> (const_cast<genfit::MeasuredStateOnPlane*> (&(kfi->getFittedState(true))));
			const genfit::MeasuredStateOnPlane* temp_state = &(kfi->getFittedState(true));
			gf_state = std::shared_ptr <genfit::MeasuredStateOnPlane> (new genfit::MeasuredStateOnPlane(*temp_state));
		} catch (...) {
			if (verbosity > 1)
				LogWarning("Exrapolation failed!");
		}
		if (!gf_state) {
			if (verbosity > 1)
				LogWarning("Exrapolation failed!");
			continue;
		}
		genfit::MeasuredStateOnPlane temp;
#ifdef _DEBUG_
	  cout << __LINE__ << " Trying to extrapolate"<<id<<endl;
#endif
	  float pathlength = -phgf_track->extrapolateToPoint(temp,vertex_position,id);

		std::shared_ptr<SvtxTrackState> state = std::shared_ptr<SvtxTrackState> (new SvtxTrackState_v1(pathlength));
		state->set_x(gf_state->getPos().x());
		state->set_y(gf_state->getPos().y());
		state->set_z(gf_state->getPos().z());

		state->set_px(gf_state->getMom().x());
		state->set_py(gf_state->getMom().y());
		state->set_pz(gf_state->getMom().z());

		//gf_state->getCov().Print();

		for (int i = 0; i < 6; i++) {
			for (int j = i; j < 6; j++) {
				state->set_error(i, j, gf_state->get6DCov()[i][j]);
			}
		}

		out_track->insert_state(state.get());

		#ifdef _DEBUG_
		cout
		<<__LINE__
		<<": " << id
		<<": " << pathlength <<" => "
		<<sqrt(state->get_x()*state->get_x() + state->get_y()*state->get_y())
		<<endl;
#endif
	}

	return out_track;
}

/*
 * Fill SvtxVertexMap from GFRaveVertexes and Tracks
 */
bool PHG4TrackKalmanFitter::FillSvtxVertexMap(
		const std::vector<genfit::GFRaveVertex*>& rave_vertices,
		const std::vector<genfit::Track*>& gf_tracks) {

	if(_over_write_svtxvertexmap){
		_vertexmap->clear();
	}

	for (unsigned int ivtx = 0; ivtx < rave_vertices.size(); ++ivtx) {
		genfit::GFRaveVertex* rave_vtx = rave_vertices[ivtx];

		if (!rave_vtx) {
			cerr << PHWHERE << endl;
			return false;
		}

		std::shared_ptr<SvtxVertex> svtx_vtx(new SvtxVertex_v1());

		svtx_vtx->set_chisq(rave_vtx->getChi2());
		svtx_vtx->set_ndof(rave_vtx->getNdf());
		svtx_vtx->set_position(0, rave_vtx->getPos().X());
		svtx_vtx->set_position(1, rave_vtx->getPos().Y());
		svtx_vtx->set_position(2, rave_vtx->getPos().Z());

		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				svtx_vtx->set_error(i, j, rave_vtx->getCov()[i][j]);

		for (unsigned int i = 0; i < rave_vtx->getNTracks(); i++) {
			//TODO Assume id's are sync'ed between _trackmap_refit and gf_tracks, need to change?
			const genfit::Track* rave_track =
					rave_vtx->getParameters(i)->getTrack();
			for (unsigned int j = 0; j < gf_tracks.size(); j++) {
				if (rave_track == gf_tracks[j])
					svtx_vtx->insert_track(j);
			}
		}

		if (_over_write_svtxvertexmap) {
			if (_vertexmap) {
				_vertexmap->insert(svtx_vtx.get());
			} else {
				LogError("!_vertexmap");
			}
		} else {
			if (_vertexmap_refit) {
				_vertexmap_refit->insert(svtx_vtx.get());
			} else {
				LogError("!_vertexmap_refit");
			}
		}

//		if (verbosity >= 2) {
//			cout << PHWHERE << endl;
//			svtx_vtx->Print();
//			_vertexmap_refit->Print();
//		}

		//delete svtx_vtx;
	} //loop over RAVE vertices

	return true;
}

//bool PHG4TrackKalmanFitter::pos_cov_uvn_to_rz(const TVector3 u, const TVector3 v,
//		const TVector3 n, const TMatrixF pos_in, const TMatrixF cov_in,
//		TMatrixF& pos_out, TMatrixF& cov_out) const {
//
//	if(pos_in.GetNcols() != 1 || pos_in.GetNrows() != 3) {
//		if(verbosity > 0) LogWarning("pos_in.GetNcols() != 1 || pos_in.GetNrows() != 3");
//		return false;
//	}
//
//	if(cov_in.GetNcols() != 3 || cov_in.GetNrows() != 3) {
//		if(verbosity > 0) LogWarning("cov_in.GetNcols() != 3 || cov_in.GetNrows() != 3");
//		return false;
//	}
//
//	TVector3 up = TVector3(0., 0., 1.).Cross(n);
//	if(up.Mag() < 0.00001){
//		if(verbosity > 0) LogWarning("n is parallel to z");
//		return false;
//	}
//
//	TMatrixF R(3, 3);
//	TMatrixF R_inv(3,3);
//	TMatrixF R_inv_T(3,3);
//
//	try {
//		TMatrixF ROT1(3, 3);
//		TMatrixF ROT2(3, 3);
//		TMatrixF ROT3(3, 3);
//
//		// rotate n along z to xz plane
//		float phi = -TMath::ATan2(n.Y(), n.X());
//		ROT1[0][0] = cos(phi);
//		ROT1[0][1] = -sin(phi);
//		ROT1[0][2] = 0;
//		ROT1[1][0] = sin(phi);
//		ROT1[1][1] = cos(phi);
//		ROT1[1][2] = 0;
//		ROT1[2][0] = 0;
//		ROT1[2][1] = 0;
//		ROT1[2][2] = 1;
//
//		// rotate n along y to z
//		TVector3 n1(n);
//		n1.RotateZ(phi);
//		float theta = -TMath::ATan2(n1.X(), n1.Z());
//		ROT2[0][0] = cos(theta);
//		ROT2[0][1] = 0;
//		ROT2[0][2] = sin(theta);
//		ROT2[1][0] = 0;
//		ROT2[1][1] = 1;
//		ROT2[1][2] = 0;
//		ROT2[2][0] = -sin(theta);
//		ROT2[2][1] = 0;
//		ROT2[2][2] = cos(theta);
//
//		// rotate u along z to x
//		TVector3 u2(u);
//		u2.RotateZ(phi);
//		u2.RotateY(theta);
//		float phip = -TMath::ATan2(u2.Y(), u2.X());
//		phip -= -TMath::ATan2(up.Y(), up.X());
//		ROT3[0][0] = cos(phip);
//		ROT3[0][1] = -sin(phip);
//		ROT3[0][2] = 0;
//		ROT3[1][0] = sin(phip);
//		ROT3[1][1] = cos(phip);
//		ROT3[1][2] = 0;
//		ROT3[2][0] = 0;
//		ROT3[2][1] = 0;
//		ROT3[2][2] = 1;
//
//		// R: rotation from u,v,n to (z X n), v', z
//		R = ROT3 * ROT2 * ROT1;
//		R_inv = R.Invert();
//		R_inv_T.Transpose(R_inv);
//
//	} catch (...) {
//		if (verbosity > 0)
//			LogWarning("Can't get rotation matrix");
//
//		return false;
//	}
//
//	pos_out.ResizeTo(3, 1);
//	cov_out.ResizeTo(3, 3);
//
//	pos_out = R_inv * pos_in;
//	cov_out = R_inv * cov_in * R_inv_T;
//
//	return true;
//}

bool PHG4TrackKalmanFitter::pos_cov_uvn_to_rz(const TVector3& u, const TVector3& v,
		const TVector3& n, const TMatrixF& pos_in, const TMatrixF& cov_in,
		TMatrixF& pos_out, TMatrixF& cov_out) const {

	if(pos_in.GetNcols() != 1 || pos_in.GetNrows() != 3) {
		if(verbosity > 0) LogWarning("pos_in.GetNcols() != 1 || pos_in.GetNrows() != 3");
		return false;
	}

	if(cov_in.GetNcols() != 3 || cov_in.GetNrows() != 3) {
		if(verbosity > 0) LogWarning("cov_in.GetNcols() != 3 || cov_in.GetNrows() != 3");
		return false;
	}

	TVector3 Z_uvn(u.Z(),v.Z(),n.Z());
	TVector3 up_uvn = TVector3(0., 0., 1.).Cross(Z_uvn); // n_uvn X Z_uvn

	if(up_uvn.Mag() < 0.00001){
		if(verbosity > 0) LogWarning("n is parallel to z");
		return false;
	}


	// R: rotation from u,v,n to n X Z, nX(nXZ), n
	TMatrixF R(3, 3);
	TMatrixF R_T(3,3);

	try {
		// rotate u along z to up
		float phi = -TMath::ATan2(up_uvn.Y(), up_uvn.X());
		R[0][0] = cos(phi);
		R[0][1] = -sin(phi);
		R[0][2] = 0;
		R[1][0] = sin(phi);
		R[1][1] = cos(phi);
		R[1][2] = 0;
		R[2][0] = 0;
		R[2][1] = 0;
		R[2][2] = 1;

		R_T.Transpose(R);
	} catch (...) {
		if (verbosity > 0)
			LogWarning("Can't get rotation matrix");

		return false;
	}

	pos_out.ResizeTo(3, 1);
	cov_out.ResizeTo(3, 3);

	pos_out = R * pos_in;
	cov_out = R * cov_in * R_T;

	return true;
}

bool PHG4TrackKalmanFitter::get_vertex_error_uvn(const TVector3& u,
		const TVector3& v, const TVector3& n, const TMatrixF& cov_in,
		TMatrixF& cov_out) const {

	/*!
	 * Get matrix that rotates frame (u,v,n) to (x,y,z)
	 * or the matrix that rotates vector defined in (x,y,z) to defined (u,v,n)
	 */

	TMatrixF R = get_rotation_matrix(u, v, n);
//
//	LogDebug("PHG4TrackKalmanFitter::get_vertex_error_uvn::R = ");
//	R.Print();
//	cout<<"R.Determinant() = "<<R.Determinant()<<"\n";

	if(!(abs(R.Determinant()-1)<0.01)) {
		if (verbosity > 0)
			LogWarning("!(abs(R.Determinant()-1)<0.0001)");
		return false;
	}

	if (R.GetNcols() != 3 || R.GetNrows() != 3) {
		if (verbosity > 0)
			LogWarning("R.GetNcols() != 3 || R.GetNrows() != 3");
		return false;
	}

	if (cov_in.GetNcols() != 3 || cov_in.GetNrows() != 3) {
		if (verbosity > 0)
			LogWarning("cov_in.GetNcols() != 3 || cov_in.GetNrows() != 3");
		return false;
	}

	TMatrixF R_T(3,3);

	R_T.Transpose(R);

	cov_out.ResizeTo(3, 3);

	cov_out = R * cov_in * R_T;

	return true;
}


bool PHG4TrackKalmanFitter::pos_cov_XYZ_to_RZ(
		const TVector3& n, const TMatrixF& pos_in, const TMatrixF& cov_in,
		TMatrixF& pos_out, TMatrixF& cov_out) const {

	if(pos_in.GetNcols() != 1 || pos_in.GetNrows() != 3) {
		if(verbosity > 0) LogWarning("pos_in.GetNcols() != 1 || pos_in.GetNrows() != 3");
		return false;
	}

	if(cov_in.GetNcols() != 3 || cov_in.GetNrows() != 3) {
		if(verbosity > 0) LogWarning("cov_in.GetNcols() != 3 || cov_in.GetNrows() != 3");
		return false;
	}

	TVector3 r = n.Cross(TVector3(0.,0.,1.));

	if(r.Mag() < 0.00001){
		if(verbosity > 0) LogWarning("n is parallel to z");
		return false;
	}


	// R: rotation from u,v,n to n X Z, nX(nXZ), n
	// these do not scale, they're pure rotation
	TMatrixF R(3, 3);
	TMatrixF R_T(3,3);

	try {
		// rotate u along z to up
		float phi = -TMath::ATan2(r.Y(), r.X());
		R[0][0] = cos(phi);
		R[0][1] = -sin(phi);
		R[0][2] = 0;
		R[1][0] = sin(phi);
		R[1][1] = cos(phi);
		R[1][2] = 0;
		R[2][0] = 0;
		R[2][1] = 0;
		R[2][2] = 1;

		R_T.Transpose(R);
	} catch (...) {
		if (verbosity > 0)
			LogWarning("Can't get rotation matrix");

		return false;
	}

	pos_out.ResizeTo(3, 1);
	cov_out.ResizeTo(3, 3);

	pos_out = R * pos_in;
	cov_out = R * cov_in * R_T;

	return true;
}

/*!
 * Get 3D Rotation Matrix that rotates frame (x,y,z) to (x',y',z')
 * Default rotate local to global, or rotate vector in global to local representation
 */
TMatrixF PHG4TrackKalmanFitter::get_rotation_matrix(const TVector3 x,
		const TVector3 y, const TVector3 z, const TVector3 xp, const TVector3 yp,
		const TVector3 zp) const {

	TMatrixF R(3, 3);

	TVector3 xu = x.Unit();
	TVector3 yu = y.Unit();
	TVector3 zu = z.Unit();

	const float max_diff = 0.01;

	if(!(
			abs(xu*yu) < max_diff and
			abs(xu*zu) < max_diff and
			abs(yu*zu) < max_diff
			)) {
		if (verbosity > 0)
			LogWarning("input frame error!");
		return R;
	}

	TVector3 xpu = xp.Unit();
	TVector3 ypu = yp.Unit();
	TVector3 zpu = zp.Unit();

	if(!(
			abs(xpu*ypu) < max_diff and
			abs(xpu*zpu) < max_diff and
			abs(ypu*zpu) < max_diff
			)) {
		if (verbosity > 0)
			LogWarning("output frame error!");
		return R;
	}

	/*!
	 * Decompose x',y',z' in x,y,z and call them u,v,n
	 * Then the question will be rotate the standard X,Y,Z to u,v,n
	 */

	TVector3 u(xpu.Dot(xu),xpu.Dot(yu),xpu.Dot(zu));
	TVector3 v(ypu.Dot(xu),ypu.Dot(yu),ypu.Dot(zu));
	TVector3 n(zpu.Dot(xu),zpu.Dot(yu),zpu.Dot(zu));


	try {
		std::shared_ptr<TRotation> rotation(new TRotation());
		//TRotation *rotation = new TRotation();

		//! Rotation that rotate standard (X, Y, Z) to (u, v, n)
		rotation->RotateAxes(u, v, n);

		R[0][0] = rotation->XX();
		R[0][1] = rotation->XY();
		R[0][2] = rotation->XZ();
		R[1][0] = rotation->YX();
		R[1][1] = rotation->YY();
		R[1][2] = rotation->YZ();
		R[2][0] = rotation->ZX();
		R[2][1] = rotation->ZY();
		R[2][2] = rotation->ZZ();
//
//		LogDebug("PHG4TrackKalmanFitter::get_rotation_matrix: TRotation:");
//		R.Print();
//		cout<<"R.Determinant() = "<<R.Determinant()<<"\n";

		//delete rotation;

//		TMatrixF ROT1(3, 3);
//		TMatrixF ROT2(3, 3);
//		TMatrixF ROT3(3, 3);
//
//		// rotate n along z to xz plane
//		float phi = -TMath::ATan2(n.Y(), n.X());
//		ROT1[0][0] = cos(phi);
//		ROT1[0][1] = -sin(phi);
//		ROT1[0][2] = 0;
//		ROT1[1][0] = sin(phi);
//		ROT1[1][1] = cos(phi);
//		ROT1[1][2] = 0;
//		ROT1[2][0] = 0;
//		ROT1[2][1] = 0;
//		ROT1[2][2] = 1;
//
//		// rotate n along y to z
//		TVector3 n1(n);
//		n1.RotateZ(phi);
//		float theta = -TMath::ATan2(n1.X(), n1.Z());
//		ROT2[0][0] = cos(theta);
//		ROT2[0][1] = 0;
//		ROT2[0][2] = sin(theta);
//		ROT2[1][0] = 0;
//		ROT2[1][1] = 1;
//		ROT2[1][2] = 0;
//		ROT2[2][0] = -sin(theta);
//		ROT2[2][1] = 0;
//		ROT2[2][2] = cos(theta);
//
//		// rotate u along z to x
//		TVector3 u2(u);
//		u2.RotateZ(phi);
//		u2.RotateY(theta);
//		float phip = -TMath::ATan2(u2.Y(), u2.X());
//		ROT3[0][0] = cos(phip);
//		ROT3[0][1] = -sin(phip);
//		ROT3[0][2] = 0;
//		ROT3[1][0] = sin(phip);
//		ROT3[1][1] = cos(phip);
//		ROT3[1][2] = 0;
//		ROT3[2][0] = 0;
//		ROT3[2][1] = 0;
//		ROT3[2][2] = 1;
//
//		// R: rotation from u,v,n to (z X n), v', z
//		R = ROT3 * ROT2 * ROT1;
//
//		R.Invert();
//		LogDebug("PHG4TrackKalmanFitter::get_rotation_matrix: Home Brew:");
//		R.Print();
//		cout<<"R.Determinant() = "<<R.Determinant()<<"\n";

	} catch (...) {
		if (verbosity > 0)
			LogWarning("Can't get rotation matrix");

		return R;
	}

	return R;
}
