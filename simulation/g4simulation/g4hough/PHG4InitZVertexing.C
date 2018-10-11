
#include "PHG4InitZVertexing.h"
#include "CellularAutomaton_v1.h"


// g4hough includes
#include "SvtxVertexMap.h"
#include "SvtxVertexMap_v1.h"
#include "SvtxVertex.h"
#include "SvtxVertex_v1.h"
#include "SvtxTrackMap.h"
#include "SvtxTrackMap_v1.h"
#include "SvtxTrack.h"
#include "SvtxTrack_v1.h"
#include "SvtxTrackState.h"
#include "SvtxClusterMap.h"
#include "SvtxCluster.h"
#include "SvtxHit_v1.h"
#include "SvtxHitMap.h"

// sPHENIX Geant4 includes
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>

#include <g4bbc/BbcVertexMap.h>
#include <g4bbc/BbcVertex.h>

// sPHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phgeom/PHGeomUtility.h>

// sGeant4 includes
#include <Geant4/G4MagneticField.hh>
#include <Geant4/G4TransportationManager.hh>
#include <Geant4/G4FieldManager.hh>

// standard includes
#include <TNtuple.h>
#include <cmath>
#include <float.h>
#include <iostream>
#include <fstream>
#include <bitset>
#include <tuple>
#include <map>


#define LogDebug(exp)		std::cout<<"DEBUG: "  <<__FILE__<<": "<<__LINE__<<": "<< exp
#define LogError(exp)		std::cout<<"ERROR: "  <<__FILE__<<": "<<__LINE__<<": "<< exp
#define LogWarning(exp)	std::cout<<"WARNING: "<<__FILE__<<": "<<__LINE__<<": "<< exp

//#define _DEBUG_
//#define _HOUGHTRANSF_
#define _MULTIVTX_
#define _TRIPLETS_

using namespace std;
using namespace Eigen;


PHG4InitZVertexing::PHG4InitZVertexing(unsigned int nlayers,
                                       unsigned int min_nlayers,
                                       const string& name)
    : SubsysReco(name),
	  _t_output_io(nullptr),
	  _seeding_layer(),
      _nlayers(nlayers),
      _min_nlayers(min_nlayers),
      _ca_nlayers(nlayers),
      _ca_min_nlayers(min_nlayers),
      _ca_chi2(2.),
      _radii(),
      _material(),
      _user_material(),
      _mag_field(1.4),
      _remove_hits(true),
      _use_max_kappa(false),
      fill_multi_zvtx(true),
      _min_pt(0.2),
      _min_d(-1.0),
      _max_d(1.0),
      _min_dzdl(-0.9),
      _max_dzdl(0.9),
      _min_z0(-14.0),
      _max_z0(+14.0),
      _cut_on_dca(false),
      _dcaxy_cut(0.2),
      _dcaz_cut(0.05),
      _pt_rescale(1.0),
      _ca_phi_cut(M_PI/36.),
      _ca_z_cut(20.),
      _z_cut(1.),
      _mult_onebin(1.),
      _mult_twobins(1.),
      _mult_threebins(1.),
      _min_zvtx_tracks(1),
      bin(0),
      ik(0),
      ip(0),
      id(0),
      il(0),
      iz(0),
      nkappa(0),
      nphi(0),
      nd(0),
      ndzdl(0),
      nz0(0),
      cluster_id(0),
      _layer_ilayer_map(),
      _temp_tracks(),
      _tracks(),
      _track_states(),
      _track_errors(),
      _track_covars(),
      _vertex(),
      _vertex_list(),
      _bbc_vertexes(NULL),
      _clustermap(NULL),
      _trackmap(NULL),
      _vertexmap(NULL),
      _vertex_finder(),
	_hough_space(NULL),
	_hough_funcs(NULL),
	_mode(0),
	_ntp_zvtx_by_event(NULL),
	_ntp_zvtx_by_track(NULL),
	_z0_dzdl(NULL),
	_kappa_phi(NULL),
	_d_phi(NULL),
	_kappa_d_phi(NULL),
	_ofile(NULL),
	_ofile2(NULL),
	  _nlayers_all(67),
	  _layer_ilayer_map_all(),
	  _radii_all(),
	  zooms_vec(),
	  hits_used(),
	  hits_map(),
	  bins_map(),
	  bins_map_prev(),
	  bins_map_cur(),
	  bins_map_sel(),
	  kappa_phi_map(),
	  d_phi_map(),
	  kappa_d_phi_map(),
	  hit(),
	  hitpos3d(),
	  nzooms(0),
	  separate_helicity(false),
	  helicity(0),
	  n_vtx_tracks(0)
	  {
	_event = 0;
	zooms_vec.clear();
}

int PHG4InitZVertexing::Init(PHCompositeNode* topNode) {

#ifdef _DEBUG_
    _ofile = new TFile("z0_dzdl_kappa_phi_d.root","recreate");
#endif

    if (!_use_max_kappa){
    _max_kappa = pt_to_kappa(_min_pt);
    cout<<"kappa max "<<_max_kappa<<endl;
    }

    set_nzooms();
    _hough_space = new HelixHoughSpace_v1();
    for (unsigned int izoom =0; izoom<nzooms ; ++izoom)
    _hough_space->add_one_zoom(zooms_vec[izoom]);

    _hough_space->set_kappa_min(0);
    _hough_space->set_kappa_max(_max_kappa);
    _hough_space->set_phi_min(0);
    _hough_space->set_phi_max(2.*M_PI);
    _hough_space->set_d_min(_min_d);
    _hough_space->set_d_max(_max_d);
    _hough_space->set_dzdl_min(_min_dzdl);// 0 if separate helicity/chargge
    _hough_space->set_dzdl_max(_max_dzdl);// 0.9 default
    _hough_space->set_z0_min(_min_z0);// -14 for non-zero vertex 
    _hough_space->set_z0_max(_max_z0);

    _hough_funcs = new HelixHoughFuncs_v1();
    _hough_funcs->set_hough_space(_hough_space);

#ifdef _DEBUG_
//    _hough_space->print_zoom_profile();
//    _hough_space->print_para_range();
    unsigned int n_z0_bins= _hough_space->get_n_z0_bins(0);
    unsigned int n_dzdl_bins= _hough_space->get_n_dzdl_bins(0);
    unsigned int n_kappa_bins = _hough_space->get_n_kappa_bins(0);
    unsigned int n_d_bins = _hough_space->get_n_d_bins(0); 
    unsigned int n_phi_bins = _hough_space->get_n_phi_bins(0);
  
    _z0_dzdl = new TH2D("z0_dzdl","z0_dzdl",n_z0_bins,_hough_space->get_z0_min(),_hough_space->get_z0_max(),n_dzdl_bins,_hough_space->get_dzdl_min(),_hough_space->get_dzdl_max());  
//    _z0_dzdl = new TH2D("z0_dzdl","z0_dzdl",n_z0_bins*_hough_space->get_n_z0_bins(1),_hough_space->get_z0_min(),_hough_space->get_z0_max(),n_dzdl_bins*_hough_space->get_n_dzdl_bins(1),_hough_space->get_dzdl_min(),_hough_space->get_dzdl_max());
    _kappa_phi = new TH2D("kappa_phi","kappa_phi",n_kappa_bins,_hough_space->get_kappa_min(),_hough_space->get_kappa_max(),n_phi_bins,_hough_space->get_phi_min(),_hough_space->get_phi_max());

//    _kappa_phi = new TH2D("kappa_phi","kappa_phi",n_kappa_bins*_hough_space->get_n_kappa_bins(1),_hough_space->get_kappa_min(),_hough_space->get_kappa_max(),n_phi_bins*_hough_space->get_n_phi_bins(1),_hough_space->get_phi_min(),_hough_space->get_phi_max());
    _d_phi = new TH2D("d_phi","d_phi",n_d_bins,_hough_space->get_d_min(),_hough_space->get_d_max(),n_phi_bins,_hough_space->get_phi_min(),_hough_space->get_phi_max());
//    _d_phi = new TH2D("d_phi","d_phi",n_d_bins*_hough_space->get_n_d_bins(1),_hough_space->get_d_min(),_hough_space->get_d_max(),n_phi_bins*_hough_space->get_n_phi_bins(1),_hough_space->get_phi_min(),_hough_space->get_phi_max());

    _kappa_d_phi = new TH3D("kappa_d_phi","kappa_d_phi",n_kappa_bins,_hough_space->get_kappa_min(),_hough_space->get_kappa_max(),n_d_bins,_hough_space->get_d_min(),_hough_space->get_d_max(),n_phi_bins,_hough_space->get_phi_min(),_hough_space->get_phi_max());
#endif
	return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4InitZVertexing::InitRun(PHCompositeNode* topNode) {

	int code = Fun4AllReturnCodes::ABORTRUN;

	code = create_nodes(topNode);
	if(code != Fun4AllReturnCodes::EVENT_OK)
		return code;

	code = initialize_geometry(topNode);
	if(code != Fun4AllReturnCodes::EVENT_OK)
		return code;

#ifdef _MULTIVTX_
    _ofile2 = new TFile(_fname.c_str(),"recreate");
    _ntp_zvtx_by_event = new TNtuple("ntp_zvtx_by_event","all vertices found event-by-event","event:zvtx");
    _ntp_zvtx_by_track = new TNtuple("ntp_zvtx_by_track","track-by-track (zvtx + z0) distribution", "event:zvtx");
#endif

	_t_output_io = new PHTimer("_t_output_io");
	_t_output_io->stop();

	if (Verbosity() > 0) {
		cout
				<< "====================== PHG4InitZVertexing::InitRun() ======================"
				<< endl;
		cout << " Magnetic field set to: " << _mag_field << " Tesla" << endl;
		cout << " Number of tracking layers: " << _nlayers << endl;
		for (unsigned int i = 0; i < _nlayers; ++i) {
			cout << "   Tracking layer #" << i << " " << "radius = "
					<< _radii[i] << " cm, " << "material = " << _material[i]
					<< endl;
		}
		cout << " Required hits: " << _min_nlayers << endl;
		cout << " Minimum pT: " << _min_pt << endl;
		cout << " Maximum DCA: " << boolalpha << _cut_on_dca << noboolalpha
				<< endl;
		if (_cut_on_dca) {
			cout << "   Maximum DCA cut: " << _dcaxy_cut << endl;
		}
		cout << "   Maximum DCAZ cut: " << _dcaz_cut << endl;
		cout << " Momentum rescale factor: " << _pt_rescale << endl;
		cout
				<< "==========================================================================="
				<< endl;
	}


	return code;
}

int PHG4InitZVertexing::process_event(PHCompositeNode *topNode) 
{

	if (Verbosity() > 0)
		cout << "PHG4InitZVertexing::process_event -- entered" << endl;

	// start fresh
	for(unsigned int i=0; i<_tracks.size(); ++i) _tracks[i].reset();
	_tracks.clear();
	_track_states.clear();
	_track_errors.clear();
	_track_covars.clear();
	_vertex.clear();
	_vertex.assign(3, 0.0);
	_vertex_list.clear();

	hits_map.clear();
	hits_used.clear();

	//-----------------------------------
	// Get Objects off of the Node Tree
	//-----------------------------------

	get_nodes(topNode);

        int code = translate_input(topNode);
        if (code != Fun4AllReturnCodes::EVENT_OK)
        return code;

	// First iteration

#ifdef _HOUGHTRANSF_
        reset_zooms();
        add_zoom(1,18,1,16,16);// zoom 0
        add_zoom(1,3,1,3,2);// zoom 1
        add_zoom(1,3,1,3,2);// zoom 2
        separate_helicity = true;
#endif

        Init(topNode);

	bool zvtx_found = false;
        _ca_nlayers = 3; // 3(vertexing)/4(track seeding) for p+p

	// p+p 
//	unsigned int nseq =4;		//		 65(0.5)       103(0.5)
//      unsigned int nattempt =4;
	// N+N
	unsigned int nseq =  1;
	unsigned int nattempt = 1;  //		         59(1.5) 63(1.5)   

	unsigned int iseq=0;
	while(1)
	{// iseq 

		if (iseq==nseq) break;
		for(unsigned int i=0; i<_tracks.size(); ++i) _tracks[i].reset();
                _tracks.clear();
                _track_states.clear();
                _track_errors.clear();
                _track_covars.clear();
                reset_hits_used();

        	float shift_dx = -_vertex[0];
        	float shift_dy = -_vertex[1];
        	float shift_dz = -_vertex[2];
        	shift_coordinate_system(shift_dx,shift_dy,shift_dz);

		for (unsigned int iattempt =0; iattempt<nattempt; ++iattempt ){
			cout<<iattempt << " th attempt "<<endl;
			helicity = 1;

#ifdef _TRIPLETS_ 
			_temp_tracks.clear();

			build_triplets_to_Track3D(_temp_tracks, true); // false : backward
#endif

#ifdef _HOUGHTRANSF_

			unsigned int zoomlevel = 0;			
			zoomlevel = 0;
			set_nbins(zoomlevel);

#ifdef _DEBUG_
//        cout<<"InitZVertexing:: nkappa " <<nkappa<<" nphi "<<nphi<<" nd "<<nd<<" ndzdl "<<ndzdl<<" nz0 " <<nz0<<endl;
#endif

			for(unsigned int i=0; i<_temp_tracks.size(); ++i) _temp_tracks[i].reset();
			_temp_tracks.clear();
                        for (std::map<unsigned int,HelixHoughBin* >::iterator it=bins_map.begin(); it!=bins_map.end(); ++it){ delete it->second; bins_map.erase(it);}
                        for (std::map<unsigned int,HelixHoughBin* >::iterator it=bins_map_prev.begin(); it!=bins_map_prev.end(); ++it){ delete it->second; bins_map_prev.erase(it);}
                        for (std::map<unsigned int,HelixHoughBin* >::iterator it=bins_map_cur.begin(); it!=bins_map_cur.end(); ++it){ delete it->second; bins_map_cur.erase(it);}
                        for (std::map<unsigned int,HelixHoughBin* >::iterator it=bins_map_sel.begin(); it!=bins_map_sel.end(); ++it){ delete it->second; bins_map_sel.erase(it);}
                        bins_map.clear();
                        bins_map_prev.clear();
                        bins_map_cur.clear();
                        bins_map_sel.clear();

			initialize_houghbin();
        		vote_z_init(zoomlevel); 
			find_track_candidates_z_init(zoomlevel);

			vote_xy(zoomlevel);
			find_track_candidates_xy(zoomlevel); 
			prune_xy(zoomlevel);// on prev

			bins_map_prev.swap(bins_map_sel);
			prune_z(zoomlevel); // on sel
			bins_map_sel.swap(bins_map_prev);

			for (zoomlevel =1; zoomlevel<2; ++zoomlevel){ // depends on the total number of tracks thrown in
			vote_z(zoomlevel); // on prev
			find_track_candidates_z(zoomlevel);
			prune_z(zoomlevel);
			vote_xy(zoomlevel);
			find_track_candidates_xy(zoomlevel);
			prune_xy(zoomlevel);
			} //zoomlevel


			// high pt, primary vertex originating tracks
			bins_to_Track3D(_temp_tracks,1,zoomlevel);// 0(z) 1(xy) : create Track3D objects, save cluster_id
#endif

			code = 0;
			code = cellular_automaton_zvtx_init(_temp_tracks); // remove bad hits, save fit parameters && turn off used hits
			if (code != Fun4AllReturnCodes::EVENT_OK)
			{
			cout<<"CellularAutomaton failed. "<<endl;
			exit(1);
			}


		}// iattempt

		code = 0;
		code = fit_vertex();
		cout<<"seq "<<iseq<<" : vertex_z = "<< _vertex[2] <<", shift_z = "<<shift_dz<<endl;
		if (iseq==0)
		{
			if (code==-1) zvtx_found = false;
			else {
			zvtx_found = true;
			shift_coordinate_system(-shift_dx,-shift_dy,-shift_dz);
			break;
			}
			++iseq;
			continue;
//		cout<<"Errors in fitting vertex. "<<endl;
//		exit(1);	
		}else if (iseq>0 && code<0 && !zvtx_found){
			++_ca_nlayers;
			_ca_chi2 +=0.5;

                        shift_coordinate_system(-shift_dx,-shift_dy,-shift_dz);
                        ++iseq;
			if(zvtx_found)
                        cout<< "z-vertex not fitted. "<<endl;
			else
			cout<< "z-vertex not found. "<<endl;
                        continue;
		}else {
		if (iseq<nseq) cout<<"z-vertex fitted. "<<endl;
		shift_coordinate_system(-shift_dx,-shift_dy,-shift_dz);
		break;
		}
	}//iseq


        code = export_output();
        if (code != Fun4AllReturnCodes::EVENT_OK)
        return code;

	++_event;

//        BbcVertex* vertex = _bbc_vertexes->begin()->second;
//	cout<< "true z-vertex  " << vertex->get_z()<<endl;

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4InitZVertexing::End(PHCompositeNode *topNode) {
#ifdef _DEBUG_
	_ofile->cd();
	_z0_dzdl->Write();
	_kappa_phi->Write();
	_d_phi->Write();
	_kappa_d_phi->Write();
	_ofile->Close();
#endif

#ifdef _MULTIVTX_
        _ofile2->cd();
        _ntp_zvtx_by_event->Write();
	_ntp_zvtx_by_track->Write();
        _ofile2->Close();
#endif

	delete _t_output_io;
	if (_hough_space != NULL) delete _hough_space;
	if (_hough_funcs != NULL) delete _hough_funcs;
	if (ca != NULL) delete ca;
	_temp_tracks.clear();

	return Fun4AllReturnCodes::EVENT_OK;
}


void PHG4InitZVertexing::set_material(int layer, float value) {
	_user_material[layer] = value;
}

void PHG4InitZVertexing::add_zoom(unsigned int n_kappa, unsigned int n_phi, unsigned int n_d,unsigned int n_dzdl, unsigned int n_z0) {
        std::vector<unsigned int> zoom {n_kappa, n_phi, n_d, n_dzdl, n_z0};
	zooms_vec.push_back(zoom);
}

int PHG4InitZVertexing::create_nodes(PHCompositeNode* topNode) {
	// create nodes...
	PHNodeIterator iter(topNode);

	PHCompositeNode* dstNode = static_cast<PHCompositeNode*>(iter.findFirst(
			"PHCompositeNode", "DST"));
	if (!dstNode) {
		cerr << PHWHERE << "DST Node missing, doing nothing." << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}
	PHNodeIterator iter_dst(dstNode);

	// Create the SVTX node
	PHCompositeNode* tb_node =
			dynamic_cast<PHCompositeNode*>(iter_dst.findFirst("PHCompositeNode",
					"SVTX"));
	if (!tb_node) {
		tb_node = new PHCompositeNode("SVTX");
		dstNode->addNode(tb_node);
		if (Verbosity() > 0)
			cout << "SVTX node added" << endl;
	}

	_trackmap = new SvtxTrackMap_v1;
	PHIODataNode<PHObject>* tracks_node = new PHIODataNode<PHObject>(_trackmap,
			"SvtxTrackMap", "PHObject");
	tb_node->addNode(tracks_node);
	if (Verbosity() > 0)
		cout << "Svtx/SvtxTrackMap node added" << endl;

	_vertexmap = new SvtxVertexMap_v1;
	PHIODataNode<PHObject>* vertexes_node = new PHIODataNode<PHObject>(
			_vertexmap, "SvtxVertexMap", "PHObject");
	tb_node->addNode(vertexes_node);
	if (Verbosity() > 0)
		cout << "Svtx/SvtxVertexMap node added" << endl;


	return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4InitZVertexing::initialize_geometry(PHCompositeNode *topNode) {

	//---------------------------------------------------------
	// Grab Run-Dependent Detector Geometry and Configure Hough
	//---------------------------------------------------------

	PHG4CylinderCellGeomContainer* cellgeos = findNode::getClass<
			PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
	PHG4CylinderGeomContainer* laddergeos = findNode::getClass<
			PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_SILICON_TRACKER");
	PHG4CylinderGeomContainer* mapsladdergeos = findNode::getClass<
			PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MAPS");


	_nlayers = _seeding_layer.size();


	_radii.assign(_nlayers, 0.0);
	map<float, int> radius_layer_map;

	_radii_all.assign(_nlayers_all, 0.0);

	if (cellgeos) {
		PHG4CylinderCellGeomContainer::ConstRange layerrange =
				cellgeos->get_begin_end();
		for (PHG4CylinderCellGeomContainer::ConstIterator layeriter =
				layerrange.first; layeriter != layerrange.second; ++layeriter) {
			radius_layer_map.insert(
					make_pair(layeriter->second->get_radius(),
							layeriter->second->get_layer()));
		}
	}

	if (laddergeos) {
		PHG4CylinderGeomContainer::ConstRange layerrange =
				laddergeos->get_begin_end();
		for (PHG4CylinderGeomContainer::ConstIterator layeriter =
				layerrange.first; layeriter != layerrange.second; ++layeriter) {
			radius_layer_map.insert(
				make_pair(layeriter->second->get_radius(),
				layeriter->second->get_layer()));
		}
	}

	if (mapsladdergeos) {
		PHG4CylinderGeomContainer::ConstRange layerrange =
				mapsladdergeos->get_begin_end();
		for (PHG4CylinderGeomContainer::ConstIterator layeriter =
				layerrange.first; layeriter != layerrange.second; ++layeriter) {
			radius_layer_map.insert(
				make_pair(layeriter->second->get_radius(),
				layeriter->second->get_layer()));
		}
	}


	// now that the layer ids are sorted by radius, I can create a storage
	// index, ilayer, that is 0..N-1 and sorted by radius

	int ilayer = 0;
	for (map<float, int>::iterator iter = radius_layer_map.begin();
			iter != radius_layer_map.end(); ++iter) {

		_layer_ilayer_map_all.insert(make_pair(iter->second, _layer_ilayer_map_all.size()));

		if (std::find(_seeding_layer.begin(), _seeding_layer.end(),
			iter->second) != _seeding_layer.end()) {
			_layer_ilayer_map.insert(make_pair(iter->second, ilayer));
			++ilayer;
		}
	}


	// now we extract the information from the cellgeos first
	if (cellgeos) {
		PHG4CylinderCellGeomContainer::ConstRange begin_end =
			cellgeos->get_begin_end();
		PHG4CylinderCellGeomContainer::ConstIterator miter = begin_end.first;
		for (; miter != begin_end.second; ++miter) {
			PHG4CylinderCellGeom *geo = miter->second;

			//if(cellgeo->get_layer() > (int) _radii.size() ) continue;

//			if (Verbosity() >= 2)
//			cellgeo->identify();

			//TODO
			_radii_all[_layer_ilayer_map_all[geo->get_layer()]] =
			geo->get_radius() + 0.5*geo->get_thickness();

			if (_layer_ilayer_map.find(geo->get_layer())
				!= _layer_ilayer_map.end()) {
			_radii[_layer_ilayer_map[geo->get_layer()]] =
			geo->get_radius();
			}
		}
	}

	if (laddergeos) {
		PHG4CylinderGeomContainer::ConstRange begin_end =
			laddergeos->get_begin_end();
		PHG4CylinderGeomContainer::ConstIterator miter = begin_end.first;
		for (; miter != begin_end.second; ++miter) {
			PHG4CylinderGeom *geo = miter->second;

			//if(geo->get_layer() > (int) _radii.size() ) continue;

//			if (Verbosity() >= 2)
//			geo->identify();

			_radii_all[_layer_ilayer_map_all[geo->get_layer()]] =
			geo->get_radius() + 0.5*geo->get_thickness();

			if (_layer_ilayer_map.find(geo->get_layer())
				!= _layer_ilayer_map.end()) {
				_radii[_layer_ilayer_map[geo->get_layer()]] = geo->get_radius();
			}
		}
	}

	if (mapsladdergeos) {
		PHG4CylinderGeomContainer::ConstRange begin_end =
			mapsladdergeos->get_begin_end();
		PHG4CylinderGeomContainer::ConstIterator miter = begin_end.first;
		for (; miter != begin_end.second; ++miter) {
			PHG4CylinderGeom *geo = miter->second;

			//if(geo->get_layer() > (int) _radii.size() ) continue;

//			if (Verbosity() >= 2)
//				geo->identify();

			_radii_all[_layer_ilayer_map_all[geo->get_layer()]] =
				geo->get_radius() + 0.5*geo->get_thickness();

			if (_layer_ilayer_map.find(geo->get_layer())
				!= _layer_ilayer_map.end()) {
				_radii[_layer_ilayer_map[geo->get_layer()]] = geo->get_radius();
			}
		}
	}

	// set material on each layer

	_material.assign(_radii.size(), 0.03);

	map<int, float>::iterator mat_it;
	for (map<int, float>::iterator iter = _user_material.begin();
			iter != _user_material.end(); ++iter) {
		if (_layer_ilayer_map.find(iter->first) != _layer_ilayer_map.end()) {
			_material[_layer_ilayer_map[iter->first]] = iter->second;
		}
	}

	// initialize the pattern recogition tools


	return Fun4AllReturnCodes::EVENT_OK;
}


int PHG4InitZVertexing::get_nodes(PHCompositeNode* topNode) {

	//---------------------------------
	// Get Objects off of the Node Tree
	//---------------------------------

	// Pull the reconstructed track information off the node tree...
	_bbc_vertexes = findNode::getClass<BbcVertexMap>(topNode, "BbcVertexMap");

	_clustermap = findNode::getClass<SvtxClusterMap>(topNode, "SvtxClusterMap");
	if (!_clustermap) {
		cerr << PHWHERE << " ERROR: Can't find node SvtxClusterMap" << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	// Pull the reconstructed track information off the node tree...
	_trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
	if (!_trackmap) {
		cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap." << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	// Pull the reconstructed track information off the node tree...
	_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
	if (!_vertexmap) {
		cerr << PHWHERE << " ERROR: Can't find SvtxVertexMap." << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4InitZVertexing::translate_input(PHCompositeNode* topNode) {

	
        for (SvtxClusterMap::Iter iter = _clustermap->begin();
                iter != _clustermap->end(); ++iter) 
	{

        	SvtxCluster* cluster = iter->second;

                unsigned int ilayer = UINT_MAX;
                std::map<int, unsigned int>::const_iterator it = _layer_ilayer_map.find(cluster->get_layer());
                if(it != _layer_ilayer_map.end())
                        ilayer = it->second;
                if(ilayer >= _nlayers) continue;

                Cluster3D hit3d;

                hit3d.set_id(cluster->get_id());
                hit3d.set_layer(ilayer);

//		cout<<" adding cluster, id =  "<< cluster->get_id()<<", layer = "<<cluster->get_layer()<<endl;
                hit3d.set_x(cluster->get_x());
                hit3d.set_y(cluster->get_y());
                hit3d.set_z(cluster->get_z());
                for (int i = 0; i < 3; ++i) {
                        for (int j = i; j < 3; ++j) {
                                hit3d.set_error(i, j, cluster->get_error(i, j));

	                //hit3d.set_size(i, j, cluster->get_size(i, j)); // original
        	        hit3d.set_size(i, j, cluster->get_error(i, j)*sqrt(12.)); // yuhw 2017-05-08
                        }
                }
                hits_map.insert(std::pair<unsigned int, Cluster3D>(hit3d.get_id(),hit3d));
		hits_used.insert(std::pair<unsigned int, bool>(hit3d.get_id(),false));
        }
	
        if (Verbosity() > 10) {
        cout << "-------------------------------------------------------------------"
             << endl;
        cout << "PHG4InitZVertexing::process_event has the following input clusters:"
             << endl;

	cout << "n init clusters = " << hits_map.size() << endl;

		for (std::map<unsigned int,Cluster3D>::iterator it=hits_map.begin(); 
			it!=hits_map.end(); 
			++it)
		{
                	Cluster3D hit = it->second;
			hit.print();
		}
	   

        cout << "-------------------------------------------------------------------"
             << endl;
        }

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4InitZVertexing::export_output(){

	_vertexmap->clear();
	_trackmap->clear();
        if (_tracks.empty())
                return Fun4AllReturnCodes::EVENT_OK;

	std::vector<SvtxVertex_v1> svtx_vertex_list;
	unsigned int nvertex = _vertex_list.size();
	for (unsigned int vid = 0; vid < nvertex; ++vid ){

        SvtxVertex_v1 vertex;
//	vertex.set_id(vid);
        vertex.set_t0(0.0);
        for (int i = 0; i < 3; ++i)
        vertex.set_position(i, _vertex_list[vid][i]);
        vertex.set_chisq(0.0);
        vertex.set_ndof(0);
        vertex.set_error(0, 0, 0.0);
        vertex.set_error(0, 1, 0.0);
        vertex.set_error(0, 2, 0.0);
        vertex.set_error(1, 0, 0.0);
        vertex.set_error(1, 1, 0.0);
        vertex.set_error(1, 2, 0.0);
        vertex.set_error(2, 0, 0.0);
        vertex.set_error(2, 1, 0.0);
        vertex.set_error(2, 2, 0.0);

	svtx_vertex_list.push_back(vertex);
	}

	cout<<"vertex list "<<endl;
        // at this point we should already have an initial pt and pz guess...
        // need to translate this into the PHG4Track object...

        std::vector<Cluster3D> track_hits;
	track_hits.clear();
        int clusterID;

        for (unsigned int itrack = 0; itrack < _tracks.size(); ++itrack) {
		if (_tracks.at(itrack).vertex_id >= nvertex) continue;
                SvtxTrack_v1 track;
                track.set_id(itrack);
                track_hits.clear();
                track_hits = _tracks.at(itrack).hits;
//		cout<<"hits vector assigned " <<endl;
                for (unsigned int ihit = 0; ihit < track_hits.size(); ++ihit) {
                        if ((track_hits.at(ihit).get_id()) >= _clustermap->size()) {
                                continue;
                        }
                        SvtxCluster* cluster = _clustermap->get(
                                        track_hits.at(ihit).get_id());
                        clusterID = cluster->get_id();

#ifdef _DEBUG_
                        cout
                        <<__LINE__
                        <<": itrack: " << itrack
                        <<": nhits: " << track_hits.size()
                        <<": clusterID: " << clusterID
                        <<": layer: " << cluster->get_layer()
                        <<endl;
#endif

                        track.insert_cluster(clusterID);
                }
//		cout<<"hits added to a track "<<endl;
                float kappa = _tracks.at(itrack).kappa;
                float d = _tracks.at(itrack).d;
                float phi = _tracks.at(itrack).phi;
                float dzdl = _tracks.at(itrack).dzdl;
                float z0 = _tracks.at(itrack).z0;

                float pT = kappa_to_pt(kappa);

                float x_center = cos(phi) * (d + 1 / kappa); // x coordinate of circle center
                float y_center = sin(phi) * (d + 1 / kappa);  // y    "      "     " "

                // find helicity from cross product sign
                short int helicity;
		if ((track_hits[0].get_x()- x_center)*(track_hits[track_hits.size()-1].get_y()- y_center)
		- (track_hits[0].get_y()- y_center)*(track_hits[track_hits.size()-1].get_x()- x_center)> 0)
		{
                        helicity = 1;
                } else {
                        helicity = -1;
                }

                float pZ = 0;
                if (dzdl != 1) {
                        pZ = pT * dzdl / sqrt(1.0 - dzdl * dzdl);
                }
                int ndf = 2 * _tracks.at(itrack).hits.size() - 5;
//              track.set_chisq(_track_errors[itrack]);
                track.set_ndf(ndf);
                track.set_px(pT * cos(phi - helicity * M_PI / 2));
                track.set_py(pT * sin(phi - helicity * M_PI / 2));
                track.set_pz(pZ);

                track.set_dca2d(d);
//                track.set_dca2d_error(sqrt(_track_covars[itrack](1, 1)));

//		cout<<"track parameters set "<<endl;
                if (_mag_field > 0) {
                        track.set_charge(helicity);
                } else {
                        track.set_charge(-1.0 * helicity);
                }


  /*              Eigen::Matrix<float, 6, 6> euclidean_cov =
                Eigen::Matrix<float, 6, 6>::Zero(6, 6);
                convertHelixCovarianceToEuclideanCovariance(_mag_field, phi, d, kappa,
                                z0, dzdl, _track_covars[itrack], euclidean_cov);
*/
                for (unsigned int row = 0; row < 6; ++row) {
                        for (unsigned int col = 0; col < 6; ++col) {
  //                      track.set_error(row, col, euclidean_cov(row, col));
                        }
                }
//		cout<<"set track "<<endl;
		unsigned int vid = _tracks.at(itrack).vertex_id;
//		cout<<"vertex_id "<<vid<<endl;
                track.set_x(svtx_vertex_list[vid].get_x() + d * cos(phi));
                track.set_y(svtx_vertex_list[vid].get_y() + d * sin(phi));
                track.set_z(svtx_vertex_list[vid].get_z() + z0);

//		cout<<"set vertex "<<endl;
                _trackmap->insert(&track);
                svtx_vertex_list[vid].insert_track(track.get_id());

                if (Verbosity() > 5) {
                        cout << "track " << itrack << " quality = " << track.get_quality()
                        << endl;
                        cout << "px = " << track.get_px() << " py = " << track.get_py()
                        << " pz = " << track.get_pz() << endl;
                }
        }  // track loop

//	cout<<"track loop"<<endl;
	for (unsigned int vid = 0; vid < _vertex_list.size(); ++vid ){
        SvtxVertex *vtxptr = _vertexmap->insert(&svtx_vertex_list[vid]);
        if (Verbosity() > 5) vtxptr->identify();
	}

        hits_map.clear();

	for(unsigned int i=0; i<_tracks.size(); ++i) _tracks[i].reset();
        _tracks.clear();
        _track_errors.clear();
        _track_covars.clear();
        _vertex.clear();
        _vertex.assign(3, 0.0);
	_vertex_list.clear();	

        return Fun4AllReturnCodes::EVENT_OK;

}

float PHG4InitZVertexing::kappa_to_pt(float kappa) {
	return _pt_rescale * _mag_field / 333.6 / kappa;
}

float PHG4InitZVertexing::pt_to_kappa(float pt) {
	return _pt_rescale * _mag_field / 333.6 / pt;
}

void PHG4InitZVertexing::set_max_kappa(float kappa_max){
	_max_kappa = kappa_max;
	_use_max_kappa = true;
}

void PHG4InitZVertexing::set_nbins(unsigned int izoom){
        nkappa = _hough_space->get_n_kappa_bins(izoom);
        nphi = _hough_space->get_n_phi_bins(izoom);
        nd = _hough_space->get_n_d_bins(izoom);
        ndzdl = _hough_space->get_n_dzdl_bins(izoom);
        nz0 = _hough_space->get_n_z0_bins(izoom);
}

void PHG4InitZVertexing::initialize_houghbin(){
        for (ik=0; ik<nkappa; ++ik){
          for (ip=0; ip<nphi; ++ip) {
            for (id=0;id<nd; ++id) {
              for (il=0; il<ndzdl; ++il) {
                for (iz=0; iz<nz0; ++iz) {

                unsigned int bins[5]={ik,ip,id,il,iz};
//              cout<<"InitZVertexing:: izoom "<<izoom<<" il " <<il<<" iz "<<iz<<endl;
                bin = _hough_space->get_bin(0,bins);
//              cout<<"InitZVertexing:: bin " << bin<<endl;
                HelixHoughBin* _bin = new HelixHoughBin_v1(bin);
                _bin->set_zoomlevel(0);
                _bin->set_hough_space(_hough_space);
                _bin->init();
                if (bin == 0) _bin->identify();
                bins_map.insert(std::pair<unsigned int, HelixHoughBin*>(bin,_bin));
//              cout<<"InitZVertexing:: lbin "<<_bin->get_dzdl_bin(0)<<" zbin "<<_bin->get_z0_bin(0)<<endl;
                                }
              }
            }
          }
        }

        cout<<"bins_map.size " <<bins_map.size()<<endl;
}

void PHG4InitZVertexing::vote_z_init(unsigned int zoomlevel){

        float hitpos3d[3];
        unsigned int cluster_id=-999;
        std::vector<float> kappa_phi_d_ranges;
        float kappa_min = _hough_space->get_kappa_min();
        float kappa_max = _hough_space->get_kappa_max();
        float phi_min = _hough_space->get_phi_min();
        float phi_max = _hough_space->get_phi_max();
        float d_min = _hough_space->get_d_min();
        float d_max = _hough_space->get_d_max();
        kappa_phi_d_ranges.push_back(kappa_min);
        kappa_phi_d_ranges.push_back(kappa_max);
        kappa_phi_d_ranges.push_back(phi_min);
        kappa_phi_d_ranges.push_back(phi_max);
        kappa_phi_d_ranges.push_back(d_min);
        kappa_phi_d_ranges.push_back(d_max);

	float dzdl_min = _hough_space->get_dzdl_min();
	float dzdl_max = _hough_space->get_dzdl_max();
        int icluster=0;//test

//	cout<<"kappa min " <<kappa_min<<" kappa max "<<kappa_max<<" phi_min " <<phi_min<<endl;
        for (std::map<unsigned int,Cluster3D>::iterator it=hits_map.begin();
                        it!=hits_map.end();
                        ++it)
        {

 		cluster_id = it->first;	
		if (cluster_id != cluster_id) continue;
		bool used = false;
		auto hitused = hits_used.find(cluster_id);
		if (hitused != hits_used.end())
		used = hits_used.find(cluster_id)->second;
		if(used) continue;
		Cluster3D hit = it->second;
		hitpos3d[0] = hit.get_x();
		hitpos3d[1] = hit.get_y();
		hitpos3d[2] = hit.get_z();

//              cout<<"x "<< hitpos3d[0]<<" y "<<hitpos3d[1]<<" z "<<hitpos3d[2]<<endl;

        	for (iz=0; iz < _hough_space->get_n_z0_bins(zoomlevel); ++iz ) {
                std::vector<float> z0_range;
		float z0_min = _hough_space->get_z0_min()+iz*_hough_space->get_z0_bin_size(zoomlevel);
		float z0_max = _hough_space->get_z0_min()+(iz+1)*_hough_space->get_z0_bin_size(zoomlevel);
                z0_range.push_back(z0_min);
		z0_range.push_back(z0_max);
                  float dzdl_range[2] = {4,5};//test
                  //compute min_dzdl & max_dzdl for given iz 
                  _hough_funcs->set_current_zoom(zoomlevel);

#ifdef _DEBUG_
//                cout<<"set zoom "<<endl;
//                cout<<"i "<<i<<", iz "<<iz<<endl;
#endif

                  _hough_funcs->calculate_dzdl_range(hitpos3d,z0_range, kappa_phi_d_ranges, dzdl_range);
		  unsigned int dzdl_bin_range[2];
			dzdl_bin_range[0] = _hough_space->get_dzdl_bin(0,dzdl_range[0]);
			dzdl_bin_range[1] = _hough_space->get_dzdl_bin(0,dzdl_range[1]);
                        if (dzdl_range[0] > dzdl_max) continue;
                        else if (dzdl_range[0] <dzdl_min) dzdl_bin_range[0]=0;
                        if (dzdl_range[1] < dzdl_min) continue;
                        else if (dzdl_range[1] >dzdl_max) dzdl_bin_range[1]=_hough_space->get_n_dzdl_bins(zoomlevel)-1;
#ifdef _DEBUG_
                        cout<<"dzdl_min " <<dzdl_range[0]<<" dzdl_max "<<dzdl_range[1]<<endl;
#endif
                        for (il =dzdl_bin_range[0]; il<dzdl_bin_range[1]+1; ++il) {
                        	unsigned int zbins[5] = {0,0,0,il,iz};
                          	bin = _hough_space->get_bin(zoomlevel,zbins);
#ifdef _DEBUG_
                        	cout<<"bin "<<bin<<endl;
#endif
                          	auto search = bins_map.find(bin);
			  	if (search != bins_map.end())
			  	{
                          	bins_map.find(bin)->second->add_cluster_ID(cluster_id);
			 	}
#ifdef _DEBUG_
                        cout<<"count "<<bins_map.find(bin)->second->get_count()<<endl; 
#endif
                	} // bins_map_prev.clear(); // for zoomlevel >0
		z0_range.clear();
        	} //iz

	kappa_phi_d_ranges.clear();
        hitpos3d[0]=-999; hitpos3d[1]= -999; hitpos3d[2]=-999;
        ++icluster;//test
        }
        cout<<"total number of clusters "<<icluster<<endl;

}


void PHG4InitZVertexing::find_track_candidates_z_init(unsigned int zoomlevel){

          for (std::map<unsigned int,HelixHoughBin* >::iterator it=bins_map_sel.begin(); it!=bins_map_sel.end(); ++it){ bins_map_sel.erase(it);}
	bins_map_sel.clear();

	unsigned int max_i=-1;
	unsigned int max_count = 0;

        for(unsigned int i =0; i<_hough_space->get_n_z0_bins(0); ++i){
        for(unsigned int j =0; j<_hough_space->get_n_dzdl_bins(0); ++j ){
                unsigned int bins[5]={0,0,0,j,i};
                bin = _hough_space->get_bin(0,bins);
                unsigned int count = bins_map.find(bin)->second->get_count();
                if (count <1) continue;
		if (count>max_count) {
			max_i=i;
			max_count=count;
		}
#ifdef _DEBUG_
                cout<<"bin "<<bin<<endl;
                cout<<"z0 bin "<<i<<" dzdl bin "<<j<< " count "<< count<<endl;
                cout<<"z0 bin "<< bins_map.find(bin)->second->get_z0_bin(0)<<" dzdl bin "<<bins_map.find(bin)->second->get_dzdl_bin(0)<<endl;

                _z0_dzdl->Fill(bins_map.find(bin)->second->get_z0_center(zoomlevel),bins_map.find(bin)->second->get_dzdl_center(zoomlevel), count);
#endif
//                if (count>= _min_nlayers)
//		{
//                bins_map_sel.insert(std::pair<unsigned int, HelixHoughBin*>(bin,bins_map.find(bin)->second) );
//              cout<<"inserting selected bin.. "<<endl;
		
//		print cluster ids
//              for(HelixHoughBin::ClusterIter iter = bins_map.find(bin)->second->begin_clusters();
//                      iter != bins_map.find(bin)->second->end_clusters();
//                      iter++){
//	      cout<<"z0 bin "<< bins_map.find(bin)->second->get_z0_bin(0)<<" dzdl bin "<<bins_map.find(bin)->second->get_dzdl_bin(0)<<endl;
//              cout<< "cluster_id "<<*iter<<endl;
//              }       
//                }                                                                                                                       
          }
        }

	for (unsigned int k =0; k<_hough_space->get_n_dzdl_bins(0); ++k){
                unsigned int bins[5]={0,0,0,k,max_i};
                bin = _hough_space->get_bin(0,bins);
                unsigned int count = bins_map.find(bin)->second->get_count();
                if (count <1) continue;

		if (count>= _min_nlayers)
                {
                bins_map_sel.insert(std::pair<unsigned int, HelixHoughBin*>(bin,bins_map.find(bin)->second) );
		}

	}
        cout<<"bins_map_sel.size at zoom "<<zoomlevel << " : (find_track_candidates_z_init) " <<bins_map_sel.size()<<endl;

}


void PHG4InitZVertexing::vote_z(unsigned int zoomlevel){

            set_nbins(zoomlevel); // at next zoom, bin number changes

            float hitpos3d[3];
            unsigned int cluster_id=-999;
            for (std::map<unsigned int,HelixHoughBin*>::iterator it=bins_map_prev.begin(); it!=bins_map_prev.end(); ++it)
            {
                bin = it->first;
                HelixHoughBin* houghbin = it->second;

#ifdef _DEBUG_
                unsigned int izprev = houghbin->get_z0_bin(zoomlevel-1);
                unsigned int ilprev = houghbin->get_dzdl_bin(zoomlevel-1);
                unsigned int ipprev = houghbin->get_phi_bin(zoomlevel-1);

//                unsigned int ikprev = houghbin->get_kappa_bin(zoomlevel-1);
//                unsigned int idprev = houghbin->get_d_bin(0);
//                bool fillhisto = (izprev==5) && (ilprev==11)&& (ipprev==21);
//                if (fillhisto) cout<<"bin "<<bin<<" ik " <<ikprev<<" ip "<<ipprev<<" id "<<idprev<<endl;
#endif
                for(HelixHoughBin::ClusterIter iter = bins_map_prev.find(bin)->second->begin_clusters();
                        iter != bins_map_prev.find(bin)->second->end_clusters();
                        ++iter){
                cluster_id = *iter;

                Cluster3D hit = hits_map.find(cluster_id)->second;
                hitpos3d[0] = hit.get_x();
                hitpos3d[1] = hit.get_y();
                hitpos3d[2] = hit.get_z();
#ifdef _DEBUG_
//                if (fillhisto) cout<<"cluster_id "<<cluster_id<<endl;
#endif

                std::vector<float> kappa_phi_d_ranges;
                float kappa_min = houghbin->get_kappa_center(zoomlevel-1)-0.5*_hough_space->get_kappa_bin_size(zoomlevel-1);
                float kappa_max = houghbin->get_kappa_center(zoomlevel-1)+0.5*_hough_space->get_kappa_bin_size(zoomlevel-1);
                float phi_min = houghbin->get_phi_center(zoomlevel-1)-0.5*_hough_space->get_phi_bin_size(zoomlevel-1);
                float phi_max = houghbin->get_phi_center(zoomlevel-1)+0.5*_hough_space->get_phi_bin_size(zoomlevel-1);
                float d_min = houghbin->get_d_center(zoomlevel-1)-0.5*_hough_space->get_d_bin_size(zoomlevel-1);
                float d_max = houghbin->get_d_center(zoomlevel-1)+0.5*_hough_space->get_d_bin_size(zoomlevel-1);
                kappa_phi_d_ranges.push_back(kappa_min);
                kappa_phi_d_ranges.push_back(kappa_max);
                kappa_phi_d_ranges.push_back(phi_min);
                kappa_phi_d_ranges.push_back(phi_max);
                kappa_phi_d_ranges.push_back(d_min);
                kappa_phi_d_ranges.push_back(d_max);
                float dzdl_min = houghbin->get_dzdl_center(zoomlevel-1)-0.5*_hough_space->get_dzdl_bin_size(zoomlevel-1);
                float dzdl_max = houghbin->get_dzdl_center(zoomlevel-1)+0.5*_hough_space->get_dzdl_bin_size(zoomlevel-1);
#ifdef _DEBUG_
//                if (fillhisto){
//                cout<<"kappa_min "<<kappa_min <<" kappa_max " <<kappa_max<<" phi_min "<<phi_min<<" phi_max "<<phi_max<<endl;
//                cout<<"d_min "<<d_min<< " d_max "<<d_max<<" dzdl_min "<<dzdl_min<<" dzdl_max "<<dzdl_max<<endl;
//                }
#endif

                for (iz=0; iz < _hough_space->get_n_z0_bins(zoomlevel); ++iz ) {
                std::vector<float> z0_range;
                float z0_min = houghbin->get_z0_center(zoomlevel-1)-0.5*_hough_space->get_z0_bin_size(zoomlevel-1)+iz*_hough_space->get_z0_bin_size(zoomlevel);
                float z0_max = z0_min + _hough_space->get_z0_bin_size(zoomlevel);;
                z0_range.push_back(z0_min);
                z0_range.push_back(z0_max);
                float dzdl_range[2];
#ifdef _DEBUG_
//                if (fillhisto) cout<<"z0 min "<<z0_min<<" z0 max "<<z0_max<<endl;
#endif
                unsigned int dzdl_bin_range[2];
                //compute min_dzdl & max_dzdl for given iz 
                _hough_funcs->set_current_zoom(zoomlevel);
                _hough_funcs->calculate_dzdl_range(hitpos3d,z0_range, kappa_phi_d_ranges, dzdl_range);
                dzdl_bin_range[0] = _hough_space->get_dzdl_bin(zoomlevel,dzdl_range[0]);
                dzdl_bin_range[1] = _hough_space->get_dzdl_bin(zoomlevel,dzdl_range[1]);

                        if (dzdl_range[0] > dzdl_max) continue;
                        else if (dzdl_range[0] <dzdl_min) dzdl_bin_range[0]=0;
                        if (dzdl_range[1] < dzdl_min) continue;
                        else if (dzdl_range[1] >dzdl_max) dzdl_bin_range[1]=_hough_space->get_n_dzdl_bins(zoomlevel)-1;
#ifdef _DEBUG_
//                       if (fillhisto)cout<<" dzdl_bin_range[0] " <<dzdl_bin_range[0]<<" dzdl_bin_range[1] "<<dzdl_bin_range[1]<<endl;
#endif
                        for (il =dzdl_bin_range[0]; il<dzdl_bin_range[1]+1; ++il) {
                                unsigned int zbins[5] = {0,0,0,il,iz};
                                unsigned int curbin = _hough_space->get_bin(zoomlevel,zbins);
                                houghbin->set_bin(zoomlevel,curbin);
                                houghbin->set_bins(zoomlevel,curbin);
                                unsigned int curgbin = houghbin->get_global_bin(zoomlevel);
#ifdef _DEBUG_
//                                if (fillhisto) cout<<"global bin "<<curgbin<<" z0 bin "<<houghbin->get_z0_bin(zoomlevel)<<" dzdl bin "<<houghbin->get_dzdl_bin(zoomlevel) <<endl;
#endif
                                auto search = bins_map_cur.find(curgbin);
                                if(search != bins_map_cur.end())
                                {
                                bins_map_cur.find(curgbin)->second->add_cluster_ID(cluster_id);
                                }
                                else
                                {
                                HelixHoughBin* cur_bin = houghbin->Clone();
				cur_bin->clear_clusters();
                                cur_bin->set_zoomlevel(zoomlevel);
                                cur_bin->set_bin(zoomlevel,curbin);
                                cur_bin->set_bins(zoomlevel,curbin);
                                cur_bin->add_cluster_ID(cluster_id);
                                bins_map_cur.insert(std::pair<unsigned int,HelixHoughBin*>(curgbin,cur_bin));
#ifdef _DEBUG_
//        			if (fillhisto) cout<<"global bin "<<cur_bin->get_global_bin(zoomlevel)<<" z0 bin "<<cur_bin->get_z0_bin(zoomlevel)<<" dzdl bin "<<cur_bin->get_dzdl_bin(zoomlevel) <<endl;
#endif
                                }
#ifdef _DEBUG_
        if (fillhisto){
//        cout<<" il "<<il<<" iz "<<iz<<" il from bin"<< bins_map_cur.find(curgbin)->second->get_dzdl_bin(zoomlevel)<<" iz from bin "<< bins_map_cur.find(curgbin)->second->get_z0_bin(zoomlevel)<<endl;

        unsigned int count = bins_map_cur.find(curgbin)->second->get_count();
//        cout<<"count "<<count<<endl;
//        cout<<"z0 "<<bins_map_cur.find(curgbin)->second->get_z0_center(zoomlevel)<<" dzdl "<< bins_map_cur.find(curgbin)->second->get_dzdl_center(zoomlevel)<<endl;
        _z0_dzdl->Fill(bins_map_cur.find(curgbin)->second->get_z0_center(zoomlevel),bins_map_cur.find(curgbin)->second->get_dzdl_center(zoomlevel), count);
        }
#endif
                        } // il
                } //iz
        hitpos3d[0]=-999; hitpos3d[1]= -999; hitpos3d[2]=-999;
	kappa_phi_d_ranges.clear();
                }//clusters
            } //bins_map_prev
            for (std::map<unsigned int,HelixHoughBin* >::iterator it=bins_map_prev.begin(); it!=bins_map_prev.end(); ++it){ bins_map_prev.erase(it);}
	    bins_map_prev.clear();

	    cout<<"bins_map_cur.size at zoom "<<zoomlevel << " (vote_z) : " <<bins_map_cur.size()<<endl;
}


void PHG4InitZVertexing::find_track_candidates_z(unsigned int zoomlevel){
        for (std::map<unsigned int,HelixHoughBin* >::iterator it=bins_map_sel.begin(); it!=bins_map_sel.end(); ++it){ bins_map_sel.erase(it);}
	bins_map_sel.clear();


	for (std::map<unsigned int,HelixHoughBin*>::iterator it = bins_map_cur.begin(); it != bins_map_cur.end(); ++it)
	{
                bin = it->first;
                HelixHoughBin* houghbin = it->second;
                unsigned int count = houghbin->get_count();
                if (count <1) continue;

                il = houghbin->get_dzdl_bin(zoomlevel);
                iz = houghbin->get_z0_bin(zoomlevel);

#ifdef _DEBUG_
//                ik = houghbin->get_kappa_bin(zoomlevel);
//                ip = houghbin->get_phi_bin(zoomlevel);
//                id = houghbin->get_d_bin(zoomlevel);

//                cout<<"bin "<<bin<<endl;
//                cout<<"kappa bin "<<ik<<" phi bin "<<ip<< " d bin "<<id <<" il "<<il<<" iz "<<iz<< " count "<< count<<endl;
#endif
                if (count>= _min_nlayers)
                {
                bins_map_sel.insert(std::pair<unsigned int, HelixHoughBin*>(bin,bins_map_cur.find(bin)->second) );
#ifdef _DEBUG_
		cout<<"il "<<il<<" iz "<<iz<< " count "<< count<<endl;
#endif
		}

	}

        for (std::map<unsigned int,HelixHoughBin* >::iterator it=bins_map_cur.begin(); it!=bins_map_cur.end(); ++it){ bins_map_cur.erase(it);}
	bins_map_cur.clear();

        cout<<"bins_map_sel.size at zoom "<< zoomlevel<<" (find_track_candidates_z) : " <<bins_map_sel.size()<<endl;
}

void PHG4InitZVertexing::vote_xy(unsigned int zoomlevel){
	
	int count = -1;
#ifdef _DEBUG_
	unsigned int ilfill= -1;
	unsigned int izfill= -1;	
#endif
	set_nbins(zoomlevel);
        // hand select bins over to bins_map_cur
        kappa_phi_map.clear();
        d_phi_map.clear();
        unsigned int phi_bin_range[2];
        for (std::map<unsigned int,HelixHoughBin*>::iterator it=bins_map_sel.begin(); it!=bins_map_sel.end(); ++it)
        {
        bin = it->first;
        HelixHoughBin* houghbin = it->second;
        il = houghbin->get_dzdl_bin(zoomlevel);
        iz = houghbin->get_z0_bin(zoomlevel);
	++count;

#ifdef _DEBUG_
	if (count==0){  ilfill = il; izfill = iz;}

        unsigned int ikprev =0;
	unsigned int ipprev =0;
	unsigned int idprev =0;
	if (zoomlevel>0){
	ikprev = houghbin->get_kappa_bin(zoomlevel-1);
        ipprev = houghbin->get_phi_bin(zoomlevel-1);
        idprev = houghbin->get_d_bin(zoomlevel-1);
	}
	
        cout<<"bin "<< bin<<endl;
        cout<<"il " <<il<<" iz "<<iz<<" ikprev "<<ikprev<<" ipprev "<<ipprev<<" idprev "<<idprev<<" count"<<houghbin->get_count()<<endl;
#endif
        
	float phi_min;
        if (zoomlevel==0) phi_min = _hough_space->get_phi_min();
        else
        phi_min = houghbin->get_phi_center(zoomlevel-1)-0.5*_hough_space->get_phi_bin_size(zoomlevel-1);
        float phi_max;
        if (zoomlevel==0) phi_max = _hough_space->get_phi_max();
        else
        phi_max = phi_min + _hough_space->get_phi_bin_size(zoomlevel-1);

#ifdef _DEBUG_
//	cout<<"phi_min "<<phi_min <<" phi_max "<<phi_max<<endl;
#endif
                for (HelixHoughBin::ClusterIter iter = houghbin->begin_clusters();
                        iter != houghbin->end_clusters(); ++iter){
                        cluster_id = *iter;

			Cluster3D hit = hits_map.find(cluster_id)->second;
                        hitpos3d[0] = hit.get_x();
                        hitpos3d[1] = hit.get_y();
#ifdef _DEBUG_
                      cout<<"cluster id "<<cluster_id<<" x "<< hitpos3d[0]<<" y "<<hitpos3d[1]<<endl;
#endif
                        for (ik=0; ik<nkappa; ++ik){
                        for (id=0; id<nd; ++id){
                        std::vector<float> kappa_d_ranges;

			float kappa_min;
			if (zoomlevel == 0) kappa_min = _hough_space->get_kappa_min();
			else
			kappa_min= houghbin->get_kappa_center(zoomlevel-1)-0.5*_hough_space->get_kappa_bin_size(zoomlevel-1);
			kappa_min += ik*_hough_space->get_kappa_bin_size(zoomlevel);
			float kappa_max = kappa_min + _hough_space->get_kappa_bin_size(zoomlevel);
			float d_min;
			if (zoomlevel==0) d_min = _hough_space->get_d_min();
			else d_min = houghbin->get_d_center(zoomlevel-1)-0.5*_hough_space->get_d_bin_size(zoomlevel-1);
			d_min += id*_hough_space->get_d_bin_size(zoomlevel);
			float d_max = d_min + _hough_space->get_d_bin_size(zoomlevel);

                        kappa_d_ranges.push_back(kappa_min);
			kappa_d_ranges.push_back(kappa_max);
			kappa_d_ranges.push_back(d_min);
                        kappa_d_ranges.push_back(d_max);
                        _hough_funcs->set_current_zoom(zoomlevel);

			if (!separate_helicity){
			float phi_r_range[2];
			float phi_l_range[2];
                        unsigned int phi_r_bin_range[2];
                        unsigned int phi_l_bin_range[2];
                        //cout<<"before calc phi "<<endl;
                        _hough_funcs->calculate_phi_range(hitpos3d,kappa_d_ranges, phi_r_range,phi_l_range);
//                      cout<<"phi_l_min "<<phi_l_range[0]<<" phi_l_max "<<phi_l_range[1]<<endl;
//                      cout<<"phi_r_min "<<phi_r_range[1]<<" phi_r_max "<<phi_r_range[1]<<endl;
	 		phi_r_bin_range[0] = _hough_space->get_phi_bin(zoomlevel,phi_r_range[0]);	
			phi_r_bin_range[1] = _hough_space->get_phi_bin(zoomlevel,phi_r_range[1]);
			phi_l_bin_range[0] = _hough_space->get_phi_bin(zoomlevel,phi_l_range[0]);
			phi_l_bin_range[1] = _hough_space->get_phi_bin(zoomlevel,phi_l_range[1]);
//		cout<<"phi_r_bin_min "<<phi_r_bin_range[0]<<" phi_r_bin_max "<<phi_r_bin_range[1]<<endl;

			bool phi_r_out_of_range=false;
			bool phi_l_out_of_range=false;
                        if (phi_r_range[0] > phi_max) phi_r_out_of_range=true;
                        else if (phi_r_range[0] <phi_min) phi_r_bin_range[0]=0;
                        if (phi_r_range[1] < phi_min) phi_r_out_of_range=true;
                        else if (phi_r_range[1] >phi_max) phi_r_bin_range[1]=_hough_space->get_n_phi_bins(zoomlevel)-1;

                        if (phi_l_range[0] > phi_max) phi_r_out_of_range=true;
                        else if (phi_r_range[0] <phi_min) phi_r_bin_range[0]=0;
                        if (phi_r_range[1] < phi_min) phi_r_out_of_range=true;
                        else if (phi_r_range[1] >phi_max) phi_r_bin_range[1]=_hough_space->get_n_phi_bins(zoomlevel)-1;

			for (int i = 0; i<2;++i){// helicity
				//cout<<"i "<<i<<endl;
				if ((i==0&&phi_r_out_of_range) || (i==1&&phi_l_out_of_range)) continue;
                                switch (i){
                                case 0:
                                phi_bin_range[0]=phi_r_bin_range[0];
                                phi_bin_range[1]=phi_r_bin_range[1];
                                break;
                                case 1:
                                phi_bin_range[0]=phi_l_bin_range[0];
                                phi_bin_range[1]=phi_l_bin_range[1];
                                break;
                                }
//				cout<<"phi_bin_range[0] " <<phi_bin_range[0]<<" phi_bin_range[1] "<<phi_bin_range[1]<<endl;

                                for (ip= phi_bin_range[0]; ip<phi_bin_range[1]+1;++ip){
                                unsigned int xyzbins[5]={ik,ip,id,il,iz};
                                unsigned int curbin = _hough_space->get_bin(zoomlevel,xyzbins);
                                houghbin->set_bin(zoomlevel,curbin);
                                houghbin->set_bins(zoomlevel,curbin);
                                unsigned int curgbin = houghbin->get_global_bin(zoomlevel);
//                              cout<<" curgbin "<<curgbin<<endl;
                                auto search = bins_map_cur.find(curgbin);
                                if(search != bins_map_cur.end())
                                {
                                bins_map_cur.find(curgbin)->second->add_cluster_ID(cluster_id);
                                }
                                else
                                {
                                HelixHoughBin* cur_bin = houghbin->Clone();
				cur_bin->clear_clusters();
                                cur_bin->set_zoomlevel(zoomlevel);
                                cur_bin->set_bin(zoomlevel,curbin);
                                cur_bin->set_bins(zoomlevel,curbin);
                                cur_bin->add_cluster_ID(cluster_id);
                                bins_map_cur.insert(std::pair<unsigned int,HelixHoughBin*>(curgbin,cur_bin));
//                              cout<<" il "<<il<<" iz "<<iz<<" il from bin"<< bins_map_cur.find(curbin)->second->get_dzdl_bin(0)<<" iz from bin "<< bins_map_cur.find(curbin)->second->get_z0_bin(0)<<endl;
                                }

				auto search0 = kappa_d_phi_map.find(curgbin);
				if(search0 != kappa_d_phi_map.end())
				{
				kappa_d_phi_map.find(curgbin)->second = kappa_d_phi_map.find(curgbin)->second + 1;
				}
				else
				{
				kappa_d_phi_map.insert(std::pair<unsigned int,unsigned int>(curgbin,1));
				}

                                unsigned int kpbins[5]={ik,ip,0,il,iz};
                                unsigned int kpbin = _hough_space->get_bin(zoomlevel,kpbins);
                                auto search1 = kappa_phi_map.find(kpbin);
                                if(search1 != kappa_phi_map.end())
                                {
                                kappa_phi_map.find(kpbin)->second = kappa_phi_map.find(kpbin)->second +1;
                                }
                                else
                                {
                                kappa_phi_map.insert(std::pair<unsigned int,unsigned int>(kpbin,1));
                                }

                                unsigned int dpbins[5]={0,ip,id,il,iz};
                                unsigned int dpbin = _hough_space->get_bin(zoomlevel,dpbins);
                                auto search2 = d_phi_map.find(dpbin);
                                if (search2 != d_phi_map.end())
                                {
                                d_phi_map.find(dpbin)->second = d_phi_map.find(dpbin)->second + 1;
                                }
                                else
                                {
                                d_phi_map.insert(std::pair<unsigned int,unsigned int>(dpbin,1));
                                }

#ifdef _DEBUG_
 			if (il ==ilfill && iz == izfill){
//                      for(HelixHoughBin::ClusterIter iter = bins_map.find(bin)->second->begin_clusters();
//                         iter != bins_map.find(bin)->second->end_clusters();
//                         ++iter){
//                         cout<< "cluster_id "<<*iter<<endl;
//                      }                                                          

			cout << "il "<<il<<" iz "<<iz<<endl;
                        _kappa_phi->Fill(bins_map_cur.find(curgbin)->second->get_kappa_center(zoomlevel),bins_map_cur.find(curgbin)->second->get_phi_center(zoomlevel),1);
                        _d_phi->Fill(bins_map_cur.find(curgbin)->second->get_d_center(zoomlevel),bins_map_cur.find(curgbin)->second->get_phi_center(zoomlevel),1);
			_kappa_d_phi->Fill(bins_map_cur.find(curgbin)->second->get_kappa_center(zoomlevel),bins_map_cur.find(curgbin)->second->get_d_center(zoomlevel),bins_map_cur.find(curgbin)->second->get_phi_center(zoomlevel),1);
                        }
#endif
                                }//ip                   
			}// helicity

			}// if not separating helicities 
			else 
			{
                        float phi_range[2];
			float phi_prev_range[2];
			float phi_next_range[2];
                        unsigned int phi_bin_range[2];
                        //cout<<"before calc phi "<<endl;

                        if (id==0)
                        _hough_funcs->calculate_phi_range(hitpos3d,kappa_d_ranges, helicity, phi_range,phi_next_range);
			else 
			_hough_funcs->calculate_phi_range(hitpos3d,kappa_d_ranges, helicity, phi_range,phi_prev_range,phi_next_range);
		
			phi_prev_range[0]= phi_next_range[0];
                        phi_prev_range[1]= phi_next_range[1];	
			
                      //cout<<"phi_min "<<phi_range[0]<<" phi_max "<<phi_range[1]<<endl;
                        phi_bin_range[0] = _hough_space->get_phi_bin(zoomlevel,phi_range[0]);
                        phi_bin_range[1] = _hough_space->get_phi_bin(zoomlevel,phi_range[1]);
              	//cout<<"phi_bin_min "<<phi_bin_range[0]<<" phi_bin_max "<<phi_bin_range[1]<<endl;

                        bool phi_out_of_range=false;
                        if (phi_range[0] > phi_max) phi_out_of_range=true;
                        else if (phi_range[0] <phi_min) phi_bin_range[0]=0;
                        if (phi_range[1] < phi_min) phi_out_of_range=true;
                        else if (phi_range[1] >phi_max) phi_bin_range[1]=_hough_space->get_n_phi_bins(zoomlevel)-1;

			if (phi_out_of_range) continue;

                                for (ip= phi_bin_range[0]; ip<phi_bin_range[1]+1;++ip){
                                unsigned int xyzbins[5]={ik,ip,id,il,iz};
                                unsigned int curbin = _hough_space->get_bin(zoomlevel,xyzbins);
                                houghbin->set_bin(zoomlevel,curbin);
                                houghbin->set_bins(zoomlevel,curbin);
                                unsigned int curgbin = houghbin->get_global_bin(zoomlevel);
#ifdef _DEBUG_
                                cout<<" curgbin "<<curgbin<<endl;
#endif
                                auto search = bins_map_cur.find(curgbin);
                                if(search != bins_map_cur.end())
                                {
                                bins_map_cur.find(curgbin)->second->add_cluster_ID(cluster_id);
                                }
                                else
                                {
                                HelixHoughBin* cur_bin = houghbin->Clone();
                                cur_bin->clear_clusters();
                                cur_bin->set_zoomlevel(zoomlevel);
                                cur_bin->set_bin(zoomlevel,curbin);
                                cur_bin->set_bins(zoomlevel,curbin);
                                cur_bin->add_cluster_ID(cluster_id);
                                bins_map_cur.insert(std::pair<unsigned int,HelixHoughBin*>(curgbin,cur_bin));
#ifdef _DEBUG_
                                cout<<" il "<<il<<" iz "<<iz<<" il from bin"<< bins_map_cur.find(curbin)->second->get_dzdl_bin(0)<<" iz from bin "<< bins_map_cur.find(curbin)->second->get_z0_bin(0)<<endl;
#endif
                                }

#ifdef _DEBUG_
//                        if (il ==ilfill && iz == izfill){

//                      for(HelixHoughBin::ClusterIter iter = bins_map.find(bin)->second->begin_clusters();
//                         iter != bins_map.find(bin)->second->end_clusters();
//                         ++iter){
//                         cout<< "cluster_id "<<*iter<<endl;
//                      }
//                        cout << "il "<<il<<" iz "<<iz<<endl;

//                        _kappa_phi->Fill(bins_map_cur.find(curgbin)->second->get_kappa_center(zoomlevel),bins_map_cur.find(curgbin)->second->get_phi_center(zoomlevel),1);
//                        _d_phi->Fill(bins_map_cur.find(curgbin)->second->get_d_center(zoomlevel),bins_map_cur.find(curgbin)->second->get_phi_center(zoomlevel),1);
//                        _kappa_d_phi->Fill(bins_map_cur.find(curgbin)->second->get_kappa_center(zoomlevel),bins_map_cur.find(curgbin)->second->get_d_center(zoomlevel),bins_map_cur.find(curgbin)->second->get_phi_center(zoomlevel),1);
//			}
#endif
				}//ip

			}// if separating helicities
			
			kappa_d_ranges.clear();
                        }//id
                        }//ik bins
                }//cluster

        }//binsmap_sel
        cout<<"bins_map_cur.size at zoom "<<zoomlevel <<" (vote_xy) : " <<bins_map_cur.size()<<endl;
          for (std::map<unsigned int,HelixHoughBin* >::iterator it=bins_map_sel.begin(); it!=bins_map_sel.end(); ++it){ bins_map_sel.erase(it);}
	bins_map_sel.clear();

}

void PHG4InitZVertexing::find_track_candidates_xy(unsigned int zoomlevel){

        for (std::map<unsigned int,HelixHoughBin* >::iterator it=bins_map_prev.begin(); it!=bins_map_prev.end(); ++it){ bins_map_prev.erase(it);}
	bins_map_prev.clear();

        for (std::map<unsigned int,HelixHoughBin*>::iterator it=bins_map_cur.begin(); it!=bins_map_cur.end(); ++it)
        {
                bin = it->first;
                HelixHoughBin* houghbin = it->second;
                unsigned int count = houghbin->get_count();
                if (count <1) continue;

                il = houghbin->get_dzdl_bin(zoomlevel);
                iz = houghbin->get_z0_bin(zoomlevel);
                ik = houghbin->get_kappa_bin(zoomlevel);
                ip = houghbin->get_phi_bin(zoomlevel);
                id = houghbin->get_d_bin(zoomlevel);

#ifdef _DEBUG_
                cout<<"bin "<<bin<<endl;
                cout<<"kappa bin "<<ik<<" phi bin "<<ip<< " d bin "<<id <<" il "<<il<<" iz "<<iz<< " count "<< count<<endl;
#endif
		if (count >= _min_nlayers 
//&& count >= cnt_kappa_up_d_phi && count >= cnt_kappa_dw_d_phi && count >= cnt_kappa_d_dw_phi && count >= cnt_kappa_d_up_phi && cnt_kappa_d_phi >= cnt_kappa_d_phi_up && count >= cnt_kappa_d_phi_dw
  		)
		{
                bins_map_prev.insert(std::pair<unsigned int, HelixHoughBin*>(bin,bins_map_cur.find(bin)->second) );
#ifdef _DEBUG_
                cout<<"inserting selected bin.. "<<endl;
                cout<<"bin "<<bin<<endl;
                cout<<"kappa bin "<<ik<<" phi bin "<<ip<< " d bin "<<id <<" il "<<il<<" iz "<<iz<< " count "<< count<<endl;
#endif
		}
		

        }

        cout<<"bins_map_prev.size at zoom "<<zoomlevel <<" (find_track_candidates_xy) " <<bins_map_prev.size()<<endl;
          for (std::map<unsigned int,HelixHoughBin* >::iterator it=bins_map_cur.begin(); it!=bins_map_cur.end(); ++it){ bins_map_cur.erase(it);}
	bins_map_cur.clear();

}

void PHG4InitZVertexing::prune_z(unsigned int zoomlevel){
	set_nbins(zoomlevel);
	for (std::map<unsigned int,HelixHoughBin*>::iterator it=bins_map_sel.begin(); it!=bins_map_sel.end(); ++it){
		HelixHoughBin* houghbin = it->second;
		unsigned int nclust = houghbin->get_count();
                unsigned int globalbin = houghbin->get_global_bin(zoomlevel);
#ifdef _DEBUG_
		cout<<"global bin "<<houghbin->get_global_bin(zoomlevel)<<endl;

                if (zoomlevel>0)
                cout<<"ikprev "<<houghbin->get_kappa_bin(zoomlevel-1)<<" ipprev "<<houghbin->get_phi_bin(zoomlevel-1)<<" idprev "<< houghbin->get_d_bin(zoomlevel-1) <<" ilprev "<< houghbin->get_dzdl_bin(zoomlevel-1)<<" izprev "<<houghbin->get_z0_bin(zoomlevel-1)<<endl; 
		cout<<"il "<<houghbin->get_dzdl_bin(zoomlevel)<<" iz "<<houghbin->get_z0_bin(zoomlevel)<<" count "<<houghbin->get_count()<<endl;
#endif
		unsigned int var[3] = {1<<3,1<<4,(1<<3)+(1<<4)};
		unsigned int sign[3][4] ={0<<3,1<<3,0,0,0<<4,1<<4,0,0,(0<<3)+(0<<4),(0<<3)+(1<<4),(1<<3)+(0<<4),(1<<3)+(1<<4)};
		for (unsigned int lz =0; lz<3; ++lz){
		for (unsigned int ss = 0; ss<4; ++ss ){
			if (lz<2 && ss>1) continue;
			unsigned int neighborbin = houghbin->get_neighbors_global_bin(zoomlevel,var[lz],sign[lz][ss]);
#ifdef _DEBUG_
			cout<<"var "<<var[lz]<<" ss "<<ss<<" neighbor_global bin "<<neighborbin<<endl;
#endif
			if (globalbin==neighborbin) continue;
			auto search = bins_map_sel.find(neighborbin);
			if (search != bins_map_sel.end()){
#ifdef _DEBUG_
			cout<<"neighbor bin found"<<endl;
#endif
			unsigned int nclust_neighbor = bins_map_sel.find(neighborbin)->second->get_count();
			if(nclust >= nclust_neighbor && nclust <= 1.1*_nlayers) {
#ifdef _DEBUG_
			cout<<"merging "<<endl;
#endif
			for(HelixHoughBin::ClusterIter iter = bins_map_sel.find(neighborbin)->second->begin_clusters();
                        iter != bins_map_sel.find(neighborbin)->second->end_clusters();
                        ++iter){
			houghbin->add_cluster_ID(*iter);
			}
			}//nclust
			}//search
			bins_map_sel.erase(neighborbin);
		}//ss
		}//lz

//		}//i

	}

        cout<<"bins_map_sel.size at zoom " <<zoomlevel<<" (prune_z) : " <<bins_map_sel.size()<<endl;	
}

void PHG4InitZVertexing::prune_xy(unsigned int zoomlevel){

        set_nbins(zoomlevel);
        cout<<"bins_map_prev.size at zoom " <<zoomlevel<<" (pre prune_xy) : " <<bins_map_prev.size()<<endl;
        for (std::map<unsigned int,HelixHoughBin*>::iterator it=bins_map_prev.begin(); it!=bins_map_prev.end(); ++it){
                HelixHoughBin* houghbin = it->second;
                unsigned int nclust = houghbin->get_count();
		unsigned int globalbin = houghbin->get_global_bin(zoomlevel);
#ifdef _DEBUG_
		cout<<"global bin "<<globalbin<<endl;

		if (zoomlevel>0) 
		cout<<"ikprev "<<houghbin->get_kappa_bin(zoomlevel-1)<<" ipprev "<<houghbin->get_phi_bin(zoomlevel-1)<<" idprev "<< houghbin->get_d_bin(zoomlevel-1) <<" ilprev "<< houghbin->get_dzdl_bin(zoomlevel-1)<<" izprev "<<houghbin->get_z0_bin(zoomlevel-1)<<endl;
		cout<<"ik "<<houghbin->get_kappa_bin(zoomlevel)<<" ip "<<houghbin->get_phi_bin(zoomlevel)<<" id "<<houghbin->get_d_bin(zoomlevel)<<" il "<<houghbin->get_dzdl_bin(zoomlevel)<<" iz "<<houghbin->get_z0_bin(zoomlevel)<<" count "<<houghbin->get_count()<<endl;
#endif
		// k, p, d, kd, kp, dp
                unsigned int var[6] = {1,1<<1,1<<2,(1)+(1<<1),(1)+(1<<2),(1<<1)+(1<<2)};
		unsigned int sign[6][4]={0,1,0,0,
					 0<<1,1<<1,0,0,
					 0<<2,1<<2,0,0,
					 0,1<<1,1,(1)+(1<<1),
					 0,1<<2,1, (1)+(1<<2),
					 0,1<<2,1<<1,(1<<1)+(1<<2)};
                for (unsigned int kpd =0; kpd<6; ++kpd){
                for (unsigned int ss = 0; ss<4; ++ss ){
			if (kpd<3 && ss>1) continue;
                        unsigned int neighborbin = houghbin->get_neighbors_global_bin(zoomlevel,var[kpd],sign[kpd][ss]);
//                      cout<<"var "<<var[kpd]<<" ss "<<ss<<" neighbor_global bin "<<neighborbin<<endl;
			if (globalbin==neighborbin) continue;
                        auto search = bins_map_prev.find(neighborbin);
                        if (search != bins_map_prev.end()){
//                      cout<<"neighbor bin found"<<endl;
                        unsigned int nclust_neighbor = bins_map_prev.find(neighborbin)->second->get_count();
                        if(nclust >= nclust_neighbor && nclust <= 1.1*_nlayers) {
#ifdef _DEBUG_
                      cout<<"merging "<<endl;
#endif
                        for(HelixHoughBin::ClusterIter iter = bins_map_prev.find(neighborbin)->second->begin_clusters();
                        iter != bins_map_prev.find(neighborbin)->second->end_clusters();
                        ++iter){
                        houghbin->add_cluster_ID(*iter);
                        }
                        }//nclust
                        }//search
                        bins_map_prev.erase(neighborbin);
                }//ss
                }//lz

//              }//i

        }

        cout<<"bins_map_prev.size at zoom " <<zoomlevel<<" (prune_xy) : " <<bins_map_prev.size()<<endl;

}

void PHG4InitZVertexing::reset_hits_used(){

	for (std::map<unsigned int,bool>::iterator it = hits_used.begin(); 
		it != hits_used.end(); ++it){
		it->second = false;
	}
}

void PHG4InitZVertexing::bins_to_Track3D(std::vector<Track3D>& new_tracks, int imap, unsigned int zoomlevel){

	unsigned int icluster = 0;
	switch (imap){
		case 0: // not implemented yet		 
		for (std::map<unsigned int,HelixHoughBin*>::iterator it = bins_map_sel.begin(); it != bins_map_sel.end(); ++it){
			HelixHoughBin *houghbin = it->second;
			Track3D track;
			icluster = 0;
			for(HelixHoughBin::ClusterIter iter = houghbin->begin_clusters();
                        	iter != houghbin->end_clusters();
                        	++iter){
                	//cout<< "cluster_id "<<*iter<<endl;    
                	Cluster3D hit = hits_map.find(*iter)->second;
               		track.hits.push_back(hit);
			track.hits[icluster].set_id(*iter);
			++icluster;
			}
			_tracks.push_back(track);
		}
		break;

		case 1:
		cout<<"bins_map_prev selected, total bins : "<< bins_map_prev.size()<<endl;
		for (std::map<unsigned int,HelixHoughBin*>::iterator it = bins_map_prev.begin(); it != bins_map_prev.end(); ++it){
			HelixHoughBin *houghbin = it->second;	
			Track3D track;
			icluster = 0;
			track.hits.assign(houghbin->get_count(),Cluster3D());
#ifdef _DEBUG_
//			cout<<"bin "<< houghbin->get_global_bin(zoomlevel);
//			cout<<" : start loop over clusters "<<endl;
#endif
			for(HelixHoughBin::ClusterIter iter = houghbin->begin_clusters();
				iter != houghbin->end_clusters();
				++iter){
#ifdef _DEBUG_
//			cout <<"cluster_id "<< *iter<<endl;
#endif
			auto search =  hits_map.find(*iter);
			if (search != hits_map.end())
			{
			Cluster3D hit = hits_map.find(*iter)->second;
			track.hits[icluster] = hit;
			track.hits[icluster].set_id(*iter);
			track.kappa = houghbin->get_kappa_center(zoomlevel);
			track.phi = houghbin->get_phi_center(zoomlevel);
			track.d = 0.;
			track.dzdl = houghbin->get_dzdl_center(zoomlevel);
			track.z0 = houghbin->get_z0_center(zoomlevel);
			++icluster;
			}
			}
//			_tracks.push_back(track);// test
			new_tracks.push_back(track);
		}
		break;
	}
	cout<< "Number of tracks added : (to be used for initial vertexing or track seeding)"<<new_tracks.size()<<endl;

}

int PHG4InitZVertexing::build_triplets_to_Track3D(std::vector<Track3D>& new_tracks, bool forward){
 
	unsigned int layer0 = 0;
	unsigned int layer1 = 1;
	unsigned int layer2 = 2;
	if (!forward) {
		layer0 = _nlayers -1;
		layer1 = _nlayers -2;
		layer2 = _nlayers -3;
	}

	std::vector<unsigned int> layer_sorted_0;
	std::vector<unsigned int> layer_sorted_1;
	std::vector<unsigned int> layer_sorted_2;

        for (std::map<unsigned int,Cluster3D>::iterator it=hits_map.begin();
                        it!=hits_map.end();
                        ++it)
        {
                cluster_id = it->first;
                if (cluster_id != cluster_id) continue;
                bool used = false;
                auto hitused = hits_used.find(cluster_id);
                if (hitused != hits_used.end())
                used = hits_used.find(cluster_id)->second;
                if(used) continue;

                Cluster3D hit = it->second;
                unsigned int layer = hit.get_layer();
		unsigned int hitid =  hit.get_id();
		if (hitid != hitid) continue;
		//cout<<"layer "<< layer<<" hitid "<<hitid<<endl;
		if (layer == layer0 ) layer_sorted_0.push_back(hitid);                
		else if (layer == layer1) layer_sorted_1.push_back(hitid);
		else if (layer == layer2) layer_sorted_2.push_back(hitid);
		else continue;
        }

	cout<<"layer 0 " <<layer_sorted_0.size()<<endl;
	cout<<"layer 1 " <<layer_sorted_1.size()<<endl;
	cout<<"layer 2 " <<layer_sorted_2.size()<<endl;
		

	std::set<unsigned int> clusters;
	bool fill_track = false;
	unsigned int cluster_layer0 = 0;
	for (std::vector<unsigned int>::iterator it0 = layer_sorted_0.begin() ; it0 != layer_sorted_0.end(); ++it0)
	{
#ifdef _DEBUG_
		cout<<"icluster layer 0 "<<cluster_layer0<<endl;
#endif
		++cluster_layer0;
		clusters.clear();
		fill_track = false;
//		cout<<"cluster_id for hit 0 "<< *it0<<endl;
 		auto search0 = hits_map.find(*it0);
		if (search0 == hits_map.end()) continue;
 		Cluster3D hit0 = hits_map.find(*it0)->second;
		float phi_h0 = atan2(hit0.get_y(),hit0.get_x());
		phi_h0 = shift_phi_range(phi_h0);
		float z_h0 = hit0.get_z();

                // one track candidate for each hit on the first layer
                Track3D track;
		clusters.insert(*it0);

		for (std::vector<unsigned int>::iterator it1 = layer_sorted_1.begin(); it1 != layer_sorted_1.end(); ++it1 )
		{
			Cluster3D hit1 = hits_map.find(*it1)->second;
	                float phi_h1 = atan2(hit1.get_y(),hit1.get_x());
			phi_h1 = shift_phi_range(phi_h1);
			float z_h1 = hit1.get_z();

                        float phi_h0h1 = phi_h1-phi_h0;
                        if (phi_h1< M_PI/2. && phi_h0 > 3*M_PI/2.) phi_h0h1 += 2.*M_PI;
                        else if (phi_h1> 3*M_PI/2. && phi_h0 <M_PI/2.) phi_h0h1 -= 2.*M_PI;

			if (fabs(phi_h0h1)> _ca_phi_cut || fabs(z_h1-z_h0) > _z_cut) continue;
			

			for(std::vector<unsigned int>::iterator it2 = layer_sorted_2.begin(); it2 != layer_sorted_2.end(); ++it2)
			{
				Cluster3D hit2 = hits_map.find(*it2)->second;
				float phi_h2 = atan2(hit2.get_y(),hit2.get_x());
				phi_h2 = shift_phi_range(phi_h2);
				float phi_h1h2 = phi_h2-phi_h1;
				if (phi_h2< M_PI/2. && phi_h1 > 3*M_PI/2.) phi_h1h2 += 2.*M_PI;
				else if (phi_h2> 3*M_PI/2. && phi_h1 <M_PI/2.) phi_h1h2 -= 2.*M_PI;
				 

				float z_h2 = hit2.get_z();
				if (fabs(phi_h1h2) < _ca_phi_cut && (z_h1-z_h0)*(z_h2-z_h1)>0 && fabs(z_h2-z_h1)<_z_cut ){
				
					clusters.insert(*it1);
					clusters.insert(*it2);
					fill_track = true;
				} 
			} // layer 2
		} // layer 1

		if (fill_track)	{
//			cout<<" fill_track "<<endl;
			unsigned int nclusters =0 ;
			track.hits.assign(clusters.size(),Cluster3D());
			for (std::set<unsigned int>::iterator it=clusters.begin(); it!=clusters.end(); ++it ){
			//cout<<"cluster id  "<<*it<<endl;
			Cluster3D hit = hits_map.find(*it)->second;
			track.hits[nclusters] = hit;
			track.hits[nclusters].set_id(*it);
			++nclusters;
			}
			new_tracks.push_back(track);
		}
	}// layer 0

	
	cout<<"number of triplets "<<new_tracks.size()<<endl;	
	return 1;	
}


int PHG4InitZVertexing::cellular_automaton_zvtx_init(std::vector<Track3D>& candidate_tracks){

	cout<<"Entering cellular autumaton : processing "<< candidate_tracks.size()<<" tracks. "<<endl;
	ca =	new CellularAutomaton_v1(candidate_tracks,_radii,_material);
	ca->set_hough_space(_hough_space);
	ca->set_mag_field(_mag_field);
	ca->set_pt_rescale(_pt_rescale);
	ca->set_remove_hits(true);
//	ca->set_propagate_forward(false);// need to implement triplet in forward propagation
	ca->set_propagate_forward(true);

//	ca->set_mode(0);
        ca->set_triplet_mode(true); // triplet
	ca->set_seeding_mode(false);
        ca->set_hits_map(hits_map);

	ca->set_remove_inner_hits(true);
	ca->set_n_layers(_ca_nlayers);
	ca->set_required_layers(_ca_nlayers);
	ca->set_ca_chi2_layer(2.0);// 1.0
	ca->set_ca_chi2(_ca_chi2);
	ca->set_ca_phi_cut(_ca_phi_cut);
	ca->set_ca_z_cut(_ca_z_cut);
	ca->set_ca_dcaxy_cut(_dcaxy_cut);
	int code = ca->run(_tracks, _track_states, hits_used); // push_back output tracks into _tracks
	ca->Reset();
	for(unsigned int i=0; i<candidate_tracks.size(); ++i) candidate_tracks[i].reset();
	candidate_tracks.clear();
	if (!code)
	return 0;
	
	return Fun4AllReturnCodes::EVENT_OK;

}

/*
int PHG4InitZVertexing::cellular_automaton_zvtx_second(std::vector<Track3D>& candidate_tracks){


        cout<<"Entering cellular autumaton zvtx_second : processing "<< candidate_tracks.size()<<" tracks. "<<endl;
        ca =    new CellularAutomaton_v1(candidate_tracks,_radii,_material);
        ca->set_hough_space(_hough_space);
        ca->set_mag_field(_mag_field);
        ca->set_pt_rescale(_pt_rescale);
      ca->set_n_layers(_ca_nlayers);
      ca->set_required_layers(_ca_min_nlayers);
        ca->set_ca_chi2(_ca_chi2);
	ca->set_ca_chi2_layer(2.);
        ca->set_propagate_forward(true);
	ca->set_mode(0);
	ca->set_remove_hits(true);
	ca->set_remove_inner_hits(true);
	ca->set_require_inner_hits(false);//>1
        int code = ca->run(_tracks,_track_states, hits_used); // push_back output tracks into _tracks
	ca->Reset();
	for(unsigned int i=0; i<candidate_tracks.size(); ++i) candidate_tracks[i].reset();
	candidate_tracks.clear();
        if (!code) return 0;
        return Fun4AllReturnCodes::EVENT_OK;

}

int PHG4InitZVertexing::cellular_automaton_zvtx_third(std::vector<Track3D>& candidate_tracks){

        cout<<"Entering cellular autumaton zvtx_third : processing "<< candidate_tracks.size()<<" tracks. "<<endl;
        ca =    new CellularAutomaton_v1(candidate_tracks,_radii,_material);
        ca->set_hough_space(_hough_space);
        ca->set_mag_field(_mag_field);
        ca->set_pt_rescale(_pt_rescale);
        ca->set_n_layers(_ca_nlayers);
        ca->set_required_layers(_ca_min_nlayers);
        ca->set_ca_chi2(_ca_chi2);
        ca->set_ca_chi2_layer(2.);
        ca->set_propagate_forward(true);
        ca->set_remove_hits(true);
        ca->set_remove_inner_hits(true);
        int code = ca->run(_tracks,_track_states, hits_used); // push_back output tracks into _tracks
        ca->Reset();
	for(unsigned int i=0; i<candidate_tracks.size(); ++i) candidate_tracks[i].reset();
        candidate_tracks.clear();
        if (!code) return 0;
        return Fun4AllReturnCodes::EVENT_OK;

}

*/

int PHG4InitZVertexing::fit_vertex(){
	cout<<"tracks size " << _tracks.size()<<endl;

	if (_tracks.empty()) return -1;
	std::vector<Track3D> vtx_tracks;
	vtx_tracks.clear();
	_tracks.swap(vtx_tracks);

	unsigned int nzvtx = vtx_tracks.size();

        // range (-32, 32)
        unsigned int nzbins= 2*1280;
	float binsize = (_max_z0-_min_z0)/nzbins;
        float zvalues[2*1280];
        unsigned int zcounts[2*1280];
        for (unsigned int i = 0; i<nzbins; ++i)
        {
        zvalues[i] = _min_z0 + (i+0.5)*binsize;
	zcounts[i] = 0;
        }

	for (unsigned int i = 0; i < nzvtx; ++i) 
	{
        double zvtx = vtx_tracks[i].z0;
        if (zvtx != zvtx || zvtx<_min_z0 || zvtx>_max_z0) continue;
	unsigned int zbin = (zvtx-_min_z0)/binsize;	
        cout<<"zvtx " <<zvtx <<" zbin "<< zbin  <<endl;
	++zcounts[zbin];
	}
	
	std::vector<float> zvertices; // peak bin position
	std::vector<float> zmeans; // averaged z for tracks within 500um window
	std::vector<float> zsums;
	std::vector<float> zsigmas;
	std::vector<unsigned int> nzvtxes;

	for (unsigned int j=2; j<(nzbins-3); ++j)
	{
		cout<<"z countz "<<j<<" "<<zcounts[j]<<endl;
//		bool threebinspeak = _mult_threebins*(zcounts[j-2]) < (zcounts[j-1]+zcounts[j]+zcounts[j+1]) 
//				&& _mult_threebins*(zcounts[j+2]) < (zcounts[j-1]+zcounts[j]+zcounts[j+1]);
		bool twobinspeak =  _mult_twobins*(zcounts[j-1]+zcounts[j-2])<(zcounts[j]+zcounts[j+1]) 
				&& _mult_twobins*(zcounts[j+2]+zcounts[j+3])< (zcounts[j]+zcounts[j+1]);
		bool onebinpeak = _mult_onebin*(zcounts[j-1])<zcounts[j] && _mult_onebin*(zcounts[j+1]< zcounts[j]);
		
		
		if ((zcounts[j]>=_min_zvtx_tracks && onebinpeak) || ( (zcounts[j]+zcounts[j+1])>= _min_zvtx_tracks && twobinspeak) 
			/*|| ((zcounts[j-1]+zcounts[j]+zcounts[j+1]) > _min_zvtx_tracks && threebinspeak)*/ ){
		  float zvertex=-999.;
//		  if (threebinspeak) zvertex = (zvalues[j-1]*zcounts[j-1] + zvalues[j]*zcounts[j] + zvalues[j+1]*zcounts[j+1])
//						/(zcounts[j-1]+zcounts[j]+zcounts[j+1]);

		  if (onebinpeak) zvertex = zvalues[j];
		  if (twobinspeak) zvertex = (zvalues[j]*zcounts[j] + zvalues[j+1]*zcounts[j+1])/(zcounts[j]+zcounts[j+1]);  
		zvertices.push_back(zvertex);
		nzvtxes.push_back(0);
		zmeans.push_back(0.0);
		zsums.push_back(0.0);
		zsigmas.push_back(0.0);
		cout<<"vertex "<< zvertex<<endl;	
		}
	}

	cout<<"number of candidate z-vertices : "<< zvertices.size()<<endl;
	for (unsigned int j = 0; j <zvertices.size(); ++j)
	cout<<"primary vertex position : "<<zvertices[j]<<endl;

	unsigned int izvtx= 999;
    	for (unsigned int i = 0; i < nzvtx; ++i) {
	double zvtx = vtx_tracks[i].z0;
	if (zvtx != zvtx) continue;
		bool pass_dcaz= false;
		for (unsigned int k = 0; k<zvertices.size(); ++k) {
		bool new_vtx = fabs(zvtx-zvertices[k]) < _dcaz_cut;
		if (new_vtx) izvtx=k;
		pass_dcaz = pass_dcaz || new_vtx;
		}

	if (!pass_dcaz || izvtx>99) continue;
	++nzvtxes[izvtx];
	zsums[izvtx] += zvtx;
	}

	for (unsigned int iz = 0; iz<zvertices.size();++iz ){
	if (nzvtxes[iz]==0) continue;
	zmeans[iz] = zsums[iz]/nzvtxes[iz];
	cout<<"zmean for vertex "<< iz <<" = "<<zmeans[iz]<<endl;
	zvertices[iz] = zmeans[iz];
	zsums[iz]=0.;
        	if (fill_multi_zvtx){
        	float ntp_data[2];
        	ntp_data[0] = _event;
        	ntp_data[1] = zmeans[iz];
		_ntp_zvtx_by_event->Fill(ntp_data);
		}
	}

	for (unsigned int j = 0; j < nzvtx; ++j){
		double zvtx = vtx_tracks[j].z0;
		if (zvtx != zvtx) continue;

#ifdef _MULTIVTX_
		if (fill_multi_zvtx){
		float ntp_data[2];
		ntp_data[0] = _event;
		ntp_data[1] = vtx_tracks[j].z0;
//		ntp_data[1] = vtx_tracks[j].d;
//		cout<<"event "<<ntp_data[0]<<" "<< "vtx "<<ntp_data[1]<<endl;
        	_ntp_zvtx_by_track->Fill(ntp_data);
		}
#endif
                bool pass_dcaz= false;
                for (unsigned int k = 0; k<zvertices.size(); ++k) {
                bool new_vtx = fabs(zvtx-zvertices[k]) < _dcaz_cut;
                if (new_vtx) izvtx=k;
                pass_dcaz = pass_dcaz || new_vtx;
                }

		if (!pass_dcaz || izvtx>99) continue;
		zsums[izvtx] += pow(fabs(zmeans[izvtx] - vtx_tracks[j].z0),2);
	}

	std::vector<float> _multi_vtx;
	std::vector<std::vector<Track3D> > _multi_vtx_tracks;
        std::vector<std::vector<double> > _multi_vtx_track_errors;     ///< working array of track chisq
        std::vector<std::vector<Eigen::Matrix<float, 5, 5> > > _multi_vtx_track_covars; ///< working array of track covariances


	for (unsigned int i = 0; i < zvertices.size(); ++i)
	{
		if (nzvtxes[i]==0) continue;
		_multi_vtx.push_back(zvertices[i]);
		std::vector<Track3D> one_vtx_tracks;
		_multi_vtx_tracks.push_back(one_vtx_tracks);
		std::vector<double> one_track_errors;
		_multi_vtx_track_errors.push_back(one_track_errors);
		std::vector<Eigen::Matrix<float, 5, 5> > one_track_covars;
		_multi_vtx_track_covars.push_back(one_track_covars);
		zsigmas[i] = sqrt(zsums[i]/(nzvtxes[i]-1));
		cout<<"zsigma for vertex "<<i <<" = "<<zsigmas[i]<<endl;
	}

	unsigned int nzvtx_final = 0;
	for (unsigned int k = 0; k < nzvtx; ++k){

		float zvtx = vtx_tracks[k].z0;
                bool pass_dcaz= false;
                for (unsigned int iz = 0; iz<_multi_vtx.size(); ++iz) {
                bool new_vtx = fabs(zvtx-_multi_vtx[iz]) < _dcaz_cut;
                if (new_vtx) izvtx=iz;
                pass_dcaz = pass_dcaz || new_vtx;
                }

	        if (!pass_dcaz || izvtx >= _multi_vtx.size()) continue;
//		cout<<"adding a track to vtx "<<izvtx<<endl;
//		if (fabs(vtx_tracks[k].z0-zmeans[k])> _dcaz_cut/*10*zsigmz*/)continue;
		++nzvtx_final;

		_multi_vtx_tracks[izvtx].push_back(vtx_tracks[k]);
		_multi_vtx_track_covars[izvtx].push_back(_track_states[k].C);
		_multi_vtx_track_errors[izvtx].push_back(_track_states[k].chi2);

    	}

	cout<<"start fitting vertex.. "<<endl;
	for (unsigned int i = 0; i<_multi_vtx.size(); ++i)
	{
		if (_multi_vtx_tracks[i].size()==0) continue;
		_vertex[2] = _multi_vtx[i];
    		if (Verbosity() > 0) {
      		cout << " seed track vertex pre-fit: "
           		<< _vertex[0] << " "
        		<< _vertex[1] << " "
        		<< _vertex[2] << endl;
    		}

		_vertex_finder.findVertex( _multi_vtx_tracks[i], _multi_vtx_track_covars[i], _vertex, 0.10, true);
		_vertex_finder.findVertex( _multi_vtx_tracks[i], _multi_vtx_track_covars[i], _vertex, 0.02, true);
		_vertex_finder.findVertex( _multi_vtx_tracks[i], _multi_vtx_track_covars[i], _vertex, 0.005, true);

		n_vtx_tracks = _multi_vtx_tracks[i].size();
		cout<<"number of fitted tracks for vertex "<<i<< " : "<<n_vtx_tracks <<endl;

  		if (Verbosity() > 0) {
    			cout << " seed track vertex post-fit: "
         		<< _vertex[0] << " " << _vertex[1] << " " << _vertex[2] << endl;
  		}
	
		// add vertex_id to Track3D
//		std::vector<float> pre_vtx(3,0.0);

		_vertex_list.push_back(_vertex);

	}// loop over vertices

	cout<<"final vertices : "<<endl;
	for (unsigned int i = 0; i<_vertex_list.size(); ++i){
		cout<<i<<" "<<_vertex_list[i][2]<<endl;
	}


        for (unsigned int k = 0; k < nzvtx; ++k){

                float zvtx = vtx_tracks[k].z0;
                bool pass_dcaz= false;
                for (unsigned int iz = 0; iz < _vertex_list.size(); ++iz) {
                bool new_vtx = fabs(zvtx - _vertex_list[iz][2]) < _dcaz_cut;
                if (new_vtx) izvtx=iz;
                pass_dcaz = pass_dcaz || new_vtx;
                }


		if (!pass_dcaz || izvtx >= _vertex_list.size()) continue;
//              if (pass_dcaz && izvtx <9999) 
                ++nzvtx_final;

                Track3D vtx_track = vtx_tracks[k];
		vtx_track.set_vertex_id(izvtx);
//		cout<<"vertex_id "<< izvtx<<endl;
                _tracks.push_back(vtx_track);
                _track_covars.push_back(_track_states[k].C);
                _track_errors.push_back(_track_states[k].chi2 );
		
		izvtx=0;
        }

	for(unsigned int i=0; i<vtx_tracks.size(); ++i) vtx_tracks[i].reset();
	vtx_tracks.clear();
        _multi_vtx_tracks.clear();
        _multi_vtx_track_covars.clear();
        _multi_vtx_track_errors.clear();

 	return Fun4AllReturnCodes::EVENT_OK;

} 

void PHG4InitZVertexing::convertHelixCovarianceToEuclideanCovariance(float B,
		float phi, float d, float kappa, float z0, float dzdl,
		Eigen::Matrix<float, 5, 5> const& input,
		Eigen::Matrix<float, 6, 6>& output) {

	Eigen::Matrix<float, 6, 5> J = Eigen::Matrix<float, 6, 5>::Zero(6, 5);

	// phi,d,nu,z0,dzdl
	// -->
	// x,y,z,px,py,pz

	float nu = sqrt(kappa);
	float dk_dnu = 2 * nu;

	float cosphi = cos(phi);
	float sinphi = sin(phi);

	J(0, 0) = -sinphi * d;
	J(0, 1) = cosphi;
	J(1, 0) = d * cosphi;
	J(1, 1) = sinphi;
	J(2, 3) = 1;

	float pt = 0.003 * B / kappa;
	float dpt_dk = -0.003 * B / (kappa * kappa);

	J(3, 0) = -cosphi * pt;
	J(3, 2) = -sinphi * dpt_dk * dk_dnu;
	J(4, 0) = -sinphi * pt;
	J(4, 2) = cosphi * dpt_dk * dk_dnu;

	float alpha = 1. / (1. - dzdl * dzdl);
	float alpha_half = pow(alpha, 0.5);
	float alpha_3_half = alpha * alpha_half;

	J(5, 2) = dpt_dk * dzdl * alpha_half * dk_dnu;
	J(5, 4) = pt * (alpha_half + dzdl * dzdl * alpha_3_half) * dk_dnu;

	output = J * input * (J.transpose());
}

void PHG4InitZVertexing::shift_coordinate_system(double dx, double dy, double dz) {

        for (std::map<unsigned int,Cluster3D>::iterator it=hits_map.begin();
                        it!=hits_map.end();
                        ++it)
        {
		Cluster3D hit = it->second;
                hit.set_x(hit.get_x() + dx);
		hit.set_y(hit.get_y() + dy);
		hit.set_z(hit.get_z() + dz);
		it->second = hit;

	}

        for (unsigned int itrack = 0; itrack < _tracks.size(); ++itrack) {
                for (unsigned int ihit = 0; ihit < _tracks[itrack].hits.size(); ++ihit) {
                        _tracks[itrack].hits[ihit].set_x(_tracks[itrack].hits[ihit].get_x() + dx);
                        _tracks[itrack].hits[ihit].set_y(_tracks[itrack].hits[ihit].get_y() + dy);
                        _tracks[itrack].hits[ihit].set_z(_tracks[itrack].hits[ihit].get_z() + dz);
                }
        }

        _vertex[0] += dx;
        _vertex[1] += dy;
        _vertex[2] += dz;

        return;
}

float PHG4InitZVertexing::shift_phi_range(float _phi){

	if (_phi < 0.) _phi += 2.*M_PI;
	return _phi;
}

