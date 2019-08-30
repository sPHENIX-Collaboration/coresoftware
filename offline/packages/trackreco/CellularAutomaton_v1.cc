#include "CellularAutomaton_v1.h"

#include "HelixHoughSpace.h"    // for HelixHoughSpace
#include "HelixKalmanFilter.h"  // for HelixKalmanFilter


#include <phool/phool.h>        // for PHWHERE

#include <Eigen/Core>

#include <algorithm>            // for sort
#include <cfloat>
#include <cmath>
#include <cstdlib>             // for exit, NULL
#include <iostream>
#include <memory>               // for allocator_traits<>::value_type
#include <sys/time.h>
#include <utility>              // for swap, pair, make_pair

using namespace std;
using namespace Eigen;

//#define _DEBUG_
// if _FULL_TEST_ is not chosen, only triplets are translated into SimpleTrack3D
#define _FULL_TEST_


CellularAutomaton_v1::CellularAutomaton_v1(std::vector<SimpleTrack3D>& input_tracks, std::vector<float>& detector_radii, std::vector<float>& detector_materials)
	:
	_hough_space(nullptr),
	_kalman(nullptr),
	in_tracks(std::vector<SimpleTrack3D>()),
	ca_tracks(std::vector<SimpleTrack3D>()),
	ca_track_states(std::vector<HelixKalmanState>()),
	temp_combo(std::vector<unsigned int>()),
	combos(std::set<std::vector<unsigned int> >()),
	layer_sorted(std::vector<std::vector<SimpleHit3D> >()),
	nlayers(10),
	rlayers(8),
	allowed_missing_inner_hits(0),
	ca_cos_ang_cut(0.985),
	ca_chi2_cut(2.0),
	ca_chi2_layer_cut(2.0),
	ca_phi_cut(M_PI/144.),
	ca_z_cut(20),
	ca_dcaxy_cut(0.02),
	fast_chi2_cut_max(FLT_MAX),
	fast_chi2_cut_par0(12.),
	fast_chi2_cut_par1(0.),
	_pt_rescale(1.0),
	_mag_field(1.4),
	_detector_radii(std::vector<float>()),
	_detector_scatter(std::vector<float>()),
	_detector_materials(std::vector<float>()),
	_integrated_scatter(std::vector<float>()),
	_hits_used(std::map<unsigned int, bool>()),
	_hits_map(std::map<unsigned int,SimpleHit3D>()),
	CAtime(0.),
	KALtime(0.),
	forward(false),
	remove_hits(false),
        remove_inner_hits(false),
        require_inner_hits(false),
	triplet_mode(true),
        seeding_mode(false),
        verbose(0)
{
	set_input_tracks(input_tracks);
	set_detector_radii(detector_radii);
	set_detector_materials(detector_materials);

}

void CellularAutomaton_v1::Reset(){

        delete _hough_space;
        delete _kalman;
        in_tracks.clear();
        ca_tracks.clear();
        ca_track_states.clear();
        
        temp_combo.clear();
	
        _detector_radii.clear();
        _detector_scatter.clear();
        _detector_materials.clear();
        _integrated_scatter.clear();

}


int CellularAutomaton_v1::run(std::vector<SimpleTrack3D>& output_tracks, std::vector<HelixKalmanState>& output_track_states, std::map<unsigned int, bool>& hits_used)
{
	if (remove_hits){
//	cout<< "hits_used size "<< hits_used.size()<<endl;
//	cout<< "_hits_used size "<<_hits_used.size()<<endl; 
	_hits_used.swap(hits_used);
	}

	int code = 0;
	if(verbose > 0) cout<< "CellularAutomaton:: initializing..."<<endl;
	code = init();
	if(verbose > 0) cout<<code<<endl;
	if (!code) 
        {
          cout << PHWHERE << "::Error - Initialization failed. " << endl;
          exit(1);
        }
	code = 0;
	if(verbose > 0) cout<<"CellularAutomaton:: processing tracks... "<<endl;
	code = process_tracks();
	if (!code)		
        {
          cout << PHWHERE << "::Error - Processing tracks failed. " << endl;
          exit(1);
        }
	code = 0;
	if(verbose > 0) cout<<"CellularAutomaton:: outputting ca tracks..."<<endl;
	code = get_ca_tracks(output_tracks, output_track_states);
	if (!code)
        {
          cout << PHWHERE << "::Error - Outputting tracks failed. " << endl;
          exit(1);
        }

	if(remove_hits){
	_hits_used.swap(hits_used);		
	if(verbose > 0)  cout<< "hits_used size "<< hits_used.size()<<endl;
	if(verbose > 0) cout<< "_hits_used size "<<_hits_used.size()<<endl;
	}

	for (unsigned int i = 0; i<in_tracks.size(); ++i) in_tracks[i].reset();
	in_tracks.clear();
	for (unsigned int i = 0; i<ca_tracks.size(); ++i) ca_tracks[i].reset();
	ca_tracks.clear();
	ca_track_states.clear();
	layer_sorted.clear();
	return 1;
}

int CellularAutomaton_v1::init()
{
        if (!_hough_space)
        {
          cout << PHWHERE << "::Error - Hough Space is not set. " << endl;
          exit(1);
        }
/*
        if (!in_tracks)
        {
          cout << PHWHERE << "::Error - Input tracks are not set" << endl;
          exit(1);
        }
*/
        if (!_detector_radii.size())
        {
          cout << PHWHERE << "::Error - Detector radii are not set" << endl;
          exit(1);
        }

        if (!_detector_materials.size())
        {
          cout << PHWHERE << "::Error - Detector materials are not set" << endl;
          exit(1);
        }

	temp_combo.clear();
	temp_combo.assign(nlayers,0);

	combos.clear();

        std::vector<SimpleHit3D> one_layer;
        layer_sorted.clear();
        layer_sorted.assign(nlayers, one_layer);

	ca_tracks.clear();
	ca_track_states.clear();
	set_cylinder_kalman();
	
	return 1;
}


void CellularAutomaton_v1::set_hough_space(HelixHoughSpace* hough_space) {

  _hough_space = dynamic_cast<HelixHoughSpace*> (hough_space->CloneMe());
  assert(_hough_space);
}

void CellularAutomaton_v1::set_mag_field(float mag_field) {
	_mag_field = mag_field;
}

void CellularAutomaton_v1::set_pt_rescale(float pt_rescale){
	_pt_rescale = pt_rescale;
}

void CellularAutomaton_v1::set_detector_radii(std::vector<float>& radii)
{
  for (unsigned int i = 0; i < radii.size(); ++i) {
    _detector_radii.push_back(radii[i]);
  }
}

void CellularAutomaton_v1::set_detector_materials(std::vector<float>& materials)
{
  for (unsigned int i = 0; i < materials.size(); ++i) {
    _detector_scatter.push_back(1.41421356237309515 * 0.0136 *
                               sqrt(3. * materials[i]));
    _detector_materials.push_back(3. * materials[i]);
  }

  _integrated_scatter.assign(_detector_scatter.size(), 0.);
  float total_scatter_2 = 0.;
  for (unsigned int l = 0; l < _detector_scatter.size(); ++l) {
    total_scatter_2 += _detector_scatter[l] * _detector_scatter[l];
    _integrated_scatter[l] = sqrt(total_scatter_2);
  }
}

void CellularAutomaton_v1::set_cylinder_kalman(){
  	_kalman =
      	new HelixKalmanFilter(_detector_radii, _detector_materials, _mag_field);
}

void CellularAutomaton_v1::set_input_tracks(std::vector<SimpleTrack3D>& input_tracks) 
{
/*
	for (unsigned int i = 0; i < input_tracks.size(); ++i)
	{
	in_tracks.push_back(input_tracks[i]);
	}

*/
	in_tracks = input_tracks;
	//cout<<"Setting input tracks : size = " << in_tracks.size()<<endl; 

}

int CellularAutomaton_v1::get_ca_tracks(std::vector<SimpleTrack3D>& output_tracks, std::vector<HelixKalmanState>& output_track_states)
{
	// push back new ca processed tracks into _tracks
	   
	if (ca_tracks.size() != ca_track_states.size()) 
	return 0;

	for (unsigned int i = 0; i <ca_tracks.size(); ++i)
	{
	output_tracks.push_back(ca_tracks[i]);
	output_track_states.push_back(ca_track_states[i]);
	}

	//cout<<"newly added ca tracks : "<< ca_tracks.size() <<" tracks. "<<endl;
	return 1;
}


int CellularAutomaton_v1::process_tracks()
{

	for (unsigned int i = 0; i < in_tracks.size(); ++i) 
	{ // loop over input tracks
	  if(verbose > 0) cout<<"track candidate "<<i<<endl;

	  if (triplet_mode)
	  {
		process_single_triplet(in_tracks[i]);
/*
		SimpleTrack3D itrack = in_tracks[i]; // hit triplets with track parameters from HT  
		int code = 0;
		unsigned int missing_layers = 0;
		bool last_layer = false;
		code = triplet_to_segment(itrack); // from hit triplets to a segment with estimated kappa & dzdl 
 		for (unsigned int n=3; n < nlayers; ++n )  // extend segments{	
			if (n== nlayers) last_layer = true;
			code = process_single_track(itrack,n, last_layer); // search layer n
			if (code==0) ++missing_layers;
			if (missing_layers > (nlayers-rlayers)) break;
		}
*/
	  }
	  else
	  {
		process_single_track(in_tracks[i]);
	  }
	}

	return 1;

}


int CellularAutomaton_v1::process_single_triplet(SimpleTrack3D& track){ // track : from hough transform

        std::vector<TrackSegment> segments1;
        std::vector<TrackSegment> segments2;
	
        std::vector<TrackSegment>* cur_seg = &segments1;
        std::vector<TrackSegment>* next_seg = &segments2;
        unsigned int cur_seg_size = 0;
        unsigned int next_seg_size = 0;

	std::vector<TrackSegment> complete_segments;

	// from hit triplets to a segment with estimated kappa & dzdl
//      triplet_to_segments(track, cur_segments); // from hit triplets to a segment with estimated kappa & dzdl 
	for (unsigned int l = 0; l < 3; ++l) {
		layer_sorted[l].clear();
	}

	if(verbose > 0) cout<<"track.hits.size "<<track.hits.size()<<endl;
  	for (unsigned int i = 0; i < track.hits.size(); ++i) {
    		SimpleHit3D hit = track.hits[i];
    		unsigned int layer = (unsigned int) hit.get_layer();
		if(verbose > 0) cout<<"layer "<<layer<< endl;
       	 	if (!forward) layer = nlayers-layer-1;
                if (layer > (nlayers-1)) continue;
    		unsigned int min = (layer - allowed_missing_inner_hits);
    		if (allowed_missing_inner_hits > layer) {
      		min = 0;
    		}
    		for (unsigned int l = min; l <= layer; ++l) {
      		layer_sorted[l].push_back(hit);
		if(verbose > 0) cout<<"adding hit in layer "<<l<<endl;
		}
	}

//        for (unsigned int l = 0; l < nlayers; ++l) {
	for (unsigned int l = 0; l< 3; ++l){
	  if(verbose > 0) cout<<"layer_sorted["<<l<<"].size = "<< layer_sorted[l].size()<<endl;
        	if (layer_sorted[l].size() == 0) {
        	return 0;
    		}
  	}

#ifdef _FULL_TEST_
	float ca_cos_ang_cut_diff = 1. - ca_cos_ang_cut;
	float ca_cos_ang_cut_diff_inv = 1. / ca_cos_ang_cut_diff;
#endif
	float ca_sin_ang_cut = sqrt(1. - ca_cos_ang_cut * ca_cos_ang_cut);

	std::vector<float> inv_layer;
	inv_layer.assign(nlayers, 1.);
	for (unsigned int l = 3; l < nlayers; ++l) {
	inv_layer[l] = 1. / (((float)l) - 2.);
	}

	// l = 3, 4,   5,   6,   7
	//     1, 1/2, 1/3, 1/4, 1/5
	//     1, 1/3, 1/5, 1/7, 1/9         

	float x1, x2, x3;
	float y1, y2, y3;
	float z1, z2, z3;
	float dx1, dx2, dx3;
	float dy1, dy2, dy3;
	float dz1, dz2, dz3;

 	float kappa;
	float dkappa;

	float ux_mid;
 	float uy_mid;
	float ux_end;
	float uy_end;

	float dzdl_1;
	float dzdl_2;
	float ddzdl_1;
	float ddzdl_2;

#ifdef _FULL_TEST_
	float cur_kappa;
	float cur_dkappa;
	float cur_ux;
	float cur_uy;
	float cur_chi2;
	float chi2;
#endif
	unsigned int hit1;
	unsigned int hit2;
	unsigned int hit3;

	TrackSegment temp_segment;
	temp_segment.hits.assign(nlayers, 0);

  	for (unsigned int i = 0; i < layer_sorted[0].size(); ++i) {
		for (unsigned int j = 0; j < layer_sorted[1].size(); ++j) {
			for (unsigned int k = 0; k < layer_sorted[2].size() ; ++k) {

			unsigned int layer0 = layer_sorted[0][i].get_layer();
			unsigned int layer1 = layer_sorted[1][j].get_layer();
			unsigned int layer2 = layer_sorted[2][k].get_layer();
			if (!forward)
			{
			layer0 = nlayers-layer0-1;
			layer1 = nlayers-layer1-1;
			layer2 = nlayers-layer2-1;
			}

        		x1 = layer_sorted[0][i].get_x();
        		y1 = layer_sorted[0][i].get_y();
       		 	z1 = layer_sorted[0][i].get_z();

        		dx1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[0][i].get_size(0,0));
        		dy1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[0][i].get_size(1,1));
        		dz1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[0][i].get_size(2,2));

        		x2 = layer_sorted[1][j].get_x();
        		y2 = layer_sorted[1][j].get_y();
        		z2 = layer_sorted[1][j].get_z();

        		dx2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[1][j].get_size(0,0));
        		dy2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[1][j].get_size(1,1));
        		dz2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[1][j].get_size(2,2));

        		x3 = layer_sorted[2][k].get_x();
        		y3 = layer_sorted[2][k].get_y();
        		z3 = layer_sorted[2][k].get_z();
        		dx3 = 0.5*sqrt(12.0)*sqrt(layer_sorted[2][k].get_size(0,0));
        		dy3 = 0.5*sqrt(12.0)*sqrt(layer_sorted[2][k].get_size(1,1));
        		dz3 = 0.5*sqrt(12.0)*sqrt(layer_sorted[2][k].get_size(2,2));

        		hit1 = i;
       	 		hit2 = j;
        		hit3 = k;

        		calculate_kappa_tangents(x1, y1, z1, x2, y2, z2, x3, y3,
                                z3, dx1, dy1, dz1, dx2, dy2, dz2,
                                dx3, dy3, dz3, kappa, dkappa,
                                ux_mid, uy_mid, ux_end, uy_end,
                                dzdl_1, dzdl_2, ddzdl_1, ddzdl_2);

#ifdef _DEBUG_
        cout<<"Triplet : "<<endl;
        cout<<"kappa "<<kappa<< " dkappa "<<dkappa<<" ux_mid "<<ux_mid<<" uy_mid "<<" ux_end "<<ux_end<<
        " uy_end " <<uy_end<<" dzdl_1 "<<dzdl_1<<" dzdl_2 "<<dzdl_2<<" ddzdl_1 "<<ddzdl_1<<" ddzdl_2 "<<ddzdl_2<<endl;
#endif

          		temp_segment.chi2 = pow(
              		(dzdl_1 - dzdl_2) /
              		(ddzdl_1 + ddzdl_2 + fabs(dzdl_1 * ca_sin_ang_cut)),2);

         		if (temp_segment.chi2 > ca_chi2_layer_cut) continue;
          		temp_segment.ux = ux_end;
          		temp_segment.uy = uy_end;
          		temp_segment.kappa = kappa;
//          if (temp_segment.kappa > _hough_space->get_kappa_max()) continue;
                      	temp_segment.dkappa = dkappa;

          		temp_segment.hits[0] = hit1;
          		temp_segment.hits[1] = hit2;
          		temp_segment.hits[2] = hit3;
          		temp_segment.n_hits = 3;


          		unsigned int outer_layer =
                	layer_sorted[2][temp_segment.hits[2]].get_layer();
          		if (!forward) outer_layer = nlayers-outer_layer-1;
          		// make sure we have required number of layers with hits
                       	if ((outer_layer - 2) > allowed_missing_inner_hits) continue;
          		// finish up if required number of layers is reached,
                       	if ((nlayers - 3) <= allowed_missing_inner_hits) {
                	complete_segments.push_back(temp_segment);
          		}
          		if (next_seg->size() == next_seg_size) { // first new segment
                		next_seg->push_back(temp_segment);
                		next_seg_size += 1;
          		} else { // next new segments
               			(*next_seg)[next_seg_size] = temp_segment;
                		next_seg_size += 1;
          		}

			}
		}
	}
	
  	swap(cur_seg, next_seg);
  	swap(cur_seg_size, next_seg_size);


  	//cout<<"number of complete segments : " << complete_segments.size()<<endl;
  	//cout<<"number of current segments : "<< cur_seg_size<<endl;

  	// copy complete segments over to current segments
       	for (unsigned int i = 0; i < complete_segments.size(); ++i) {
    		if (cur_seg->size() == cur_seg_size) {
      		cur_seg->push_back(complete_segments[i]);
      		++cur_seg_size;
    		} else {
     	 	(*cur_seg)[cur_seg_size] = complete_segments[i];
      		++cur_seg_size;
    		}
  	}

  	std::set<unsigned int> comp1;
  	std::set<unsigned int> comp2;

//  	cout<<"number of segments generated "<< cur_seg_size<<endl;
  	if (cur_seg_size==0  || cur_seg_size>10000) return 1;
  	for (unsigned int i = cur_seg_size-1; i>0; --i){
    		if ((*cur_seg)[i].n_hits==0) continue;
    		comp1.clear();
    		temp_combo.assign((*cur_seg)[i].n_hits, 0);
    		for (unsigned int l = 0; l < (*cur_seg)[i].n_hits; ++l) {
    		temp_combo[l] = layer_sorted[l][(*cur_seg)[i].hits[l]].get_id();
    		comp1.insert(temp_combo[l]);
    		}

    		sort(temp_combo.begin(), temp_combo.end());
    		set<vector<unsigned int> >::iterator it = combos.find(temp_combo);
    		if (it != combos.end()) {
     	 	(*cur_seg)[i].n_hits = 0.;
    		}
    		if (combos.size() > 10000) {
      		combos.clear();
    		}
    		combos.insert(temp_combo);

    		for (unsigned int j = i-1; j>=0; --j){
      			comp2.clear();
      			for (unsigned int m = 0; m < (*cur_seg)[j].n_hits; ++m){
      				comp2.insert(layer_sorted[m][(*cur_seg)[j].hits[m]].get_id());
      			};
      			for (std::set<unsigned int>::iterator it=comp1.begin(); it!=comp1.end(); ++it){
      			auto it2 = comp2.find(*it);
      			if (it2 != comp2.end()) comp2.erase(*it2);
      			}
      			if (comp2.empty()) {
        		(*cur_seg)[j].n_hits = 0;
      			}
      			if (j==0) break;
    		}
  	}

	comp1.clear();
	comp2.clear();

  	unsigned int nsegs = cur_seg_size;
  	for (unsigned int i = 0; i<cur_seg_size; i++) {
  		if ((*cur_seg)[i].n_hits ==0) --nsegs;
  	}
  	//cout<<"number of segments from triplets to be processed "<< nsegs<<endl;


#ifdef _FULL_TEST_

  	unsigned int which_seg;
	// seperate counting from triplets
//	unsigned int allowed_missing = nlayers - rlayers;
////	cout<<"allowed missing "<< allowed_missing<<endl;

	std::map<unsigned int, unsigned int> missing_layers_map; // segment_id, missing_layers
	std::map<unsigned int, unsigned int> missing_layers_map_next;
	missing_layers_map.clear();
	missing_layers_map_next.clear();
	for (unsigned int n = 0 ; n < cur_seg_size; ++n) missing_layers_map.insert(make_pair(n,0));
	unsigned int added_next_segments = 0;

 	// Extend segments, start from searching the 4th layer
 	for (unsigned int l=3; l < nlayers; ++l ) {
////		cout<<"nlayers "<<nlayers<<" l "<<l<<endl;
		next_seg_size = 0;
		next_seg->clear();
		// Loop over current segments
		for (unsigned int i = 0; i < cur_seg_size; ++i) {
			if ((*cur_seg)[i].n_hits ==0) continue;
			added_next_segments = 0;
			auto search = missing_layers_map.find(i);
                        unsigned int missing_layers = search->second;
////			cout<<"segment "<< i <<", missing layer "<<missing_layers<<endl;
//                next_segments(cur_segments, next_segments, n ); // search layer n
			// loop over all hits on layer n
			// if multiple legitimate hits in next layer, duplicate current segment, copy over missing layer and assign new segment_id
			//// ** drop layer sorted keep good hits in layer_sorted, hits[l] : id 
			// keep cluster ids in hits[l] for good hits of a segment
			// for a layer without a good hit, store with hits[l] = 999

			if ( (l-2) < 3){
		        x1 = layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_x();
		        y1 = layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_y();
        		z1 = layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_z();
		        dx1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_size(0,0));
                        dy1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_size(1,1));
                        dz1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_size(2,2));

			}else {
			// get positions from hits_map
			auto searchl2 = _hits_map.find((*cur_seg)[i].hits[l-2]);
			SimpleHit3D clusterl2 = searchl2->second;
			x1 = clusterl2.get_x();
			y1 = clusterl2.get_y();
			z1 = clusterl2.get_z();
			dx1 = 0.5*sqrt(12.0)*sqrt(clusterl2.get_size(0,0));
			dy1 = 0.5*sqrt(12.0)*sqrt(clusterl2.get_size(1,1));
			dz1 = 0.5*sqrt(12.0)*sqrt(clusterl2.get_size(2,2));
			}

			if  ( (l-1) < 3) {
        		x2 = layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_x();
        		y2 = layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_y();
        		z2 = layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_z();
                        dx2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_size(0,0));
                        dy2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_size(1,1));
                        dz2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_size(2,2));
			} else {
			// get positions from hits_map
			auto searchl1 = _hits_map.find((*cur_seg)[i].hits[l-1]);
			SimpleHit3D clusterl1 = searchl1->second;
			x2 = clusterl1.get_x();
			y2 = clusterl1.get_y();
			z2 = clusterl1.get_z(); 
                        dx2 = 0.5*sqrt(12.0)*sqrt(clusterl1.get_size(0,0));
                        dy2 = 0.5*sqrt(12.0)*sqrt(clusterl1.get_size(1,1));
                        dz2 = 0.5*sqrt(12.0)*sqrt(clusterl1.get_size(2,2));
			}
	

                        // convert segment to track to process it through kalman filter or just call fit_track

			SimpleTrack3D init_track;
			bool fit_layer = (l >= 6);
			if (fit_layer){
                                //init_track.hits.assign(l, SimpleHit3D());
                                init_track.hits.assign((*cur_seg)[i].n_hits, SimpleHit3D());

                		for (unsigned int ll = 0; ll < (*cur_seg)[i].n_hits; ++ll) {
                			if (ll<3){
                			init_track.hits[ll] = layer_sorted[ll][(*cur_seg)[i].hits[ll]];
                			}else {
                			auto search = _hits_map.find((*cur_seg)[i].hits[ll]);
                			SimpleHit3D cluster = search->second;
                			init_track.hits[ll] = cluster;
                			}
		                }
			}

                        	for (std::map<unsigned int,SimpleHit3D>::iterator jt = _hits_map.begin();
                                                jt!= _hits_map.end();
                                                ++jt) {

                        	SimpleHit3D hit3d = jt->second;
                        	unsigned int layer = hit3d.get_layer();
                        	if (layer != l) continue;
				which_seg = i;
                                hit1 = hit3d.get_id();
				
                                x3 = hit3d.get_x();
                                y3 = hit3d.get_y();
                                z3 = hit3d.get_z();

                                float phi_prev =shift_phi_range(atan2(y2,x2));
                                float phi_cur = shift_phi_range(atan2(y3,x3));
				float phi_diff = phi_cur-phi_prev;
				if (phi_cur< M_PI/2. && phi_prev > 3*M_PI/2.) phi_diff += 2.*M_PI;
				else if (phi_cur>3*M_PI/2 && phi_prev<M_PI/2.) phi_diff -= 2.*M_PI;
				if (!seeding_mode){
                                if ((fabs(phi_diff)> ca_phi_cut || abs(z3-z2)> ca_z_cut)) continue;
				} else {
				if ((fabs(phi_diff)> ca_phi_cut || abs(z3-z2)> ca_z_cut) && cur_seg_size!=1) continue;
				}

//                		auto search = _hits_used.find(hit1);
//                		if(search != _hits_used.end() && search->second ) continue;

				if (fit_layer)
				{
				// fit init_track to get kappa to compare with new kappa
				// copy init_track over to temp_track and add a hit in kalman filter
				SimpleTrack3D temp_track;
			        temp_track.hits.assign(init_track.hits.size()+1, SimpleHit3D());

      				for (unsigned int ll = 0; ll < init_track.hits.size(); ++ll) {
        			temp_track.hits[ll] = init_track.hits[ll];
      				}
				temp_track.hits[init_track.hits.size()] = hit3d;

				// track fitting instead of computing from triplets

 				float temp_chi2 = temp_track.fit_track();
//				cout<<"chi2 from fit_track "<<init_chi2 <<endl;			
				if (temp_chi2 != temp_chi2) continue;
		                if (temp_track.kappa != temp_track.kappa )  continue;
                		if (temp_track.z0 != temp_track.z0) continue;

		                HelixKalmanState state;
		                state.phi = temp_track.phi;
		                if (state.phi < 0.) {
		                state.phi += 2. * M_PI;
		                }
		                state.d = temp_track.d;
		                state.kappa = temp_track.kappa;
		                state.nu = sqrt(state.kappa);
		                state.z0 = temp_track.z0;
		                state.dzdl = temp_track.dzdl;
		                state.C = Matrix<float, 5, 5>::Zero(5, 5);
		                state.C(0, 0) = pow(0.01, 2.);
		                state.C(1, 1) = pow(0.01, 2.);
		                state.C(2, 2) = pow(0.01 * state.nu, 2.);
		                state.C(3, 3) = pow(0.05, 2.);
		                state.C(4, 4) = pow(0.05, 2.);
		                state.chi2 = 0.;
		                state.position = 0;
		                state.x_int = 0.;
		                state.y_int = 0.;
		                state.z_int = 0.;
			
				// place holder for kalman filter
/*
	                	unsigned int nfits = 0;
                		for (unsigned int h = 0; h < temp_track.hits.size(); ++h) {
                		_kalman->addHit(temp_track.hits[h], state);
                		nfits += 1;
                		cout<<"nfits "<<nfits<<endl;
                		}
                		cout<<"z0 after kalman "<<state.z0<<endl;

                		// fudge factor for non-gaussian hit sizes
                                state.C *= 3.;
                		state.chi2 *= 6.;
                		// kappa cut here drives both efficiency and ghost track rates down at the same time
                                if (!(temp_track.kappa == temp_track.kappa) ) {
                		continue;
                		}
                		if (!(state.chi2 == state.chi2)) {
                		continue;
                		}

*/
              			if (state.chi2 / (2. * ((float)(temp_track.hits.size())) - 5.) < ca_chi2_cut /* 10. */) {
					// translate temp_track into temp_segment (only hit info is saved) and save in next segments
					for (unsigned int ll = 0; ll < l; ++ll) {
                                        temp_segment.hits[ll] = (*cur_seg)[which_seg].hits[ll];
                                        }
                                        temp_segment.hits[l] = hit1;
					temp_segment.n_hits = l + 1;

                                        if (next_seg->size() == next_seg_size) { // first new segment
                                                next_seg->push_back(temp_segment);
                                                missing_layers_map_next.insert(make_pair(next_seg->size()-1, missing_layers));
                                                next_seg_size += 1;
                                        } else { // next new segments
                                        #ifdef _DEBUG_
                                        cout<<"Next segment size inconsistent "<<endl;
                                        #endif
                                                (*next_seg)[next_seg_size] = temp_segment;
                                                missing_layers_map_next.insert(make_pair(next_seg->size()-1, missing_layers));
                                                next_seg_size += 1;
                                        }
                                        ++added_next_segments;
#ifdef _DEBUG_
                                        cout<<"segment "<< which_seg << " added segment "<<added_next_segments<<endl;
#endif
                		
                		}// chi2 from track fitting and/or kalman filter : looser cuts than computing kappa from triplets
	
				}
				else
				{
	//        		x3 = layer_sorted[l][j].get_x();
	//        		y3 = layer_sorted[l][j].get_y();
	//        		z3 = layer_sorted[l][j].get_z();

	//        		dx1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_size(0,0));
	//        		dy1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_size(1,1));
	//        		dz1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_size(2,2));
	//        		dx2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_size(0,0));
	//        		dy2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_size(1,1));
	//        		dz2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_size(2,2));

	//        		dx3 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l][j].get_size(0,0));
	//        		dy3 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l][j].get_size(1,1));
	//        		dz3 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l][j].get_size(2,2));
        	                dx3 = 0.5*sqrt(12.0)*sqrt(hit3d.get_size(0,0));
                	        dy3 = 0.5*sqrt(12.0)*sqrt(hit3d.get_size(1,1));
                        	dz3 = 0.5*sqrt(12.0)*sqrt(hit3d.get_size(2,2));

       		 		cur_kappa = (*cur_seg)[i].kappa;
        			cur_dkappa = (*cur_seg)[i].dkappa;
        			cur_ux = (*cur_seg)[i].ux;
        			cur_uy = (*cur_seg)[i].uy;
        			cur_chi2 = (*cur_seg)[i].chi2;

				chi2=9999;
	        		calculate_kappa_tangents(
        	      			x1, y1, z1, x2, y2, z2, x3, y3, z3,
              				dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3,
              				kappa, dkappa, ux_mid, uy_mid, ux_end, uy_end,
              				dzdl_1, dzdl_2, ddzdl_1, ddzdl_2,
              				ca_sin_ang_cut, ca_cos_ang_cut_diff_inv,
              				cur_kappa, cur_dkappa, cur_ux, cur_uy, cur_chi2, chi2);
#ifdef _DEBUG_
        			cout<<"Extended layers for segment "<<which_seg<<endl;
        			cout<<"kappa "<<kappa<< " dkappa "<<dkappa
        			<<" ux_mid "<<ux_mid<<" uy_mid "<<" ux_end "<<ux_end<<" uy_end " <<uy_end
        			<<" dzdl_1 "<<dzdl_1<<" dzdl_2 "<<dzdl_2<<" ddzdl_1 "<<ddzdl_1<<" ddzdl_2 "<<ddzdl_2
        			<<" cur_chi2 "<<cur_chi2<<" chi2 "<<chi2<<" chi2*inv_layer "<<l <<" "<<chi2*inv_layer[l]
        			<<" chi2_cut " <<ca_chi2_layer_cut<<endl;
#endif


			// if segment with good chi2, store it in next segment
        			if (chi2 * inv_layer[l] < ca_chi2_layer_cut ) {
              				temp_segment.chi2 = chi2;
              				temp_segment.ux = ux_end;
              				temp_segment.uy = uy_end;
              				temp_segment.kappa = kappa;
					if (seeding_mode){
              					if (temp_segment.kappa > _hough_space->get_kappa_max() &&cur_seg_size!=1) continue;
					} else {
						if (temp_segment.kappa > _hough_space->get_kappa_max()) continue;
					} 
					
                        		temp_segment.dkappa = dkappa;
              				for (unsigned int ll = 0; ll < l; ++ll) {
                			temp_segment.hits[ll] = (*cur_seg)[which_seg].hits[ll];
              				}
              				temp_segment.hits[l] = hit1;
              				//unsigned int outer_layer =
                  			//	layer_sorted[l][temp_segment.hits[l]].get_layer();
              				//if (!forward) outer_layer = nlayers-outer_layer-1;
              				temp_segment.n_hits = l + 1;
              				// finish up if required number of layers is reached
                        		//if ((nlayers - (l + 1)) <= allowed_missing) {
                			//	complete_segments.push_back(temp_segment);
              				//}
              				// make sure we have required number of layers with hits
                        		//if ((outer_layer - l) > allowed_missing) {
                			//continue;
              				//}

		                	// Discard segment if too many missing layers                 		
              				if (next_seg->size() == next_seg_size) { // first new segment
                				next_seg->push_back(temp_segment);
                                        	missing_layers_map_next.insert(make_pair(next_seg->size()-1, missing_layers));
                				next_seg_size += 1;
              				} else { // next new segments
					#ifdef _DEBUG_
					cout<<"Next segment size inconsistent "<<endl;
					#endif
                				(*next_seg)[next_seg_size] = temp_segment;
				        	missing_layers_map_next.insert(make_pair(next_seg->size()-1, missing_layers));
                				next_seg_size += 1;
              				}
					++added_next_segments;
#ifdef _DEBUG_					
					cout<<"segment "<< which_seg << " added segment "<<added_next_segments<<endl;
#endif
         			} // chi2 cut from segment building method -> change to switch - case block
				}

			}// clusters on layer l		
			
/*
			if (added_next_segments==0){
			 ++missing_layers;
			if (missing_layers <= allowed_missing){
                                if (next_seg->size() == next_seg_size) { // first new segment
                                        next_seg->push_back((*cur_seg)[i]);
					missing_layers_map_next.insert(make_pair(next_seg->size()-1, missing_layers));
                                        next_seg_size += 1;
                                } else { // next new segments
                                        (*next_seg)[next_seg_size] = (*cur_seg)[i];
					missing_layers_map_next.insert(make_pair(next_seg->size()-1, missing_layers));
                                        next_seg_size += 1;
                                }
			}
			}
*/


        	} // current segments
	    	swap(cur_seg, next_seg);
    		swap(cur_seg_size, next_seg_size);
		missing_layers_map.swap(missing_layers_map_next);
		// Current segments now hold extended segments up to layer l, and missing layers info is in missing_layers_map 

	}// next segment


#endif

	SimpleTrack3D temp_track;
	temp_track.hits.assign(nlayers, SimpleHit3D());

	std::vector<SimpleTrack3D> best_track;
	std::vector<HelixKalmanState> best_track_state;
	float best_chi2 = 9999;
	for (unsigned int i = 0; i< cur_seg_size; ++i) {

		if ((*cur_seg)[i].n_hits ==0) continue;

#ifdef _DEBUG_
		cout<<"segment " <<i <<endl;
#endif
      		temp_track.hits.assign((*cur_seg)[i].n_hits, SimpleHit3D());
//		float seg_kappa = (*cur_seg)[i].kappa;

      		for (unsigned int l = 0; l < (*cur_seg)[i].n_hits; ++l) {
		if (l<3){
        	temp_track.hits[l] = layer_sorted[l][(*cur_seg)[i].hits[l]];
		}else {
		auto search = _hits_map.find((*cur_seg)[i].hits[l]);
		SimpleHit3D cluster = search->second;
		temp_track.hits[l] = cluster;
		}
      		}

      		float init_chi2 = temp_track.fit_track();
		if(verbose > 0)  cout<<"chi2 from fit_track "<<init_chi2 <<" kappa "<< temp_track.kappa <<endl;

#ifdef _DEBUG_
          	cout	<<" kappa " <<temp_track.kappa <<" phi "<<temp_track.phi<<" d "<<temp_track.d
          		<<" z0 "<<temp_track.z0<<" dzdl "<<temp_track.dzdl<< endl;
#endif
		
		if (seeding_mode){
			if (temp_track.kappa != temp_track.kappa && cur_seg_size!=1) continue;
			if (temp_track.z0 != temp_track.z0 && cur_seg_size != 1) continue;
		} else {
		        if (temp_track.kappa != temp_track.kappa) continue;
                        if (temp_track.z0 != temp_track.z0) continue;
		}

    		HelixKalmanState state;
    		state.phi = temp_track.phi;
    		if (state.phi < 0.) {
      		state.phi += 2. * M_PI;
    		}
    		state.d = temp_track.d;
   		state.kappa = temp_track.kappa;
    		state.nu = sqrt(state.kappa);
    		state.z0 = temp_track.z0;
    		state.dzdl = temp_track.dzdl;
    		state.C = Matrix<float, 5, 5>::Zero(5, 5);
    		state.C(0, 0) = pow(0.01, 2.);
    		state.C(1, 1) = pow(0.01, 2.);
    		state.C(2, 2) = pow(0.01 * state.nu, 2.);
    		state.C(3, 3) = pow(0.05, 2.);
    		state.C(4, 4) = pow(0.05, 2.);
    		state.chi2 = 0.;
    		state.position = 0;
    		state.x_int = 0.;
    		state.y_int = 0.;
    		state.z_int = 0.;

    		unsigned int nfits = 0;
    		for (unsigned int h = 0; h < temp_track.hits.size(); ++h) {
      		_kalman->addHit(temp_track.hits[h], state);
      		nfits += 1;
#ifdef _DEBUG_
		cout<<"nfits "<<nfits<<endl;
#endif
    		}
#ifdef _DEBUG_
    		cout<<"z0 after kalman "<<state.z0<<endl;
#endif
   		// fudge factor for non-gaussian hit sizes
           	state.C *= 3.;
    		state.chi2 *= 6.;
		// kappa cut *here* drives both efficiency and ghost track rates down at the same time
		if (seeding_mode){
    		if (!(temp_track.kappa == temp_track.kappa) && (cur_seg_size !=1)) continue;
    		if (!(state.chi2 == state.chi2) && (cur_seg_size !=1)) continue;
		} else {
                if (!(temp_track.kappa == temp_track.kappa)) continue;
                if (!(state.chi2 == state.chi2)) continue;
		}

/*
      		if (fabs(temp_track.d) > ca_dcaxy_cut) continue;
      		if (fabs(temp_track.z0) > dca_cut) continue;
*/

	// no chi2 cut for tests on triplets

		//cout<<"state.chi2 from kalman "<<state.chi2<<endl;
		if (seeding_mode){
    			if (state.chi2 / (2. * ((float)(temp_track.hits.size())) - 5.) > ca_chi2_cut && cur_seg_size !=1) continue;
		} else {
			if (state.chi2 / (2. * ((float)(temp_track.hits.size())) - 5.) > ca_chi2_cut) continue;
		}

    		if (best_chi2 > state.chi2 || (seeding_mode && cur_seg_size==1)){
      			if (!best_track.empty()){
      			best_track.pop_back();
      			best_track_state.pop_back();
      			}
      			best_track.push_back(temp_track);
     			best_track_state.push_back(state);
    		}
  	}

	if (best_track.empty()) return 1;

	if(verbose > 0) cout<<"best_track.size "<<best_track.size()<<endl;
  	ca_tracks.push_back(best_track.back());
  	ca_track_states.push_back(best_track_state.back());
	if(verbose > 0) cout <<"ca track added, chi2 =  "<< (best_track_state.back().chi2)/(2. * ((float)(temp_track.hits.size())) - 5.) <<" z0 = "<<best_track_state.back().z0<<endl;


	//    if ((remove_hits == true) && (state.chi2 < chi2_removal_cut) &&
      	//        (temp_track.hits.size() >= n_removal_hits)) {
        temp_track = best_track.back();

        if  (remove_hits){
        for (unsigned int i = 0; i < temp_track.hits.size(); ++i) {
                if (!remove_inner_hits && temp_track.hits[i].get_layer()<3) continue;
                auto search = _hits_used.find(temp_track.hits[i].get_id());
                if(search != _hits_used.end())
                {
                _hits_used.find(temp_track.hits[i].get_id())->second = true;
                }
        }
        }

	segments1.clear();
	segments2.clear();
	return 1;
}

int CellularAutomaton_v1::process_single_track(SimpleTrack3D& track)
{

        std::vector<TrackSegment> segments1;
        std::vector<TrackSegment> segments2;


  std::vector<TrackSegment>* cur_seg = &segments1;
  std::vector<TrackSegment>* next_seg = &segments2;
  unsigned int cur_seg_size = 0;
  unsigned int next_seg_size = 0;

  std::vector<TrackSegment> complete_segments;

  unsigned int allowed_missing = nlayers - rlayers;
  cout<<"allowed missing "<< allowed_missing<<endl;


  	// l is not actual layer number, it is ith layer
  	for (unsigned int l = 0; l < nlayers; ++l) {
    	layer_sorted[l].clear();
  	}
 
//  cout<<"track.hits.size "<<track.hits.size()<<endl;
  for (unsigned int i = 0; i < track.hits.size(); ++i) {
    SimpleHit3D hit = track.hits[i];
    unsigned int layer = (unsigned int) hit.get_layer();
    if (layer > (nlayers-1)) continue;
    if (!forward) layer = nlayers-layer-1;
    unsigned int min = (layer - allowed_missing);
    if (allowed_missing > layer) {
      min = 0;
    }
    for (unsigned int l = min; l <= layer; ++l) {
      layer_sorted[l].push_back(hit);
    }

  }

  for (unsigned int l = 0; l < nlayers; ++l) {
//    cout<<"layer_sorted["<<l<<"].size = "<< layer_sorted[l].size()<<endl;
    if (layer_sorted[l].size() == 0) {
      return 0;
    }
  }

  timeval t1, t2;
  double time1 = 0.;
  double time2 = 0.;

  gettimeofday(&t1, nullptr);

  float ca_cos_ang_cut_diff = 1. - ca_cos_ang_cut;
  float ca_cos_ang_cut_diff_inv = 1. / ca_cos_ang_cut_diff;
  float ca_sin_ang_cut = sqrt(1. - ca_cos_ang_cut * ca_cos_ang_cut);

  std::vector<float> inv_layer;
  inv_layer.assign(nlayers, 1.);
  for (unsigned int l = 3; l < nlayers; ++l) {
    inv_layer[l] = 1. / (((float)l) - 2.);
  }

  float x1, x2, x3;
  float y1, y2, y3;
  float z1, z2, z3;
  float dx1, dx2, dx3;
  float dy1, dy2, dy3;
  float dz1, dz2, dz3;

  float kappa;
  float dkappa;

  float ux_mid;
  float uy_mid;
  float ux_end;
  float uy_end;

  float dzdl_1;
  float dzdl_2;
  float ddzdl_1;
  float ddzdl_2;


  float cur_kappa;
  float cur_dkappa;
  float cur_ux;
  float cur_uy;
  float cur_chi2;
  float chi2;

  unsigned int hit1;
  unsigned int hit2;
  unsigned int hit3;

  TrackSegment temp_segment;
  temp_segment.hits.assign(nlayers, 0);

  for (unsigned int i = 0; i < layer_sorted[0].size(); ++i) {
    for (unsigned int j = 0; j < layer_sorted[1].size(); ++j) {
      for (unsigned int k = 0; k < layer_sorted[2].size() ; ++k) {

	unsigned int layer0 = layer_sorted[0][i].get_layer(); 
	unsigned int layer1 = layer_sorted[1][j].get_layer();
	unsigned int layer2 = layer_sorted[2][k].get_layer();
	if (!forward) 
	{
	layer0 = nlayers-layer0-1;
	layer1 = nlayers-layer1-1;
	layer2 = nlayers-layer2-1;
	}
        if ((layer0 >= layer1) || (layer1 >= layer2)) {
          continue;
        }

        x1 = layer_sorted[0][i].get_x();
        y1 = layer_sorted[0][i].get_y();
        z1 = layer_sorted[0][i].get_z();

	// sigma ?= half pictch
        dx1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[0][i].get_size(0,0));
        dy1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[0][i].get_size(1,1));
        dz1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[0][i].get_size(2,2));

        x2 = layer_sorted[1][j].get_x();
        y2 = layer_sorted[1][j].get_y();
        z2 = layer_sorted[1][j].get_z();

        dx2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[1][j].get_size(0,0));
        dy2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[1][j].get_size(1,1));
        dz2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[1][j].get_size(2,2));

        x3 = layer_sorted[2][k].get_x();
        y3 = layer_sorted[2][k].get_y();
        z3 = layer_sorted[2][k].get_z();
        dx3 = 0.5*sqrt(12.0)*sqrt(layer_sorted[2][k].get_size(0,0));
        dy3 = 0.5*sqrt(12.0)*sqrt(layer_sorted[2][k].get_size(1,1));
        dz3 = 0.5*sqrt(12.0)*sqrt(layer_sorted[2][k].get_size(2,2));

	// layer of hit
        hit1 = i;
        hit2 = j;
        hit3 = k;

	calculate_kappa_tangents(x1, y1, z1, x2, y2, z2, x3, y3,
				z3, dx1, dy1, dz1, dx2, dy2, dz2,
				dx3, dy3, dz3, kappa, dkappa,
				ux_mid, uy_mid, ux_end, uy_end,
				dzdl_1, dzdl_2, ddzdl_1, ddzdl_2);

#ifdef _DEBUG_
	cout<<"Triplet : "<<endl;
	cout<<"kappa "<<kappa<< " dkappa "<<dkappa<<" ux_mid "<<ux_mid<<" uy_mid "<<" ux_end "<<ux_end<<
	" uy_end " <<uy_end<<" dzdl_1 "<<dzdl_1<<" dzdl_2 "<<dzdl_2<<" ddzdl_1 "<<ddzdl_1<<" ddzdl_2 "<<ddzdl_2<<endl;
#endif

          temp_segment.chi2 = pow(
              (dzdl_1 - dzdl_2) /
              (ddzdl_1 + ddzdl_2 + fabs(dzdl_1 * ca_sin_ang_cut)),2);
          if (temp_segment.chi2 > ca_chi2_layer_cut) continue;
          temp_segment.ux = ux_end;
          temp_segment.uy = uy_end;
          temp_segment.kappa = kappa;
//          if (temp_segment.kappa > _hough_space->get_kappa_max()) continue;
          temp_segment.dkappa = dkappa;
          temp_segment.hits[0] = hit1;
          temp_segment.hits[1] = hit2;
          temp_segment.hits[2] = hit3;
          temp_segment.n_hits = 3;
          unsigned int outer_layer =
          	layer_sorted[2][temp_segment.hits[2]].get_layer();
	  if (!forward) outer_layer = nlayers-outer_layer-1;
	  // make sure we have required number of layers with hits
          if ((outer_layer - 2) > allowed_missing) continue;
	  // finish up if required number of layers is reached,
          if ((nlayers - 3) <= allowed_missing) {
          	complete_segments.push_back(temp_segment);
          }
          if (next_seg->size() == next_seg_size) { // first new segment
          	next_seg->push_back(temp_segment);
          	next_seg_size += 1;
          } else { // next new segments
          	(*next_seg)[next_seg_size] = temp_segment;
          	next_seg_size += 1;
          }
	
      }
    }
  }

  swap(cur_seg, next_seg);
  swap(cur_seg_size, next_seg_size);

  // cout<<"number of segments from first 3 layers : " << cur_seg_size<<endl;
  unsigned int which_seg;

  // add hits to segments layer-by-layer, cutting out bad segments
  for (unsigned int l = 3; l < nlayers; ++l) {
	if (l == (nlayers - 1)) {
//	ca_chi2_cut_layer = 0.25* ca_chi2_cut;// 2.*0.25 = 0.8 less loose cut after adding all clusters 
	}
    	next_seg_size = 0;
    for (unsigned int i = 0; i < cur_seg_size; ++i) {
      for (unsigned int j = 0; j < layer_sorted[l].size(); ++j) {

	unsigned int layer0 = layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_layer();
	unsigned int layer1 = layer_sorted[l][j].get_layer();
	if (!forward) 
	{
		layer0 = nlayers-layer0-1;
		layer1 = nlayers-layer1-1;
	}
	if (layer0 >= layer1) continue;

        x1 = layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_x();
       	y1 = layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_y();
        z1 = layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_z();
        x2 = layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_x();
        y2 = layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_y();
        z2 = layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_z();
        x3 = layer_sorted[l][j].get_x();
        y3 = layer_sorted[l][j].get_y();
        z3 = layer_sorted[l][j].get_z();

        dx1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_size(0,0));
        dy1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_size(1,1));
        dz1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_size(2,2));
        dx2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_size(0,0));
        dy2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_size(1,1));
        dz2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_size(2,2));
        dx3 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l][j].get_size(0,0));
        dy3 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l][j].get_size(1,1));
        dz3 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l][j].get_size(2,2));

        cur_kappa = (*cur_seg)[i].kappa;
        cur_dkappa = (*cur_seg)[i].dkappa;
        cur_ux = (*cur_seg)[i].ux;
        cur_uy = (*cur_seg)[i].uy;
        cur_chi2 = (*cur_seg)[i].chi2;

        which_seg = i;
        hit1 = j;

	calculate_kappa_tangents(
              x1, y1, z1, x2, y2, z2, x3, y3, z3, 
	      dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3, 
	      kappa, dkappa, ux_mid, uy_mid, ux_end, uy_end, 
	      dzdl_1, dzdl_2, ddzdl_1, ddzdl_2, 
	      ca_sin_ang_cut, ca_cos_ang_cut_diff_inv,
              cur_kappa, cur_dkappa, cur_ux, cur_uy, cur_chi2, chi2);

#ifdef _DEBUG_
	cout<<"Extended layers for segment "<<which_seg<<endl;
	cout<<"kappa "<<kappa<< " dkappa "<<dkappa
	<<" ux_mid "<<ux_mid<<" uy_mid "<<" ux_end "<<ux_end<<" uy_end " <<uy_end
	<<" dzdl_1 "<<dzdl_1<<" dzdl_2 "<<dzdl_2<<" ddzdl_1 "<<ddzdl_1<<" ddzdl_2 "<<ddzdl_2
	<<" cur_chi2 "<<cur_chi2<<" chi2 "<<chi2<<" chi2*inv_layer "<<l <<" "<<chi2*inv_layer[l]
	<<" chi2_cut" <<ca_chi2_layer_cut<<endl;
#endif
	if (chi2 * inv_layer[l] < ca_chi2_layer_cut) {
              temp_segment.chi2 = chi2;
              temp_segment.ux = ux_end;
              temp_segment.uy = uy_end;
              temp_segment.kappa = kappa;
              //if (temp_segment.kappa > _hough_space->get_kappa_max()) continue;
              temp_segment.dkappa = dkappa;
              for (unsigned int ll = 0; ll < l; ++ll) {
                temp_segment.hits[ll] = (*cur_seg)[which_seg].hits[ll];
              }
              temp_segment.hits[l] = hit1;
              unsigned int outer_layer =
                  layer_sorted[l][temp_segment.hits[l]].get_layer();
	      if (!forward) outer_layer = nlayers-outer_layer-1;
              temp_segment.n_hits = l + 1;
	      // finish up if required number of layers is reached
              if ((nlayers - (l + 1)) <= allowed_missing) {
                complete_segments.push_back(temp_segment);
              }
	      // make sure we have required number of layers with hits
              if ((outer_layer - l) > allowed_missing) {
                continue;
              }
              if (next_seg->size() == next_seg_size) { // first new segment
                next_seg->push_back(temp_segment);
                next_seg_size += 1;
              } else { // next new segments
                (*next_seg)[next_seg_size] = temp_segment;
                next_seg_size += 1;
              }
         }
      }// j 
    }// i

    swap(cur_seg, next_seg);
    swap(cur_seg_size, next_seg_size);

  }//l = 3

  //cout<<"number of complete segments : " << complete_segments.size()<<endl;
  //cout<<"number of current segments : "<< cur_seg_size<<endl;

  // copy complete segments over to current segments
  for (unsigned int i = 0; i < complete_segments.size(); ++i) {
    if (cur_seg->size() == cur_seg_size) {
      cur_seg->push_back(complete_segments[i]);
      ++cur_seg_size;
    } else {
      (*cur_seg)[cur_seg_size] = complete_segments[i];
      ++cur_seg_size;
    }
  }

  gettimeofday(&t2, nullptr);
  time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
  time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
  CAtime += (time2 - time1);

  std::set<unsigned int> comp1;
  std::set<unsigned int> comp2;

  //cout<<"number of segments generated "<< cur_seg_size<<endl;
  if (cur_seg_size==0  || cur_seg_size>10000) return 1;
  for (unsigned int i = cur_seg_size-1; i>0; --i){
    if ((*cur_seg)[i].n_hits==0) continue;
    comp1.clear();
    temp_combo.assign((*cur_seg)[i].n_hits, 0);
    for (unsigned int l = 0; l < (*cur_seg)[i].n_hits; ++l) {
    temp_combo[l] = layer_sorted[l][(*cur_seg)[i].hits[l]].get_id();
    comp1.insert(temp_combo[l]);
    }

    sort(temp_combo.begin(), temp_combo.end());
    set<vector<unsigned int> >::iterator it = combos.find(temp_combo);
    if (it != combos.end()) {
      (*cur_seg)[i].n_hits = 0.;
    }

//    if (combos.size() > 10000) {
    if (combos.size() > 100000) { 
      combos.clear();
    }
    combos.insert(temp_combo);

    for (unsigned int j = i-1; j>=0; --j){
      comp2.clear();
      for (unsigned int m = 0; m < (*cur_seg)[j].n_hits; ++m){
      comp2.insert(layer_sorted[m][(*cur_seg)[j].hits[m]].get_id());
      };
      for (std::set<unsigned int>::iterator it=comp1.begin(); it!=comp1.end(); ++it){
      auto it2 = comp2.find(*it);
      if (it2 != comp2.end()) comp2.erase(*it2);
      }	
      if (comp2.empty()) {
	(*cur_seg)[j].n_hits = 0;
      }
      if (j==0) break;
    }
  }

  unsigned int nsegs = cur_seg_size;
  for (unsigned int i = 0; i<cur_seg_size; i++) {
  if ((*cur_seg)[i].n_hits ==0) --nsegs;
  }
  //cout<<"number of segments to be processed "<< nsegs<<endl;

  SimpleTrack3D temp_track;
  temp_track.hits.assign(nlayers, SimpleHit3D());

  std::vector<SimpleTrack3D> best_track;
  std::vector<HelixKalmanState> best_track_state;
  float best_chi2 = 9999;
  for (unsigned int i = 0; i< cur_seg_size; ++i) {

    if ((*cur_seg)[i].n_hits ==0) continue;

#ifdef _DEBUG_
      cout<<"segment " <<i <<endl;
#endif
      temp_track.hits.assign((*cur_seg)[i].n_hits, SimpleHit3D());

      for (unsigned int l = 0; l < (*cur_seg)[i].n_hits; ++l) {
        temp_track.hits[l] = layer_sorted[l][(*cur_seg)[i].hits[l]];
      }

      unsigned int ninner_hits =0;
      for (unsigned int n = 0; n < temp_track.hits.size(); ++n){
        if (temp_track.hits[n].get_layer()<3)
	++ninner_hits;
      }
     
      if(require_inner_hits && ninner_hits<2) continue;

      gettimeofday(&t1, nullptr);

      float init_chi2 = temp_track.fit_track();
#ifdef _DEBUG_
      cout<<"chi2 from fit_track "<<init_chi2
  	  <<" kappa " <<temp_track.kappa <<" phi "<<temp_track.phi<<" d "<<temp_track.d
	  <<" z0 "<<temp_track.z0<<" dzdl "<<temp_track.dzdl<< endl;
#endif

    // not being used for the time being
    if (init_chi2 > fast_chi2_cut_max) { 
      if (init_chi2 > fast_chi2_cut_par0 +
                          fast_chi2_cut_par1 / kappa_to_pt(temp_track.kappa)) {
        gettimeofday(&t2, nullptr);
        time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
        time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
        KALtime += (time2 - time1);
        continue;
      }
    }

    HelixKalmanState state;
    state.phi = temp_track.phi;
    if (state.phi < 0.) {
      state.phi += 2. * M_PI;
    }
    state.d = temp_track.d;
    state.kappa = temp_track.kappa;
    state.nu = sqrt(state.kappa);
    state.z0 = temp_track.z0;
    state.dzdl = temp_track.dzdl;
    state.C = Matrix<float, 5, 5>::Zero(5, 5);
    state.C(0, 0) = pow(0.01, 2.);
    state.C(1, 1) = pow(0.01, 2.);
    state.C(2, 2) = pow(0.01 * state.nu, 2.);
    state.C(3, 3) = pow(0.05, 2.);
    state.C(4, 4) = pow(0.05, 2.);
    state.chi2 = 0.;
    state.position = 0;
    state.x_int = 0.;
    state.y_int = 0.;
    state.z_int = 0.;
 
    unsigned int nfits = 0;
    for (unsigned int h = 0; h < temp_track.hits.size(); ++h) {
      _kalman->addHit(temp_track.hits[h], state);
      nfits += 1;
    }
    if(verbose > 0) cout<<"z0 after kalman "<<state.z0<<endl;

    // fudge factor for non-gaussian hit sizes
    state.C *= 3.;
    state.chi2 *= 6.;

    gettimeofday(&t2, nullptr);
    time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
    time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
    KALtime += (time2 - time1);

    if (!(temp_track.kappa == temp_track.kappa)) {
      continue;
    }
/*
    if (temp_track.kappa > _hough_space->get_kappa_max()) {
      continue;
    }
*/
    if (!(state.chi2 == state.chi2)) {
      continue;
    }
    if (state.chi2 / (2. * ((float)(temp_track.hits.size())) - 5.) > ca_chi2_cut) {
      continue;
    }
/*
    if (cut_on_dca == true) {
      if (fabs(temp_track.d) > dca_cut) {
        continue;
      }
      if (fabs(temp_track.z0) > dca_cut) {
        continue;
      }
    }
*/
    if (best_chi2 > state.chi2){
      if (!best_track.empty()){
      best_track.pop_back();
      best_track_state.pop_back();
      }
      best_track.push_back(temp_track);
      best_track_state.push_back(state);
    }
  }
  
  if (best_track.empty()) return 1;

  ca_tracks.push_back(best_track.back());
  ca_track_states.push_back(best_track_state.back());
  if(verbose > 0) cout <<"ca track added, chi2 =  "<< (best_track_state.back().chi2)/(2. * ((float)(temp_track.hits.size())) - 5.) <<" z0 = "<<best_track_state.back().z0<<endl;


//    if ((remove_hits == true) && (state.chi2 < chi2_removal_cut) &&
//        (temp_track.hits.size() >= n_removal_hits)) {
	temp_track = best_track.back();

	if  (remove_hits){
	for (unsigned int i = 0; i < temp_track.hits.size(); ++i) {
		if (!remove_inner_hits && temp_track.hits[i].get_layer()<3) continue;
		auto search = _hits_used.find(temp_track.hits[i].get_id());
		if(search != _hits_used.end())
        	{
        	_hits_used.find(temp_track.hits[i].get_id())->second = true;
        	}
	}
	}

        segments1.clear();
        segments2.clear();


  return 1;
}

int CellularAutomaton_v1::calculate_kappa_tangents(
			float x1, float y1, float z1, float x2, float y2, float z2, 
			float x3, float y3, float z3, 
			float dx1, float dy1, float dz1, float dx2, float dy2, float dz2,
                        float dx3, float dy3, float dz3, 
			float& kappa, float& dkappa,
                        float& ux_mid, float& uy_mid, float& ux_end, float& uy_end,
                        float& dzdl_1, float& dzdl_2, float& ddzdl_1, float& ddzdl_2)
{

	float D12 = sqrt(pow(x2-x1,2)+pow(y2-y1,2));
	float D23 = sqrt(pow(x3-x2,2)+pow(y3-y2,2));
	float D13 = sqrt(pow(x3-x1,2)+pow(y3-y1,2));
	kappa = 1./(D12*D23*D13);
	float num = (D12+D23+D13)*(D23+D13-D12)*(D12+D13-D23)*(D12+D23-D13);
	if (num<0) num = 0;
	num = sqrt(num);	
	kappa *=num;
	
	float kappa_inv = 1/kappa;
	float D12_inv = 1./D12;
	float D23_inv = 1./D23;
	float D13_inv = 1./D13;

	float dr1 = sqrt(pow(dx1,2)+pow(dy1,2));
	float dr2 = sqrt(pow(dx2,2)+pow(dy2,2));
	float dr3 = sqrt(pow(dx3,2)+pow(dy3,2));

	float dk1 = (dr1+dr2)* D12_inv*D12_inv;
	float dk2 = (dr2+dr3)* D23_inv*D23_inv;
	dkappa = dk1+dk2;

	float ux12 = (x2-x1)*D12_inv;
	float uy12 = (y2-y1)*D12_inv;
	float ux23 = (x3-x2)*D23_inv;
	float uy23 = (y3-y2)*D23_inv;
	float ux13 = (x3-x1)*D13_inv;
	float uy13 = (y3-y1)*D13_inv;

	// cos(alpha) = cos(alpha12 - alpha13)
	// sin(alpha) = sin(alpha12 - alpha13)
	float cosalpha = ux12*ux13 + uy12*uy13;
	float sinalpha = uy12*ux13 - ux12*uy13;
	// alpha23 + alpha  
	ux_mid = ux23 * cosalpha - uy23 * sinalpha;
	uy_mid = ux23 * sinalpha + uy23 * cosalpha;
	// alpha23 - alpha
	ux_end = ux23 * cosalpha + uy23 * sinalpha;
	uy_end = uy23 * cosalpha - ux23 * sinalpha;

	// dzdl = dz/sqrt(ds^2 + dz^2)
	float ds23 = 2.*kappa_inv*atan(sinalpha/(1.+ sqrt(1.-pow(sinalpha,2))));
	if (kappa<=0)  ds23 = D23;

	float dz23 = z3 - z2;
	dzdl_2 =  dz23/sqrt(pow(ds23,2) + pow(dz23,2));
	ddzdl_2 = (dz2 + dz3)*D23_inv;

	// sin(alpha) = sin(alpha13 -alpha23)
	sinalpha = ux13 *uy23 - ux23 *  uy13;
	float ds12 = 2.*kappa_inv*atan(sinalpha/(1.+sqrt(1.-pow(sinalpha,2))));
	if (kappa<=0) ds12 = D12;

	float dz12 = z2 - z1;
	dzdl_1 = dz12/sqrt(pow(ds12,2) + pow(dz12,2));
	ddzdl_1 = (dz1 + dz2) * D12_inv;

	return 1;

}

int CellularAutomaton_v1::calculate_kappa_tangents(
                        float x1, float y1, float z1, float x2, float y2, float z2,
                        float x3, float y3, float z3,
                        float dx1, float dy1, float dz1, float dx2, float dy2, float dz2,
                        float dx3, float dy3, float dz3,
                        float& kappa, float& dkappa,
                        float& ux_mid, float& uy_mid, float& ux_end, float& uy_end,
                        float& dzdl_1, float& dzdl_2, float& ddzdl_1, float& ddzdl_2,
			float ca_sin_ang_cut, float ca_cos_ang_cut_diff_inv,
              		float cur_kappa, float cur_dkappa, float cur_ux, float cur_uy, 
			float cur_chi2, float& chi2)
{


        float D12 = sqrt(pow(x2-x1,2)+pow(y2-y1,2));
        float D23 = sqrt(pow(x3-x2,2)+pow(y3-y2,2));
        float D13 = sqrt(pow(x3-x1,2)+pow(y3-y1,2));
        kappa = 1./(D12*D23*D13);
        float num = (D12+D23+D13)*(D23+D13-D12)*(D12+D13-D23)*(D12+D23-D13);
        if (num<0) num = 0;
        num = sqrt(num);
        kappa *=num;

        float kappa_inv = 1/kappa;
        float D12_inv = 1./D12;
        float D23_inv = 1./D23;
        float D13_inv = 1./D13;

        float dr1 = sqrt(pow(dx1,2)+pow(dy1,2));
        float dr2 = sqrt(pow(dx2,2)+pow(dy2,2));
        float dr3 = sqrt(pow(dx3,2)+pow(dy3,2));

        float dk1 = (dr1+dr2)* D12_inv*D12_inv;
        float dk2 = (dr2+dr3)* D23_inv*D23_inv;
        dkappa = dk1+dk2;

        float ux12 = (x2-x1)*D12_inv;
        float uy12 = (y2-y1)*D12_inv;
        float ux23 = (x3-x2)*D23_inv;
        float uy23 = (y3-y2)*D23_inv;
        float ux13 = (x3-x1)*D13_inv;
        float uy13 = (y3-y1)*D13_inv;

        // cos(alpha) = cos(alpha12 - alpha13)
        // sin(alpha) = sin(alpha12 - alpha13)
        float cosalpha = ux12*ux13 + uy12*uy13;
        float sinalpha = uy12*ux13 - ux12*uy13;
        // alpha23 + alpha  
        ux_mid = ux23 * cosalpha - uy23 * sinalpha;
        uy_mid = ux23 * sinalpha + uy23 * cosalpha;
        // alpha23 - alpha
        ux_end = ux23 * cosalpha + uy23 * sinalpha;
        uy_end = uy23 * cosalpha - ux23 * sinalpha;

        // dzdl = dz/sqrt(ds^2 + dz^2)
        float ds23 = 2.*kappa_inv*atan(sinalpha/(1.+ sqrt(1.-pow(sinalpha,2))));
        if (kappa<=0)  ds23 = D23;

        float dz23 = z3 - z2;
        dzdl_2 =  dz23/sqrt(pow(ds23,2) + pow(dz23,2));
        ddzdl_2 = (dz2 + dz3)*D23_inv;

        // sin(alpha) = sin(alpha13 -alpha23)
        sinalpha = ux13 *uy23 - ux23 *  uy13;
        float ds12 = 2.*kappa_inv*atan(sinalpha/(1.+sqrt(1.-pow(sinalpha,2))));
        if (kappa<=0) ds12 = D12;

        float dz12 = z2 - z1;
        dzdl_1 = dz12/sqrt(pow(ds12,2) + pow(dz12,2));
        ddzdl_1 = (dz1 + dz2) * D12_inv;

	float kappa_diff = cur_kappa - kappa;
	float n_dk = cur_dkappa	+ dkappa + ca_sin_ang_cut * kappa;
	float chi2_kappa = pow(kappa_diff,2)/pow(n_dk,2);

	float cos_scatter = cur_ux * ux_mid + cur_uy * uy_mid;
	float chi2_ang = pow((1-cos_scatter)*ca_cos_ang_cut_diff_inv,2);

	float sin_scatter = dzdl_1 * ca_sin_ang_cut;
	float chi2_dzdl = 0.5*pow((dzdl_1-dzdl_2)/(ddzdl_1+ddzdl_2+fabs(sin_scatter)),2);

	chi2 = cur_chi2 + chi2_ang + chi2_kappa + chi2_dzdl;
	return 1;
}


float CellularAutomaton_v1::kappa_to_pt(float kappa) {
        return _pt_rescale * _mag_field / 333.6 / kappa;
}

float CellularAutomaton_v1::shift_phi_range(float _phi){

        if (_phi < 0.) _phi += 2.*M_PI;
        return _phi;
}


