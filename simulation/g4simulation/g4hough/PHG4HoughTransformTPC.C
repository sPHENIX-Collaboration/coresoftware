#include "PHG4HoughTransformTPC.h"

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

// PHENIX Geant4 includes
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4main/PHG4InEvent.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>

// PHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

// Geant4 includes
#include <Geant4/G4MagneticField.hh>
#include <Geant4/G4TransportationManager.hh>
#include <Geant4/G4FieldManager.hh>

// Helix Hough includes
#include <SimpleHit3D.h>
#include <SimpleTrack3D.h>
#include <HelixResolution.h>
#include <HelixRange.h>
#include <HelixHough.h>
#include <VertexFinder.h>

// ROOT includes
#include <TH1D.h>

// standard includes
#include <cmath>
#include <iostream>

#include "SimpleTrack.h"
#include "TTree.h"
#include "TFile.h"

using findNode::getClass;
using namespace std;
using namespace Eigen;


static void convertHelixCovarianceToEuclideanCovariance( float B, float phi, float d, float kappa, float z0, float dzdl, Eigen::Matrix<float,5,5> const& input, Eigen::Matrix<float,6,6>& output )
{
  Matrix<float,6,5> J = Matrix<float,6,5>::Zero(6,5);
  // phi,d,nu,z0,dzdl
  // -->
  // x,y,z,px,py,pz
  float nu = sqrt(kappa);
  float dk_dnu = 2*nu;

  float cosphi = cos(phi);
  float sinphi = sin(phi);

  J( 0, 0 ) = -sinphi*d;
  J( 0, 1 ) = cosphi;
  J( 1, 0 ) = d*cosphi;
  J( 1, 1 ) = sinphi;
  J( 2, 3 ) = 1;

  float pt = 0.003*B/kappa;
  float dpt_dk = -0.003*B/(kappa*kappa);

  J( 3, 0 ) = -cosphi*pt;
  J( 3, 2 ) = -sinphi*dpt_dk*dk_dnu;
  J( 4, 0 ) = -sinphi*pt;
  J( 4, 2 ) = cosphi*dpt_dk*dk_dnu;

  float alpha = 1./(1. - dzdl*dzdl);
  float alpha_half = pow( alpha, 0.5 );
  float alpha_3_half = alpha*alpha_half;

  J( 5, 2 ) = dpt_dk*dzdl*alpha_half*dk_dnu;
  J( 5, 4 ) = pt*( alpha_half + dzdl*dzdl*alpha_3_half )*dk_dnu;

  output = J*input*(J.transpose());
}


static inline double sign(double x)
{
  return ((double)(x > 0.)) - ((double)(x < 0.));
}


void PHG4HoughTransformTPC::projectToRadius(const SvtxTrack* track,
					    double B,
					    double radius,
					    std::vector<double>& intersection) {
  intersection.clear();
  intersection.assign(3,NAN);

  // start from the inner most state vector
  const SvtxTrackState* state = track->get_state(0.0);
  projectToRadius(state,track->get_charge(),B,radius,intersection);
  
  // iterate once to see if there is a state vector closer to the intersection
  if (track->size_states() == 1) return;

  const SvtxTrackState* closest = NULL;  
  float min_dist = FLT_MAX;
  for (SvtxTrack::ConstStateIter iter = track->begin_states();
       iter != track->end_states();
       ++iter) {
    const SvtxTrackState* candidate = iter->second;
    float dist = sqrt(pow(candidate->get_x()-intersection[0],2) +
		      pow(candidate->get_y()-intersection[1],2) +
		      pow(candidate->get_z()-intersection[2],2));

    if (dist < min_dist) {
      closest = candidate;
      min_dist = dist;
    }
  }

  // if we just got back the previous case, bail
  if (closest->get_pathlength() == 0.0) return;

  // recompute using the closer state vector
  projectToRadius(closest,track->get_charge(),B,radius,intersection);

  return;
}

void PHG4HoughTransformTPC::projectToRadius(const SvtxTrackState* state,
					    int charge,
					    double B,
					    double radius,
					    std::vector<double>& intersection) {
  intersection.clear();
  intersection.assign(3,NAN);
  
  // find 2d intersections in x,y plane
  std::set<std::vector<double> > intersections;
  if (B != 0.0) {
    // magentic field present, project track as a circle leaving the state position
    
    // compute the center of rotation and the helix parameters
    // x(u) = r*cos(q*u+cphi) + cx
    // y(u) = r*sin(q*u+cphi) + cy
    // z(u) = b*u + posz;
    
    double cr = state->get_pt() * 333.6 / B;                                  // radius of curvature
    double cx = state->get_x() - (state->get_py()*cr)/charge/state->get_pt(); // center of rotation, x
    double cy = (state->get_px()*cr)/charge/state->get_pt() + state->get_y(); // center of rotation, y
    double cphi = atan2(state->get_y()-cy,state->get_x()-cx);                 // phase of state position
    double b = state->get_pz()/state->get_pt()*cr;                            // pitch along z

    if (!circle_circle_intersections(0.0,0.0,radius,
				    cx,cy,cr,
				    &intersections)) {
      return;
    }

    if (intersections.empty()) return;
    
    // downselect solutions based on track direction
    // we want the solution where the state vector would exit the cylinder
    // this can be determined by the direction that the track circulates in

    // rotate the px,py to the postion of the solution
    // then ask if the dot product of the displacement vector between the solution
    // and the cylinder center with the rotated momentum vector is positive
    std::set<std::vector<double> >::iterator remove_iter = intersections.end();
    double intersection_z = 0.0;
    for (std::set<std::vector<double> >::iterator iter = intersections.begin();
	 iter != intersections.end();
	 ++iter) {
      double x = iter->at(0);
      double y = iter->at(1);

      // find the azimuthal rotation about the center of rotation between the state vector and the solution

      // displacement between center of rotation and solution
      double dx = x - cx;
      double dy = y - cy;
      double dphi = atan2(dy,dx);
       
      // displacement between center of rotation and state position
      double state_dx = state->get_x() - cx;
      double state_dy = state->get_y() - cy;
      double state_dphi = atan2(state_dy,state_dx);

      // relative rotation angle
      double rotphi = (dphi-state_dphi);
        
      // rotated momentum at the solution
      double rotpx = cos(rotphi)*state->get_px() - sin(rotphi)*state->get_py();
      double rotpy = sin(rotphi)*state->get_px() + cos(rotphi)*state->get_py();
    
      // assumes cylinder is centered at 0,0
      double dot = rotpx*x + rotpy*y;

      // our solution will have a momentum vector leaving the cylinder surface
      if (dot >= 0.0) {
	// find the solution for z
	double u = (dphi - cphi)/charge;
	
	// look only along the projection (not backward)
	if (u > 2.0*M_PI) {
	  u = u - 2.0*M_PI;
	} else if (u < 0.0) {
	  u = u + 2.0*M_PI;
	}

	intersection_z = b*u+state->get_z();      
      } else {
	remove_iter = iter;
      }
    }

    if (remove_iter != intersections.end()) {
      intersections.erase(remove_iter);
    }

    if (intersections.empty()) return;
    
    intersection[0] = intersections.begin()->at(0);
    intersection[1] = intersections.begin()->at(1);
    intersection[2] = intersection_z;

    return;
    
  } else {
    // no magnetic field project track as a line
    
    circle_line_intersections(0.0,0.0,radius,
     			      state->get_x(),state->get_y(),state->get_px(),state->get_py(),
     			      &intersections);

    if (intersections.empty()) return;
    
    // downselect solutions based on track direction
    // we want the solution where the state vector would point outward
    // since the track doesn't bend this will be the solution where
    // the dot product of the displacement between the solution and the cylinder center
    // and the momentum vector is positive
    std::set<std::vector<double> >::iterator remove_iter = intersections.end();
    double intersection_z = 0.0;
    for (std::set<std::vector<double> >::iterator iter = intersections.begin();
      	 iter != intersections.end();
      	 ++iter) {
       double x = iter->at(0);
       double y = iter->at(1);

      // assumes cylinder is centered at 0,0
      double dot = state->get_px()*x + state->get_py()*y;
      if (dot >= 0.0) {

	// x(u) = px*u + x1
	// y(u) = py*u + y1
	// z(u) = pz*u + z1
	
	double u = NAN;
	if (state->get_px() != 0) {
	  u = (intersection[0] - state->get_x()) / state->get_px();
	} else if (state->get_py() != 0) {
	  u = (intersection[1] - state->get_y()) / state->get_py();      
	}

	intersection_z = state->get_pz()*u+state->get_z();    
      } else {
	remove_iter = iter;
      }
    }

    if (remove_iter != intersections.end()) {
      intersections.erase(remove_iter);
    }

    if (intersections.empty()) return;

    intersection[0] = intersections.begin()->at(0);
    intersection[1] = intersections.begin()->at(1);
    intersection[2] = intersection_z;

    return;
  }
  
  return;
}

float PHG4HoughTransformTPC::kappaToPt(float kappa) {  
  return _pt_rescale * _magField / 333.6 / kappa;
}

float PHG4HoughTransformTPC::ptToKappa(float pt) {  
  return _pt_rescale * _magField / 333.6 / pt;
}

PHG4HoughTransformTPC::PHG4HoughTransformTPC(unsigned int seed_layers, unsigned int req_seed, const string &name) :
  SubsysReco(name),
  _timer(PHTimeServer::get()->insert_new("PHG4HoughTransformTPC")),
  _timer_initial_hough(PHTimeServer::get()->insert_new("PHG4HoughTransformTPC::track finding")),
  _min_pT(0.2), 
  _min_pT_init(0.2), 
  _seed_layers(seed_layers), 
  _req_seed(req_seed), 
  _reject_ghosts(true), 
  _remove_hits(true), 
  _use_cell_size(false),
  _max_cluster_error(3.0),
  _bin_scale(0.8), 
  _z_bin_scale(0.8), 
  _cut_on_dca(false), 
  _dca_cut(0.1),
  _dcaz_cut(0.2),
  _write_reco_tree(false)
{
  _magField = 1.5; // Tesla
  _use_vertex = true;
  _chi2_cut_init = 4.0;
  _chi2_cut_fast_par0 = 16.0;
  _chi2_cut_fast_par1 = 0.0;
  _chi2_cut_fast_max = FLT_MAX;
  _chi2_cut_full = 4.0;
  _ca_chi2_cut = 4.0;
  _cos_angle_cut = 0.985;

  _beta = 1;
  _lambda = 1;

  _pt_rescale = 1.0;
  
  _vote_error_scale.assign(_seed_layers, 1.0);
  _fit_error_scale.assign(_seed_layers, 1.0/sqrt(12.));

  _layer_ilayer_map.clear();
}

int PHG4HoughTransformTPC::Init(PHCompositeNode *topNode)
{
  if(_write_reco_tree == true)
  {
    _reco_tree = new TTree("reco_events", "a tree of SimpleRecoEvent");
    _recoevent = new SimpleRecoEvent();
    _reco_tree->Branch("tracks", "SimpleRecoEvent", &_recoevent);
  }


  return Fun4AllReturnCodes::EVENT_OK;
}


int PHG4HoughTransformTPC::InitRun(PHCompositeNode *topNode)
{
  int code = CreateNodes(topNode);

  if (verbosity > 0) {
    cout << "====================== PHG4HoughTransformTPC::InitRun() ======================" << endl;
    cout << " CVS Version: $Id: PHG4HoughTransformTPC.C,v 1.101 2015/04/21 23:47:09 pinkenbu Exp $" << endl;
    cout << " Magnetic field set to: " << _magField << " Tesla" << endl;
    cout << " Number of tracking layers: " << _nlayers << endl;
    for (int i=0; i<_nlayers; ++i) {
      cout << "   Tracking layer #" << i << " "
	   << "radius = " << _radii[i] << " cm, "
	   << "material = " << _material[i]
	   << endl;
      cout << "   Tracking layer #" << i << " "
	   << "vote error scale = " << _vote_error_scale[i] << ", "
	   << "fit error scale = " << _fit_error_scale[i]
	   << endl;
    }
    cout << " Required hits: " << _req_seed << endl;
    cout << " Minimum pT: " << _min_pT << endl;
    cout << " Fast fit chisq cut min(par0+par1/pt,max): min( "
	 << _chi2_cut_fast_par0 << " + " << _chi2_cut_fast_par1 << " / pt, "
	 << _chi2_cut_fast_max << " )" << endl;
    cout << " Maximum chisq (kalman fit): " << _chi2_cut_full << endl;
    cout << " Cell automaton chisq: " << _ca_chi2_cut << endl;
    cout << " Cos Angle Cut: " << _cos_angle_cut << endl;
    cout << " Ghost rejection: " << boolalpha << _reject_ghosts << noboolalpha << endl;
    cout << " Hit removal: " << boolalpha << _remove_hits << noboolalpha << endl;
    cout << " Use cell size in place of cluster sizes: " << boolalpha << _use_cell_size << noboolalpha << endl;
    if (!_use_cell_size) cout << " Max cluster size error = " << _max_cluster_error << endl;
    cout << " Maximum DCA: " << boolalpha << _cut_on_dca << noboolalpha << endl;
    if (_cut_on_dca) {
      cout << "   Maximum DCA cut: " << _dca_cut << endl;
    }
    cout << "   Maximum DCAZ cut: " << _dcaz_cut << endl;
    cout << " Phi bin scale: " << _bin_scale << endl;
    cout << " Z bin scale: " << _z_bin_scale << endl;
    cout << " Produce an initial vertex for tracking: " << boolalpha << _use_vertex << noboolalpha << endl;
    if (_use_vertex) {
      cout << "   Initial vertex minimum pT: " << _min_pT_init << endl;
      cout << "   Initial vertex maximum chisq: " << _chi2_cut_init << endl;
    }
    cout << " Momentum rescale factor: " << _pt_rescale << endl; 
    cout << "===========================================================================" << endl;
  }

  return code;
}

int PHG4HoughTransformTPC::process_event(PHCompositeNode *topNode)
{
  _timer.get()->restart();
  if(_write_reco_tree==true){ _recoevent->tracks.clear();}

  if(verbosity > 0) cout << "PHG4HoughTransformTPC::process_event -- entered" << endl;

  // moving clearing to the beginning of event or we will have
  // return bugs from early exits!
  _clusters_init.clear();
  _clusters.clear();
  _tracks.clear();
  
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------
  
  GetNodes(topNode);
  
  // Translate into Helix_Hough objects
  //-----------------------------------
  //wrap_clusters_timer.get()->restart();

  for (SvtxClusterMap::Iter iter = _g4clusters->begin();
       iter != _g4clusters->end();
       ++iter) {
    SvtxCluster* cluster = iter->second;

    //cluster->identify();
    
    float phi = atan2(cluster->get_position(1),cluster->get_position(0));
    unsigned int ilayer = _layer_ilayer_map[cluster->get_layer()];
    
    float xy_error=0.;float z_error=0.;
    if (_use_cell_size) {
      xy_error = _smear_xy_layer[ilayer] * _vote_error_scale[ilayer];
      z_error  = _smear_z_layer[ilayer] * _vote_error_scale[ilayer];
      
    }
    else {
      if( cluster->get_phi_size() <= _max_cluster_error*_smear_xy_layer[ilayer] ){xy_error = cluster->get_phi_size() * _vote_error_scale[ilayer];}
      else{xy_error = _max_cluster_error*_smear_xy_layer[ilayer] * _vote_error_scale[ilayer];}
      if(cluster->get_z_size() <= _max_cluster_error*_smear_z_layer[ilayer]){z_error  = cluster->get_z_size() * _vote_error_scale[ilayer];}
      else{z_error  = _max_cluster_error*_smear_z_layer[ilayer] * _vote_error_scale[ilayer];}
    }

    vector<SimpleHit3D>* which_vec = &_clusters;
    if (ilayer<_seed_layers) {which_vec=&_clusters_init;}

    //SimpleHit3D(float xx, float dxx, float yy, float dyy, float zz, float dzz, unsigned int ind, int lyr=-1)
    SimpleHit3D hit3d(cluster->get_x(),fabs(xy_error*sin(phi)),
		      cluster->get_y(),fabs(xy_error*cos(phi)),
		      cluster->get_z(),z_error,
		      cluster->get_id(),ilayer);

    // copy covariance over
    for (int i=0; i<3; ++i) {
      for (int j=i; j<3; ++j) {
	hit3d.set_error(i,j,cluster->get_error(i,j));
      }
    }

    which_vec->push_back(hit3d);
  }

  if (verbosity > 20) {
    cout << "-------------------------------------------------------------------" << endl;
    cout << "PHG4HoughTransformTPC::process_event has the following input clusters:" << endl;

    if (!_clusters_init.empty()) {
      for (unsigned int i = 0; i < _clusters_init.size(); ++i) {
	cout << "n init clusters = "<<_clusters_init.size() << endl;
	_clusters_init[i].print();
      }
    } else {
      for (unsigned int i = 0; i < _clusters.size(); ++i) {
	cout << "n clusters = "<<_clusters.size() << endl;
	_clusters[i].print();
      }
    }
    
    cout << "-------------------------------------------------------------------" << endl;
  }
  
  //------------------------------------
  // Perform the initial zvertex finding
  //------------------------------------

  

  // Grab some initial tracks for initial z-vertex finding
  _tracks.clear();

  _vertex.clear();
  _vertex.push_back(0.0); // x guess
  _vertex.push_back(0.0); // y guess
  _vertex.push_back(0.0); // z guess

  if(_use_vertex) {

    if(verbosity > 0) cout << "PHG4HoughTransformTPC::process_event -- initial vertex finding..." << endl;
    
    // find maxtracks tracks
    // unsigned int maxtracks = 100;
    // _tracker->findHelices(_clusters_init, _req_seed, _max_hits_init, _tracks, maxtracks);
    _tracker->findHelices(_clusters_init, _req_seed, _max_hits_init, _tracks);
    // _tracker_vertex.at(0)->findHelices(_clusters_init, _req_seed, _max_hits_init, _tracks, maxtracks);
    // _tracker->setRemoveHits(_remove_hits);

    cout<<"found "<<_tracks.size()<<" initial tracks"<<endl;
    
    if(_tracks.size() == 0){return Fun4AllReturnCodes::EVENT_OK;}
    else if(_tracks.size() == 1)
    {
      _vertex[0] = cos(_tracks[0].phi) * _tracks[0].d;
      _vertex[1] = sin(_tracks[0].phi) * _tracks[0].d;
      _vertex[2] = _tracks[0].z0;
    }
    else
    {
      vector<vector<double> > pTmap;
      for(unsigned int i=0;i<_tracks.size();++i)
      {
        if(_tracks[i].kappa == 0.0){continue;}
        double pT = kappaToPt(_tracks[i].kappa);
        pTmap.push_back(vector<double>());
        pTmap.back().push_back(pT);
        pTmap.back().push_back((double)i);
      }
      sort(pTmap.begin(), pTmap.end());
      vector<SimpleTrack3D> vtxtracks;
      unsigned int maxvtxtracks=100;
      if(_tracks.size() < maxvtxtracks)
      {
        vtxtracks = _tracks;
      }
      else
      {
        for(unsigned int i=0;i<maxvtxtracks;++i)
        {
          vtxtracks.push_back(_tracks[ (int)(pTmap[pTmap.size()-1-i][1]) ]);
        }
      }
      
      
      vector<double> zvertices(3,0.);
      vector<float> temp_vertex(3,0.);
      vector<unsigned int> vtracks(3,0);
      for(unsigned int iter = 0;iter < 3; ++iter)
      {
        temp_vertex[2] = 0.;
        
        TH1D z0_hist("z0_hist","z0_hist", 20, -10., 10.);
        for(unsigned int i=0;i<vtxtracks.size();++i)
        {
          z0_hist.Fill(vtxtracks[i].z0);
        }
        temp_vertex[2] = z0_hist.GetBinCenter( z0_hist.GetMaximumBin() );
        
        _vertexFinder.findVertex(vtxtracks, temp_vertex, 3., true);
        _vertexFinder.findVertex(vtxtracks, temp_vertex, 0.1, true);
        _vertexFinder.findVertex(vtxtracks, temp_vertex, 0.02, false);
        
        
        vector<SimpleTrack3D> ttracks;
        for(unsigned int t=0;t<vtxtracks.size();++t)
        {
          if( fabs(vtxtracks[t].z0 - temp_vertex[2]) < 0.1 ){vtracks[iter] += 1;}
          else{ttracks.push_back(vtxtracks[t]);}
        }
        vtxtracks = ttracks;
        zvertices[iter] = temp_vertex[2];
      }
      _vertex[2] = zvertices[0];
      unsigned int zbest = 0;
      for(unsigned int iter = 1;iter < 3; ++iter)
      {
        if(vtracks[iter] > vtracks[zbest])
        {
          _vertex[2] = zvertices[iter];
          zbest = iter;
        }
      }
    }
    
    
    
    
    
    
    if(verbosity > 0) cout << "PHG4HoughTransformTPC::process_event -- found initial vertex : " << _vertex[0] << " " << _vertex[1] << " " << _vertex[2] << endl;
    
    _tracks.clear();
    
    // shift the vertex to the origin
    for(unsigned int ht=0;ht<_clusters_init.size();++ht)
    {
      _clusters_init[ht].x -= _vertex[0];
      _clusters_init[ht].y -= _vertex[1];
      _clusters_init[ht].z -= _vertex[2];
    }
    for(unsigned int ht=0;ht<_clusters.size();++ht)
    {
      _clusters[ht].x -= _vertex[0];
      _clusters[ht].y -= _vertex[1];
      _clusters[ht].z -= _vertex[2];
    }
    
    

    
  }  // if(_use_vertex)
  
  //----------------------------------
  // Preform the track finding
  //----------------------------------
  _tracker->clear();
  _tracks.clear();
  _timer_initial_hough.get()->restart();
  _tracker->findHelices(_clusters_init, _min_hits_init, _max_hits_init, _tracks);

  _timer_initial_hough.get()->stop();
  
  

  if(verbosity > 0)
  {
    cout << "PHG4HoughTransformTPC::process_event -- full track finding pass found: " << _tracks.size() << " tracks" << endl;
  }    
   
  //----------------------------
  // Re-center event on detector
  //----------------------------

  if(verbosity > 0) cout << "PHG4HoughTransformTPC::process_event -- recentering event on detector..." << endl;
  vector<double> chi_squareds;
  for(unsigned int tt=0;tt<_tracks.size();tt++)
  {
    // move the hits in the track back to their original position                
    for(unsigned int hh=0;hh<_tracks[tt].hits.size();hh++)
    {
      _tracks[tt].hits[hh].x = _tracks[tt].hits[hh].x + _vertex[0];
      _tracks[tt].hits[hh].y = _tracks[tt].hits[hh].y + _vertex[1];
      _tracks[tt].hits[hh].z = _tracks[tt].hits[hh].z + _vertex[2];
      // _tracks[tt].z0 += _vertex[2];
    }
    chi_squareds.push_back(_tracker->getKalmanStates()[tt].chi2);}

  if(verbosity > 0)
  {
    cout << "PHG4HoughTransformTPC::process_event -- final track count: " << _tracks.size() << endl;
  }

  //---------------------------
  // Final vertex determination
  //---------------------------
  
  // final best guess of the primary vertex position here...
  if(verbosity > 0)
  {
    cout<< "PHG4HoughTransformTPC::process_event -- calculating final vertex" << endl;
  }
  
  // sort the tracks by pT
  vector<vector<double> > pTmap;
  for(unsigned int i=0;i<_tracks.size();++i)
  {
    double pT = kappaToPt(_tracks[i].kappa);
    pTmap.push_back(vector<double>());
    pTmap.back().push_back(pT);
    pTmap.back().push_back((double)i);
  }
  sort(pTmap.begin(), pTmap.end());
  vector<SimpleTrack3D> vtxtracks;
  vector<Matrix<float,5,5> > vtxcovariances;
  unsigned int maxvtxtracks=100;
  if(_tracks.size() < maxvtxtracks){vtxtracks = _tracks;}
  else
  {
    for(unsigned int i=0;i<maxvtxtracks;++i)
    {
      vtxtracks.push_back(_tracks[ (int)(pTmap[pTmap.size()-1-i][1]) ]);
      vtxcovariances.push_back( (_tracker->getKalmanStates())[ (int)(pTmap[pTmap.size()-1-i][1]) ].C );
    }
  }
  
  double vx = _vertex[0];
  double vy = _vertex[1];
  double vz = _vertex[2];
  
  _vertex[0] = 0.;
  _vertex[1] = 0.;
  _vertex[2] = 0.;
  
  _vertexFinder.findVertex(vtxtracks, vtxcovariances, _vertex, 0.3, false);
  _vertexFinder.findVertex(vtxtracks, vtxcovariances, _vertex, 0.1, false);
  _vertexFinder.findVertex(vtxtracks, vtxcovariances, _vertex, 0.02, false);
  _vertexFinder.findVertex(vtxtracks, vtxcovariances, _vertex, 0.005, false);
  
  _vertex[0] += vx;
  _vertex[1] += vy;
  _vertex[2] += vz;
  
  if(verbosity > 0)
  {
    cout << "PHG4HoughTransformTPC::process_event -- final vertex: " << _vertex[0] << " " << _vertex[1] << " " << _vertex[2] << endl;
  }

  //--------------------------------
  // Translate back into PHG4 objects
  //--------------------------------

  if(verbosity > 0)
  {
    cout << "PHG4HoughTransformTPC::process_event -- producing PHG4Track objects..." << endl;
  }

  SvtxVertex_v1 vertex;
  vertex.set_t0(0.0);
  for (int i=0;i<3;++i) vertex.set_position(i,_vertex[i]);
  vertex.set_chisq(0.0);
  vertex.set_ndof(0); 
  vertex.set_error(0,0,0.0);
  vertex.set_error(0,1,0.0);
  vertex.set_error(0,2,0.0);
  vertex.set_error(1,0,0.0);
  vertex.set_error(1,1,0.0);
  vertex.set_error(1,2,0.0);
  vertex.set_error(2,0,0.0);
  vertex.set_error(2,1,0.0);
  vertex.set_error(2,2,0.0);
  
  // copy out the reconstructed vertex position
  //_g4tracks->setVertex(_vertex[0],_vertex[1],_vertex[2]);
  //_g4tracks->setVertexError(0.0,0.0,0.0);
 
  // at this point we should already have an initial pt and pz guess...
  // need to translate this into the PHG4Track object...

  vector<SimpleHit3D> track_hits;
  int clusterID;
  int clusterLayer;

  for(unsigned int itrack=0; itrack<_tracks.size();itrack++)
  {
    SvtxTrack_v1 track;
    track.set_id(itrack);
    track_hits.clear();
    track_hits = _tracks.at(itrack).hits;
    
    for(unsigned int ihit = 0; ihit<track_hits.size();ihit++)
    {
      //      dEdx1=0;
      //      dEdx2=0;
      if( (track_hits.at(ihit).index) >= _g4clusters->size()){continue;}
      SvtxCluster *cluster = _g4clusters->get(track_hits.at(ihit).index);
      clusterID = cluster->get_id();
      clusterLayer = cluster->get_layer();
      if( (clusterLayer < (int)_seed_layers) && (clusterLayer >= 0) )
      {
        track.insert_cluster(clusterID);
      }
    }
    float kappa = _tracks.at(itrack).kappa;
    float d = _tracks.at(itrack).d;
    float phi = _tracks.at(itrack).phi;

    float dzdl = _tracks.at(itrack).dzdl;
    float z0 = _tracks.at(itrack).z0;
    
    // track.set_helix_phi(phi);
    // track.set_helix_kappa(kappa);
    // track.set_helix_d(d);
    // track.set_helix_z0(z0);
    // track.set_helix_dzdl(dzdl);
    
    float pT = kappaToPt(kappa);

    float x_center = cos(phi)*(d+1/kappa); // x coordinate of circle center
    float y_center = sin(phi)*(d+1/kappa); // y    "      "     "      "

    // find helicity from cross product sign
    short int helicity;
    if((track_hits[0].x-x_center)*(track_hits[track_hits.size()-1].y-y_center) -
       (track_hits[0].y-y_center)*(track_hits[track_hits.size()-1].x-x_center) > 0)
    {
      helicity = 1;
    }
    else
    { 
      helicity = -1;
    }
    float pZ = 0;
    if(dzdl != 1)
    {
      pZ = pT * dzdl / sqrt(1.0 - dzdl*dzdl);
    }
    int ndf = 2*_tracks.at(itrack).hits.size() - 5;
    track.set_chisq(chi_squareds[itrack]);
    track.set_ndf(ndf);
    track.set_px( pT*cos(phi-helicity*M_PI/2) );
    track.set_py( pT*sin(phi-helicity*M_PI/2) );
    track.set_pz( pZ );

    track.set_dca2d( d );
    track.set_dca2d_error(sqrt(_tracker->getKalmanStates()[itrack].C(1,1)));  

    if(_write_reco_tree==true)
    {
      _recoevent->tracks.push_back( SimpleRecoTrack() );
      _recoevent->tracks.back().px = pT*cos(phi-helicity*M_PI/2);
      _recoevent->tracks.back().py = pT*sin(phi-helicity*M_PI/2);
      _recoevent->tracks.back().pz = pZ;
      _recoevent->tracks.back().d = d;
      _recoevent->tracks.back().z0 = z0;
      _recoevent->tracks.back().quality = chi_squareds[itrack]/((float)ndf);
      _recoevent->tracks.back().charge = (-1*helicity);
    }
    

    if(_magField > 0)
    {
      track.set_charge( helicity );
    }
    else
    {
      track.set_charge( -1.0*helicity );
    }

    Matrix<float,6,6> euclidean_cov = Matrix<float,6,6>::Zero(6,6);
    convertHelixCovarianceToEuclideanCovariance( _magField, phi, d, kappa, z0, dzdl, _tracker->getKalmanStates()[itrack].C, euclidean_cov );
    
    for(unsigned int row=0;row<6;++row)
    {
      for(unsigned int col=0;col<6;++col)
      {
	track.set_error(row,col,euclidean_cov(row,col));
      }
    }

    track.set_x( vertex.get_x() + d*cos(phi) );
    track.set_y( vertex.get_y() + d*sin(phi) );
    track.set_z( vertex.get_z() + z0 );
    
    _g4tracks->insert(&track);
    vertex.insert_track(track.get_id());

    if (verbosity > 5) {
      cout << "track " << itrack << " quality = "
           << track.get_quality() << endl;
      cout << "px = " << track.get_px()
           << " py = " << track.get_py()
           << " pz = " << track.get_pz() << endl;
    }
  } // track loop

  SvtxVertex *vtxptr = _g4vertexes->insert(&vertex);
  if (verbosity > 5) vtxptr->identify();
  
  if(verbosity > 0)
  {
    cout << "PHG4HoughTransformTPC::process_event -- leaving process_event" << endl;
  }

  if(_write_reco_tree==true){ _reco_tree->Fill(); }

  _timer.get()->stop();
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4HoughTransformTPC::End(PHCompositeNode *topNode) {
  for (unsigned int i = 0; i < _tracker_vertex.size(); ++i) {
    delete _tracker_vertex[i];
  }
  delete _tracker;

  if(_write_reco_tree==true)
  {
    TFile recofile( "recotracks.root", "recreate" );
    recofile.cd();
    _reco_tree->Write();
    recofile.Close();
  }
  
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4HoughTransformTPC::InitializeGeometry(PHCompositeNode *topNode) {

  //---------------------------------------------------------
  // Grab Run-Dependent Detector Geometry and Configure Hough
  //---------------------------------------------------------

  bool default_geo = false;
  
  PHG4CylinderCellGeomContainer* cellgeos = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode,"CYLINDERCELLGEOM_SVTX");
  PHG4CylinderGeomContainer* laddergeos = findNode::getClass<PHG4CylinderGeomContainer>(topNode,"CYLINDERGEOM_SILICON_TRACKER");
  //PHG4CylinderCellContainer* cells = findNode::getClass<PHG4CylinderCellContainer>(topNode,"G4CELL_SVTX");
  
  if (cellgeos||laddergeos) {
    unsigned int ncelllayers = 0;
    if (cellgeos) ncelllayers += cellgeos->get_NLayers();
    unsigned int nladderlayers = 0;
    if (laddergeos) nladderlayers += laddergeos->get_NLayers();
    _nlayers = ncelllayers + nladderlayers;
    default_geo = false;
  } else {
    cerr << PHWHERE
	 << "Neither CYLINDERCELLGEOM_SVTX nor CYLINDERGEOM_SILICON_TRACKER available, reverting to a default geometry"
	 << std::endl;    
    _nlayers = 6;
    default_geo = true;
  }

  //=================================================//                                                                                                      
  //  Initializing HelixHough objects                //                                                                                                      
  //=================================================//                                                                                            

  _radii.assign(_nlayers, 0.0);
  _smear_xy_layer.assign(_nlayers, 0.0);
  _smear_z_layer.assign(_nlayers, 0.0);
  float sqrt_12 = sqrt(12.);
	
  if (default_geo) {

    // default geometry
    _radii[0] = 2.5;
    _radii[1] = 5.0;
    _radii[2] = 10.0;
    _radii[3] = 14.0;
    _radii[4] = 40.0;
    _radii[5] = 60.0;
    
    _smear_xy_layer[0] = (50.0e-4/sqrt_12);
    _smear_z_layer[0] = (425.0e-4/sqrt_12);
    _smear_xy_layer[1] = (50.0e-4/sqrt_12);
    _smear_z_layer[1] = (425.0e-4/sqrt_12);
    _smear_xy_layer[2] = (80.0e-4/sqrt_12);
    _smear_z_layer[2] = (1000.0e-4/sqrt_12);
    _smear_xy_layer[3] = (80.0e-4/sqrt_12);
    _smear_z_layer[3] = (1000.0e-4/sqrt_12);
    
    for(int il=4; il<_nlayers; ++il) {
      _smear_xy_layer[il] = (80.0e-4/sqrt_12);
      _smear_z_layer[il] = (30000.0e-4/sqrt_12);
    }

    _layer_ilayer_map.clear();
    for (int ilayer = 0; ilayer < _nlayers; ++ilayer) {
      _layer_ilayer_map.insert(make_pair(ilayer,ilayer));
    }
    
  } else {

    // Since the G4 layers don't necessarily correspond to the
    // silicon layers, and don't necessarily start from zero (argh),
    // we create our own layers numbers that are consecutive
    // starting from zero.

    // Now that we have two kinds of layers, I won't know in principle
    // which type is in what order, so I figure that out now...
    
    map<float,int> radius_layer_map;

    if (cellgeos) {
      PHG4CylinderCellGeomContainer::ConstRange layerrange = cellgeos->get_begin_end();
      for(PHG4CylinderCellGeomContainer::ConstIterator layeriter = layerrange.first;
	  layeriter != layerrange.second;
	  ++layeriter) {
	radius_layer_map.insert( make_pair(layeriter->second->get_radius(),
					   layeriter->second->get_layer()) );
      }
    }

    if (laddergeos) {
      PHG4CylinderGeomContainer::ConstRange layerrange = laddergeos->get_begin_end();
      for(PHG4CylinderGeomContainer::ConstIterator layeriter = layerrange.first;
	  layeriter != layerrange.second;
	  ++layeriter) {
	radius_layer_map.insert( make_pair(layeriter->second->get_radius(),
					   layeriter->second->get_layer()) );
      }
    }

    // now that the layer ids are sorted by radius, I can create a storage
    // index, ilayer, that is 0..N-1 and sorted by radius
    
    int ilayer = 0;
    for(map<float,int>::iterator iter = radius_layer_map.begin();
	iter != radius_layer_map.end();
	++iter) {
      _layer_ilayer_map.insert( make_pair(iter->second,ilayer) );
      ++ilayer;
    }   

    // now we extract the information from the cellgeos first
    if (cellgeos) {    
      PHG4CylinderCellGeomContainer::ConstRange begin_end = cellgeos->get_begin_end();
      PHG4CylinderCellGeomContainer::ConstIterator miter = begin_end.first;
      for( ; miter != begin_end.second; miter++) {
	PHG4CylinderCellGeom *cellgeo = miter->second;
      
	if (verbosity > 1) cellgeo->identify();

	_radii[_layer_ilayer_map[cellgeo->get_layer()]] = cellgeo->get_radius();      
	_smear_xy_layer[_layer_ilayer_map[cellgeo->get_layer()]] = cellgeo->get_radius()*cellgeo->get_phistep();
	_smear_z_layer[_layer_ilayer_map[cellgeo->get_layer()]] = cellgeo->get_zstep();     
      }
    }

    if (laddergeos) {    
      PHG4CylinderGeomContainer::ConstRange begin_end = laddergeos->get_begin_end();
      PHG4CylinderGeomContainer::ConstIterator miter = begin_end.first;
      for( ; miter != begin_end.second; miter++) {
	PHG4CylinderGeom *geo = miter->second;
	
	if (verbosity > 1) geo->identify();
	
	_radii[_layer_ilayer_map[geo->get_layer()]] = geo->get_radius();      
	_smear_xy_layer[_layer_ilayer_map[geo->get_layer()]] = geo->get_strip_y_spacing();
	_smear_z_layer[_layer_ilayer_map[geo->get_layer()]] = geo->get_strip_z_spacing();     
      }
    }
  }  

  // set material on each layer
  
  _material.assign(_radii.size(), 0.03);

  map<int, float>::iterator mat_it;
  for (map<int, float>::iterator iter = _user_material.begin();
       iter != _user_material.end();
       ++iter) {
    _material[_layer_ilayer_map[iter->first]] = iter->second;
  }

  float kappa_max = ptToKappa(_min_pT);

  HelixRange top_range( 0.0, 2.*M_PI,
		       -0.2, 0.2,
			0.0, kappa_max,
		       -0.9, 0.9,
		       -1.0*_dcaz_cut, 1.0*_dcaz_cut);
  if (!_use_vertex) {
    top_range.min_z0 = -10.;
    top_range.max_z0 = 10.;
  }

  top_range.min_z0 = -10.;
  top_range.max_z0 = 10.;
  
  vector<unsigned int> onezoom(5,0);
  vector<vector<unsigned int> > zoomprofile;
  zoomprofile.assign(8,onezoom);
  zoomprofile[0][0] = 16;
  zoomprofile[0][1] = 1;
  zoomprofile[0][2] = 4;
  zoomprofile[0][3] = 8;
  zoomprofile[0][4] = 4;
  
  zoomprofile[1][0] = 16;
  zoomprofile[1][1] = 1;
  zoomprofile[1][2] = 4;
  zoomprofile[1][3] = 4;
  zoomprofile[1][4] = 4;
  
  zoomprofile[2][0] = 4;
  zoomprofile[2][1] = 2;
  zoomprofile[2][2] = 3;
  zoomprofile[2][3] = 2;
  zoomprofile[2][4] = 4;
  
  for (unsigned int i = 3; i <= 7; ++i) {
    zoomprofile[i][0] = 5;
    zoomprofile[i][1] = 2;
    zoomprofile[i][2] = 5;
    zoomprofile[i][3] = 5;
    zoomprofile[i][4] = 5;
  }
    
  _tracker = new sPHENIXTracker(zoomprofile, 3, top_range, _material, _radii, _magField);
  _tracker->setIterateClustering(true);
  _tracker->setNLayers(_seed_layers);
  _tracker->requireLayers(_req_seed);
  _max_hits_init = _seed_layers*4;
  if(_seed_layers >= 10){_max_hits_init = _seed_layers*1;}
  _min_hits_init = _req_seed;
  _tracker->setClusterStartBin(1);
  // if(_seed_layers < 10){ _tracker->setClusterStartBin(1); }
  // else{ _tracker->setClusterStartBin(10); }
  _tracker->setRejectGhosts(_reject_ghosts);
  _tracker->setFastChi2Cut(_chi2_cut_fast_par0,
			   _chi2_cut_fast_par1,
			   _chi2_cut_fast_max);
  _tracker->setChi2Cut(_chi2_cut_full);
  _tracker->setChi2RemovalCut(_chi2_cut_full*0.5);
  _tracker->setCellularAutomatonChi2Cut(_ca_chi2_cut);
  _tracker->setPrintTimings(false);
  if(verbosity > 3){_tracker->setPrintTimings(true);}
  _tracker->setVerbosity(verbosity);
  _tracker->setCutOnDca(_cut_on_dca);
  _tracker->setDcaCut(_dca_cut);
  _tracker->setSmoothBack(false);
  _tracker->setBinScale(_bin_scale);
  _tracker->setZBinScale(_z_bin_scale);
  _tracker->setRemoveHits(_remove_hits);
  _tracker->setSeparateByHelicity(true);
  _tracker->setMaxHitsPairs(0);
  _tracker->setCosAngleCut(_cos_angle_cut);
  
  
  vector<vector<unsigned int> > zoomprofile_init;
  zoomprofile_init.assign(4,onezoom);
  for(unsigned int i=0;i<=1;++i)
  {
    zoomprofile_init[i][0] = 8;
    zoomprofile_init[i][1] = 1;
    zoomprofile_init[i][2] = 3;
    zoomprofile_init[i][3] = 4;
    zoomprofile_init[i][4] = 4;
  }
  for(unsigned int i=2;i<=3;++i)
  {
    zoomprofile_init[i][0] = 8;
    zoomprofile_init[i][1] = 1;
    zoomprofile_init[i][2] = 2;
    zoomprofile_init[i][3] = 2;
    zoomprofile_init[i][4] = 2;
  }
  vector<HelixRange> top_range_init;
  unsigned int nphi = 1;
  unsigned int nz0 = 5;
  double phimin = 0.;
  double phi_step = 2.0*M_PI/((double)nphi);
  float kappa_max_init = ptToKappa(_min_pT_init);
  for(unsigned int i=0;i<nphi;++i)
  {
    double z0min = -10.;
    double z0_step = 20./((double)nz0);
    for(unsigned int j=0;j<nz0;++j)
    {
      top_range_init.push_back(HelixRange(phimin, phimin+phi_step,   -0.2, 0.2,   0.0, kappa_max_init,   -0.9, 0.9,   z0min, z0min+z0_step));
      _tracker_vertex.push_back( new sPHENIXTracker(zoomprofile_init, 1, top_range_init.back(), _material, _radii, _magField) );

        (_tracker_vertex.back())->setIterateClustering(true);
        (_tracker_vertex.back())->setNLayers(_seed_layers);
        (_tracker_vertex.back())->requireLayers(_req_seed);
        _max_hits_init = _seed_layers*4;
        if(_seed_layers >= 10){_max_hits_init = _seed_layers*2;}
        _min_hits_init = _req_seed;
        if(_seed_layers < 10){ (_tracker_vertex.back())->setClusterStartBin(1); }
        else{ (_tracker_vertex.back())->setClusterStartBin(10); }
        (_tracker_vertex.back())->setRejectGhosts(_reject_ghosts);
        (_tracker_vertex.back())->setFastChi2Cut(_chi2_cut_fast_par0,
               _chi2_cut_fast_par1,
               _chi2_cut_fast_max);
        (_tracker_vertex.back())->setChi2Cut(_chi2_cut_full);
        (_tracker_vertex.back())->setChi2RemovalCut(_chi2_cut_full*0.5);
        (_tracker_vertex.back())->setCellularAutomatonChi2Cut(_ca_chi2_cut);
        (_tracker_vertex.back())->setPrintTimings(false);
        if(verbosity > 3){(_tracker_vertex.back())->setPrintTimings(true);}
        (_tracker_vertex.back())->setVerbosity(verbosity);
        (_tracker_vertex.back())->setCutOnDca(_cut_on_dca);
        (_tracker_vertex.back())->setDcaCut(_dca_cut);
        (_tracker_vertex.back())->setSmoothBack(false);
        (_tracker_vertex.back())->setBinScale(_bin_scale);
        (_tracker_vertex.back())->setZBinScale(_z_bin_scale);
        (_tracker_vertex.back())->setRemoveHits(_remove_hits);
        (_tracker_vertex.back())->setSeparateByHelicity(true);
        (_tracker_vertex.back())->setMaxHitsPairs(0);
        (_tracker_vertex.back())->setCosAngleCut(_cos_angle_cut);
      z0min += z0_step;
    }
    phimin += phi_step;
  }
  
  for(unsigned int ilayer = 0; ilayer < _fit_error_scale.size(); ++ilayer) {
    float scale1 = _fit_error_scale[ilayer];
    float scale2 = _vote_error_scale[ilayer];
    float scale = scale1/scale2;
    _tracker->setHitErrorScale(ilayer, scale);
    for(unsigned int j = 0; j < _tracker_vertex.size(); ++j) {
      _tracker_vertex[j]->setHitErrorScale(ilayer, scale);
    }
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4HoughTransformTPC::set_material(int layer, float value)
{
  _user_material[layer] = value;
}

int PHG4HoughTransformTPC::CreateNodes(PHCompositeNode *topNode)
{
  // create nodes...
  PHNodeIterator iter(topNode);
  
  PHCompositeNode *dstNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if(!dstNode)
  {
    cerr << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
      
  // Create the SVTX node
  PHCompositeNode* tb_node = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "SVTX"));
  if (!tb_node) 
  {
    tb_node = new PHCompositeNode("SVTX");
    dstNode->addNode(tb_node);
    if (verbosity>0) cout << "SVTX node added" << endl;
  }
 	
  _g4tracks = new SvtxTrackMap_v1;
  PHIODataNode<PHObject>* tracks_node = new PHIODataNode<PHObject>(_g4tracks,"SvtxTrackMap","PHObject");
  tb_node->addNode(tracks_node);
  if (verbosity>0) cout << "Svtx/SvtxTrackMap node added" << endl;

  _g4vertexes = new SvtxVertexMap_v1;
  PHIODataNode<PHObject>* vertexes_node = new PHIODataNode<PHObject>(_g4vertexes,"SvtxVertexMap","PHObject");
  tb_node->addNode(vertexes_node);
  if (verbosity>0) cout << "Svtx/SvtxVertexMap node added" << endl;
  
  //PHG4CylinderGeomContainer* geoms = getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_SVTX");
  //if(!geoms) 
  //{
  //  cerr << PHWHERE << " ERROR: Can't find CYLINDERGEOM_SVTX Node." << endl;
  //  return Fun4AllReturnCodes::ABORTEVENT;
  //}

  return InitializeGeometry(topNode);
}

int PHG4HoughTransformTPC::GetNodes(PHCompositeNode *topNode)
{
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------
  
  _ghitlist = getClass<PHG4HitContainer>(topNode,"G4HIT_SVTX");
  if(!_ghitlist) 
  {
    cerr << PHWHERE << " ERROR: Can't find node PHG4HitContainer" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  
  _g4clusters = getClass<SvtxClusterMap>(topNode,"SvtxClusterMap");
  if(!_g4clusters) 
  {
    cerr << PHWHERE << " ERROR: Can't find node SvtxClusterMap" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Pull the reconstructed track information off the node tree...
  _g4tracks = getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if(!_g4tracks) 
  {
    cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Pull the reconstructed track information off the node tree...
  _g4vertexes = getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if(!_g4vertexes) 
  {
    cerr << PHWHERE << " ERROR: Can't find SvtxVertexMap." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}


bool PHG4HoughTransformTPC::circle_line_intersections(double x0, double y0, double r0,
						      double x1, double y1, double vx1, double vy1,
						      std::set<std::vector<double> >* points) {
  // P0: center of rotation
  // P1: point on line
  // P2: second point on line
  // P3: intersections

  // dr: distance between P1 & P2
  // delta: discriminant on number of solutions
  
  points->clear();
  
  double x2 = x1 + vx1;
  double y2 = y1 + vy1;

  double dr = sqrt(pow(vx1,2)+pow(vy1,2));
  double det = x1*y2-x2*y1;

  double delta = pow(r0,2)*pow(dr,2) - pow(det,2);
  if (delta < 0) return false;

  double sgn_vy1 = 1.0;
  if (vy1 < 0.0) sgn_vy1 = -1.0;
  
  double x3 = (det*vy1 + sgn_vy1*vx1*sqrt(pow(r0,2)*pow(dr,2)-pow(det,2))) / pow(dr,2);
  double y3 = (-1.0*det*vx1 + fabs(vy1)*sqrt(pow(r0,2)*pow(dr,2)-pow(det,2))) / pow(dr,2);

  std::vector<double> p3;
  p3.push_back(x3);
  p3.push_back(y3);
  points->insert(p3);

  x3 = (det*vy1 - sgn_vy1*vx1*sqrt(pow(r0,2)*pow(dr,2)-pow(det,2))) / pow(dr,2);
  y3 = (-1.0*det*vx1 - fabs(vy1)*sqrt(pow(r0,2)*pow(dr,2)-pow(det,2))) / pow(dr,2);

  p3[0] = x3;
  p3[1] = y3;
  points->insert(p3);

  return true;
}
  
bool PHG4HoughTransformTPC::circle_circle_intersections(double x0, double y0, double r0,
							double x1, double y1, double r1,
							std::set<std::vector<double> >* points) {
  // P0: center of rotation on first circle
  // P1: center of rotation of second circle
  // P2: point between P0 & P1 on radical line
  // P3: intersection points

  // d: distance between P0 and P1
  // a: distance between P0 and P2
  // h: distance from P2 and P3s

  points->clear();
  
  // distance between two circle centers
  double d = sqrt(pow(x0-x1,2)+pow(y0-y1,2));

  // handle error conditions
  if (fabs(r0+r1) < d) return false; // no solution
  if (fabs(r0-r1) > d) return false; // no solution
  if (d == 0 && r0 == r1) return false; // infinite solutions

  // compute distances to intersection points
  double a = (pow(r0,2)-pow(r1,2)+pow(d,2)) / (2*d);
  double h = sqrt(pow(r0,2)-pow(a,2));

  // compute P2
  double x2 = x0 + a * (x1 - x0) / d;
  double y2 = y0 + a * (y1 - y0) / d;
  
  // compute intersection, p3
  double x3 = x2 + h * (y1 - y0) / d;
  double y3 = y2 - h * (x1 - x0) / d;

  std::vector<double> p3;
  p3.push_back(x3);
  p3.push_back(y3);
  points->insert(p3);
  
  // second intersection (if different than first)
  x3 = x2 - h * (y1 - y0) / d;
  y3 = y2 + h * (x1 - x0) / d;

  p3[0] = x3;
  p3[1] = y3;
  points->insert(p3);

  return points;  
}


