#include "PHG4HoughTransformTPC.h"

#include "SvtxTrackState.h"

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
#include <g4main/PHG4VtxPoint.h>
#include <g4bbc/BbcVertexMap.h>
#include <g4bbc/BbcVertex.h>

// PHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

// Helix Hough includes
#include <HelixHough/HelixHough.h>
#include <HelixHough/HelixResolution.h>
#include <HelixHough/HelixRange.h>
#include <HelixHough/SimpleHit3D.h>
#include <HelixHough/SimpleTrack3D.h>
#include <HelixHough/sPHENIXTrackerTPC.h>
#include <HelixHough/VertexFinder.h>

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


PHG4HoughTransformTPC::PHG4HoughTransformTPC(unsigned int seed_layers, unsigned int req_seed, const string &name) :
  SubsysReco(name),
  _nlayers(0),
  _seed_layers(seed_layers), 
  _req_seed(req_seed), 
  _radii(),
  _material(),
  _user_material(),
  _magField(1.5), ///< in Tesla
  _use_bbc(true),
  _reject_ghosts(true), 
  _remove_hits(true), 
  _use_cell_size(false),
  _max_cluster_error(3.0),
  _min_pt(0.2),
  _min_z0(-14.0),
  _max_z0(+14.0),
  _max_d(1.0),
  _cut_on_dca(false),
  _dcaxy_cut(0.1),
  _dcaz_cut(0.2),
  _chi2_cut_fast_par0(16.0),
  _chi2_cut_fast_par1(0.0),
  _chi2_cut_fast_max(FLT_MAX),
  _chi2_cut_init(3.0),
  _chi2_cut_full(2.0),
  _ca_chi2_cut(4.0),
  _cos_angle_cut(0.985),
  _bin_scale(0.8), 
  _z_bin_scale(0.8), 
  _min_combo_hits(req_seed),
  _max_combo_hits(3/2*seed_layers),
  _min_vtx_hits(3),
  _max_vtx_hits(14),
  _pt_rescale(1.0),
  _fit_error_scale(_seed_layers,1.0/sqrt(12.0)),
  _vote_error_scale(_seed_layers,1.0),
  _layer_ilayer_map(),
  _clusters_init(),
  _clusters(),
  _tracks(),
  _track_errors(),
  _track_covars(),
  _vertex(),
  _tracker(NULL),
  _tracker_etap_seed(NULL),
  _tracker_etam_seed(NULL),
  _vertexFinder(),
  _bbc_vertexes(NULL),
  _g4clusters(NULL),
  _g4tracks(NULL),
  _g4vertexes(NULL),
  _timer(PHTimeServer::get()->insert_new("PHG4HoughTransformTPC")),
  _timer_initial_hough(PHTimeServer::get()->insert_new("PHG4HoughTransformTPC::track finding")),
  _write_reco_tree(false),
  _reco_tree(NULL),
  _recoevent(NULL)
{  
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

  if (Verbosity() > 0) {
    cout << "====================== PHG4HoughTransformTPC::InitRun() ======================" << endl;
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
    cout << " Minimum pt: " << _min_pt << endl;
    cout << " Fast fit chisq cut min(par0+par1/pt,max): min( "
	 << _chi2_cut_fast_par0 << " + " << _chi2_cut_fast_par1 << " / pt, "
	 << _chi2_cut_fast_max << " )" << endl;
    cout << " Initial vertex maximum chisq: " << _chi2_cut_init << endl;
    cout << " Maximum chisq (kalman fit): " << _chi2_cut_full << endl;

    cout << " Cell automaton chisq: " << _ca_chi2_cut << endl;
    cout << " Cos Angle Cut: " << _cos_angle_cut << endl;
    cout << " Ghost rejection: " << boolalpha << _reject_ghosts << noboolalpha << endl;
    cout << " Hit removal: " << boolalpha << _remove_hits << noboolalpha << endl;
    cout << " Use cell size in place of cluster sizes: " << boolalpha << _use_cell_size << noboolalpha << endl;
    if (!_use_cell_size) cout << " Max cluster size error = " << _max_cluster_error << endl;
    cout << " Maximum DCA: " << boolalpha << _cut_on_dca << noboolalpha << endl;
    if (_cut_on_dca) {
      cout << "   Maximum DCA cut: " << _dcaxy_cut << endl;
    }
    cout << "   Maximum DCAZ cut: " << _dcaz_cut << endl;
    cout << " Phi bin scale: " << _bin_scale << endl;
    cout << " Z bin scale: " << _z_bin_scale << endl;
    cout << " Momentum rescale factor: " << _pt_rescale << endl; 
    cout << "===========================================================================" << endl;
  }

  return code;
}

int PHG4HoughTransformTPC::process_event(PHCompositeNode *topNode)
{
  if(Verbosity()>1000) std::cout << "PHG4HoughTransformTPC::Process_Event" << std::endl;
  _timer.get()->restart();
  if(_write_reco_tree==true){ _recoevent->tracks.clear();}

  if(Verbosity() > 0) cout << "PHG4HoughTransformTPC::process_event -- entered" << endl;

  _clusters_init.clear();
  _clusters.clear();
  _tracks.clear();
  _track_errors.clear();
  _vertex.clear();
  _vertex.assign(3,0.0); 
 
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------
  
  GetNodes(topNode);
  
  //----------------------------------
  // Translate into Helix_Hough objects
  //-----------------------------------

  int code = translate_input();
  if (code != Fun4AllReturnCodes::EVENT_OK) return code;  

  //------------------------------------
  // Find an initial z-vertex position
  //------------------------------------

  if (Verbosity()>0 || _use_bbc){
    fast_vertex_from_bbc();
    if(!_use_bbc)  _vertex[2]=0.0;
  }

  code = initial_zvertex_finding();
  if (code != Fun4AllReturnCodes::EVENT_OK) return code;
  
  //---------------------------------------------------------
  // Preform the track finding and final vertex determination
  //---------------------------------------------------------

  code = full_tracking_and_vertexing();
  if (code != Fun4AllReturnCodes::EVENT_OK) return code;

  //--------------------------------
  // Translate back into SVTX objects
  //--------------------------------

  code = export_output(); 
  if (code != Fun4AllReturnCodes::EVENT_OK) return code;

  if(_write_reco_tree==true) {
    _reco_tree->Fill(); 
  }

  _timer.get()->stop();
  if(Verbosity()>1000) std::cout << "PHG4HoughTransformTPC::Process_Event DONE" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4HoughTransformTPC::End(PHCompositeNode *topNode) {
  
  delete _tracker_etap_seed; _tracker_etap_seed = NULL;
  delete _tracker_etam_seed; _tracker_etam_seed = NULL;
  delete _tracker; _tracker = NULL;

  if(_write_reco_tree==true)
  {
    TFile recofile( "recotracks.root", "recreate" );
    recofile.cd();
    _reco_tree->Write();
    recofile.Close();
  }
  
  
  return Fun4AllReturnCodes::EVENT_OK;
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
  PHNodeIterator iter_dst(dstNode);

  // Create the SVTX node
  PHCompositeNode* tb_node = dynamic_cast<PHCompositeNode*>(iter_dst.findFirst("PHCompositeNode", "SVTX"));
  if (!tb_node)
  {
  tb_node = new PHCompositeNode("SVTX");
  dstNode->addNode(tb_node);
  if (Verbosity()>0) cout << "SVTX node added" << endl;
  }

  _g4tracks = new SvtxTrackMap_v1;
  PHIODataNode<PHObject>* tracks_node = new PHIODataNode<PHObject>(_g4tracks,"SvtxTrackMap","PHObject");
  tb_node->addNode(tracks_node);
  if (Verbosity()>0) cout << "Svtx/SvtxTrackMap node added" << endl;
  
  _g4vertexes = new SvtxVertexMap_v1;
  PHIODataNode<PHObject>* vertexes_node = new PHIODataNode<PHObject>(_g4vertexes,"SvtxVertexMap","PHObject");
  tb_node->addNode(vertexes_node);
  if (Verbosity()>0) cout << "Svtx/SvtxVertexMap node added" << endl;

  return InitializeGeometry(topNode);
}

int PHG4HoughTransformTPC::InitializeGeometry(PHCompositeNode *topNode) {

  //---------------------------------------------------------
  // Grab Run-Dependent Detector Geometry and Configure Hough
  //---------------------------------------------------------

  bool default_geo = false;
  
  PHG4CylinderCellGeomContainer* cellgeos = 
      findNode::getClass<PHG4CylinderCellGeomContainer>(
          topNode,"CYLINDERCELLGEOM_SVTX");
  PHG4CylinderGeomContainer* laddergeos = 
      findNode::getClass<PHG4CylinderGeomContainer>(
          topNode,"CYLINDERGEOM_SILICON_TRACKER");
  PHG4CylinderGeomContainer* mapsladdergeos = findNode::getClass<PHG4CylinderGeomContainer>(
	  topNode,"CYLINDERGEOM_MAPS");

 if (cellgeos||laddergeos||mapsladdergeos) {
    unsigned int ncelllayers = 0;
    if (cellgeos) ncelllayers += cellgeos->get_NLayers();
    unsigned int nladderlayers = 0;
    if (laddergeos) nladderlayers += laddergeos->get_NLayers();
    unsigned int nmapsladderlayers = 0;
    if (mapsladdergeos) nmapsladderlayers += mapsladdergeos->get_NLayers();
    _nlayers = ncelllayers + nladderlayers+nmapsladderlayers;
    default_geo = false;
  } else {
    cerr << PHWHERE
	 << "None of CYLINDERCELLGEOM_SVTX or CYLINDERGEOM_MAPS or CYLINDERGEOM_SILICON_TRACKER available, reverting to a default geometry"
	 << std::endl;    
    _nlayers = 6;
    default_geo = true;
 }

  //=================================================
  //  Initializing HelixHough objects
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

    if (mapsladdergeos) {
      PHG4CylinderGeomContainer::ConstRange layerrange = mapsladdergeos->get_begin_end();
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
      
	if (Verbosity() > 1) cellgeo->identify();

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
	
	if (Verbosity() > 1) geo->identify();
	
	_radii[_layer_ilayer_map[geo->get_layer()]] = geo->get_radius();      
	_smear_xy_layer[_layer_ilayer_map[geo->get_layer()]] = geo->get_strip_y_spacing();
	_smear_z_layer[_layer_ilayer_map[geo->get_layer()]] = geo->get_strip_z_spacing();     
      }
    }

  if (mapsladdergeos) {    
      PHG4CylinderGeomContainer::ConstRange begin_end = mapsladdergeos->get_begin_end();
      PHG4CylinderGeomContainer::ConstIterator miter = begin_end.first;
      for( ; miter != begin_end.second; miter++) {
	PHG4CylinderGeom *geo = miter->second;
	
	if (Verbosity() > 1) geo->identify();
	
	_radii[_layer_ilayer_map[geo->get_layer()]] = geo->get_radius();      
	_smear_xy_layer[_layer_ilayer_map[geo->get_layer()]] = geo->get_pixel_x();
	_smear_z_layer[_layer_ilayer_map[geo->get_layer()]] = geo->get_pixel_z();     
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

  // initialize the pattern recognition tools
  setup_seed_tracker_objects();
  setup_tracker_object();

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4HoughTransformTPC::setup_seed_tracker_objects() {

  float kappa_max = ptToKappa(_min_pt);

  std::vector<unsigned int> onezoom(5,0);
  std::vector<vector<unsigned int> > zoomprofile;
  zoomprofile.assign(5,onezoom);
  zoomprofile[0][0] = 16;
  zoomprofile[0][1] = 1;
  zoomprofile[0][2] = 8;
  zoomprofile[0][3] = 8;
  zoomprofile[0][4] = 2;

  zoomprofile[1][0] = 8;
  zoomprofile[1][1] = 1;
  zoomprofile[1][2] = 4;
  zoomprofile[1][3] = 4;
  zoomprofile[1][4] = 2;

  zoomprofile[2][0] = 4;
  zoomprofile[2][1] = 3;
  zoomprofile[2][2] = 2;
  zoomprofile[2][3] = 2;
  zoomprofile[2][4] = 3;

  for (unsigned int i = 3; i <= 4; ++i) {
    zoomprofile[i][0] = 3;
    zoomprofile[i][1] = 3;
    zoomprofile[i][2] = 3;
    zoomprofile[i][3] = 3;
    zoomprofile[i][4] = 3;
  }


  HelixRange pos_range( 0.0, 2.*M_PI,      // center of rotation azimuthal angles
                        -0.1, 0.1,   // 2d dca range
                        0.0, kappa_max,    // curvature range
                        0.0, 0.9,          // dzdl range
                        _min_z0, _max_z0); // dca_z range

  _tracker_etap_seed = new sPHENIXTrackerTPC(zoomprofile, 1, pos_range, _material, _radii, _magField);
  _tracker_etap_seed->setEndLayer(3); // layers numbered from 0 towards outside
  _tracker_etap_seed->requireLayers(4);
  _tracker_etap_seed->setClusterStartBin(1);
  _tracker_etap_seed->setRejectGhosts(_reject_ghosts);
  _tracker_etap_seed->setFastChi2Cut(_chi2_cut_fast_par0,
                                     _chi2_cut_fast_par1,
                                     _chi2_cut_fast_max);
  _tracker_etap_seed->setChi2Cut(_chi2_cut_init);
  _tracker_etap_seed->setChi2RemovalCut(_chi2_cut_init*0.5);
  _tracker_etap_seed->setCellularAutomatonChi2Cut(_ca_chi2_cut);
  _tracker_etap_seed->setPrintTimings(false);
  _tracker_etap_seed->setCutOnDca(false);
  _tracker_etap_seed->setSmoothBack(true);
  _tracker_etap_seed->setBinScale(_bin_scale);
  _tracker_etap_seed->setZBinScale(_z_bin_scale);
  _tracker_etap_seed->setRemoveHits(_remove_hits);
  _tracker_etap_seed->setSeparateByHelicity(true);
  _tracker_etap_seed->setMaxHitsPairs(0);
  _tracker_etap_seed->setCosAngleCut(_cos_angle_cut);
  _tracker_etap_seed->setRequirePixels(true);

  for(unsigned int ilayer = 0; ilayer < _fit_error_scale.size(); ++ilayer) {
    float scale1 = _fit_error_scale[ilayer];
    float scale2 = _vote_error_scale[ilayer];
    float scale = scale1/scale2;
    _tracker_etap_seed->setHitErrorScale(ilayer, scale);
  }

  HelixRange neg_range( 0.0, 2.*M_PI,      // center of rotation azimuthal angles
                        -0.1, 0.1,   // 2d dca range
                        0.0, kappa_max,    // curvature range
                        -0.9, 0.0,         // dzdl range
                        _min_z0, _max_z0); // dca_z range

  _tracker_etam_seed = new sPHENIXTrackerTPC(zoomprofile, 1, neg_range, _material, _radii, _magField);
  _tracker_etam_seed->setEndLayer(3);
  _tracker_etam_seed->requireLayers(4);
  _tracker_etam_seed->setClusterStartBin(1);
  _tracker_etam_seed->setRejectGhosts(_reject_ghosts);
  _tracker_etam_seed->setFastChi2Cut(_chi2_cut_fast_par0,
                                     _chi2_cut_fast_par1,
                                     _chi2_cut_fast_max);
  _tracker_etam_seed->setChi2Cut(_chi2_cut_init);
  _tracker_etam_seed->setChi2RemovalCut(_chi2_cut_init*0.5);
  _tracker_etam_seed->setCellularAutomatonChi2Cut(_ca_chi2_cut);
  _tracker_etam_seed->setPrintTimings(false);
  _tracker_etam_seed->setCutOnDca(false);
  _tracker_etam_seed->setSmoothBack(true);
  _tracker_etam_seed->setBinScale(_bin_scale);
  _tracker_etam_seed->setZBinScale(_z_bin_scale);
  _tracker_etam_seed->setRemoveHits(_remove_hits);
  _tracker_etam_seed->setSeparateByHelicity(true);
  _tracker_etam_seed->setMaxHitsPairs(0);
  _tracker_etam_seed->setCosAngleCut(_cos_angle_cut);
  _tracker_etam_seed->setRequirePixels(true);

  for(unsigned int ilayer = 0; ilayer < _fit_error_scale.size(); ++ilayer) {
    float scale1 = _fit_error_scale[ilayer];
    float scale2 = _vote_error_scale[ilayer];
    float scale = scale1/scale2;
    _tracker_etam_seed->setHitErrorScale(ilayer, scale);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}


int PHG4HoughTransformTPC::setup_tracker_object() {

  float kappa_max = ptToKappa(_min_pt);

  HelixRange top_range( 0.0, 2.*M_PI,             // center of rotation azimuthal angles
                        -0.1, 0.1,  // 2d dca range
                        0.0, kappa_max,           // curvature range
                        -0.9, 0.9,                // dzdl range
                        -0.1, 0.1);   // dca_z range

    vector<unsigned int> onezoom(5,0);
    vector<vector<unsigned int> > zoomprofile;
    zoomprofile.assign(5,onezoom);
    zoomprofile[0][0] = 8;
    zoomprofile[0][1] = 1;
    zoomprofile[0][2] = 8;
    zoomprofile[0][3] = 8;
    zoomprofile[0][4] = 1;
    
    zoomprofile[1][0] = 8;
    zoomprofile[1][1] = 1;
    zoomprofile[1][2] = 8;
    zoomprofile[1][3] = 8;
    zoomprofile[1][4] = 1;
    
    zoomprofile[2][0] = 8;
    zoomprofile[2][1] = 1;
    zoomprofile[2][2] = 8;
    zoomprofile[2][3] = 8;
    zoomprofile[2][4] = 1;

    zoomprofile[3][0] = 8;
    zoomprofile[3][1] = 3;
    zoomprofile[3][2] = 8;
    zoomprofile[3][3] = 2;
    zoomprofile[3][4] = 3;

    zoomprofile[4][0] = 8;
    zoomprofile[4][1] = 3;
    zoomprofile[4][2] = 8;
    zoomprofile[4][3] = 2;
    zoomprofile[4][4] = 3;

  _tracker = new sPHENIXTrackerTPC(zoomprofile, 3, top_range, _material, _radii, _magField);
  _tracker->setIterateClustering(true); 
  _tracker->setNLayers(_seed_layers);
  _tracker->requireLayers(_req_seed);
  if(_seed_layers >= 10){_max_combo_hits = _seed_layers;}
  _tracker->setClusterStartBin(3);
  _tracker->setRejectGhosts(_reject_ghosts);
  _tracker->setFastChi2Cut(_chi2_cut_fast_par0,
                           _chi2_cut_fast_par1,
                           _chi2_cut_fast_max);
  _tracker->setChi2Cut(_chi2_cut_full);
  _tracker->setChi2RemovalCut(_chi2_cut_full*0.5);
  _tracker->setCellularAutomatonChi2Cut(_ca_chi2_cut);
  _tracker->setPrintTimings(false);
  _tracker->setCutOnDca(_cut_on_dca);
  _tracker->setDcaCut(_dcaxy_cut);
  _tracker->setSmoothBack(true);
  _tracker->setBinScale(_bin_scale);
  _tracker->setZBinScale(_z_bin_scale);
  _tracker->setRemoveHits(_remove_hits);
  _tracker->setSeparateByHelicity(true);
  _tracker->setMaxHitsPairs(0);
  _tracker->setCosAngleCut(_cos_angle_cut);
  _tracker->setRequirePixels(true);

  for(unsigned int ilayer = 0; ilayer < _fit_error_scale.size(); ++ilayer) {
    float scale1 = _fit_error_scale[ilayer];
    float scale2 = _vote_error_scale[ilayer];
    float scale = scale1/scale2;
    _tracker->setHitErrorScale(ilayer, scale);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4HoughTransformTPC::GetNodes(PHCompositeNode *topNode)
{
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  _bbc_vertexes = getClass<BbcVertexMap>(topNode, "BbcVertexMap");
  
  
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

int PHG4HoughTransformTPC::translate_input(){

  for (SvtxClusterMap::Iter iter = _g4clusters->begin();
       iter != _g4clusters->end();
       ++iter) {
    SvtxCluster* cluster = iter->second;

//    float phi = atan2(cluster->get_position(1),cluster->get_position(0));
    unsigned int ilayer = _layer_ilayer_map[cluster->get_layer()];

//    float xy_error=0.;float z_error=0.;
//    if (_use_cell_size) {
//      xy_error = _smear_xy_layer[ilayer] * _vote_error_scale[ilayer];
//      z_error  = _smear_z_layer[ilayer] * _vote_error_scale[ilayer];

//    }
//    else {
//      if( cluster->get_phi_size() <= _max_cluster_error*_smear_xy_layer[ilayer] ){xy_error = cluster->get_phi_size() * _vote_error_scale[ilayer];}
//      else{xy_error = _max_cluster_error*_smear_xy_layer[ilayer] * _vote_error_scale[ilayer];}
//      if(cluster->get_z_size() <= _max_cluster_error*_smear_z_layer[ilayer]){z_error  = cluster->get_z_size() * _vote_error_scale[ilayer];}
//      else{z_error  = _max_cluster_error*_smear_z_layer[ilayer] * _vote_error_scale[ilayer];}
//    }

    // cout<<"layer : "<<ilayer<<" , xy_error : "<<xy_error<<" , z_error : "<<z_error<<endl;

//    if(ilayer < 3)
//    {
//      xy_error = 0.002;
//      z_error = 0.002;
//    }
//    else
//    {
//      xy_error = 0.011;
//      z_error = 0.03;
//    }

    vector<SimpleHit3D>* which_vec = &_clusters;
    if (ilayer<_seed_layers) {which_vec=&_clusters_init;}

    SimpleHit3D hit3d;

    hit3d.set_id(cluster->get_id());
    hit3d.set_layer(ilayer);

    hit3d.set_x(cluster->get_x());
    hit3d.set_y(cluster->get_y());
    hit3d.set_z(cluster->get_z());

//    hit3d.set_ex(fabs(xy_error*sin(phi)));
//    hit3d.set_ey(fabs(xy_error*cos(phi)));
//    hit3d.set_ez(z_error);

    // copy covariance over
    for (int i=0; i<3; ++i) {
      for (int j=i; j<3; ++j) {
        hit3d.set_error(i,j,cluster->get_error(i,j));
        hit3d.set_size(i,j,cluster->get_size(i,j));
      }
    }

    which_vec->push_back(hit3d);
  }

  if (Verbosity() > 20) {
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

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4HoughTransformTPC::initial_zvertex_finding() {

  // fast vertex guessing uses two tracker objects
  // one looks for postive going eta tracks, the other for negative going tracks
  // it asks for just a handful of each, but searches
  // over the broadest possible phase space for the vertex origin

  // limit each window to no more than N tracks
  unsigned int maxtracks = 250;

  _tracks.clear();
  _track_errors.clear();
  _track_covars.clear();

  // loop over initial tracking windows
  std::vector<SimpleTrack3D> newtracks;

  if(Verbosity()) std::cout <<" positive eta seeded track finding.. "<< std::endl;
  _tracker_etap_seed->clear();
  _tracker_etap_seed->findHelices(_clusters_init, _min_vtx_hits, _max_vtx_hits,
                                  newtracks, maxtracks);

  if (Verbosity()) cout << " seed track finding count: " << newtracks.size() << endl;
  for (unsigned int t = 0; t < newtracks.size(); ++t) {
    if(newtracks[t].z0 < _min_z0 || newtracks[t].z0 > _max_z0 
				|| (newtracks[t].z0)!=newtracks[t].z0){
      if (Verbosity()) cout << " z0 out of bound: " << newtracks[t].z0 << endl;
      continue;
    }
    if (Verbosity()) cout << " z0: " << newtracks[t].z0 << endl;
    _tracks.push_back(newtracks[t]);
    _track_covars.push_back( (_tracker_etap_seed->getKalmanStates())[t].C );
  }

  _tracker_etap_seed->clear();
  newtracks.clear();

  if(Verbosity()) std::cout <<" negative eta seeded track finding.. "<< std::endl;
  _tracker_etam_seed->clear();
  _tracker_etam_seed->findHelices(_clusters_init, _min_vtx_hits, _max_vtx_hits,
                                  newtracks, maxtracks);

  if (Verbosity()) cout << " seed track finding count: " << newtracks.size() << endl;
  for (unsigned int t = 0; t < newtracks.size(); ++t) {
    if(newtracks[t].z0 < _min_z0 || newtracks[t].z0 > _max_z0 
				|| (newtracks[t].z0)!=newtracks[t].z0){
      if (Verbosity()) cout << " z0 out of bound: " << newtracks[t].z0 << endl;
      continue;
    }
    if (Verbosity()) cout << " z0: " << newtracks[t].z0 << endl;
    _tracks.push_back(newtracks[t]);
    _track_covars.push_back( (_tracker_etam_seed->getKalmanStates())[t].C );
  }

  _tracker_etam_seed->clear();
  newtracks.clear();

  _vertex.clear();
  _vertex.assign(3,0.0);


  if (Verbosity()) cout << " seed track used count: " << _tracks.size() << endl;

  if (!_tracks.empty()) {

    // --- compute seed vertex from initial track finding --------------------

    double zsum = 0.0;
    double zmse = 0.0;

    for (unsigned int i = 0; i < _tracks.size(); ++i) {
	zsum += _tracks[i].z0;
    }

    _vertex[2] = zsum / _tracks.size();

    for(unsigned int i = 0; i< _tracks.size(); ++i) {	
	zmse += pow(_vertex[2]-_tracks[i].z0, 2);
    }
    zmse /= (_tracks.size()-1);
    zmse = sqrt(zmse);

    if (Verbosity() > 0) {
	cout << " seed track vertex pre-fit: "
           << _vertex[0] << " "
           << _vertex[1] << " "
           << _vertex[2] << "+-" << zmse << endl;
    }

    // start with the average position and converge from there
    _vertexFinder.findVertex( _tracks, _track_covars, _vertex, 3.00, true);
    _vertexFinder.findVertex( _tracks, _track_covars, _vertex, 0.10, true);
    _vertexFinder.findVertex( _tracks, _track_covars, _vertex, 0.02, true);
  }

  if (Verbosity() > 0) {
    cout << " seed track vertex post-fit: "
         << _vertex[0] << " " << _vertex[1] << " " << _vertex[2] << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4HoughTransformTPC::full_tracking_and_vertexing() {

  float shift_dx = -_vertex[0];
  float shift_dy = -_vertex[1];
  float shift_dz = -_vertex[2];

  // shift to initial vertex position  
  shift_coordinate_system(shift_dx,shift_dy,shift_dz);

  // reset track storage and tracker
  _tracks.clear();
  _track_errors.clear();
  _track_covars.clear();

  _tracker->clear();
  _timer_initial_hough.get()->restart();

  // final track finding
  if(Verbosity()) std::cout <<" final track finding.. "<< std::endl;
  _tracker->findHelices(_clusters_init, _min_combo_hits, _max_combo_hits, _tracks);

  if (Verbosity() > 0){
    cout << " final track count, 1st pass: " << _tracks.size() << endl;
  }


  vector<SimpleHit3D> remaining_hits;
  vector<int> used( _clusters_init.size(), 0 );
  for(unsigned int tr=0;tr<_tracks.size();++tr)
  {
    for(unsigned int h=0;h<_tracks[tr].hits.size();++h)
    {
      used[_tracks[tr].hits[h].get_id()] = 1;
    }
  }
  for(unsigned int i=0;i<_clusters_init.size();++i)
  {
    if( used[i] == 0 )
    {
      remaining_hits.push_back( _clusters_init[i] );
    }
  }
  vector<SimpleTrack3D> append_tracks;
  _tracker->findHelices(remaining_hits, _min_combo_hits, _max_combo_hits, append_tracks);
  for(unsigned int tr=0;tr<append_tracks.size();++tr)
  {
    _tracks.push_back(append_tracks[tr]);
  }
  
  _timer_initial_hough.get()->stop();

  for (unsigned int tt = 0; tt < _tracks.size(); ++tt) {
    _track_covars.push_back( (_tracker->getKalmanStates())[tt].C );
    _track_errors.push_back( _tracker->getKalmanStates()[tt].chi2 );
  }

  // we will need the tracker object below to refit the tracks... so we won't
  // reset it here
  if (Verbosity() > 0){ 
    cout << " final track count, 2nd pass: " << _tracks.size() << endl;
    cout<< "PHG4HoughTransformTPC::process_event -- calculating final vertex" << endl;
  }

  if (!_tracks.empty()) {

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
    vector<SimpleTrack3D> vtx_tracks;
    vector<Matrix<float,5,5> > vtx_covars;
    unsigned int maxtracks=100;
    if(_tracks.size() < maxtracks){vtx_tracks = _tracks;}
    else
    {
      for(unsigned int i=0;i<maxtracks;++i)
      {
        vtx_tracks.push_back(_tracks[ (int)(pTmap[pTmap.size()-1-i][1]) ]);
        vtx_covars.push_back( (_tracker->getKalmanStates())[ (int)(pTmap[pTmap.size()-1-i][1]) ].C );
      }
    }

    if (Verbosity() > 0) {
      cout << " final vertex pre-fit: "
           << _vertex[0] - shift_dx << " "
           << _vertex[1] - shift_dy << " "
           << _vertex[2] - shift_dz << endl;
    }

    _vertexFinder.findVertex( vtx_tracks, vtx_covars, _vertex, 0.300, false );
    _vertexFinder.findVertex( vtx_tracks, vtx_covars, _vertex, 0.100, false );
    _vertexFinder.findVertex( vtx_tracks, vtx_covars, _vertex, 0.020, false );
    _vertexFinder.findVertex( vtx_tracks, vtx_covars, _vertex, 0.005, false );

    if (Verbosity() > 0) {
      cout << " final vertex post-fit: "
           << _vertex[0] - shift_dx << " "
           << _vertex[1] - shift_dy << " "
           << _vertex[2] - shift_dz << endl;
    }
  }

  // shift back to global coordinates
  shift_coordinate_system(-shift_dx,-shift_dy,-shift_dz);

  if(Verbosity() > 0)
  {
    cout << "PHG4HoughTransformTPC::process_event -- final vertex: " 
	 << _vertex[0] << " " << _vertex[1] << " " << _vertex[2] << endl;
  }

/*
  // we still need to update the track fields for DCA and PCA
  // we can do that from the final vertex position
  shift_dx = -_vertex[0];
  shift_dy = -_vertex[1];
  shift_dz = -_vertex[2];

  // shift to precision final vertex
  shift_coordinate_system(shift_dx,shift_dy,shift_dz);

  // recompute track fits to fill dca and pca + error fields
  std::vector<SimpleTrack3D> refit_tracks;
  std::vector<double> refit_errors;
  std::vector<Eigen::Matrix<float,5,5> > refit_covars;

  _tracker->finalize(_tracks,refit_tracks);

  for (unsigned int tt = 0; tt < refit_tracks.size(); ++tt) {
    refit_errors.push_back( _tracker->getKalmanStates()[tt].chi2);
    refit_covars.push_back( _tracker->getKalmanStates()[tt].C );
  }

  // okay now we are done with the tracker
  _tracker->clear();

  _tracks = refit_tracks;
  _track_errors = refit_errors;
  _track_covars = refit_covars;
  // shift back to global coordinates
  shift_coordinate_system(-shift_dx,-shift_dy,-shift_dz);
*/
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4HoughTransformTPC::export_output() {

  if(Verbosity() > 0)
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
      if( (track_hits.at(ihit).get_id()) >= _g4clusters->size()){continue;}
      SvtxCluster *cluster = _g4clusters->get(track_hits.at(ihit).get_id());
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
    
    float pT = kappaToPt(kappa);

    float x_center = cos(phi)*(d+1/kappa); // x coordinate of circle center
    float y_center = sin(phi)*(d+1/kappa); // y    "      "     "      "
  
    // find helicity from cross product sign
    short int helicity;
    if ((track_hits[0].get_x()-x_center) * 
        (track_hits[track_hits.size()-1].get_y()-y_center) -
        (track_hits[0].get_y()-y_center) * 
        (track_hits[track_hits.size()-1].get_x()-x_center) > 0){
      helicity = 1;
    } else {
      helicity = -1;
    }
    float pZ = 0;
    if(dzdl != 1)
    {
      pZ = pT * dzdl / sqrt(1.0 - dzdl*dzdl);
    }
    int ndf = 2*_tracks.at(itrack).hits.size() - 5;
    track.set_chisq(_track_errors[itrack]);
    track.set_ndf(ndf);
    track.set_px( pT*cos(phi-helicity*M_PI/2) );
    track.set_py( pT*sin(phi-helicity*M_PI/2) );
    track.set_pz( pZ );

    track.set_dca2d( d );
    track.set_dca2d_error(sqrt(_track_covars[itrack](1,1)));

    if(_write_reco_tree==true) {
      _recoevent->tracks.push_back( SimpleRecoTrack() );
      _recoevent->tracks.back().px = pT*cos(phi-helicity*M_PI/2);
      _recoevent->tracks.back().py = pT*sin(phi-helicity*M_PI/2);
      _recoevent->tracks.back().pz = pZ;
      _recoevent->tracks.back().d = d;
      _recoevent->tracks.back().z0 = z0;
      _recoevent->tracks.back().quality = _track_errors[itrack]/((float)ndf);
      _recoevent->tracks.back().charge = (-1*helicity);
    }

    if(_magField > 0) {
      track.set_charge( helicity );
    } else {
      track.set_charge( -1.0*helicity );
    }

    Eigen::Matrix<float,6,6> euclidean_cov = Matrix<float,6,6>::Zero(6,6);
    convertHelixCovarianceToEuclideanCovariance( 
        _magField, phi, d, kappa, z0, dzdl, 
        _track_covars[itrack], euclidean_cov );

    for(unsigned int row=0;row<6;++row) {
      for(unsigned int col=0;col<6;++col) {
        track.set_error(row,col,euclidean_cov(row,col));
      }
    }

    track.set_x(vertex.get_x() + d * cos(phi));
    track.set_y(vertex.get_y() + d * sin(phi));
    track.set_z(vertex.get_z() + z0);

    _g4tracks->insert(&track);
    vertex.insert_track(track.get_id());

    if (Verbosity() > 5) {
      cout << "track " << itrack << " quality = " << track.get_quality()
           << endl;
      cout << "px = " << track.get_px() << " py = " << track.get_py()
           << " pz = " << track.get_pz() << endl;
    }
  }  // track loop

  SvtxVertex *vtxptr = _g4vertexes->insert(&vertex);
  if (Verbosity() > 5) vtxptr->identify();

  if (Verbosity() > 0) {
    cout << "PHG4HoughTransform::process_event -- leaving process_event"
         << endl;
  }

  _clusters.clear();
  _tracks.clear();
  _track_errors.clear();
  _track_covars.clear();
  _vertex.clear();
  _vertex.assign(3,0.0);

  return Fun4AllReturnCodes::EVENT_OK;
}

float PHG4HoughTransformTPC::kappaToPt(float kappa) {
  return _pt_rescale * _magField / 333.6 / kappa;
}

float PHG4HoughTransformTPC::ptToKappa(float pt) {
  return _pt_rescale * _magField / 333.6 / pt;
}

void PHG4HoughTransformTPC::convertHelixCovarianceToEuclideanCovariance( 
     float B, float phi, float d, float kappa, float z0, float dzdl, 
     Eigen::Matrix<float,5,5> const& input, 
     Eigen::Matrix<float,6,6>& output ){

  Eigen::Matrix<float,6,5> J = Matrix<float,6,5>::Zero(6,5);

  // phi,d,nu,z0,dzdl
  //--> 
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

void PHG4HoughTransformTPC::shift_coordinate_system(double dx,
                                                 double dy,
                                                 double dz) {

  for (unsigned int ht = 0; ht < _clusters_init.size(); ++ht) {
    _clusters_init[ht].set_x( _clusters_init[ht].get_x() + dx);
    _clusters_init[ht].set_y( _clusters_init[ht].get_y() + dy);
    _clusters_init[ht].set_z( _clusters_init[ht].get_z() + dz);
  }

  for (unsigned int tt = 0; tt < _tracks.size(); ++tt) {
    for (unsigned int hh = 0; hh < _tracks[tt].hits.size(); ++hh) {
      _tracks[tt].hits[hh].set_x( _tracks[tt].hits[hh].get_x() + dx);
      _tracks[tt].hits[hh].set_y( _tracks[tt].hits[hh].get_y() + dy);
      _tracks[tt].hits[hh].set_z( _tracks[tt].hits[hh].get_z() + dz);
    }
  }

  _vertex[0] += dx;
  _vertex[1] += dy;
  _vertex[2] += dz;

  return;
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

int PHG4HoughTransformTPC::fast_vertex_from_bbc() {
  
  // fail over to bbc vertex if no tracks were found...
  if (_bbc_vertexes) {

    BbcVertex* vertex = _bbc_vertexes->begin()->second;

    if (vertex) {
	
      _vertex[0] = 0.0;
      _vertex[1] = 0.0;
      _vertex[2] = vertex->get_z();
     
      if (Verbosity()) cout << " initial bbc vertex guess: "
			  << _vertex[0] << " "
			  << _vertex[1] << " "
			  << _vertex[2] << endl;      
    }  
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}
