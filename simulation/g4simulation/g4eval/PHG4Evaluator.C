
#include "PHG4Evaluator.h"

#include "EvalLinks.h"
#include "EvalLinksV1.h"

// PHENIX includes
#include <phool/PHCompositeNode.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/getClass.h>

// PHG4 includes
#include <g4hough/SvtxVertexMap.h>
#include <g4hough/SvtxVertex.h>
#include <g4hough/SvtxTrackMap.h>
#include <g4hough/SvtxTrack.h>
#include <g4hough/SvtxClusterMap.h>
#include <g4hough/SvtxCluster.h>
#include <g4hough/SvtxHitMap.h>
#include <g4hough/SvtxHit.h>

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4detectors/PHG4CylinderCell.h>
#include <g4detectors/PHG4CylinderCellDefs.h>

// HelixHough includes (just for helix fitting algorithms)
#include <SimpleTrack3D.h>
#include <SimpleHit3D.h>
#include <sPHENIXTracker.h>

// ROOT includes
#include <TNtuple.h>

// standard includes
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

template <class T>
inline T sqr(T x)
{
  T result;
  result = x*x;
  return result;
}

PHG4Evaluator::PHG4Evaluator(const string &name, const string &filename) :
  SubsysReco("PHG4Evaluator"),
  _timer(PHTimeServer::get()->insert_new("PHG4Evaluator")),
  _nlayers(6),
  _min_hits(6), // inclusive limit
  _max_hits(7), // exclusive limit
  _truth_info_container(0) {

  _filename = filename;
  _internal_timer.push_back(PHTimeServer::get()->insert_new("PHG4Evaluator get nodes"));
  _internal_timer.push_back(PHTimeServer::get()->insert_new("PHG4Evaluator printInputInfo()"));
  _internal_timer.push_back(PHTimeServer::get()->insert_new("PHG4Evaluator fillGtrackObjects()"));
  _internal_timer.push_back(PHTimeServer::get()->insert_new("PHG4Evaluator fillClusterToG4HitMap()"));
  _internal_timer.push_back(PHTimeServer::get()->insert_new("PHG4Evaluator fillTrackToGtrackMap()"));
  _internal_timer.push_back(PHTimeServer::get()->insert_new("PHG4Evaluator fillTrackFindingFields()"));
  _internal_timer.push_back(PHTimeServer::get()->insert_new("PHG4Evaluator fitGtrackProducedClusters()"));
  _internal_timer.push_back(PHTimeServer::get()->insert_new("PHG4Evaluator fillOutputNtuples()"));
  _internal_timer.push_back(PHTimeServer::get()->insert_new("PHG4Evaluator printOutputInfo()"));
}

int PHG4Evaluator::Init(PHCompositeNode *topNode)
{
  _ievent = 0;
  _trackingWasRun = true;

  _gtrack_list.clear();
  _g4hit_gtrack_map.clear();

  _cluster_g4hit_mmap.clear();
  _cell_g4hit_mmap.clear();
  _track_gtrack_map.clear();

  _g4hit_cluster_mmap.clear();
  _g4hit_cell_mmap.clear();
  _gtrack_track_mmap.clear();

  _cluster_allg4hits_mmap.clear();
  _allg4hits_cluster_mmap.clear();
  _track_cluster_mmap.clear();

  _gtrack_track_match.clear();
  _track_gtrack_match.clear();

  _track_purity_map.clear();

  _g4hitList.clear();
  
  _tfile = new TFile(_filename.c_str(), "RECREATE");

  _ntp_event = new TNtuple("ntp_event","event-wise ntuple",
                           "event:vx:vy:vz:"
                           "gvx:gvy:gvz:ntracks:ngtracks:"
                           "nclusters:ng4hits:"
			   "hit_occupancy_layer0:"
			   "hit_occupancy_layer1:"
			   "hit_occupancy_layer2:"
			   "hit_occupancy_layer3:"
			   "hit_occupancy_layer4:"
			   "hit_occupancy_layer5");
  
  _ntp_g4hit = new TNtuple("ntp_g4hit","g4hit-wise ntuple",
                           "event:g4hitID:x:y:z:edep:"
                           "layer:nclusters:gtrackID:gflavor:"
                           "px:py:pz:vx:vy:vz:"
                           "fpx:fpy:fpz:fx:fy:fz:last:"
			   "embed:primary:match");

  _ntp_cell = new TNtuple("ntp_cell","cell-wise ntuple",
			  "event:cellID:e:layer:"
			  "g4hitID:gedep:gx:gy:gz:"
			  "gtrackID:gflavor:"
			  "gpx:gpy:gpz:gvx:gvy:gvz:"
			  "gfpx:gfpy:gfpz:gfx:gfy:gfz:glast:"
			  "gembed:gprimary");

  _ntp_cluster = new TNtuple("ntp_cluster","cluster-wise ntuple",
                             "event:hitID:x:y:z:"
                             "e:adc:layer:size:phisize:"
                             "zsize:g4hitID:gx:"
                             "gy:gz:gtrackID:gflavor:"
                             "gpx:gpy:gpz:gvx:gvy:gvz:"
                             "gfpx:gfpy:gfpz:gfx:gfy:gfz:glast:"
                             "gembed:gprimary:nhits:purity");

  _ntp_gtrack  = new TNtuple("ntp_gtrack","gtrack-wise ntuple",
                             "event:gtrackID:flavor:ng4hit0:ng4hit1:"
                             "ng4hit2:ng4hit3:px:py:pz:"
                             "ntracks:vx:vy:vz:"
                             "fpx:fpy:fpz:fx:fy:fz:last:"
                             "embed:primary:chisq:chisqv:match:bestpurity:bestdpp");

  _ntp_track = new TNtuple("ntp_track","track-wise ntuple",
                           "event:trackID:charge:quality:chisq:"
                           "chisqv:ndf:primary:nhits:layers:"
                           "dedx1:dedx2:dca:dca2d:dca2dsigma:"
                           "px:py:pz:"
			   "presdphi:presdeta:prese3x3:prese:"
			   "cemcdphi:cemcdeta:cemce3x3:cemce:"
			   "hcalindphi:hcalindeta:hcaline3x3:hcaline:"
			   "hcaloutdphi:hcaloutdeta:hcaloute3x3:hcaloute:"
			   "gtrackID:gflavor:gpx:gpy:"
                           "gpz:gvx:gvy:gvz:"
                           "gfpx:gfpy:gfpz:gfx:gfy:gfz:glast:"
                           "gembed:gprimary:purity:pcax:pcay:pcaz:"
                           "phi:d:kappa:z0:dzdl");

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4Evaluator::InitRun(PHCompositeNode *topNode) {

  //--------------------------
  // Add Cluster to G4Hit Node
  //--------------------------

  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode 
    = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode","DST"));
  if (!dstNode) {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
    
  // Create the SVTX_EVAL node if required
  PHCompositeNode* svxNode 
    = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode","SVTX_EVAL"));
  if (!svxNode) {
    svxNode = new PHCompositeNode("SVTX_EVAL");
    dstNode->addNode(svxNode);
  }

  if (findNode::getClass<PHG4CylinderCellContainer>(topNode,"G4CELL_SVTX")) {
    EvalLinks *links = findNode::getClass<EvalLinks>(topNode,"SvtxClusterMap_G4HIT_SVTX_Links");
    if (!links) {
      links = new EvalLinksV1("SvtxClusterMap","G4HIT_SVTX","edep");
      PHIODataNode<PHObject> *linksNode =
	new PHIODataNode<PHObject>(links, "SvtxClusterMap_G4HIT_SVTX_Links", "PHObject");
      svxNode->addNode(linksNode);
    }
  }

  if (findNode::getClass<PHG4CylinderCellContainer>(topNode,"G4CELL_SILICON_TRACKER")) {
    EvalLinks *links = findNode::getClass<EvalLinks>(topNode,"SvtxClusterMap_G4HIT_SILICON_TRACKER_Links");
    if (!links) {
      links = new EvalLinksV1("SvtxClusterMap","G4HIT_SILICON_TRACKER","edep");
      PHIODataNode<PHObject> *linksNode =
	new PHIODataNode<PHObject>(links, "SvtxClusterMap_G4HIT_SILICON_TRACKER_Links", "PHObject");
      svxNode->addNode(linksNode);
    }
  }

  if (findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap")) {
    EvalLinks *links = findNode::getClass<EvalLinks>(topNode,"SvtxTrackMap_G4TruthInfo_Links");
    if (!links) {
      links = new EvalLinksV1("SvtxTrackMap","G4TruthInfo","ng4hits");
      PHIODataNode<PHObject> *linksNode =
  	new PHIODataNode<PHObject>(links, "SvtxTrackMap_G4TruthInfo_Links", "PHObject");
      svxNode->addNode(linksNode);
    }
  }

  if (findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap")) {
    EvalLinks *links = findNode::getClass<EvalLinks>(topNode,"G4TruthInfo_SvtxTrackMap_Links");
    if (!links) {
      links = new EvalLinksV1("G4TruthInfo","SvtxTrackMap","nclusters");
      PHIODataNode<PHObject> *linksNode =
  	new PHIODataNode<PHObject>(links, "G4TruthInfo_SvtxTrackMap_Links", "PHObject");
      svxNode->addNode(linksNode);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
  
int PHG4Evaluator::process_event(PHCompositeNode *topNode)
{
  _timer.get()->restart();
  _internal_timer[0].get()->restart();

  if ((verbosity > 0)&&(_ievent%100==0)) {
    cout << "PHG4Evaluator::process_event - Event = " << _ievent << endl;
  }

  //------------------------
  // clear out ancestry maps
  //------------------------
  
  _gtrack_list.clear();
  _g4hit_gtrack_map.clear();

  _cluster_g4hit_mmap.clear();
  _cell_g4hit_mmap.clear();
  _track_gtrack_map.clear();

  _g4hit_cluster_mmap.clear();
  _g4hit_cell_mmap.clear();
  _gtrack_track_mmap.clear();

  _cluster_allg4hits_mmap.clear();
  _allg4hits_cluster_mmap.clear();
  _track_cluster_mmap.clear();

  _gtrack_track_match.clear();
  _track_gtrack_match.clear();

  _track_purity_map.clear();
  
  _g4hitList.clear();

  _cluster_g4hit_svtx_links = NULL;
  _cluster_g4hit_silicon_tracker_links = NULL;

  //------------------------
  // fill the Layer Type Map
  //------------------------

  // pulls the geos and cells
  int code = fillLayerTypeMap(topNode);
  if (code != Fun4AllReturnCodes::EVENT_OK) return code;

  // pull the g4hits
  code = fillGhitList(topNode);
  if (code != Fun4AllReturnCodes::EVENT_OK) return code;

  // pull the clusters
  _clusterList = findNode::getClass<SvtxClusterMap>(topNode,"SvtxClusterMap");
  if (!_clusterList) {
    cout << PHWHERE << " WARNING: Can't find SvtxClusterMap" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // pull the hits
  _hitList = findNode::getClass<SvtxHitMap>(topNode,"SvtxHitMap");
  if (!_hitList) {
    cout << PHWHERE << " WARNING: Can't find SvtxHitMap" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  
  // pull the tracks
  _trackList = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
  if(!_trackList) 
    {
      static bool firsttrackwarn = true;
      if (firsttrackwarn) {
	cout << PHWHERE << " WARNING: Can't find SvtxTrackMap. Tracking evaluation will not be performed." << endl;
	_trackingWasRun = false;
	firsttrackwarn = false;
      }
    }

  // pull the vertexs
  _vertexList = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMap");
  if (!_vertexList) {
    static bool firstvertexwarn = true;
    if (firstvertexwarn) {
      cout << PHWHERE << " WARNING: Can't find SvtxVertexMap. Vertexing evaluation will not be performed." << endl;
      _trackingWasRun = false;
      firstvertexwarn = false;
    }
  }

  // pull the tracks
  _truth_info_container = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if(!_truth_info_container) 
    {
      cout << PHWHERE << " WARNING: Can't find PHG4TruthInfoContainer." << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  if ((!_cellGeos)&&(!_ladderGeos)) {
    std::cout << PHWHERE << "CYLINDERCELLGEOM_SVTX and CYLINDERGEOM_SILICON_TRACKER Node missing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;     
  } else {

    _nlayers = 0;
    if (_cellGeos) _nlayers += _cellGeos->get_NLayers();
    if (_ladderGeos) _nlayers += _ladderGeos->get_NLayers();
    
    _nhits_per_layer.assign(_nlayers,0);
    _nchannels_per_layer.assign(_nlayers,0);

    int ilayer = 0;
    if (_cellGeos) {
      for (ilayer = 0; ilayer < _cellGeos->get_NLayers(); ++ilayer) {
	PHG4CylinderCellGeom *geom = _cellGeos->GetLayerCellGeom(ilayer);
	_nchannels_per_layer[ilayer] = geom->get_zbins()*geom->get_phibins();
      }
    }
    if (_ladderGeos) {
      for (int jlayer = ilayer; jlayer < (int)_nlayers; ++jlayer) {
	PHG4CylinderGeom *geom = _ladderGeos->GetLayerGeom(jlayer);
	_nchannels_per_layer[jlayer]
	  = geom->get_N_strips_per_column()
	  * geom->get_N_strip_columns()
	  * geom->get_N_sensors_in_layer();
      }
    }

    //_min_hits = _nlayers*0.75; // gives us 3,4,6 for 4,6,8 layers
    _min_hits = _nlayers;
    _max_hits = _nlayers+1;
    static bool firstevent=true;
    if (verbosity > 0 && firstevent) {
      cout << "PHG4Evaluator layers: " << _nlayers
	   << " _min_hits: " << _min_hits << " _max_hits: " << _max_hits << endl;
      firstevent = false;
    }
  }

  _internal_timer[0].get()->stop();

  //-----------------------------------
  // print what is coming into the code
  //-----------------------------------
  
  _internal_timer[1].get()->restart();
  printInputInfo();
  _internal_timer[1].get()->stop();
  
  //-------------------------------
  // fill the Gtrack storage vector
  //-------------------------------
  
  _internal_timer[2].get()->restart();
  fillGtrackObjects();
  _internal_timer[2].get()->stop();
  
  //----------------------------------
  // fill cluster to ghit associations
  //----------------------------------
  
  _internal_timer[3].get()->restart();
  fillCellToG4HitMap();
  fillClusterToG4HitMap();

  fillClusterToG4HitLinks(topNode); //---new-method-----------------------------

  _internal_timer[3].get()->stop();
  
  //---------------------------------------
  // fill track to particle id associations
  //---------------------------------------
  
  _internal_timer[4].get()->restart();
  if (_trackingWasRun) fillTrackToGtrackMap();
  //if (_trackingWasRun) fillTrackToG4TruthInfoLinks(topNode);
  //if (_trackingWasRun) fillG4TruthInfoToTrackLinks(topNode);
  _internal_timer[4].get()->stop();
  
  //------------------------
  // fill a evaluator fields
  //------------------------
  
  _internal_timer[5].get()->restart();
  if(_trackingWasRun) fillTrackPurityMap();
  _internal_timer[5].get()->stop();
  _internal_timer[6].get()->restart();
  if(_trackingWasRun) fitGtrackProducedClusters();
  _internal_timer[6].get()->stop();
  
  //---------------------------
  // fill the Evaluator NTuples
  //---------------------------
  
  _internal_timer[7].get()->restart();
  fillOutputNtuples();
  _internal_timer[7].get()->stop();
  
  //--------------------------------------------------
  // Print out the ancestry information for this event
  //--------------------------------------------------
  
  _internal_timer[8].get()->restart();
  printOutputInfo();
  _internal_timer[8].get()->stop();
  
  ++_ievent;
  _timer.get()->stop();
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4Evaluator::End(PHCompositeNode *topNode)
{
  _tfile->cd();

  _ntp_event->Write();
  _ntp_g4hit->Write();
  _ntp_cell->Write();
  _ntp_cluster->Write();
  _ntp_gtrack->Write();
  _ntp_track->Write();

  _tfile->Close();

  delete _tfile;

  if (verbosity >= 0) {
    cout << "========================= PHG4Evaluator::End() ============================" << endl;
    cout << " " << _ievent << " events of output written to: " << _filename << endl;
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

/// loop over the geometry objects and determine
/// which layers are of what type (ladder or cylinder style)
int PHG4Evaluator::fillLayerTypeMap(PHCompositeNode *topNode) {

  _layer_type_map.clear();
  
  PHG4CylinderCellGeomContainer* cylcellgeos = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode,"CYLINDERCELLGEOM_SVTX");
  if (cylcellgeos) {

    _cellGeos = cylcellgeos;
    
    PHG4CylinderCellGeomContainer::ConstRange range = cylcellgeos->get_begin_end();
    for (PHG4CylinderCellGeomContainer::ConstIterator citer = range.first;
	 citer != range.second;
	 ++citer) {

      const PHG4CylinderCellGeom *geom = citer->second;
      int layerid = geom->get_layer();

      if (_layer_type_map.find(layerid) == _layer_type_map.end()) {
	_layer_type_map.insert(make_pair(layerid,CylinderLayer));
      } else {
	cout << PHWHERE << "Error: The layer id #" << layerid << " appears multiple times in the SVTX" << endl;
	return Fun4AllReturnCodes::ABORTRUN;
      }
    }

    _cellList = findNode::getClass<PHG4CylinderCellContainer>(topNode,"G4CELL_SVTX");
    if (!_cellList) {
      cerr << PHWHERE << " ERROR: Can't find G4CELL_SVTX" << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }    
  } else {
    _cellGeos = NULL;
    _cellList = NULL;
  }

  PHG4CylinderGeomContainer* laddergeos = findNode::getClass<PHG4CylinderGeomContainer>(topNode,"CYLINDERGEOM_SILICON_TRACKER");
  if (laddergeos) {

    _ladderGeos = laddergeos;

    PHG4CylinderGeomContainer::ConstRange range = laddergeos->get_begin_end();
    for (PHG4CylinderGeomContainer::ConstIterator citer = range.first;
	 citer != range.second;
	 ++citer) {

      const PHG4CylinderGeom *geom = citer->second;
      int layerid = geom->get_layer();

      if (_layer_type_map.find(layerid) == _layer_type_map.end()) {
	_layer_type_map.insert(make_pair(layerid,LadderLayer));
      } else {
	cout << PHWHERE << "Error: The layer id #" << layerid << " appears multiple times in the SVTX" << endl;
	return Fun4AllReturnCodes::ABORTRUN;
      }
    }

    _ladderCellList = findNode::getClass<PHG4CylinderCellContainer>(topNode,"G4CELL_SILICON_TRACKER");
    if (!_ladderCellList) {
      cerr << PHWHERE << " ERROR: Can't find G4Cell_SILICON_TRACKER" << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }     
  } else {
    _ladderGeos = NULL;
    _ladderCellList = NULL;
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4Evaluator::fillGhitList(PHCompositeNode *topNode) {

  _g4hitList.clear();
  
  // pull the g4hits
  PHG4HitContainer* g4hitList1 = findNode::getClass<PHG4HitContainer>(topNode,"G4HIT_SVTX");
  PHG4HitContainer* g4hitList2 = findNode::getClass<PHG4HitContainer>(topNode,"G4HIT_SILICON_TRACKER");
  if ((!g4hitList1) && (!g4hitList2)) {
    cerr << PHWHERE << " ERROR: Can't find G4HIT_SVTX or G4HIT_SILICON_TRACKER, aborting!" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (g4hitList1) {
    PHG4HitContainer::ConstRange g4hit_range = g4hitList1->getHits();
    for(PHG4HitContainer::ConstIterator g4hit_iter = g4hit_range.first;
	g4hit_iter != g4hit_range.second;
	g4hit_iter++) {

      if (_g4hitList.find(g4hit_iter->second->get_hit_id())
	  != _g4hitList.end()) {
	cout << PHWHERE << "Error: The g4hit id #" << g4hit_iter->second->get_hit_id() << " appears multiple times in the SVTX" << endl;
	return Fun4AllReturnCodes::ABORTRUN;
      }
      
      _g4hitList[g4hit_iter->second->get_hit_id()] = g4hit_iter->second;
    }
  }
  
  if (g4hitList2) {
    PHG4HitContainer::ConstRange g4hit_range = g4hitList2->getHits();
    for(PHG4HitContainer::ConstIterator g4hit_iter = g4hit_range.first;
	g4hit_iter != g4hit_range.second;
	g4hit_iter++) {

      if (_g4hitList.find(g4hit_iter->second->get_hit_id())
	  != _g4hitList.end()) {
	cout << PHWHERE << "Error: The g4hit id #" << g4hit_iter->second->get_hit_id() << " appears multiple times in the SVTX" << endl;
	return Fun4AllReturnCodes::ABORTRUN;
      }
      
      _g4hitList[g4hit_iter->second->get_hit_id()] = g4hit_iter->second;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4Evaluator::fillGtrackObjects()
{
  if(verbosity > 1) cout << "PHG4Evaluator::fillGtrackObjects() entered" << endl;

  // fill a vector of truth particles that know which g4hits they made

  // reset the vector
  _gtrack_list.clear();
  _g4hit_gtrack_map.clear();
  _particleid_g4hitid_mmap.clear();
  
  map<int, unsigned int> id_gtrack_map;
  
  // loop over the g4hits
  for(HitMap::const_iterator g4hit_iter = _g4hitList.begin();
      g4hit_iter != _g4hitList.end();
      g4hit_iter++)
    {
      PHG4Hit *g4hit = g4hit_iter->second;
      PHG4HitDefs::keytype g4hit_id = g4hit->get_hit_id();

      int track_id = g4hit->get_trkid();
      
      _particleid_g4hitid_mmap.insert(make_pair(track_id,g4hit_id));
      
      bool isNewGtrack = true;
      bool isNewG4hit = true;
    
      map<int, unsigned int>::iterator giter = id_gtrack_map.find(track_id);
      if(giter != id_gtrack_map.end())
	{
	  unsigned int igtrack = giter->second;
	  isNewGtrack = false;
      
	  // loop to see if we have this hit already
	  for(unsigned int ig4hit = 0; ig4hit < _gtrack_list[igtrack].get_ng4hits(); ig4hit++)
	    {
	      if(g4hit == _gtrack_list[igtrack].get_g4hit(ig4hit))
		{
		  isNewG4hit = false;
		}
	    }
      
	  // if the hit isn't on the gtrack already, add it
	  if(isNewG4hit)
	    {
	      SvxGtrack *gtrack = &_gtrack_list[igtrack];
        
	      gtrack->add_g4hit(g4hit);
	    }
	}
    
      // if the g4hit doesn't make a gtrack, create a new one
      if(isNewGtrack)
	{
	  SvxGtrack gtrack;

	  PHG4Particle* particle = _truth_info_container->GetHit( track_id );

	  gtrack.set_track_id(track_id);
	  gtrack.set_particle(particle);
        
	  gtrack.set_flavor(particle->get_pid());
	  gtrack.set_px( particle->get_px() );
	  gtrack.set_py( particle->get_py() );
	  gtrack.set_pz( particle->get_pz() );
        
	  gtrack.set_embed(_truth_info_container->isEmbeded(particle->get_track_id()));

	  gtrack.set_primary(false);
	  PHG4TruthInfoContainer::Map primary_map = _truth_info_container->GetPrimaryMap();
	  for (PHG4TruthInfoContainer::ConstIterator iter = primary_map.begin(); 
	       iter != primary_map.end(); 
	       ++iter) {
	    if (iter->second->get_track_id() == particle->get_track_id() ) {
	      gtrack.set_primary(true);
	    }
	  }

	  // get the vertex
	  PHG4VtxPoint* vertex = _truth_info_container->GetVtx( particle->get_vtx_id() );

	  float vx = NAN;
	  float vy = NAN;
	  float vz = NAN;

	  if (vertex)
	    {
	      vx = vertex->get_x();
	      vy = vertex->get_y();
	      vz = vertex->get_z();
	    }	  

	  gtrack.set_vx( vx );
	  gtrack.set_vy( vy );
	  gtrack.set_vz( vz );
	  
	  gtrack.set_standalone_match(0.0);
	  
	  gtrack.add_g4hit(g4hit);	  
	  
	  _gtrack_list.push_back(gtrack);
	  id_gtrack_map[track_id] = (_gtrack_list.size() - 1);
	}
      
      // otherwise keep scanning the g4hits for new mc tracks
    } // loop over all g4hits

  // now loop over all gtracks and fill a reverse lookup map from the
  // g4hit back to the gtrack
  for(unsigned int igtrack = 0; igtrack < _gtrack_list.size(); igtrack++)
    {
      SvxGtrack *gtrack = &_gtrack_list[igtrack];
      
      unsigned int outer_hit = 0;
      unsigned int outer_hit_layer = 0;

      for(unsigned int ig4hit = 0; ig4hit < gtrack->get_ng4hits(); ig4hit++)
	{
	  PHG4Hit *g4hit = gtrack->get_g4hit(ig4hit);
	 
	  if( g4hit->get_layer() > outer_hit_layer )
	    {
	      outer_hit = ig4hit;
	      outer_hit_layer = g4hit->get_layer();
	    }
 
	  // fill map
	  _g4hit_gtrack_map.insert(make_pair(g4hit,gtrack));
	}
      
      // also fill the exit parameters from the last hit
      PHG4Hit *g4hit = gtrack->get_g4hit(outer_hit);

      gtrack->set_fx( g4hit->get_x(1) );
      gtrack->set_fy( g4hit->get_y(1) );
      gtrack->set_fz( g4hit->get_z(1) );

      gtrack->set_fpx( g4hit->get_px(1) );
      gtrack->set_fpy( g4hit->get_py(1) );
      gtrack->set_fpz( g4hit->get_pz(1) );
    }

  // get the last particle index from the truth container
  int last_particle_id = _truth_info_container->GetLastParticleIndex();
  
  // set last shower flag
  for(unsigned int igtrack = 0; igtrack < _gtrack_list.size(); igtrack++)
    {
      if(_gtrack_list[igtrack].get_track_id() == last_particle_id)
	{
	  _gtrack_list[igtrack].set_is_last(true);
	}
    }
}

void PHG4Evaluator::fillCellToG4HitMap()
{
  if(verbosity > 1) cout << "PHG4Evaluator::fillCellToG4hitMap() entered" << endl;

  // descend into the ancestry information and fill a map between
  // the reconstructed cluster and the primary contributing g4hit
  // (as g4hits are the closest monte carlo analog to the cluster)

  // reset the map for this event
  _cell_g4hit_mmap.clear();
  _g4hit_cell_mmap.clear();

  if (_cellList) {
    PHG4CylinderCellContainer::ConstRange cellrange = _cellList->getCylinderCells();
    for(PHG4CylinderCellContainer::ConstIterator celliter = cellrange.first;
	celliter != cellrange.second;
	++celliter)
      {
	PHG4CylinderCell* cell = celliter->second;
	  
	++_nhits_per_layer[cell->get_layer()];

	PHG4Hit* max_g4hit = NULL;
	float max_edep = 0.0;

	// loop over all ghits within those cells
	PHG4CylinderCell::EdepConstRange g4hits = cell->get_g4hits();
	PHG4CylinderCell::EdepConstIterator g4iter = g4hits.first;
	for (g4iter = g4hits.first; g4iter != g4hits.second; g4iter++)
	  {
	    // ask each ghit how much energy it deposited and where that energy came from
	    PHG4HitDefs::keytype g4hit_id = g4iter->first;
	      
	    HitMap::const_iterator tmpiter = _g4hitList.find(g4hit_id);
	    if(tmpiter==_g4hitList.end()) {
	      continue;
	    }
	      
	    PHG4Hit* g4hit = tmpiter->second;
	      
	    // fill our handy map so we never have to do this again for this module...
	    //_cell_allg4hits_mmap.insert( make_pair(cell, g4hit) );
	    //if (_allg4hits_cluster_mmap.find(g4hit) == _allg4hits_cluster_mmap.end()) {
	    //  _allg4hits_cluster_mmap.insert( make_pair(g4hit, cluster) );
	    //}
	    
	    if (g4hit->get_edep() > max_edep) {
	      max_g4hit = g4hit;
	      max_edep = g4hit->get_edep();
	    }
	  } // ghit loop
          
	// fill our handy map so we never have to do this again for this module...
	_cell_g4hit_mmap.insert( make_pair(cell, max_g4hit) );
	  
	// fill our handy reverse lookup multimap
	_g4hit_cell_mmap.insert( make_pair(max_g4hit, cell) );
      }
  }

  if (_ladderCellList) {
    PHG4CylinderCellContainer::ConstRange cellrange = _ladderCellList->getCylinderCells();
    for(PHG4CylinderCellContainer::ConstIterator celliter = cellrange.first;
	celliter != cellrange.second;
	++celliter)
      {
	PHG4CylinderCell* cell = celliter->second;
	  
	++_nhits_per_layer[cell->get_layer()];

	PHG4Hit* max_g4hit = NULL;
	float max_edep = 0.0;

	// loop over all ghits within those cells
	PHG4CylinderCell::EdepConstRange g4hits = cell->get_g4hits();
	PHG4CylinderCell::EdepConstIterator g4iter = g4hits.first;
	for (g4iter = g4hits.first; g4iter != g4hits.second; g4iter++)
	  {
	    // ask each ghit how much energy it deposited and where that energy came from
	    PHG4HitDefs::keytype g4hit_id = g4iter->first;
	      
	    HitMap::const_iterator tmpiter = _g4hitList.find(g4hit_id);
	    if(tmpiter==_g4hitList.end()) {
	      continue;
	    }
	      
	    PHG4Hit* g4hit = tmpiter->second;
	      
	    // fill our handy map so we never have to do this again for this module...
	    //_cell_allg4hits_mmap.insert( make_pair(cell, g4hit) );
	    //if (_allg4hits_cluster_mmap.find(g4hit) == _allg4hits_cluster_mmap.end()) {
	    //  _allg4hits_cluster_mmap.insert( make_pair(g4hit, cluster) );
	    //}
	    
	    if (g4hit->get_edep() > max_edep) {
	      max_g4hit = g4hit;
	      max_edep = g4hit->get_edep();
	    }
	  } // ghit loop
          
	// fill our handy map so we never have to do this again for this module...
	_cell_g4hit_mmap.insert( make_pair(cell, max_g4hit) );
	  
	// fill our handy reverse lookup multimap
	_g4hit_cell_mmap.insert( make_pair(max_g4hit, cell) );
      }
  }

  return;
}

/// loop over all the reco'd SvtxClusters and create the ancestry links to
/// the G4Hits from the cylinder cell layers or ladder layers.
int PHG4Evaluator::fillClusterToG4HitLinks(PHCompositeNode *topNode) {

  if (verbosity > 1) cout << "PHG4Evaluator::fillClusterToG4HitLinks() entered" << endl;

  EvalLinks *evalcyllinks = NULL;
  if (_cellList) {
    PHTypedNodeIterator<EvalLinks> evaliter(topNode);
    PHIODataNode<EvalLinks> *EvalLinksNode = evaliter.find("SvtxClusterMap_G4HIT_SVTX_Links");
    if (!EvalLinksNode) {
      cout << PHWHERE << " ERROR: Can't find SvtxClusterMap_G4HIT_SVTX_Links" << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    } else {
      evalcyllinks = (EvalLinks*)EvalLinksNode->getData();
    }
    evalcyllinks->Reset();
    evalcyllinks->set_names("SvtxClusterMap","G4HIT_SVTX","edep");
  }
  
  EvalLinks *evalladderlinks = NULL;
  if (_ladderCellList) {
    PHTypedNodeIterator<EvalLinks> evaliter(topNode);
    PHIODataNode<EvalLinks> *EvalLinksNode = evaliter.find("SvtxClusterMap_G4HIT_SILICON_TRACKER_Links");
    if (!EvalLinksNode) {
      cout << PHWHERE << " ERROR: Can't find SvtxClusterMap_G4HIT_SILICON_TRACKER_Links" << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    } else {
      evalladderlinks = (EvalLinks*)EvalLinksNode->getData();
    }
    evalladderlinks->Reset();
    evalladderlinks->set_names("SvtxClusterMap","G4HIT_SILICON_TRACKER","edep");
  }
  
  // loop over all the clusters
  for (SvtxClusterMap::Iter iter = _clusterList->begin();
       iter != _clusterList->end();
       ++iter) {
  
    SvtxCluster* cluster = &iter->second;
    unsigned int cluster_id = cluster->get_id();

    // loop over all SvtxHits
    for (SvtxCluster::HitIter hiter = cluster->begin_hits();
	 hiter != cluster->end_hits();
	 ++hiter) {
      unsigned int hit_id = *hiter;

      SvtxHit* hit = _hitList->get(hit_id);   
      PHG4CylinderCellDefs::keytype cell_id = hit->get_cellid();
      
      // get this cell...
      if ((_cellList)&&(_cellList->findCylinderCell(cell_id))) {    
	PHG4CylinderCell *cell = _cellList->findCylinderCell(cell_id);     

	// loop over all g4hits within the cell
	PHG4CylinderCell::EdepConstRange range = cell->get_g4hits();
	for (PHG4CylinderCell::EdepConstIterator iter = range.first;
	     iter != range.second;
	     ++iter) {
	  PHG4HitDefs::keytype ghit_id = iter->first;
	  float edep = iter->second;

	  // create new link
	  evalcyllinks->link(cluster_id,ghit_id,edep);
	}
      } else if ((_ladderCellList)&&(_ladderCellList->findCylinderCell(cell_id))) {
	PHG4CylinderCell *cell = _ladderCellList->findCylinderCell(cell_id);
	
	// loop over all g4hits within the cell
	PHG4CylinderCell::EdepConstRange range = cell->get_g4hits();
	for (PHG4CylinderCell::EdepConstIterator iter = range.first;
	     iter != range.second;
	     ++iter) {
	  PHG4HitDefs::keytype ghit_id = iter->first;
	  float edep = iter->second;
	  
	  // create new link
	  evalladderlinks->link(cluster_id,ghit_id,edep);
	}   
      } else {
	cout << "SVTX evalution traces to missing cell_id, corrupt ancestry error" << endl;
	exit(-1);
      } 	
    } // hit loop      
  } // cluster loop

  //evalcyllinks->identify();
  
  _cluster_g4hit_svtx_links = evalcyllinks;
  _cluster_g4hit_silicon_tracker_links = evalladderlinks;

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4Evaluator::fillClusterToG4HitMap()
{
  // reset the map for this event
  _cluster_g4hit_mmap.clear();
  _g4hit_cluster_mmap.clear();

  // fill convienence mmap
  _cluster_allg4hits_mmap.clear();
  _allg4hits_cluster_mmap.clear();

  // loop over all the clusters
  for (SvtxClusterMap::Iter iter = _clusterList->begin();
       iter != _clusterList->end();
       ++iter) {

    SvtxCluster* cluster = &iter->second;
      
    // this array will hold pointers to the possible primary choices
    // for this cluster
    vector<PHG4Hit*> cluster_contributors;
    cluster_contributors.clear();

    // right now the clusters only know about one contributing
    // PHG4Hit, so I'll skip all the largest-contributor stuff
    // (adding this stuff back in -MPM)

    PHG4Hit* max_g4hit = NULL;
    float max_edep = 0.0;

    // loop over all hit cells
    for (SvtxCluster::HitIter hiter = cluster->begin_hits();
	 hiter != cluster->end_hits();
	 ++hiter) {
      unsigned int hit_id = *hiter;

      SvtxHit* hit = _hitList->get(hit_id);
      unsigned int cell_id = hit->get_cellid();
      
      PHG4CylinderCell *cell = _cellList->findCylinderCell(cell_id);
      if (!cell) {
	cell = _ladderCellList->findCylinderCell(cell_id);
      }
      if (!cell) continue;

      // loop over the cells g4hits...
      for (std::multimap<PHG4CylinderCell*,PHG4Hit*>::iterator giter = _cell_g4hit_mmap.lower_bound(cell);
	   giter != _cell_g4hit_mmap.upper_bound(cell);
	   ++giter) {

	PHG4Hit* g4hit = giter->second;
	  
	// fill our handy map so we never have to do this again for this module...
	_cluster_allg4hits_mmap.insert( make_pair(cluster, g4hit) );
	if (_allg4hits_cluster_mmap.find(g4hit) == _allg4hits_cluster_mmap.end()) {
	  _allg4hits_cluster_mmap.insert( make_pair(g4hit, cluster) );
	}
	  
	if (g4hit->get_edep() > max_edep) {
	  max_g4hit = g4hit;
	  max_edep = g4hit->get_edep();
	}
      } // ghit loop
    } // hit loop
      
      // fill our handy map so we never have to do this again for this module...
    _cluster_g4hit_mmap.insert( make_pair(cluster, max_g4hit) );

    // fill our handy reverse lookup multimap
    _g4hit_cluster_mmap.insert( make_pair(max_g4hit, cluster) );
  } // cluster loop

  return;
}

int PHG4Evaluator::fillTrackToG4TruthInfoLinks(PHCompositeNode *topNode) {

  EvalLinks* evallinks = NULL;
  PHTypedNodeIterator<EvalLinks> evaliter(topNode);
  PHIODataNode<EvalLinks> *EvalLinksNode = evaliter.find("SvtxTrackMap_G4TruthInfo_Links");
  if (!EvalLinksNode) {
    cout << PHWHERE << " ERROR: Can't find SvtxTrackMap_G4TruthInfo_Links" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  } else {
    evallinks = (EvalLinks*)EvalLinksNode->getData();
  }
  evallinks->Reset();
  evallinks->set_names("SvtxTrackMap","G4TruthInfo","ng4hits");

  _track_particle_links = evallinks;
  
  // loop over all tracks
  for (SvtxTrackMap::Iter iter = _trackList->begin();
       iter != _trackList->end();
       ++iter) {
    SvtxTrack *track = &iter->second;
    unsigned int track_id = track->getTrackID();
    std::map<unsigned int,unsigned int> track_votes; // particle id => votes

    // loop over all the layers
    short found_hits = 0;
    for (unsigned int ilayer = 0; ilayer < 100; ++ilayer) {

      // particles can get at most one vote per layer
      std::set<unsigned int> layer_votes;
      
      // get the associated clusters for the layer
      if (track->hasCluster(ilayer)) {
	++found_hits;
	
	unsigned int cluster_id = track->getClusterID(ilayer);

	// loop over all g4hits associated to cluster
	// using the previously established eval links
	std::set<unsigned long long> g4hit_ids;
	if (_cluster_g4hit_svtx_links) g4hit_ids = _cluster_g4hit_svtx_links->right(cluster_id);
	for (std::set<unsigned long long>::iterator iiter = g4hit_ids.begin();
	     iiter != g4hit_ids.end();
	     ++iiter) {
	  unsigned long long g4hit_id = *iiter;

	  HitMap::const_iterator tmpiter = _g4hitList.find(g4hit_id);
	  PHG4Hit* g4hit = tmpiter->second;

	  unsigned int particle_id = g4hit->get_trkid();

	  layer_votes.insert(particle_id);
	}

	// loop over all g4hits associated to cluster
	// using the previously established eval links
	g4hit_ids.clear();
	if (_cluster_g4hit_silicon_tracker_links) g4hit_ids = _cluster_g4hit_silicon_tracker_links->right(cluster_id);
	for (std::set<unsigned long long>::iterator iiter = g4hit_ids.begin();
	     iiter != g4hit_ids.end();
	     ++iiter) {
	  unsigned long long g4hit_id = *iiter;

	  HitMap::const_iterator tmpiter = _g4hitList.find(g4hit_id);
	  PHG4Hit* g4hit = tmpiter->second;

	  unsigned int particle_id = g4hit->get_trkid();

	  layer_votes.insert(particle_id);
	} // associated g4hit loop
      } // layer has cluster

      // loop over all layer votes and insert into vote counter
      for (std::set<unsigned int>::iterator viter = layer_votes.begin();
	   viter != layer_votes.end();
	   ++viter) {
	unsigned int particle_id = *viter;
	if (track_votes.find(particle_id) == track_votes.end()) {
	  track_votes[particle_id] = 1;
	} else {
	  ++track_votes[particle_id];
	}
      } // end layer vote counting
      
      if (found_hits >= track->getNhits()) break; // all clusters visited
    } // layer loop

    // loop over vote map
    for (std::map<unsigned int,unsigned int>::iterator iiter = track_votes.begin();
	 iiter != track_votes.end();
	 ++iiter) {
      unsigned int particle_id = iiter->first;
      unsigned int votes       = iiter->second;     
      evallinks->link(track_id,particle_id,(float)votes);
    } // end vote recording
    
  } // track loop

  return Fun4AllReturnCodes::EVENT_OK;
}


// forward look up based on number of clusters left by truth
// particles (pattern reco testing).
int PHG4Evaluator::fillG4TruthInfoToTrackLinks(PHCompositeNode *topNode) {

  EvalLinks* evallinks = NULL;
  PHTypedNodeIterator<EvalLinks> evaliter(topNode);
  PHIODataNode<EvalLinks> *EvalLinksNode = evaliter.find("G4TruthInfo_SvtxTrackMap_Links");
  if (!EvalLinksNode) {
    cout << PHWHERE << " ERROR: Can't find G4TruthInfo_SvtxTrackMap_Links" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  } else {
    evallinks = (EvalLinks*)EvalLinksNode->getData();
  }
  evallinks->Reset();
  evallinks->set_names("G4TruthInfo","SvtxTrackMap","nclusters left by particle");

  // loop over all truth particles
  PHG4TruthInfoContainer::Map map = _truth_info_container->GetMap();
  for (PHG4TruthInfoContainer::ConstIterator citer = map.begin();
       citer != map.end();
       ++citer) {
    PHG4Particle* particle = citer->second;
    int particle_id = particle->get_track_id();
    unsigned int nclusters = 0;
      
    // loop over this particle's g4hits
    for (std::multimap<int,PHG4HitDefs::keytype>::iterator iter = _particleid_g4hitid_mmap.lower_bound(particle_id);
	 iter != _particleid_g4hitid_mmap.upper_bound(particle_id);
	 ++iter) {

      PHG4HitDefs::keytype g4hitid = iter->second;
      
      // grab the number of clusters from the g4hit
      if (_cluster_g4hit_svtx_links) {
	std::set<unsigned long long> cluster_ids = _cluster_g4hit_svtx_links->left(g4hitid);
	nclusters += cluster_ids.size();
      }

      if (_cluster_g4hit_silicon_tracker_links) {
	std::set<unsigned long long> cluster_ids = _cluster_g4hit_silicon_tracker_links->left(g4hitid);
	nclusters += cluster_ids.size();
      }      
    }    

    // the number of clusters has been established
    unsigned int leading_track_id = _track_particle_links->max_left(particle_id);
    evallinks->link(particle_id,leading_track_id,nclusters);    
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4Evaluator::fillTrackToGtrackMap()
{
  if(verbosity > 1) cout << "PHG4Evaluator::fillTrackToGtrackMap() entered" << endl;

  // the idea here is to one time only decend into the ancesrty and produce a map between
  // the reconstructed tracks and truth particles that can be used to bridge this gap
  // more easily from now on

  // some time soon we should fill an SvxMCParticle so these comparisons... but this will
  // exercise the important part of digging through the ancestry
  // >I don't understand the above comment. Truth particles are already available 
  // >in SvxGtrack. -MPM

  // ok here we go...

  _track_gtrack_map.clear();
  _gtrack_track_mmap.clear();

  _track_gtrack_match.clear();
  _gtrack_track_match.clear();

  _track_cluster_mmap.clear();
  
  //------------------------------
  // fill some convenience storage
  //------------------------------

  // fill a convenience map between cluster id and cluster pointers
  map<unsigned int, SvtxCluster*> id_cluster_map;
  // redundant can be removed now
  for (SvtxClusterMap::Iter iter = _clusterList->begin();
       iter != _clusterList->end();
       ++iter) {
    SvtxCluster* cluster = &iter->second;
    id_cluster_map[cluster->get_id()] = cluster; 
  }
  
  // fill a convenience map between track id and track pointers
  map<int, SvxGtrack*> id_gtrack_map;
  for (unsigned int igtrack = 0; igtrack < _gtrack_list.size(); igtrack++) {
    id_gtrack_map[_gtrack_list[igtrack].get_track_id()] = &(_gtrack_list[igtrack]);
  }
  
  // fill a map between the reco tracks and their clusters
  // loop over all tracks
  for (SvtxTrackMap::Iter iter = _trackList->begin();
       iter != _trackList->end();
       ++iter) {
    SvtxTrack *track = &iter->second;
      
    // loop over all the layers
    for (unsigned int ilayer = 0; ilayer < 10; ilayer++) {

      // loop over the associated clusters
      if (track->hasCluster(ilayer)) {
	  
	unsigned int cluster_id = track->getClusterID(ilayer);
	  
	map<unsigned int, SvtxCluster*>::iterator clus_iter = id_cluster_map.find(cluster_id);
	if (clus_iter != id_cluster_map.end()) {
	  SvtxCluster *cluster = clus_iter->second;
	    
	  // fill a quick lookup in case i want it later...
	  _track_cluster_mmap.insert(make_pair(track,cluster));
	}
      }
    }
  }

  //------------------------------------------------------------------
  // Determine a maximum truth contributor to each reconstructed track
  //------------------------------------------------------------------

  // loop over all reco tracks
  for (SvtxTrackMap::Iter iter = _trackList->begin();
       iter != _trackList->end();
       ++iter) {
    SvtxTrack *track = &iter->second;

    // this will hold a list of possible track ids
    // each entry votes for a truth particle
    vector<unsigned int> track_contributors;

    // loop over all truth particles
    for (unsigned int igtrack = 0; igtrack < _gtrack_list.size(); ++igtrack) {

      // loop over all the layers in track
      for(unsigned int ilayer = 0; ilayer < 10; ilayer++) {

	// loop over the track clusters
	if (track->hasCluster(ilayer)) {

	  unsigned int cluster_id = track->getClusterID(ilayer);

	  bool should_vote = false;
	    
	  map<unsigned int, SvtxCluster*>::iterator clus_iter = id_cluster_map.find(cluster_id);
	  if (clus_iter != id_cluster_map.end()) {
	    SvtxCluster *cluster = clus_iter->second;
	      
	    // loop over all of the contributing g4hits to this cluster
	    typedef multimap<SvtxCluster*,PHG4Hit*>::const_iterator mmapiter;
	    pair<mmapiter,mmapiter> hitrange = _cluster_allg4hits_mmap.equal_range(cluster);
	    for(mmapiter hititer = hitrange.first;
		hititer!=hitrange.second;
		hititer++) {
	      PHG4Hit *g4hit = hititer->second;
		    
	      if (!g4hit) continue; // this must have been a noise cluster

	      // we find this truth particle had a g4hit in this cluster
	      // we will add one and only one vote for this cluster 
	      if (g4hit->get_trkid() == _gtrack_list[igtrack].get_track_id()) {
		should_vote = true;
	      }
	    } // g4hit loop

	      // if this cluster contained a g4hit from this gtrack we vote for it
	    if (should_vote) {
	      track_contributors.push_back(_gtrack_list[igtrack].get_track_id());
	    }
	  } // cluster loop
	} // associated cluster if
      }  // layer loop
    } // truth particle loop

      // loop over the contributors looking for the largest contributor
    int max_contributor = 0;
    unsigned int max_contributor_count = 0;
    for (unsigned int icontributor = 0; icontributor < track_contributors.size(); icontributor++) {
      // get the degeneracy of this contibutor in the list
      unsigned int count = std::count(track_contributors.begin(),track_contributors.end(),track_contributors.at(icontributor));
	
      // if it exceeds the current max, reset with this one
      if (count > max_contributor_count) {
	max_contributor_count = count;
	max_contributor = track_contributors.at(icontributor);
      }
    }

    // look for this greatest contributor particle in our list of gtracks
    SvxGtrack *max_contributor_gtrack = NULL;
    map<int, SvxGtrack*>::iterator cont_iter = id_gtrack_map.find(max_contributor);
    if (cont_iter != id_gtrack_map.end()) { 
      max_contributor_gtrack = cont_iter->second;
    }

    // fill the even handier map between tracks and our gtracks...
    _track_gtrack_map.insert(make_pair(track,max_contributor_gtrack));

    // fill our handy reverse lookup multimap
    _gtrack_track_mmap.insert(make_pair(max_contributor_gtrack,track));

  } // track loop

  //------------------------------------------------------------------------
  // Determine if the deposit left by each gtrack was properly reconstructed
  //------------------------------------------------------------------------
  
  // convenience storage, this gtrack left how many clusters?
  map<SvxGtrack*, unsigned int> gtrack_nhits_map;
  
  // loop over the gtracks and look for the matching reconstructions
  for (unsigned int igtrack = 0; igtrack < _gtrack_list.size(); igtrack++) {

    SvxGtrack *gtrack = &_gtrack_list[igtrack];
    gtrack_nhits_map[gtrack] = 0;
	  
    //
    // check that this gtrack left N clusters for N g4hits
    //

    // loop over all g4hits
    bool badCount = false;
    unsigned int nclusters = 0; // produced clusters
    for (unsigned int ig4hit = 0; ig4hit < gtrack->get_ng4hits(); ++ig4hit) {
      PHG4Hit *g4hit = gtrack->get_g4hit(ig4hit);
	  
      // loop over all clusters this
      unsigned int count = _allg4hits_cluster_mmap.count(g4hit);
      nclusters += count;
	      
      // this g4hit left more than a single cluster
      if (count != 1) badCount = true;
    }

    // if a cluster was split or the cluster was missed, move on
    // to the next track
    if (nclusters != gtrack->get_ng4hits()) continue;
    if (badCount) continue;

    //
    // these truth particles at this point should be reconstructed
    //

    // set the default to a bad reconstruction
    gtrack->set_standalone_match(-1.0*gtrack->get_ng4hits());

    // loop over all reconstructed tracks that think this track was its
    // greatest truth contributor
    std::pair<std::multimap<SvxGtrack*,SvtxTrack*>::iterator,
	      std::multimap<SvxGtrack*,SvtxTrack*>::iterator> 
      range = _gtrack_track_mmap.equal_range(gtrack);
    for (std::multimap<SvxGtrack*,SvtxTrack*>::iterator it = range.first; it!=range.second; ++it) {
 
      SvtxTrack *track = it->second;
	      
      // check that the number of hits match
      unsigned int ng4hits = gtrack->get_ng4hits();
      unsigned int nseg4hits = 0;

      // loop over the associated clusters
      for (int ilayer = 0; ilayer < 10; ilayer++) {
	if (track->hasCluster(ilayer)) ++nseg4hits;
      }
	  
      if (ng4hits != nseg4hits) {
	// this gtrack produced a track with a different number of hits
	// than the number of clusters left by the truth particle
	continue; 	      
      }

      //
      // these reco tracks now have the same number of clusters
      // as the truth particle left in the detector
      //

      // now check that all the g4hits and clusters are matches
      vector<bool> hitsMatch(gtrack->get_ng4hits(), false);
	
      // loop over all g4hits from the truth particle
      for (unsigned int ig4hit = 0; ig4hit < gtrack->get_ng4hits(); ++ig4hit) {

	PHG4Hit *g4hit = gtrack->get_g4hit(ig4hit);
	  
	// loop over the reconstructed track clusters
	typedef multimap<SvtxTrack*,SvtxCluster*>::iterator mapiterS2C;
	typedef pair<mapiterS2C,mapiterS2C> maprangeS2C;
	maprangeS2C rangeS2C = _track_cluster_mmap.equal_range( track );
	for (mapiterS2C iterS2C = rangeS2C.first; iterS2C != rangeS2C.second; iterS2C++) {

	  SvtxCluster *cluster = iterS2C->second;
		
	  // is the g4hits connected to this cluster
	  typedef multimap<SvtxCluster*,PHG4Hit*>::const_iterator mmapiter;
	  pair<mmapiter,mmapiter> hitrange = _cluster_allg4hits_mmap.equal_range(cluster);
	  for(mmapiter hititer = hitrange.first;
	      hititer!=hitrange.second;
	      hititer++) {
	    if (g4hit == hititer->second) {
	      // this cluster contains content from our truth particle
	      hitsMatch[ig4hit] = true;
	    }
	  }

	} // cluster loop
      } // g4hit loop
	  
      unsigned int found_hits = 0;
      for (unsigned int ig4hit = 0; ig4hit < gtrack->get_ng4hits(); ig4hit++) {
	if (hitsMatch[ig4hit]) {
	  found_hits += 1;
	}
      }
	    
      if (found_hits > gtrack_nhits_map[gtrack]) {
	gtrack_nhits_map[gtrack] = found_hits;
      }

      // we correctly found all the truth particle hits
      // we now designate this match connection
      // we only care if the track is corretly reconstructed once
      // ghosts and additional fakes don't count against this metric 
      // (see purity for those things)
      if (found_hits == gtrack->get_ng4hits()) {
	gtrack->set_standalone_match(gtrack->get_ng4hits());
	_gtrack_track_match.insert(make_pair(gtrack,track));
	_track_gtrack_match.insert(make_pair(track,gtrack));
      }

    } // reco track loop
  } // gtrack loop
  
  return;
}

void PHG4Evaluator::fillTrackPurityMap()
{
  _track_purity_map.clear();

  if(verbosity > 1) cout << "PHG4Evaluator::fillTrackPurityMap() entered" << endl;

  for (SvtxTrackMap::Iter iter = _trackList->begin();
       iter != _trackList->end();
       ++iter) {
    SvtxTrack *track = &iter->second;
    SvxGtrack  *gtrack  = _track_gtrack_map[track];
      
    unsigned int purity = 0;

    if(gtrack)
      {
	// loop over all the tracks clusters
	typedef multimap<SvtxTrack*,SvtxCluster*>::iterator mapiterS2C;
	typedef pair<mapiterS2C,mapiterS2C> maprangeS2C;
	maprangeS2C rangeS2C = _track_cluster_mmap.equal_range( track );
	for(mapiterS2C iterS2C = rangeS2C.first; iterS2C != rangeS2C.second; iterS2C++) {

	  SvtxCluster *cluster = iterS2C->second;
          
	  // loop over all the g4hits contributing to this cluster
	  // if the primary truth contributed in anyway to the cluster
	  // count it as good
	  bool good_cluster = false;
	  typedef multimap<SvtxCluster*,PHG4Hit*>::const_iterator mmapiter;
	  pair<mmapiter,mmapiter> hitrange = _cluster_allg4hits_mmap.equal_range(cluster);
	  for(mmapiter hititer = hitrange.first;
	      hititer!=hitrange.second;
	      hititer++) {
	    PHG4Hit *g4hit = hititer->second;
	    if (g4hit->get_trkid() == gtrack->get_track_id()) good_cluster = true;
	  }
	    
	  if (good_cluster) ++purity;		  
	}
      }

    if (purity > gtrack->get_best_purity()) gtrack->set_best_purity(purity);

    float dp = sqrt(pow(gtrack->get_px()-track->get3Momentum(0),2) + 
		    pow(gtrack->get_py()-track->get3Momentum(1),2) + 
		    pow(gtrack->get_pz()-track->get3Momentum(2),2));
    float p = sqrt(pow(gtrack->get_px(),2) + 
		   pow(gtrack->get_py(),2) + 
		   pow(gtrack->get_pz(),2));
    float dpp = dp/p;

    if (dpp < gtrack->get_best_dpp()) gtrack->set_best_dpp(dpp);

    _track_purity_map.insert(make_pair(track,purity));
  }

  return;
}

void PHG4Evaluator::printInputInfo()
{
  if(verbosity > 1) cout << "PHG4Evaluator::fillClusterToG4hitMap() entered" << endl;

  if(verbosity > 3)
    {
      // event information
      cout << endl;
      cout << PHWHERE << "   INPUT FOR EVENT " << _ievent << endl;

      cout << endl;
      cout << "---PHG4HITCONTAINER-------------" << endl;
      unsigned int ig4hit=0;
      for(HitMap::const_iterator g4hit_iter = _g4hitList.begin();
	  g4hit_iter != _g4hitList.end();
	  g4hit_iter++)
	{
	  PHG4Hit *g4hit = g4hit_iter->second;
	  cout << ig4hit << " of " << _g4hitList.size();
	  cout << ": PHG4Hit: " << g4hit << endl;
	  ++ig4hit;
	}

      cout << "---PHG4SVTXSIMPLECLUSTERLIST-------------" << endl;
      unsigned int icluster = 0;
      for (SvtxClusterMap::Iter iter = _clusterList->begin();
	   iter != _clusterList->end();
	   ++iter) {
	SvtxCluster* cluster = &iter->second;
	cout << icluster << " of " << _clusterList->size();	  
	cout << ": SvtxCluster: " << cluster << endl;
	++icluster;
      }

      if(_trackingWasRun)
	{
	  cout << "---SVXTRACKLIST-------------" << endl;
	  for (SvtxTrackMap::Iter iter = _trackList->begin();
	       iter != _trackList->end();
	       ++iter) {
	    SvtxTrack *track = &iter->second;
	    cout << "SvtxTrack:" << endl;
	    track->identify(cout);
	    cout << endl;
	  }
	}

    }

  return;
}

void PHG4Evaluator::printOutputInfo()
{
  if(verbosity > 1) cout << "PHG4Evaluator::printLogInfo() entered" << endl;

  //==========================================
  // print out some useful stuff for debugging
  //==========================================

  if(verbosity > 0)
    {
      // event information
      cout << endl;
      cout << PHWHERE << "   NEW OUTPUT FOR EVENT " << _ievent << endl;
      cout << endl;

      PHG4VtxPoint *gvertex = _truth_info_container->GetPrimaryVtx( _truth_info_container->GetPrimaryVertexIndex() );
      float gvx = gvertex->get_x();
      float gvy = gvertex->get_y();
      float gvz = gvertex->get_z();

      float vx = NAN;
      float vy = NAN;
      float vz = NAN;
      if (_vertexList) {
	if (!_vertexList->empty()) {
	  SvtxVertex* vertex = &(_vertexList->begin()->second);
	
	  vx = vertex->get_x();
	  vy = vertex->get_y();
	  vz = vertex->get_z();
	}
      }

      cout << "===Vertex Reconstruction=======================" << endl;
      cout << "vtrue = (" << gvx << "," << gvy << "," << gvz << ") => vreco = (" << vx << "," << vy << "," << vz << ")" << endl;
      cout << endl;

      cout << "===Tracking Summary============================" << endl;
      unsigned int ng4hits[10] = {0};
      for(unsigned int ilayer=0; ilayer<_nlayers; ++ilayer)
	{
	  for(HitMap::const_iterator g4hit_iter = _g4hitList.begin();
	      g4hit_iter != _g4hitList.end();
	      g4hit_iter++)
	    {
	      PHG4Hit *g4hit = g4hit_iter->second;
	      if (g4hit->get_layer() == ilayer) ++ng4hits[ilayer];
	    }
	}

      unsigned int nclusters[10] = {0};
      for (SvtxClusterMap::Iter iter = _clusterList->begin();
	   iter != _clusterList->end();
	   ++iter) {

	SvtxCluster* cluster = &iter->second;
	++nclusters[cluster->get_layer()];
      }

      for(unsigned int ilayer=0; ilayer<_nlayers; ++ilayer) {
	cout << "layer " << ilayer << ": nG4hits = " << ng4hits[ilayer]
	     << " => nCells = " << _nhits_per_layer[ilayer]
	     << " => nClusters = " << nclusters[ilayer] << endl;
      }
    
      cout << "nGtracks = " << _gtrack_list.size();
      cout << " => nTracks = ";
      if(_trackingWasRun) 
	{
	  cout << _trackList->size() << endl;
	}
      else
	{
	  cout << 0 << endl;
	}

      // cluster wise information
      if(verbosity > 1)
	{
 
	  unsigned int ig4hit=0;
	  for(HitMap::const_iterator g4hit_iter = _g4hitList.begin();
	      g4hit_iter != _g4hitList.end();
	      g4hit_iter++)
	    {
	      PHG4Hit *g4hit = g4hit_iter->second;

	      cout << endl;
	      cout << "===PHG4Hit===================================" << endl;
	      cout << " PHG4Hit: " << g4hit;
  
	      typedef multimap<PHG4Hit*,SvtxCluster*>::iterator mapiter2;
	      typedef pair<mapiter2,mapiter2> maprange2;
	      maprange2 therange2 = _allg4hits_cluster_mmap.equal_range( g4hit );
	      for(mapiter2 theiter2=therange2.first; theiter2!=therange2.second; theiter2++) 
		{
		  SvtxCluster *cluster = theiter2->second;	  
		  cout << "===Created-SvtxCluster================" << endl;      
		  cout << "SvtxCluster: "; cluster->identify();
		}

	      ++ig4hit;
	    } 
	}
      
      for(unsigned igtrack = 0; igtrack < _gtrack_list.size(); igtrack++)
	{
	  SvxGtrack *gtrack = &_gtrack_list[igtrack];

	  // don't print out the non-primary tracks
	  if (!gtrack->get_primary()) continue;
	  
	  // track-wise information
	  cout << endl;

	  cout << "=== Gtrack ===================================================" << endl;
	  cout << " Gtrack id = " << gtrack->get_track_id() << endl;
	  cout << " match = " << gtrack->get_standalone_match() << endl;
	  cout << " best purity = " << gtrack->get_best_purity() << endl;
	  cout << " best dp/p = " << gtrack->get_best_dpp() << endl;

	  cout << " PHG4Particle: ";
	  gtrack->get_particle()->identify(cout);
	  cout << " ptrue = (";
	  cout.width(5); cout << gtrack->get_px();
	  cout << ",";
	  cout.width(5); cout << gtrack->get_py();
	  cout << ",";
	  cout.width(5); cout << gtrack->get_pz();
	  cout << ")" << endl;

	  cout << " vtrue = (";
	  cout.width(5); cout << gtrack->get_vx();
	  cout << ",";
	  cout.width(5); cout << gtrack->get_vy();
	  cout << ",";
	  cout.width(5); cout << gtrack->get_vz();
	  cout << ")" << endl;
	  
	  cout << " pt = " << sqrt(pow(gtrack->get_px(),2)+pow(gtrack->get_py(),2)) << endl;
	  cout << " phi = " << atan2(gtrack->get_py(),gtrack->get_px()) << endl;
	  cout << " eta = " << asinh(gtrack->get_pz()/sqrt(pow(gtrack->get_px(),2)+pow(gtrack->get_py(),2))) << endl;
	  
	  cout << " chi^2 = " << gtrack->get_chisq() << endl;
	  cout << " chi^2 w/ vertex = " << gtrack->get_chisqv() << endl;
	  cout << " embed flag = " << gtrack->get_embed() << endl;
	  cout << " primary flag = " << gtrack->get_primary() << endl;
	  cout << " ---Associated-PHG4Hits-----------------------------------------" << endl;
	  
	  for(unsigned int ig4hit = 0; ig4hit < gtrack->get_ng4hits(); ig4hit++)
	    {
	      PHG4Hit *g4hit = gtrack->get_g4hit(ig4hit);
	      
	      float x = g4hit->get_x(0);
	      float y = g4hit->get_y(0);
	      float z = g4hit->get_z(0);
	      
	      cout << " #" << g4hit->get_hit_id() << " xtrue = (";
	      cout.width(5); cout << x;
	      cout << ",";
	      cout.width(5); cout << y;
	      cout << ",";
	      cout.width(5); cout << z;
	      cout << ")";

	      typedef multimap<PHG4Hit*,SvtxCluster*>::iterator mapiter2;
	      typedef pair<mapiter2,mapiter2> maprange2;
	      maprange2 therange2 = _allg4hits_cluster_mmap.equal_range( g4hit );
	      for(mapiter2 theiter2=therange2.first; theiter2!=therange2.second; theiter2++) 
		{
		  SvtxCluster *cluster = theiter2->second;

		  float x = cluster->get_x();
		  float y = cluster->get_y();
		  float z = cluster->get_z();
		 
		  cout << " => #" << cluster->get_id(); 
		  cout << " xreco = (";
		  cout.width(5); cout << x;
		  cout << ",";
		  cout.width(5); cout << y;
		  cout << ",";
		  cout.width(5); cout << z;
		  cout << ")";
		}

	      cout << endl;
	    }

	  if(_trackingWasRun)
	    {
	      typedef multimap<SvxGtrack*,SvtxTrack*>::iterator mapiter;
	      typedef pair<mapiter,mapiter> maprange;
	      maprange therange = _gtrack_track_mmap.equal_range( gtrack );
	      for(mapiter theiter=therange.first; theiter!=therange.second; theiter++) 
		{
		  SvtxTrack *track = theiter->second;

		  float px = track->get3Momentum(0);
		  float py = track->get3Momentum(1);
		  float pz = track->get3Momentum(2);

		  cout << "===Created-SvtxTrack==========================================" << endl;
		  cout << " SvtxTrack id = " << track->getTrackID() << endl;
		  cout << " preco = (";
		  cout.width(5); cout << px;
		  cout << ",";
		  cout.width(5); cout << py;
		  cout << ",";
		  cout.width(5); cout << pz;
		  cout << ")" << endl;
		  cout << " quality = " << track->getQuality() << endl;
		  cout << " purity = " << _track_purity_map[track] << endl;

		  cout << " ---Associated-SvtxClusters-to-PHG4Hits-------------------------" << endl;
	  
		  // loop over the associated track clusters
		  typedef multimap<SvtxTrack*,SvtxCluster*>::iterator mapiterS2C;
		  typedef pair<mapiterS2C,mapiterS2C> maprangeS2C;
		  maprangeS2C rangeS2C = _track_cluster_mmap.equal_range( track );
		  for(mapiterS2C iterS2C = rangeS2C.first; iterS2C != rangeS2C.second; iterS2C++) 
		    {
		      SvtxCluster *cluster = iterS2C->second;
		  
		      float x = cluster->get_x();
		      float y = cluster->get_y();
		      float z = cluster->get_z();
			  
		      cout << " #" << cluster->get_id() << " xreco = (";
		      cout.width(5); cout << x;
		      cout << ",";
		      cout.width(5); cout << y;
		      cout << ",";
		      cout.width(5); cout << z;
		      cout << ") =>";
			  
		      typedef multimap<SvtxCluster*,PHG4Hit*>::const_iterator mmapiter;
		      pair<mmapiter,mmapiter> hitrange = _cluster_allg4hits_mmap.equal_range(cluster);
		      for(mmapiter hititer = hitrange.first;
			  hititer!=hitrange.second;
			  hititer++) {
			PHG4Hit *g4hit = hititer->second;
			
			if (g4hit->get_trkid() == gtrack->get_track_id()) {
			  
			  if (g4hit) {
			    x = g4hit->get_x(0);
			    y = g4hit->get_y(0);
			    z = g4hit->get_z(0);
			    
			    cout << " #" << g4hit->get_hit_id()
				 << " xtrue = (";
			    cout.width(5); cout << x;
			    cout << ",";
			    cout.width(5); cout << y;
			    cout << ",";
			    cout.width(5); cout << z;
			    cout << ") => Gtrack id = " << g4hit->get_trkid();
			    break; // print only the first match
			  } else {
			    cout << " noise hit";
			  }
			} 
		      }
		      cout << endl;
		    }
		}
	    }
	}
      
      cout << endl;

    } // if verbosity

  return;
}

void PHG4Evaluator::fillOutputNtuples()
{
  if(verbosity > 1) cout << "PHG4Evaluator::fillOutputNtuples() entered" << endl;

  //----------------------
  // fill the Event NTuple
  //----------------------

  float vx = NAN;
  float vy = NAN;
  float vz = NAN;
  if (_vertexList) {
    if (!_vertexList->empty()) {
      SvtxVertex* vertex = &(_vertexList->begin()->second);
      
      vx = vertex->get_x();
      vy = vertex->get_y();
      vz = vertex->get_z();
    }
  }

  
  PHG4VtxPoint *gvertex = _truth_info_container->GetPrimaryVtx( _truth_info_container->GetPrimaryVertexIndex() );
  float gvx = gvertex->get_x();
  float gvy = gvertex->get_y();
  float gvz = gvertex->get_z();

  float ngtracks = _gtrack_list.size(); 
  float ng4hits = _g4hitList.size();

  float hit_occupancy_layer0 = NAN;
  float hit_occupancy_layer1 = NAN;
  float hit_occupancy_layer2 = NAN;
  float hit_occupancy_layer3 = NAN;
  float hit_occupancy_layer4 = NAN;
  float hit_occupancy_layer5 = NAN;
  if ((_nlayers > 0)&&(_nchannels_per_layer[0] > 0)) hit_occupancy_layer0 = _nhits_per_layer[0]*1.0 / _nchannels_per_layer[0];
  if ((_nlayers > 1)&&(_nchannels_per_layer[1] > 0)) hit_occupancy_layer1 = _nhits_per_layer[1]*1.0 / _nchannels_per_layer[1];
  if ((_nlayers > 2)&&(_nchannels_per_layer[2] > 0)) hit_occupancy_layer2 = _nhits_per_layer[2]*1.0 / _nchannels_per_layer[2];
  if ((_nlayers > 3)&&(_nchannels_per_layer[3] > 0)) hit_occupancy_layer3 = _nhits_per_layer[3]*1.0 / _nchannels_per_layer[3];
  if ((_nlayers > 4)&&(_nchannels_per_layer[4] > 0)) hit_occupancy_layer4 = _nhits_per_layer[4]*1.0 / _nchannels_per_layer[4];
  if ((_nlayers > 5)&&(_nchannels_per_layer[5] > 0)) hit_occupancy_layer5 = _nhits_per_layer[5]*1.0 / _nchannels_per_layer[5];

  float nclusters = _clusterList->size();
  float ntracks = 0.0;
  if(_trackingWasRun) ntracks = _trackList->size();

  float event_data[17] = {_ievent,
                          vx,
                          vy,
                          vz,
                          gvx,
                          gvy,
                          gvz,
                          ntracks,
                          ngtracks,
                          nclusters,
                          ng4hits,
			  hit_occupancy_layer0,
			  hit_occupancy_layer1,
			  hit_occupancy_layer2,
			  hit_occupancy_layer3,
			  hit_occupancy_layer4,
			  hit_occupancy_layer5
  };

  _ntp_event->Fill(event_data);

  //---------------------
  // fill the G4hit NTuple
  //---------------------
  
  for(HitMap::const_iterator g4hit_iter = _g4hitList.begin();
      g4hit_iter != _g4hitList.end();
      g4hit_iter++)
    {
      PHG4Hit *g4hit = g4hit_iter->second;
      map<PHG4Hit*,SvxGtrack*>::iterator finditer = _g4hit_gtrack_map.find(g4hit);
      if( finditer == _g4hit_gtrack_map.end() ) continue;
      SvxGtrack *gtrack = finditer->second;
  
      float g4hitID   = g4hit->get_hit_id();
      float x         = 0.5*(g4hit->get_x(0)+g4hit->get_x(1));
      float y         = 0.5*(g4hit->get_y(0)+g4hit->get_y(1));
      float z         = 0.5*(g4hit->get_z(0)+g4hit->get_z(1));
      float edep      = g4hit->get_edep();
      float layer     = g4hit->get_layer();
      float nclusters = _g4hit_cluster_mmap.count(g4hit);
      float gtrackID  = g4hit->get_trkid();
      float flavor    = gtrack->get_flavor();
      float px        = gtrack->get_px();
      float py        = gtrack->get_py();
      float pz        = gtrack->get_pz();
      float vx        = gtrack->get_vx();
      float vy        = gtrack->get_vy();
      float vz        = gtrack->get_vz();
      float fpx       = gtrack->get_fpx();
      float fpy       = gtrack->get_fpy();
      float fpz       = gtrack->get_fpz();
      float fx        = gtrack->get_fx();
      float fy        = gtrack->get_fy();
      float fz        = gtrack->get_fz();
      float last      = gtrack->get_is_last();
      float embed     = gtrack->get_embed();
      float primary   = gtrack->get_primary();

      float match     = 0.0;

      // if the g4hit produced N hits and all N hits were in the same cluster
      bool foundMatch = false;
      unsigned int ng4hitrawhits = _g4hit_cluster_mmap.count( g4hit );
   
      if(foundMatch)
	{
	  match = (float)ng4hitrawhits;
	}
      else
	{
	  match = -1.0*(float)ng4hitrawhits;
	} 

      float g4hit_data[26] = {_ievent,
			      g4hitID,
			      x,
			      y,
			      z,
			      edep,
			      layer,
			      nclusters,
			      gtrackID,
			      flavor,
			      px,
			      py,
			      pz,
			      vx,
			      vy,
			      vz,
			      fpx,
			      fpy,
			      fpz,
			      fx,
			      fy,
			      fz,
			      last,
			      embed,
			      primary,
			      match
      };

      _ntp_g4hit->Fill(g4hit_data);
    }

  //------------------------
  // fill the Cell NTuple
  //------------------------

  if (_cellList) {
    PHG4CylinderCellContainer::ConstRange cellrange = _cellList->getCylinderCells();
    for(PHG4CylinderCellContainer::ConstIterator celliter = cellrange.first;
	celliter != cellrange.second;
	++celliter)
      {
	PHG4CylinderCell* cell = celliter->second;
	PHG4Hit *g4hit = _cell_g4hit_mmap.find(cell)->second;
	map<PHG4Hit*,SvxGtrack*>::iterator finditer = _g4hit_gtrack_map.find(g4hit);
	if( finditer == _g4hit_gtrack_map.end() ) continue;
	SvxGtrack *gtrack = finditer->second;

	// _ntp_cell = new TNtuple("ntp_cell","cell-wise ntuple",
	// 			  "event:cellID:e:layer:"
	// 			  "g4hitID:gx:gy:gz:"
	// 			  "gtrackID:gflavor:"
	// 			  "gpx:gpy:gpz:gvx:gvy:gvz:"
	// 			  "gfpx:gfpy:gfpz:gfx:gfy:gfz:glast:"
	// 			  "ng4hits:purity");

	float event = _ievent;
	float cellID = cell->get_cell_id();
	float e = cell->get_edep();
	float layer = cell->get_layer();
	  
	float g4hitID   = -9999.9;
      
	float gedep    = -9999.9;

	float gx       = -9999.9;
	float gy       = -9999.9;
	float gz       = -9999.9;
	  
	float gtrackID = -9999.9;
	float gflavor  = -9999.9;
	  
	float gpx      = -9999.9;
	float gpy      = -9999.9;
	float gpz      = -9999.9;

	float gvx      = -9999.9;
	float gvy      = -9999.9;
	float gvz      = -9999.9;
	  
	float gfpx      = -9999.9;
	float gfpy      = -9999.9;
	float gfpz      = -9999.9;
	  
	float gfx      = -9999.9;
	float gfy      = -9999.9;
	float gfz      = -9999.9;

	float glast    = -9999.9;
	float gembed   = -9999.9;
	float gprimary = -9999.9;

	if(g4hit)
	  {
	    g4hitID   = g4hit->get_hit_id();
	      
	    gedep = g4hit->get_edep();

	    gx = g4hit->get_x(0);
	    gy = g4hit->get_y(0);
	    gz = g4hit->get_z(0);

	    gtrackID = gtrack->get_track_id();
	    gflavor  = gtrack->get_flavor();
	      
	    gpx      = gtrack->get_px();
	    gpy      = gtrack->get_py();
	    gpz      = gtrack->get_pz();
	      
	    gvx      = gtrack->get_vx();
	    gvy      = gtrack->get_vy();
	    gvz      = gtrack->get_vz();
	      
	    gfpx      = gtrack->get_fpx();
	    gfpy      = gtrack->get_fpy();
	    gfpz      = gtrack->get_fpz();
	      
	    gfx      = gtrack->get_fx();
	    gfy      = gtrack->get_fy();
	    gfz      = gtrack->get_fz();
	      
	    glast      = gtrack->get_is_last();
	    gembed     = gtrack->get_embed();
	    gprimary   = gtrack->get_primary();	  
	  }      

	float cell_data[26] = {
	  event,
	  cellID,
	  e,
	  layer,
	  g4hitID,
	  gedep,
	  gx,
	  gy,
	  gz,	  
	  gtrackID,
	  gflavor,	  
	  gpx,
	  gpy,
	  gpz,	    
	  gvx,
	  gvy,
	  gvz,	    
	  gfpx,
	  gfpy,
	  gfpz,	    
	  gfx,
	  gfy,
	  gfz,	    
	  glast,
	  gembed,
	  gprimary
	};
	  
	_ntp_cell->Fill(cell_data);     
      }
  }

  if (_ladderCellList) {
    PHG4CylinderCellContainer::ConstRange cellrange = _ladderCellList->getCylinderCells();
    for(PHG4CylinderCellContainer::ConstIterator celliter = cellrange.first;
	celliter != cellrange.second;
	++celliter)
      {
	PHG4CylinderCell* cell = celliter->second;
	PHG4Hit *g4hit = _cell_g4hit_mmap.find(cell)->second;
	map<PHG4Hit*,SvxGtrack*>::iterator finditer = _g4hit_gtrack_map.find(g4hit);
	if( finditer == _g4hit_gtrack_map.end() ) continue;
	SvxGtrack *gtrack = finditer->second;

	// _ntp_cell = new TNtuple("ntp_cell","cell-wise ntuple",
	// 			  "event:cellID:e:layer:"
	// 			  "g4hitID:gx:gy:gz:"
	// 			  "gtrackID:gflavor:"
	// 			  "gpx:gpy:gpz:gvx:gvy:gvz:"
	// 			  "gfpx:gfpy:gfpz:gfx:gfy:gfz:glast:"
	// 			  "ng4hits:purity");

	float event = _ievent;
	float cellID = cell->get_cell_id();
	float e = cell->get_edep();
	float layer = cell->get_layer();
	  
	float g4hitID   = -9999.9;
      
	float gedep    = -9999.9;

	float gx       = -9999.9;
	float gy       = -9999.9;
	float gz       = -9999.9;
	  
	float gtrackID = -9999.9;
	float gflavor  = -9999.9;
	  
	float gpx      = -9999.9;
	float gpy      = -9999.9;
	float gpz      = -9999.9;

	float gvx      = -9999.9;
	float gvy      = -9999.9;
	float gvz      = -9999.9;
	  
	float gfpx      = -9999.9;
	float gfpy      = -9999.9;
	float gfpz      = -9999.9;
	  
	float gfx      = -9999.9;
	float gfy      = -9999.9;
	float gfz      = -9999.9;

	float glast    = -9999.9;
	float gembed   = -9999.9;
	float gprimary = -9999.9;

	if(g4hit)
	  {
	    g4hitID   = g4hit->get_hit_id();
	      
	    gedep = g4hit->get_edep();

	    gx = g4hit->get_x(0);
	    gy = g4hit->get_y(0);
	    gz = g4hit->get_z(0);

	    gtrackID = gtrack->get_track_id();
	    gflavor  = gtrack->get_flavor();
	      
	    gpx      = gtrack->get_px();
	    gpy      = gtrack->get_py();
	    gpz      = gtrack->get_pz();
	      
	    gvx      = gtrack->get_vx();
	    gvy      = gtrack->get_vy();
	    gvz      = gtrack->get_vz();
	      
	    gfpx      = gtrack->get_fpx();
	    gfpy      = gtrack->get_fpy();
	    gfpz      = gtrack->get_fpz();
	      
	    gfx      = gtrack->get_fx();
	    gfy      = gtrack->get_fy();
	    gfz      = gtrack->get_fz();
	      
	    glast      = gtrack->get_is_last();
	    gembed     = gtrack->get_embed();
	    gprimary   = gtrack->get_primary();	  
	  }      

	float cell_data[26] = {
	  event,
	  cellID,
	  e,
	  layer,
	  g4hitID,
	  gedep,
	  gx,
	  gy,
	  gz,	  
	  gtrackID,
	  gflavor,	  
	  gpx,
	  gpy,
	  gpz,	    
	  gvx,
	  gvy,
	  gvz,	    
	  gfpx,
	  gfpy,
	  gfpz,	    
	  gfx,
	  gfy,
	  gfz,	    
	  glast,
	  gembed,
	  gprimary
	};
	  
	_ntp_cell->Fill(cell_data);     
      }
  }

  
  //------------------------
  // fill the Cluster NTuple
  //------------------------

  for (SvtxClusterMap::Iter iter = _clusterList->begin();
       iter != _clusterList->end();
       ++iter) {
    SvtxCluster* cluster = &iter->second;
    PHG4Hit *g4hit = _cluster_g4hit_mmap.find(cluster)->second;
    map<PHG4Hit*,SvxGtrack*>::iterator finditer = _g4hit_gtrack_map.find(g4hit);
    if( finditer == _g4hit_gtrack_map.end() ) continue;
    SvxGtrack *gtrack = finditer->second;
    //cout << "G4Hit: " << g4hit << " Gtrack: " << gtrack << endl;

    float hitID    = cluster->get_id();
    float x        = cluster->get_x();
    float y        = cluster->get_y();
    float z        = cluster->get_z();
    float e        = cluster->get_e();
    float adc      = cluster->get_adc();
    float layer    = cluster->get_layer();
    float size     = cluster->size_hits();
    float phisize  = cluster->get_phi_size();
    float zsize    = cluster->get_z_size();

    float g4hitID   = -9999.9;

    float gx       = -9999.9;
    float gy       = -9999.9;
    float gz       = -9999.9;

    float gtrackID = -9999.9;
    float gflavor  = -9999.9;

    float gpx      = -9999.9;
    float gpy      = -9999.9;
    float gpz      = -9999.9;

    float gvx      = -9999.9;
    float gvy      = -9999.9;
    float gvz      = -9999.9;

    float gfpx      = -9999.9;
    float gfpy      = -9999.9;
    float gfpz      = -9999.9;

    float gfx      = -9999.9;
    float gfy      = -9999.9;
    float gfz      = -9999.9;

    float glast    = -9999.9;
    float gembed   = -9999.9;
    float gprimary = -9999.9;

    float nhits    = -9999.9;
    float purity   = -9999.9;
      
    if(g4hit)
      {
	g4hitID   = g4hit->get_hit_id();

	gx = g4hit->get_x(0);
	gy = g4hit->get_y(0);
	gz = g4hit->get_z(0);

	gtrackID = gtrack->get_track_id();
	gflavor  = gtrack->get_flavor();

	gpx      = gtrack->get_px();
	gpy      = gtrack->get_py();
	gpz      = gtrack->get_pz();

	gvx      = gtrack->get_vx();
	gvy      = gtrack->get_vy();
	gvz      = gtrack->get_vz();

	gfpx      = gtrack->get_fpx();
	gfpy      = gtrack->get_fpy();
	gfpz      = gtrack->get_fpz();

	gfx      = gtrack->get_fx();
	gfy      = gtrack->get_fy();
	gfz      = gtrack->get_fz();

	glast      = gtrack->get_is_last();
	gembed     = gtrack->get_embed();
	gprimary   = gtrack->get_primary();
      }      

    // loop over the cluster => rawhits => g4hits to calculate the purity
    nhits  = 0.0;
    purity = 0.0;

    {
      typedef multimap<SvtxCluster*,PHG4Hit*>::const_iterator mmapiter;
      pair<mmapiter,mmapiter> hitrange = _cluster_allg4hits_mmap.equal_range(cluster);
      for(mmapiter hititer = hitrange.first;
	  hititer!=hitrange.second;
	  hititer++) {
	PHG4Hit *ig4hit = hititer->second;

	++nhits;

	if(ig4hit) // skip noise hits
	  {
	    if(ig4hit->get_trkid() == gtrackID) // if this g4hit came from the principle contributor truth particle increase the purity
	      {
		purity += 1.0;
	      }
	  }
      }
    }

    float cluster_data[34] = {_ievent,
			      hitID,
			      x,
			      y,
			      z,
			      e,
			      adc,
			      layer,
			      size,
			      phisize,
			      zsize,
			      g4hitID,
			      gx,
			      gy,
			      gz,
			      gtrackID,
			      gflavor,
			      gpx,
			      gpy,
			      gpz,
			      gvx,
			      gvy,
			      gvz,
			      gfpx,
			      gfpy,
			      gfpz,
			      gfx,
			      gfy,
			      gfz,
			      glast,
			      gembed,
			      gprimary,
			      nhits,
			      purity};

    _ntp_cluster->Fill(cluster_data);
  }		  

  //------------------------
  // fill the Gtrack NTuple
  //------------------------

  for(unsigned int igtrack = 0; igtrack < _gtrack_list.size(); igtrack++)
    {
      SvxGtrack *gtrack = &_gtrack_list[igtrack];

      float gtrackID  = gtrack->get_track_id();
      float flavor    = gtrack->get_flavor();

      float ng4hit0    = 0;
      float ng4hit1    = 0;
      float ng4hit2    = 0;
      float ng4hit3    = 0;
      for(unsigned int ig4hit = 0; ig4hit < gtrack->get_ng4hits(); ig4hit++)
	{
	  PHG4Hit *g4hit = gtrack->get_g4hit(ig4hit);
	  
	  if(g4hit->get_layer() == 0) ng4hit0++;
	  else if(g4hit->get_layer() == 1) ng4hit1++;
	  else if(g4hit->get_layer() == 2) ng4hit2++;
	  else if(g4hit->get_layer() == 3) ng4hit3++;
	}

      float px        = gtrack->get_px();
      float py        = gtrack->get_py();
      float pz        = gtrack->get_pz();
      float ntracks = _gtrack_track_mmap.count(gtrack);
      float vx        = gtrack->get_vx();
      float vy        = gtrack->get_vy();
      float vz        = gtrack->get_vz();

      float fpx       = gtrack->get_fpx();
      float fpy       = gtrack->get_fpy();
      float fpz       = gtrack->get_fpz();
      float fx        = gtrack->get_fx();
      float fy        = gtrack->get_fy();
      float fz        = gtrack->get_fz();
    
      float last      = gtrack->get_is_last();
      float embed    = gtrack->get_embed();
      float primary  = gtrack->get_primary();

      float chisq       = gtrack->get_chisq();
      float chisqv      = gtrack->get_chisqv();
      float match       = gtrack->get_standalone_match();
      float best_purity = gtrack->get_best_purity();
      float best_dpp    = gtrack->get_best_dpp();

      float gtrack_data[28] = {_ievent,
			       gtrackID,
			       flavor,
			       ng4hit0,
			       ng4hit1,

			       ng4hit2,
			       ng4hit3,
			       px,
			       py,
			       pz,

			       ntracks,
			       vx,
			       vy,
			       vz,
			       fpx,

			       fpy,
			       fpz,
			       fx,
			       fy,
			       fz,

			       last,
			       embed,
			       primary,
			       chisq,
			       chisqv,

			       match,
			       best_purity,
			       best_dpp
      };

      _ntp_gtrack->Fill(gtrack_data);
    }	  
	
  //------------------------
  // fill the Track NTuple
  //------------------------
  
  if(_trackingWasRun)
    {
      for (SvtxTrackMap::Iter iter = _trackList->begin();
	   iter != _trackList->end();
	   ++iter) {
	SvtxTrack *track = &iter->second;
	SvxGtrack  *gtrack  = _track_gtrack_map[track];

	float trackID   = track->getTrackID();
	float charge    = track->getCharge();
	float quality   = track->getQuality();
	float chisq     = track->getChisq();
	float chisqv    = track->getChisqv();
	float ndf       = track->getNDF();
	float primary   = track->getPrimary();
	float nhits     = track->getNhits();

	unsigned int layers = 0x0;
	for (unsigned int i = 0; i < 10; i++){	
	  if(track->hasCluster(i)) {
	    layers |= (0x1 << i);
	  }
	}

	float dedx1     = NAN;//track->get_dEdX1();
	float dedx2     = NAN;//track->get_dEdX2();
	float dca       = track->getDCA();
	float dca2d     = track->getDCA2D();
	float dca2dsigma = track->getDCA2Dsigma();
	float px        = track->get3Momentum(0);
	float py        = track->get3Momentum(1);
	float pz        = track->get3Momentum(2);
	float pcax      = track->d * sin(track->phi);
	float pcay      = track->d * cos(track->phi);
	float pcaz      = track->z0;
	float phi	  = track->phi;
	float d         = track->d;
	float kappa     = track->kappa;
	float z0        = track->z0;
	float dzdl      = track->dzdl;

	float presdphi = track->get_cal_dphi(SvtxTrack::PRES);
	float presdeta = track->get_cal_deta(SvtxTrack::PRES);
	float prese3x3 = track->get_cal_energy_3x3(SvtxTrack::PRES);
	float prese    = track->get_cal_cluster_e(SvtxTrack::PRES);

	float cemcdphi = track->get_cal_dphi(SvtxTrack::CEMC);
	float cemcdeta = track->get_cal_deta(SvtxTrack::CEMC);
	float cemce3x3 = track->get_cal_energy_3x3(SvtxTrack::CEMC);
	float cemce    = track->get_cal_cluster_e(SvtxTrack::CEMC);

	float hcalindphi = track->get_cal_dphi(SvtxTrack::HCALIN);
	float hcalindeta = track->get_cal_deta(SvtxTrack::HCALIN);
	float hcaline3x3 = track->get_cal_energy_3x3(SvtxTrack::HCALIN);
	float hcaline    = track->get_cal_cluster_e(SvtxTrack::HCALIN);

	float hcaloutdphi = track->get_cal_dphi(SvtxTrack::HCALOUT);
	float hcaloutdeta = track->get_cal_deta(SvtxTrack::HCALOUT);
	float hcaloute3x3 = track->get_cal_energy_3x3(SvtxTrack::HCALOUT);
	float hcaloute    = track->get_cal_cluster_e(SvtxTrack::HCALOUT);

	float gtrackID  = NAN;
	float gflavor   = NAN;
	float gpx       = NAN;
	float gpy       = NAN;
	float gpz       = NAN;
	float gvx       = NAN;
	float gvy       = NAN;
	float gvz       = NAN;
	float gfpx      = NAN;
	float gfpy      = NAN;
	float gfpz      = NAN;
	float gfx       = NAN;
	float gfy       = NAN;
	float gfz       = NAN;
	float glast     = NAN;
	float gembed    = NAN;
	float gprimary  = NAN;
	float purity    = 0.0;

	if(gtrack)
	  {
	    gtrackID  = gtrack->get_track_id();
	    gflavor   = gtrack->get_flavor();
	    gpx       = gtrack->get_px();
	    gpy       = gtrack->get_py();
	    gpz       = gtrack->get_pz();
	    gvx       = gtrack->get_vx();
	    gvy       = gtrack->get_vy();
	    gvz       = gtrack->get_vz();
	    gfpx       = gtrack->get_fpx();
	    gfpy       = gtrack->get_fpy();
	    gfpz       = gtrack->get_fpz();
	    gfx       = gtrack->get_fx();
	    gfy       = gtrack->get_fy();
	    gfz       = gtrack->get_fz();
	    glast      = gtrack->get_is_last();
	    gembed     = gtrack->get_embed();
	    gprimary   = gtrack->get_primary();
	    purity = _track_purity_map[track];
	  }

	float track_data[65] = {_ievent,
				trackID,
				charge,
				quality,
				chisq,
				chisqv,
				ndf,
				primary,
				nhits,
				(float)layers,
				dedx1,
				dedx2,
				dca,
				dca2d,
				dca2dsigma,
				px,
				py,
				pz,

				presdphi,
				presdeta,
				prese3x3,
				prese,
				  
				cemcdphi,
				cemcdeta,
				cemce3x3,
				cemce,
				  
				hcalindphi,
				hcalindeta,
				hcaline3x3,
				hcaline,
				  
				hcaloutdphi, 
				hcaloutdeta, 
				hcaloute3x3,
				hcaloute,    
				  
				gtrackID,
				gflavor,
				gpx,
				gpy,
				gpz,
				gvx,
				gvy,
				gvz,
				gfpx,
				gfpy,
				gfpz,
				gfx,
				gfy,
				gfz,
				glast,
				gembed,
				gprimary,
				purity,
				pcax,
				pcay,
				pcaz,
				phi,
				d,
				kappa,
				z0,
				dzdl };
      
	_ntp_track->Fill(track_data);
      }
    }

  return;
}

void PHG4Evaluator::fitGtrackProducedClusters() {

  //---------------------------------
  // fit the Gtrack produced clusters
  //---------------------------------

  float vertex[3] = {0.0,0.0,0.0};  
  float vertex_err[3] = {0.0,0.0,0.0};  

  if (_vertexList) {
    if (!_vertexList->empty()) {
      SvtxVertex* vtx = &(_vertexList->begin()->second);
      
      vertex[0] = vtx->get_x();
      vertex[1] = vtx->get_y();
      vertex[2] = vtx->get_z();

      vertex_err[0] = sqrt(vtx->get_error(0,0));
      vertex_err[1] = sqrt(vtx->get_error(1,1));
      vertex_err[2] = sqrt(vtx->get_error(2,2));
    }
  }
  

  
  // quick&dirty detector description
  vector<float> radii; // not used in fit -- used to construct the HelixHough object only
  radii.assign(_nlayers, 0.0);

  vector<float> smear_xy_layer;
  vector<float> smear_z_layer;
  smear_xy_layer.assign(_nlayers,0);
  smear_z_layer.assign(_nlayers,0);

  // these used in the fit
  float sqrt_12 = sqrt(12.);

  int ilayer = -1;
  PHG4CylinderCellGeomContainer::ConstRange layerrange = _cellGeos->get_begin_end();
  for (PHG4CylinderCellGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second;
       ++layeriter) {
    ++ilayer;

    float pitch = (layeriter->second)->get_phistep()*(layeriter->second)->get_radius();
    float length = (layeriter->second)->get_zstep();
    float radius = (layeriter->second)->get_radius();

    smear_xy_layer[ilayer] = (pitch/sqrt_12);
    smear_z_layer[ilayer] = (length/sqrt_12);
    radii[ilayer] = radius;
  }

  for (unsigned int igtrack = 0; igtrack < _gtrack_list.size(); ++igtrack) {

    SvxGtrack *gtrack = &_gtrack_list[igtrack];

    // require one of or matched or non-matched tracks
    // as these have the one ghit->one cluster dependence that we desire
    // to examine
    if (abs(gtrack->get_standalone_match()) < _min_hits) continue;
    if (abs(gtrack->get_standalone_match()) >= _max_hits) continue;

    // translate into a simple 3d track and provide a vector of chisqs
    SimpleTrack3D track;

    for (unsigned int ighit = 0; ighit < gtrack->get_ng4hits(); ++ighit) {
      PHG4Hit *ghit = gtrack->get_g4hit(ighit);

      typedef multimap<PHG4Hit*,SvtxCluster*>::iterator mapiter2;
      typedef pair<mapiter2,mapiter2> maprange2;
      maprange2 therange2 = _g4hit_cluster_mmap.equal_range( ghit );
      for (mapiter2 theiter2=therange2.first; theiter2!=therange2.second; theiter2++) {

	SvtxCluster *cluster = theiter2->second;

	float phi = atan2(cluster->get_y(),cluster->get_x());
	float xy_error = smear_xy_layer[cluster->get_layer()]*3.4641*0.5;
	float z_error  = smear_z_layer[cluster->get_layer()]*3.4641*0.5;

	SimpleHit3D hit(cluster->get_x() - vertex[0],
			fabs(xy_error*sin(phi)),
			cluster->get_y() - vertex[1],
			fabs(xy_error*cos(phi)),
			cluster->get_z() - vertex[2],
			z_error,
			cluster->get_id(),
			cluster->get_layer());
	  
	track.hits.push_back(hit);
      }
    }

    // fit without the vertex
    float chisq = sPHENIXTracker::fitTrack(track);
    gtrack->set_chisq(chisq);
    gtrack->set_chisqv(NAN);
  }

  return;
}
