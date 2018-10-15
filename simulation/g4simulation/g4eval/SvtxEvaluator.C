#include "SvtxEvaluator.h"

#include "SvtxEvalStack.h"

#include <phool/PHCompositeNode.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <phool/PHTimeServer.h>
#include <phool/PHTimer.h>

#include <g4hough/SvtxVertexMap.h>
#include <g4hough/SvtxVertex.h>
#include <g4hough/SvtxTrackMap.h>
#include <g4hough/SvtxTrack.h>
#include <g4hough/SvtxClusterMap.h>
#include <g4hough/SvtxCluster.h>
#include <g4hough/SvtxHitMap.h>
#include <g4hough/SvtxHit.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>

#include <TFile.h>
#include <TNtuple.h>
#include <TVector3.h>

#include <iostream>
#include <set>
#include <cmath>
#include <cassert>
#include <algorithm>

using namespace std;

SvtxEvaluator::SvtxEvaluator(const string &name, const string &filename, const string &trackmapname,
		unsigned int nlayers_maps,
		unsigned int nlayers_intt,
		unsigned int nlayers_tpc) :
  SubsysReco("SvtxEvaluator"),
  _ievent(0),
  _svtxevalstack(nullptr),
  _strict(false),
  _errors(0),
  _do_vertex_eval(true),
  _do_gpoint_eval(true),
  _do_g4hit_eval(true),
  _do_hit_eval(true),
  _do_cluster_eval(true),
  _do_gtrack_eval(true),
  _do_track_eval(true),
  _do_gseed_eval(false),
  _do_track_match(true),
  _do_eval_light(true),
  _scan_for_embedded(false),
  _nlayers_maps(nlayers_maps),
  _nlayers_intt(nlayers_intt),
  _nlayers_tpc(nlayers_tpc),
  _ntp_vertex(nullptr),
  _ntp_gpoint(nullptr),
  _ntp_g4hit(nullptr),
  _ntp_hit(nullptr),
  _ntp_cluster(nullptr),
  _ntp_gtrack(nullptr),
  _ntp_track(nullptr),
  _ntp_gseed(nullptr),
  _filename(filename),
  _trackmapname(trackmapname),
  _tfile(nullptr),
  _timer(nullptr)
{}

int SvtxEvaluator::Init(PHCompositeNode *topNode) {
  
  _ievent = 0;

  _tfile = new TFile(_filename.c_str(), "RECREATE");

  if (_do_vertex_eval) _ntp_vertex = new TNtuple("ntp_vertex","vertex => max truth",
                                                 "event:vx:vy:vz:ntracks:"
                                                 "gvx:gvy:gvz:gvt:gembed:gntracks:gntracksmaps:"
                                                 "gnembed:nfromtruth:"
						 "nhittpcall:nhittpcin:nhittpcmid:nhittpcout");

  if (_do_gpoint_eval) _ntp_gpoint = new TNtuple("ntp_gpoint","g4point => best vertex",
						 "event:gvx:gvy:gvz:gvt:gntracks:gembed:"
						 "vx:vy:vz:ntracks:"
						 "nfromtruth:"
						 "nhittpcall:nhittpcin:nhittpcmid:nhittpcout");
  
  if (_do_g4hit_eval) _ntp_g4hit = new TNtuple("ntp_g4hit","g4hit => best svtxcluster",
					       "event:g4hitID:gx:gy:gz:gt:gedep:geta:gphi:"
					       "gdphi:gdz:"
					       "glayer:gtrackID:gflavor:"
					       "gpx:gpy:gpz:"
					       "gvx:gvy:gvz:"
					       "gfpx:gfpy:gfpz:gfx:gfy:gfz:"
					       "gembed:gprimary:nclusters:"
					       "clusID:x:y:z:eta:phi:e:adc:layer:size:"
					       "phisize:zsize:efromtruth:dphitru:detatru:dztru:drtru:"
					       "nhittpcall:nhittpcin:nhittpcmid:nhittpcout");

  if (_do_hit_eval) _ntp_hit = new TNtuple("ntp_hit","svtxhit => max truth",
					   "event:hitID:e:adc:layer:"
					   "cellID:ecell:phibin:zbin:phi:z:"
					   "g4hitID:gedep:gx:gy:gz:gt:"
					   "gtrackID:gflavor:"
					   "gpx:gpy:gpz:gvx:gvy:gvz:gvt:"
					   "gfpx:gfpy:gfpz:gfx:gfy:gfz:"
					   "gembed:gprimary:efromtruth:"
					   "nhittpcall:nhittpcin:nhittpcmid:nhittpcout");

  if (_do_cluster_eval) _ntp_cluster = new TNtuple("ntp_cluster","svtxcluster => max truth",
						   "event:hitID:x:y:z:r:phi:eta:ex:ey:ez:ephi:"
						   "e:adc:layer:size:phisize:"
						   "zsize:trackID:g4hitID:gx:"
						   "gy:gz:gr:gphi:geta:gt:gtrackID:gflavor:"
						   "gpx:gpy:gpz:gvx:gvy:gvz:gvt:"
						   "gfpx:gfpy:gfpz:gfx:gfy:gfz:"
						   "gembed:gprimary:efromtruth:nparticles:"
						   "nhittpcall:nhittpcin:nhittpcmid:nhittpcout");

  if (_do_gtrack_eval) _ntp_gtrack  = new TNtuple("ntp_gtrack","g4particle => best svtxtrack",
						  "event:gtrackID:gflavor:gnhits:gnmaps:gnintt:gntpc:gnlmaps:gnlintt:gnltpc:"
						  "gpx:gpy:gpz:gpt:geta:gphi:"
						  "gvx:gvy:gvz:gvt:"
						  "gfpx:gfpy:gfpz:gfx:gfy:gfz:"
						  "gembed:gprimary:"
						  "trackID:px:py:pz:pt:eta:phi:"
						  "charge:quality:chisq:ndf:nhits:layers:nmaps:nintt:ntpc:nlmaps:nlintt:nltpc:"
						  "dca2d:dca2dsigma:dca3dxy:dca3dxysigma:dca3dz:dca3dzsigma:pcax:pcay:pcaz:nfromtruth:nwrong:ntrumaps:ntruintt:ntrutpc:layersfromtruth:"
						  "nhittpcall:nhittpcin:nhittpcmid:nhittpcout");
  
  if (_do_track_eval) _ntp_track = new TNtuple("ntp_track","svtxtrack => max truth",
					       "event:trackID:px:py:pz:pt:eta:phi:charge:"
					       "quality:chisq:ndf:nhits:nmaps:nintt:ntpc:nlmaps:nlintt:nltpc:layers:"
					       "dca2d:dca2dsigma:dca3dxy:dca3dxysigma:dca3dz:dca3dzsigma:pcax:pcay:pcaz:"
					       "presdphi:presdeta:prese3x3:prese:"   
					       "cemcdphi:cemcdeta:cemce3x3:cemce:"
					       "hcalindphi:hcalindeta:hcaline3x3:hcaline:"
					       "hcaloutdphi:hcaloutdeta:hcaloute3x3:hcaloute:"
					       "gtrackID:gflavor:gnhits:gnmaps:gnintt:gntpc:gnlmaps:gnlintt:gnltpc:"
					       "gpx:gpy:gpz:gpt:geta:gphi:"
					       "gvx:gvy:gvz:gvt:"
					       "gfpx:gfpy:gfpz:gfx:gfy:gfz:"
					       "gembed:gprimary:nfromtruth:nwrong:ntrumaps:ntruintt:ntrutpc:layersfromtruth:"
					       "nhittpcall:nhittpcin:nhittpcmid:nhittpcout");

  if (_do_gseed_eval) _ntp_gseed = new TNtuple("ntp_gseed","seeds from truth",
					       "event:ntrk:gx:gy:gz:gr:geta:gphi:"
					       "glayer:"
					       "gpx:gpy:gpz:gtpt:gtphi:gteta:"
					       "gvx:gvy:gvz:"
					       "gembed:gprimary:gflav:"
					       "dphiprev:detaprev:"
					       "nhittpcall:nhittpcin:nhittpcmid:nhittpcout");

  _timer = new PHTimer("_eval_timer");
  _timer->stop();

  return Fun4AllReturnCodes::EVENT_OK;
}

int SvtxEvaluator::InitRun(PHCompositeNode *topNode) {
  return Fun4AllReturnCodes::EVENT_OK;
}
  
int SvtxEvaluator::process_event(PHCompositeNode *topNode) {
  
  if ((Verbosity() > 0)&&(_ievent%100==0)) {
    cout << "SvtxEvaluator::process_event - Event = " << _ievent << endl;
  }

  if (!_svtxevalstack) {
    _svtxevalstack = new SvtxEvalStack(topNode);
    _svtxevalstack->set_strict(_strict);
    _svtxevalstack->set_verbosity(Verbosity()+1);
  } else {
    _svtxevalstack->next_event(topNode);
  }
  
  //-----------------------------------
  // print what is coming into the code
  //-----------------------------------
  
  printInputInfo(topNode);
  
  //---------------------------
  // fill the Evaluator NTuples
  //---------------------------

  fillOutputNtuples(topNode);
  
  //--------------------------------------------------
  // Print out the ancestry information for this event
  //--------------------------------------------------
  
  printOutputInfo(topNode);
  
  ++_ievent;
  return Fun4AllReturnCodes::EVENT_OK;
}

int SvtxEvaluator::End(PHCompositeNode *topNode) {

  _tfile->cd();

  if (_ntp_vertex)  _ntp_vertex->Write();
  if (_ntp_gpoint)  _ntp_gpoint->Write();
  if (_ntp_g4hit)   _ntp_g4hit->Write();
  if (_ntp_hit)     _ntp_hit->Write();
  if (_ntp_cluster) _ntp_cluster->Write();
  if (_ntp_gtrack)  _ntp_gtrack->Write();
  if (_ntp_track)   _ntp_track->Write();
  if (_ntp_gseed)   _ntp_gseed->Write();

  _tfile->Close();

  delete _tfile;

  if (Verbosity() >  0) {
    cout << "========================= SvtxEvaluator::End() ============================" << endl;
    cout << " " << _ievent << " events of output written to: " << _filename << endl;
    cout << "===========================================================================" << endl;
  }

  _errors += _svtxevalstack->get_errors();
  
  if (Verbosity() > -1) {
    if ((_errors > 0)||(Verbosity() > 0)) {
      cout << "SvtxEvaluator::End() - Error Count: " << _errors << endl;
    }
  }
  
  delete _svtxevalstack;
  
  return Fun4AllReturnCodes::EVENT_OK;
}

void SvtxEvaluator::printInputInfo(PHCompositeNode *topNode) {
  
  if (Verbosity() > 1) cout << "SvtxEvaluator::printInputInfo() entered" << endl;

  if (Verbosity() > 3) {
    
    // event information
    cout << endl;
    cout << PHWHERE << "   INPUT FOR EVENT " << _ievent << endl;

    cout << endl;
    cout << "---PHG4HITS-------------" << endl;
    _svtxevalstack->get_truth_eval()->set_strict(_strict);
    std::set<PHG4Hit*> g4hits = _svtxevalstack->get_truth_eval()->all_truth_hits();
    unsigned int ig4hit = 0;
    for(std::set<PHG4Hit*>::iterator iter = g4hits.begin();
	iter != g4hits.end();
	++iter) {
      PHG4Hit *g4hit = *iter;
      cout << ig4hit << " of " << g4hits.size();
      cout << ": PHG4Hit: " << endl;
      g4hit->identify();
      ++ig4hit;
    }

    cout << "---SVTXCLUSTERS-------------" << endl;
    SvtxClusterMap* clustermap = findNode::getClass<SvtxClusterMap>(topNode,"SvtxClusterMap");
    if (clustermap) {
      unsigned int icluster = 0;
      for (SvtxClusterMap::Iter iter = clustermap->begin();
	   iter != clustermap->end();
	   ++iter) {
	SvtxCluster* cluster = iter->second;
	cout << icluster << " of " << clustermap->size();	  
	cout << ": SvtxCluster: " << endl;
	cluster->identify();
	++icluster;
      }
    }

    cout << "---SVXTRACKS-------------" << endl;
    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode,_trackmapname.c_str());
    if (trackmap) {
      unsigned int itrack = 0;
      for (SvtxTrackMap::Iter iter = trackmap->begin();
	   iter != trackmap->end();
	   ++iter) {
	cout << itrack << " of " << trackmap->size();
	SvtxTrack *track = iter->second;
	cout << " : SvtxTrack:" << endl;
	track->identify();
	cout << endl;
      }
    }
    
    cout << "---SVXVERTEXES-------------" << endl;
    SvtxVertexMap* vertexmap = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMap");
    if (vertexmap) {
      unsigned int ivertex = 0;
      for (SvtxVertexMap::Iter iter = vertexmap->begin();
	   iter != vertexmap->end();
	   ++iter) {
	cout << ivertex << " of " << vertexmap->size();
	SvtxVertex *vertex = iter->second;
	cout << " : SvtxVertex:" << endl;
	vertex->identify();
	cout << endl;
      }
    }
  }

  return;
}

void SvtxEvaluator::printOutputInfo(PHCompositeNode *topNode) {
  
  if (Verbosity() > 1) cout << "SvtxEvaluator::printOutputInfo() entered" << endl;

  //==========================================
  // print out some useful stuff for debugging
  //==========================================

  if (Verbosity() > 0) {
    
    SvtxTrackEval*     trackeval = _svtxevalstack->get_track_eval();
    SvtxClusterEval* clustereval = _svtxevalstack->get_cluster_eval();
    SvtxTruthEval*     trutheval = _svtxevalstack->get_truth_eval();
  
    // event information
    cout << endl;
    cout << PHWHERE << "   NEW OUTPUT FOR EVENT " << _ievent << endl;
    cout << endl;

    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
    
    PHG4VtxPoint *gvertex = truthinfo->GetPrimaryVtx( truthinfo->GetPrimaryVertexIndex() );
    float gvx = gvertex->get_x();
    float gvy = gvertex->get_y();
    float gvz = gvertex->get_z();

    float vx = NAN;
    float vy = NAN;
    float vz = NAN;

    SvtxVertexMap* vertexmap = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMap");    
    if (vertexmap) {
      if (!vertexmap->empty()) {
	SvtxVertex* vertex = (vertexmap->begin()->second);
	
	vx = vertex->get_x();
	vy = vertex->get_y();
	vz = vertex->get_z();
      }
    }

    cout << "===Vertex Reconstruction=======================" << endl;
    cout << "vtrue = (" << gvx << "," << gvy << "," << gvz << ") => vreco = (" << vx << "," << vy << "," << vz << ")" << endl;
    cout << endl;

    cout << "===Tracking Summary============================" << endl;
    unsigned int ng4hits[100] = {0};
    std::set<PHG4Hit*> g4hits = trutheval->all_truth_hits();
    for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
	 iter != g4hits.end();
	 ++iter) {
      PHG4Hit *g4hit = *iter;
      ++ng4hits[g4hit->get_layer()];
    }

    SvtxHitMap* hitmap = findNode::getClass<SvtxHitMap>(topNode,"SvtxHitMap");
    unsigned int nhits[100] = {0};
    if (hitmap) {
      for (SvtxHitMap::Iter iter = hitmap->begin();
	   iter != hitmap->end();
	   ++iter) {
	SvtxHit* hit = iter->second;
	++nhits[hit->get_layer()];
      }
    }
    SvtxClusterMap* clustermap = findNode::getClass<SvtxClusterMap>(topNode,"SvtxClusterMap");
    unsigned int nclusters[100] = {0};
    if (clustermap) {
      for (SvtxClusterMap::Iter iter = clustermap->begin();
	   iter != clustermap->end();
	   ++iter) {
	SvtxCluster* cluster = iter->second;
	++nclusters[cluster->get_layer()];
      }
    }
    for (unsigned int ilayer = 0; ilayer < _nlayers_maps + _nlayers_intt+ _nlayers_tpc; ++ilayer) {
      cout << "layer " << ilayer << ": nG4hits = " << ng4hits[ilayer]
	   << " => nHits = " << nhits[ilayer]
	   << " => nClusters = " << nclusters[ilayer] << endl;
    }

    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode,_trackmapname.c_str());
    
    cout << "nGtracks = " << std::distance(truthinfo->GetPrimaryParticleRange().first,
					  truthinfo->GetPrimaryParticleRange().second);
    cout << " => nTracks = ";
    if (trackmap) cout << trackmap->size() << endl;
    else cout << 0 << endl;

    // cluster wise information
    if (Verbosity() > 1) {
 
      for(std::set<PHG4Hit*>::iterator iter = g4hits.begin();
	  iter != g4hits.end();
	  ++iter) {
	PHG4Hit *g4hit = *iter;

	cout << endl;
        cout << "===PHG4Hit===================================" << endl;
	cout << " PHG4Hit: "; g4hit->identify();

	std::set<SvtxCluster*> clusters = clustereval->all_clusters_from(g4hit);

	for (std::set<SvtxCluster*>::iterator jter = clusters.begin();
	     jter != clusters.end();
	     ++jter) {
	  SvtxCluster *cluster = *jter;
	  cout << "===Created-SvtxCluster================" << endl;      
	  cout << "SvtxCluster: "; cluster->identify();
	}
      }

      PHG4TruthInfoContainer::ConstRange range = truthinfo->GetPrimaryParticleRange();
      for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
	   iter != range.second; 
	   ++iter) {
	
	PHG4Particle *particle = iter->second;

	// track-wise information
	cout << endl;

	cout << "=== Gtrack ===================================================" << endl;
	cout << " PHG4Particle id = " << particle->get_track_id() << endl;
	particle->identify();
	cout << " ptrue = (";
	cout.width(5); cout << particle->get_px();
	cout << ",";
	cout.width(5); cout << particle->get_py();
	cout << ",";
	cout.width(5); cout << particle->get_pz();
	cout << ")" << endl;

	cout << " vtrue = (";
	cout.width(5); cout << truthinfo->GetVtx(particle->get_vtx_id())->get_x();
	cout << ",";
	cout.width(5); cout << truthinfo->GetVtx(particle->get_vtx_id())->get_y();
	cout << ",";
	cout.width(5); cout << truthinfo->GetVtx(particle->get_vtx_id())->get_z();
	cout << ")" << endl;
	  
	cout << " pt = " << sqrt(pow(particle->get_px(),2)+pow(particle->get_py(),2)) << endl;
	cout << " phi = " << atan2(particle->get_py(),particle->get_px()) << endl;
	cout << " eta = " << asinh(particle->get_pz()/sqrt(pow(particle->get_px(),2)+pow(particle->get_py(),2))) << endl;
	  
	cout << " embed flag = " << truthinfo->isEmbeded(particle->get_track_id()) << endl;

	cout << " ---Associated-PHG4Hits-----------------------------------------" << endl;
	std::set<PHG4Hit*> g4hits = trutheval->all_truth_hits(particle);
	for(std::set<PHG4Hit*>::iterator jter = g4hits.begin();
	    jter != g4hits.end();
	    ++jter) {
	  PHG4Hit *g4hit = *jter;

	  float x = 0.5*(g4hit->get_x(0)+g4hit->get_x(1));
	  float y = 0.5*(g4hit->get_y(0)+g4hit->get_y(1));
	  float z = 0.5*(g4hit->get_z(0)+g4hit->get_z(1));
	      
	  cout << " #" << g4hit->get_hit_id() << " xtrue = (";
	  cout.width(5); cout << x;
	  cout << ",";
	  cout.width(5); cout << y;
	  cout << ",";
	  cout.width(5); cout << z;
	  cout << ")";

	  std::set<SvtxCluster*> clusters = clustereval->all_clusters_from(g4hit);
	  for (std::set<SvtxCluster*>::iterator kter = clusters.begin();
	       kter != clusters.end();
	       ++kter) {
	  
	    SvtxCluster *cluster = *kter;

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

	if (trackmap&&clustermap) {

	  std::set<SvtxTrack*> tracks = trackeval->all_tracks_from(particle);
	  for (std::set<SvtxTrack*>::iterator jter = tracks.begin();
	       jter != tracks.end();
	       ++jter) {
	  
	    SvtxTrack *track = *jter;

	    float px = track->get_px();
	    float py = track->get_py();
	    float pz = track->get_pz();

	    cout << "===Created-SvtxTrack==========================================" << endl;
	    cout << " SvtxTrack id = " << track->get_id() << endl;
	    cout << " preco = (";
	    cout.width(5); cout << px;
	    cout << ",";
	    cout.width(5); cout << py;
	    cout << ",";
	    cout.width(5); cout << pz;
	    cout << ")" << endl;
	    cout << " quality = " << track->get_quality() << endl;
	    cout << " nfromtruth = " << trackeval->get_nclusters_contribution(track,particle) << endl;

	    cout << " ---Associated-SvtxClusters-to-PHG4Hits-------------------------" << endl;    

	    for (SvtxTrack::ConstClusterIter iter = track->begin_clusters();
		 iter != track->end_clusters();
		 ++iter) {
	      unsigned int cluster_id = *iter;
	      SvtxCluster* cluster = clustermap->get(cluster_id);
	      		  
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

	      PHG4Hit* g4hit = clustereval->max_truth_hit_by_energy(cluster);
	      if ((g4hit) && (g4hit->get_trkid() == particle->get_track_id())) {
			  
		x = 0.5*(g4hit->get_x(0)+g4hit->get_x(1));
		y = 0.5*(g4hit->get_y(0)+g4hit->get_y(1));
		z = 0.5*(g4hit->get_z(0)+g4hit->get_z(1));
			    
		cout << " #" << g4hit->get_hit_id()
		     << " xtrue = (";
		cout.width(5); cout << x;
		cout << ",";
		cout.width(5); cout << y;
		cout << ",";
		cout.width(5); cout << z;
		cout << ") => Gtrack id = " << g4hit->get_trkid();
	      } else {
		cout << " noise hit";
	      }
	    }
  
	    cout << endl;
	  }
	}
      }
    }
      
    cout << endl;

  } // if Verbosity()

  return;
}

void SvtxEvaluator::fillOutputNtuples(PHCompositeNode *topNode) {

  if (Verbosity() > 0) cout << "SvtxEvaluator::fillOutputNtuples() entered" << endl;

  SvtxVertexEval*   vertexeval = _svtxevalstack->get_vertex_eval();
  SvtxTrackEval*     trackeval = _svtxevalstack->get_track_eval();
  SvtxClusterEval* clustereval = _svtxevalstack->get_cluster_eval();
  SvtxHitEval*         hiteval = _svtxevalstack->get_hit_eval();
  SvtxTruthEval*     trutheval = _svtxevalstack->get_truth_eval();

  float nhit_tpc_all = 0;
  float nhit_tpc_in  = 0;
  float nhit_tpc_mid = 0;
  float nhit_tpc_out = 0;

  SvtxHitMap* hitmap_in = findNode::getClass<SvtxHitMap>(topNode,"SvtxHitMap");
  if (hitmap_in) {
    for (SvtxHitMap::Iter iter = hitmap_in->begin();
	 iter != hitmap_in->end();
	 ++iter) {
      nhit_tpc_all++;
      SvtxHit* hit             = iter->second;
      float layer  = hit->get_layer();
      if(layer== (float)_nlayers_maps+_nlayers_intt) nhit_tpc_in++;
      if(layer== (float)_nlayers_maps+_nlayers_intt+_nlayers_tpc - 1) nhit_tpc_out++;
      if(layer== (float)_nlayers_maps+_nlayers_intt+_nlayers_tpc/2 - 1) nhit_tpc_mid++;
    }
  }


  //-----------------------
  // fill the Vertex NTuple
  //-----------------------

  if (_ntp_vertex) {
    if (Verbosity() > 0){
      cout << "Filling ntp_vertex " << endl;
      cout << "start vertex time:                "<<_timer->get_accumulated_time()/1000. << " sec" <<endl;
      _timer->restart();
    }


    SvtxVertexMap* vertexmap = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMap");
    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
    if (vertexmap && truthinfo) {

      const auto prange = truthinfo->GetPrimaryParticleRange();
      map<int, unsigned int> embedvtxid_particle_count;
      map<int, unsigned int> embedvtxid_maps_particle_count;
      map<int, unsigned int> vertex_particle_count;

      if(_do_eval_light == false){
	for (auto iter = prange.first; iter != prange.second; ++iter) // process all primary paricle
	  {
	    const int point_id = iter->second->get_vtx_id();
	    int gembed = truthinfo->isEmbededVtx(iter->second->get_vtx_id());
	    ++vertex_particle_count[point_id];
	    ++embedvtxid_particle_count[gembed];
	    PHG4Particle* g4particle = iter->second;
	    
	    if (_scan_for_embedded && gembed <=0) continue;
	    
	    std::set<SvtxCluster*> g4clusters = clustereval->all_clusters_from(g4particle);
	    unsigned int nglmaps = 0;
	    unsigned int nglintt = 0;
	    unsigned int ngltpc  = 0;
	    
	    int lmaps[_nlayers_maps+1];
	    if(_nlayers_maps>0) for(unsigned int i = 0;i<_nlayers_maps;i++) lmaps[i] = 0;
	    
	    int lintt[_nlayers_intt+1];
	    if(_nlayers_intt>0) for(unsigned int i = 0;i<_nlayers_intt;i++) lintt[i] = 0;
	    
	    int ltpc[_nlayers_tpc+1];
	    if(_nlayers_tpc>0) for(unsigned int i = 0;i<_nlayers_tpc;i++) ltpc[i] = 0;
	    
	    for(const SvtxCluster* g4cluster : g4clusters){
	      unsigned int layer = g4cluster->get_layer();
	      //cout<<__LINE__<<": " << _ievent <<": " <<gtrackID << ": " << layer <<": " <<g4cluster->get_id() <<endl;
	      if(_nlayers_maps>0&&layer<_nlayers_maps) {
		lmaps[layer] = 1;
	      }
	      
	      if(_nlayers_intt>0&&layer>=_nlayers_maps&&layer<_nlayers_maps+_nlayers_intt){
		lintt[layer-_nlayers_maps] = 1;
	      }
	      
	      if(_nlayers_tpc>0&&layer>=_nlayers_maps+_nlayers_intt && layer<_nlayers_maps+_nlayers_intt+_nlayers_tpc){
		ltpc[layer-(_nlayers_maps+_nlayers_intt)] = 1;
	      }
	    }
	    if(_nlayers_maps>0) for(unsigned int i = 0;i<_nlayers_maps;i++) nglmaps+=lmaps[i];
	    if(_nlayers_intt>0) for(unsigned int i = 0;i<_nlayers_intt;i++) nglintt+=lintt[i];
	    if(_nlayers_tpc>0)  for(unsigned int i = 0;i<_nlayers_tpc;i++)  ngltpc+=ltpc[i];
	    
	    //        float gflavor   = g4particle->get_pid();
	    float gpx       = g4particle->get_px();
	    float gpy       = g4particle->get_py();
	    float gpz       = g4particle->get_pz();
	    float gpt       = NAN;
	    float geta      = NAN;
	    
	    if(gpx!=0&&gpy!=0){
	      TVector3 gv(gpx,gpy,gpz);
	      gpt  = gv.Pt();
	      geta = gv.Eta();
	      //          gphi = gv.Phi();
	    }
	    
	    if (nglmaps==3 && fabs(geta)<1.0 && gpt>0.5)
	      ++embedvtxid_maps_particle_count[gembed];
	  }
      }
      auto vrange = truthinfo->GetPrimaryVtxRange();
      map<int, bool> embedvtxid_found;
      map<int, int> embedvtxid_vertex_id;
      map<int, PHG4VtxPoint*> embedvtxid_vertex;
      for (auto iter = vrange.first; iter != vrange.second; ++iter) // process all primary vertexes
	{
	  const int point_id = iter->first;
	  int gembed = truthinfo->isEmbededVtx(point_id);
	  if (_scan_for_embedded && gembed <= 0) continue;
	  
	  auto search = embedvtxid_found.find(gembed);
	  if (search != embedvtxid_found.end())
	    {
	      embedvtxid_vertex_id[gembed] = point_id;
	      embedvtxid_vertex[gembed] = iter->second;
	    }
	  else
	    {
	      if (vertex_particle_count[embedvtxid_vertex_id[gembed]] < vertex_particle_count[point_id])
		{
		  embedvtxid_vertex_id[gembed] = point_id;
		  embedvtxid_vertex[gembed] = iter->second;
		}
	    }
	  embedvtxid_found[gembed]=false;
	}
      
      unsigned int ngembed=0;
      for (std::map<int,bool>::iterator iter = embedvtxid_found.begin();
        iter != embedvtxid_found.end();
        ++iter)
      {
        if (iter->first >= 0 || iter->first != iter->first) continue;
        ++ngembed;
      }

      for (SvtxVertexMap::Iter iter = vertexmap->begin();
	   iter != vertexmap->end();
	   ++iter) {
	SvtxVertex* vertex = iter->second;
	PHG4VtxPoint* point = vertexeval->max_truth_point_by_ntracks(vertex);
	
	float vx         = vertex->get_x();
	float vy         = vertex->get_y();
	float vz         = vertex->get_z();
	float ntracks    = vertex->size_tracks();

        float gvx        = NAN;
        float gvy        = NAN;
        float gvz        = NAN;
        float gvt        = NAN;
        float gembed     = NAN;
        float gntracks   = truthinfo->GetNumPrimaryVertexParticles();
        float gntracksmaps = NAN;
        float gnembed    = NAN;
        float nfromtruth = NAN;
	
	if (point) {
          const int point_id = point->get_id();
          gvx        = point->get_x();
          gvy        = point->get_y();
          gvz        = point->get_z();
          gvt        = point->get_t();
          gembed     = truthinfo->isEmbededVtx(point_id);
          gntracks   = embedvtxid_particle_count[(int)gembed];
          if (embedvtxid_maps_particle_count[(int)gembed]>0 && fabs(gvt)<2000.&& fabs(gvz)<13.0)
          gntracksmaps = embedvtxid_maps_particle_count[(int)gembed];
          gnembed    = (float) ngembed;
          nfromtruth = vertexeval->get_ntracks_contribution(vertex,point);
          embedvtxid_found[(int)gembed]= true;
	}
	  
        float vertex_data[18] = {(float) _ievent,
                                 vx,
                                 vy,
                                 vz,
                                 ntracks,
                                 gvx,
                                 gvy,
                                 gvz,
                                 gvt,
                                 gembed,
                                 gntracks,
                                 gntracksmaps,
                                 gnembed,
                                 nfromtruth,
				 nhit_tpc_all,
				 nhit_tpc_in,
				 nhit_tpc_mid,
				 nhit_tpc_out
	};
	      
	_ntp_vertex->Fill(vertex_data);      
      }

      if (!_scan_for_embedded) {
      for (std::map<int,bool>::iterator iter = embedvtxid_found.begin();
        iter != embedvtxid_found.end();
        ++iter)
      {
        if (embedvtxid_found[iter->first]) continue;

        float vx         = NAN;
        float vy         = NAN;
        float vz         = NAN;
        float ntracks    = NAN;

        float gvx        = NAN;
        float gvy        = NAN;
        float gvz        = NAN;
        float gvt        = NAN;
        float gembed     = iter->first;
        float gntracks   = NAN;
        float gntracksmaps = NAN;
        float gnembed    = NAN;
        float nfromtruth = NAN;

        PHG4VtxPoint* point = embedvtxid_vertex[gembed];

        if (point){
        const int point_id = point->get_id();
        gvx        = point->get_x();
        gvy        = point->get_y();
        gvz        = point->get_z();
        gvt        = point->get_t();
        gembed     = truthinfo->isEmbededVtx(point_id);
        gntracks   = embedvtxid_particle_count[(int)gembed];
        if (embedvtxid_maps_particle_count[(int)gembed]>0 && fabs(gvt)<2000&& fabs(gvz)<13.0)
        gntracksmaps = embedvtxid_maps_particle_count[(int)gembed];
        gnembed    = (float) ngembed;
//        nfromtruth = vertexeval->get_ntracks_contribution(vertex,point);
        }

        float vertex_data[18] = {(float) _ievent,
                                 vx,
                                 vy,
                                 vz,
                                 ntracks,
                                 gvx,
                                 gvy,
                                 gvz,
                                 gvt,
                                 gembed,
                                 gntracks,
                                 gntracksmaps,
                                 gnembed,
                                 nfromtruth,
				 nhit_tpc_all,
				 nhit_tpc_in,
				 nhit_tpc_mid,
				 nhit_tpc_out
        };

        _ntp_vertex->Fill(vertex_data);

      }

      }

    }
    if(Verbosity() >= 1){
      _timer->stop();
      cout << "vertex time:                "<<_timer->get_accumulated_time()/1000. << " sec" <<endl;
    }
  }
  
  //-----------------------
  // fill the gpoint NTuple
  //-----------------------

  if (_ntp_gpoint)
  {
    if (Verbosity() > 0)
    {
      cout << "Filling ntp_gpoint " << endl;
      _timer->restart();
    }
    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

    if (truthinfo)
    {
      auto vrange = truthinfo->GetPrimaryVtxRange();
      const auto prange = truthinfo->GetPrimaryParticleRange();

      map<int,unsigned int> vertex_particle_count;
      for (auto iter = prange.first; iter != prange.second; ++iter) // process all primary paricle
      {
        ++vertex_particle_count[iter->second->get_vtx_id()];
      }

      for (auto iter = vrange.first; iter != vrange.second; ++iter) // process all primary vertexes
      {
        const int point_id = iter->first;
        PHG4VtxPoint* point = iter->second;

        //      PHG4VtxPoint* point =  truthinfo->GetPrimaryVtx(truthinfo->GetPrimaryVertexIndex());

        if (point)
        {
          SvtxVertex* vertex = vertexeval->best_vertex_from(point);

          float gvx = point->get_x();
          float gvy = point->get_y();
          float gvz = point->get_z();
          float gvt = point->get_t();
          float gntracks = vertex_particle_count[point_id];

          float gembed = truthinfo->isEmbededVtx(point_id);
          float vx = NAN;
          float vy = NAN;
          float vz = NAN;
          float ntracks = NAN;
          float nfromtruth = NAN;

          if (vertex)
          {
            vx = vertex->get_x();
            vy = vertex->get_y();
            vz = vertex->get_z();
            ntracks = vertex->size_tracks();
            nfromtruth = vertexeval->get_ntracks_contribution(vertex, point);
          }

          float gpoint_data[] = {(float) _ievent,
                                   gvx,
                                   gvy,
                                   gvz,
                                   gvt,
                                   gntracks,
                                   gembed,
                                   vx,
                                   vy,
                                   vz,
                                   ntracks,
                                   nfromtruth,
				 nhit_tpc_all,
				 nhit_tpc_in,
				 nhit_tpc_mid,
				 nhit_tpc_out};
	  
          _ntp_gpoint->Fill(gpoint_data);
        }
      }
    }
    if (Verbosity() >= 1)
    {
      _timer->stop();
      cout << "gpoint time:                " << _timer->get_accumulated_time() / 1000. << " sec" << endl;
    }
  }
  
  //---------------------
  // fill the G4hit NTuple
  //---------------------

  if (_ntp_g4hit) {
    if (Verbosity() > 0) {cout << "Filling ntp_g4hit " << endl;_timer->restart();}
    std::set<PHG4Hit*> g4hits = trutheval->all_truth_hits();
    for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
	 iter != g4hits.end();
	 ++iter) {
            
      PHG4Hit *g4hit = *iter;
      PHG4Particle *g4particle = trutheval->get_particle(g4hit);
      
      float g4hitID   = g4hit->get_hit_id();
      float gx        = g4hit->get_avg_x();
      float gy        = g4hit->get_avg_y();
      float gz        = g4hit->get_avg_z();
      TVector3 vg4(gx,gy,gz);
      float gt        = g4hit->get_avg_t();
      TVector3 vin(g4hit->get_x(0),g4hit->get_y(0),g4hit->get_z(0));
      TVector3 vout(g4hit->get_x(1),g4hit->get_y(1),g4hit->get_z(1));
      float gdphi     = vin.DeltaPhi(vout);
      float gdz       = fabs(g4hit->get_z(1) - g4hit->get_z(0));
      float gedep     = g4hit->get_edep();
      float glayer    = g4hit->get_layer();
      float gtrackID  = g4hit->get_trkid();

      float gflavor   = NAN;
      float gpx       = NAN;
      float gpy       = NAN;
      float gpz       = NAN;
      TVector3 vec(g4hit->get_avg_x(),g4hit->get_avg_y(),g4hit->get_avg_z());
      float geta      = vec.Eta();
      float gphi      = vec.Phi();
      float gvx       = NAN;
      float gvy       = NAN;
      float gvz       = NAN;

      float gembed    = NAN;
      float gprimary  = NAN;

      float gfpx      = 0.;
      float gfpy      = 0.;
      float gfpz      = 0.;
      float gfx       = 0.;
      float gfy       = 0.;
      float gfz       = 0.;

      if (g4particle) {

	if (_scan_for_embedded) {
	  if (trutheval->get_embed(g4particle) <= 0) continue;
	}
	
	gflavor   = g4particle->get_pid();
	gpx       = g4particle->get_px();
	gpy       = g4particle->get_py();
	gpz       = g4particle->get_pz();

	PHG4VtxPoint* vtx = trutheval->get_vertex(g4particle);	

	if (vtx) {
	  gvx       = vtx->get_x();
	  gvy       = vtx->get_y();
	  gvz       = vtx->get_z();
	}
	PHG4Hit* outerhit = nullptr;
	if(_do_eval_light == false)
	  outerhit = trutheval->get_outermost_truth_hit(g4particle);	

	if (outerhit) {
	  gfpx      = outerhit->get_px(1);
	  gfpy      = outerhit->get_py(1);
	  gfpz      = outerhit->get_pz(1);
	  gfx       = outerhit->get_x(1);
	  gfy       = outerhit->get_y(1);
	  gfz       = outerhit->get_z(1);
	}
	
	gembed    = trutheval->get_embed(g4particle);
	gprimary  = trutheval->is_primary(g4particle);
      } //       if (g4particle)
      
      std::set<SvtxCluster*> clusters = clustereval->all_clusters_from(g4hit);  
      float nclusters = clusters.size();

      // best cluster reco'd
      //      SvtxCluster* cluster = nullptr;//clustereval->best_cluster_from(g4hit);
      SvtxCluster* cluster = clustereval->best_cluster_from(g4hit);

      float clusID     = NAN;
      float x          = NAN;
      float y          = NAN;
      float z          = NAN;
      float eta        = NAN;
      float phi        = NAN;
      float e          = NAN;
      float adc        = NAN;
      float layer      = NAN;
      float size       = NAN;
      float phisize    = NAN;
      float zsize      = NAN;
      float efromtruth = NAN;
      float dphitru    = NAN;
      float detatru    = NAN;
      float dztru      = NAN;
      float drtru      = NAN;

      if (cluster) {
	clusID     = cluster->get_id();
	x          = cluster->get_x();
	y          = cluster->get_y();
	z          = cluster->get_z();
	TVector3 vec2(x,y,z);
	eta        = vec2.Eta();
	phi        = vec2.Phi();
	e          = cluster->get_e();
	adc        = cluster->get_adc();
	layer      = cluster->get_layer();
	size       = cluster->size_hits();
	phisize    = cluster->get_phi_size();
	zsize      = cluster->get_z_size();
	dphitru    = vec2.DeltaPhi(vg4);
	detatru    = eta - geta;
	dztru      = z - gz;
	drtru	   = vec2.DeltaR(vg4);
	if (g4particle) {
	  efromtruth = clustereval->get_energy_contribution(cluster,g4particle);
	}
      }

      float g4hit_data[50] = {(float) _ievent,
			      g4hitID,
			      gx,
			      gy,
			      gz,
			      gt,
			      gedep,
			      geta,
			      gphi,
			      gdphi,
			      gdz,
			      glayer,
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
			      gembed,
			      gprimary,
			      nclusters,
			      clusID,  
			      x,      
			      y,      
			      z,  
			      eta,
			      phi,
			      e,      
			      adc,    
			      layer,  
			      size,   
			      phisize,
			      zsize,  
			      efromtruth,
			      dphitru,
			      detatru,
			      dztru,
			      drtru,
			      nhit_tpc_all,
			      nhit_tpc_in,
			      nhit_tpc_mid,
			      nhit_tpc_out
      };

      _ntp_g4hit->Fill(g4hit_data);
    }
    if(Verbosity() >= 1){
      _timer->stop();
      cout << "g4hit time:                "<<_timer->get_accumulated_time()/1000. << " sec" <<endl;
    }
  }
  
  //--------------------
  // fill the Hit NTuple
  //--------------------

  if (_ntp_hit) {
    if (Verbosity() > 0){ cout << "Filling ntp_hit " << endl;_timer->restart();}
    // need things off of the DST...
    SvtxHitMap* hitmap = findNode::getClass<SvtxHitMap>(topNode,"SvtxHitMap");
    PHG4CylinderCellGeomContainer* geom_container =
      findNode::getClass<PHG4CylinderCellGeomContainer>(topNode,"CYLINDERCELLGEOM_SVTX");
    if (!geom_container) {
      std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
      return;
    }
    
    if (hitmap) {

      for (SvtxHitMap::Iter iter = hitmap->begin();
	   iter != hitmap->end();
	   ++iter) {

	SvtxHit* hit             = iter->second;
	PHG4Hit* g4hit           = hiteval->max_truth_hit_by_energy(hit);
	PHG4Cell* g4cell = hiteval->get_cell(hit);
	PHG4Particle* g4particle = trutheval->get_particle(g4hit);

	float event  = _ievent;
	float hitID  = hit->get_id();
	float e      = hit->get_e();
	float adc    = hit->get_adc();
	float layer  = hit->get_layer();
	float cellID = hit->get_cellid();
	float ecell  = g4cell->get_edep();
	
	int phibin   = NAN;
 	int zbin     = NAN;
 	float phi    = NAN;
 	float z      = NAN;

	if(layer>=_nlayers_maps+_nlayers_intt){
	  PHG4CylinderCellGeom *GeoLayer = geom_container->GetLayerCellGeom(layer);
	  //"cellID:ecell:phibin:zbin:phi:z"
	  phibin = PHG4CellDefs::SizeBinning::get_phibin(g4cell->get_cellid());//cell->get_binphi();
	  zbin = PHG4CellDefs::SizeBinning::get_zbin(g4cell->get_cellid());//cell->get_binz();
	  phi = GeoLayer->get_phicenter( phibin );
	  z = GeoLayer->get_zcenter( zbin );
	}

	float g4hitID  = NAN;
	float gedep    = NAN;
	float gx       = NAN;
	float gy       = NAN;
	float gz       = NAN;
	float gt       = NAN;
	float gtrackID = NAN;
	float gflavor  = NAN;
	float gpx      = NAN;
	float gpy      = NAN;
	float gpz      = NAN;
	float gvx      = NAN;
	float gvy      = NAN;
	float gvz      = NAN;
	float gvt      = NAN;
	float gfpx     = NAN;
	float gfpy     = NAN;
	float gfpz     = NAN;
	float gfx      = NAN;
	float gfy      = NAN;
	float gfz      = NAN;
	float gembed   = NAN;
	float gprimary = NAN;
      
	float efromtruth = NAN;
      
	if (g4hit) {
	  g4hitID  = g4hit->get_hit_id();
	  gedep    = g4hit->get_edep();
	  gx       = g4hit->get_avg_x();
	  gy       = g4hit->get_avg_y();
	  gz       = g4hit->get_avg_z();
	  gt       = g4hit->get_avg_t();

	  if (g4particle) {

	    if (_scan_for_embedded) {
	      if (trutheval->get_embed(g4particle) <= 0) continue;
	    }
	    
	    gtrackID = g4particle->get_track_id();
	    gflavor  = g4particle->get_pid();
	    gpx      = g4particle->get_px();
	    gpy      = g4particle->get_py();
	    gpz      = g4particle->get_pz();

	    PHG4VtxPoint* vtx = trutheval->get_vertex(g4particle);

	    if (vtx) {
	      gvx      = vtx->get_x();
	      gvy      = vtx->get_y();
	      gvz      = vtx->get_z();
	      gvt      = vtx->get_t();
	    }

	    PHG4Hit* outerhit = NULL;
	    if(_do_eval_light == false)
 	      outerhit = trutheval->get_outermost_truth_hit(g4particle);
	    if (outerhit) {
	      gfpx     = outerhit->get_px(1);
	      gfpy     = outerhit->get_py(1);
	      gfpz     = outerhit->get_pz(1);
	      gfx      = outerhit->get_x(1);
	      gfy      = outerhit->get_y(1);
	      gfz      = outerhit->get_z(1);
	    }
	    gembed   = trutheval->get_embed(g4particle);
	    gprimary = trutheval->is_primary(g4particle);
	  } //   if (g4particle){
	}      

	if (g4particle) {
	  efromtruth = hiteval->get_energy_contribution(hit,g4particle);
	}

	float hit_data[40] = {
	  event,
	  hitID,
	  e,
	  adc,
	  layer,
	  cellID,
	  ecell,
	  (float) phibin,
	  (float) zbin,
	  phi,
	  z,
	  g4hitID,
	  gedep,
	  gx,
	  gy,
	  gz,
	  gt,
	  gtrackID,
	  gflavor,
	  gpx,
	  gpy,
	  gpz,
	  gvx,
	  gvy,
	  gvz,
	  gvt,
	  gfpx,
	  gfpy,
	  gfpz,
	  gfx,
	  gfy,
	  gfz,
	  gembed,
	  gprimary,
	  efromtruth,
	  nhit_tpc_all,
	  nhit_tpc_in,
	  nhit_tpc_mid,
	  nhit_tpc_out
	};
	  
	_ntp_hit->Fill(hit_data);     
      }
    }
    if(Verbosity() >= 1){
      _timer->stop();
      cout << "hit time:                "<<_timer->get_accumulated_time()/1000. << " sec" <<endl;
    }
  }
  
  //------------------------
  // fill the Cluster NTuple
  //------------------------

  if (Verbosity() > 0){ cout << "check for ntp_cluster" << endl;_timer->restart();}

  if (_ntp_cluster && !_scan_for_embedded) {
    if (Verbosity() > 0) cout << "Filling ntp_cluster (all of them) " << endl;
    // need things off of the DST...
    SvtxClusterMap* clustermap = findNode::getClass<SvtxClusterMap>(topNode,"SvtxClusterMap");
    if (clustermap) {

      for (SvtxClusterMap::Iter iter = clustermap->begin();
	   iter != clustermap->end();
	   ++iter) {
    
	SvtxCluster* cluster     = iter->second;   
	SvtxTrack* track         = trackeval->best_track_from(cluster);
       	PHG4Hit *g4hit           = clustereval->max_truth_hit_by_energy(cluster); 
	PHG4Particle *g4particle = trutheval->get_particle(g4hit);
    
	float hitID    = cluster->get_id();
	float x        = cluster->get_x();
	float y        = cluster->get_y();
	float z        = cluster->get_z();
	TVector3 pos(x,y,z);
	float r = pos.Perp();
	float phi = pos.Phi();
	float eta = pos.Eta();

	float ex       = sqrt(cluster->get_error(0,0));
	float ey       = sqrt(cluster->get_error(1,1));
	float ez       = cluster->get_z_error();

	float ephi     = cluster->get_rphi_error();
	
	float e        = cluster->get_e();
	float adc      = cluster->get_adc();
	float layer    = cluster->get_layer();
	float size     = cluster->size_hits();
	float phisize  = cluster->get_phi_size();
	float zsize    = cluster->get_z_size();

	float trackID  = NAN;
	if (track) trackID = track->get_id();

	float g4hitID  = NAN;
	float gx       = NAN;
	float gy       = NAN;
	float gz       = NAN;
	float gr       = NAN;
	float gphi     = NAN;
	float geta     = NAN;
	float gt       = NAN;
	float gtrackID = NAN;
	float gflavor  = NAN;
	float gpx      = NAN;
	float gpy      = NAN;
	float gpz      = NAN;
	float gvx      = NAN;
	float gvy      = NAN;
	float gvz      = NAN;
	float gvt      = NAN;
	float gfpx     = NAN;
	float gfpy     = NAN;
	float gfpz     = NAN;
	float gfx      = NAN;
	float gfy      = NAN;
	float gfz      = NAN;
	float gembed   = NAN;
	float gprimary = NAN;
    
	float efromtruth = NAN;

	if (g4hit) 
	  {
	    if(layer>=_nlayers_maps+_nlayers_intt)
	      {
		// This calculates the truth cluster position for the TPC from all of the contributing g4hits, typically 2-4 for the TPC
		// Complicated, since only the part of the energy that is collected within a layer contributes to the position
		//===============================================================================
		
		PHG4CylinderCellGeomContainer* geom_container =
		  findNode::getClass<PHG4CylinderCellGeomContainer>(topNode,"CYLINDERCELLGEOM_SVTX");
		if (!geom_container) 
		  {
		    std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
		    return;
		  }

		// radii of layer boundaries
		float rbin = NAN;
		float rbout = NAN;
		PHG4CylinderCellGeom *GeoLayer = geom_container->GetLayerCellGeom(layer);
		// get layer boundaries here (for nominal layer value) for later use
		rbin = GeoLayer->get_radius() - GeoLayer->get_thickness() / 2.0; 
		rbout = GeoLayer->get_radius() + GeoLayer->get_thickness() / 2.0; 	
 
		gx = 0.0; 
		gy = 0.0; 
		gz = 0.0; 
		gt = 0.0;
		float gwt = 0.0; 

		//cout << "Eval: cluster in layer " << layer << " rbin " << rbin << " rbout " << rbout << endl;		
		std::set<PHG4Hit*> truth_hits = clustereval->all_truth_hits(cluster);	  
		for (std::set<PHG4Hit*>::iterator iter = truth_hits.begin();
		     iter != truth_hits.end();
		     ++iter) 
		  {
		    PHG4Hit* this_g4hit = *iter;
		    
		    float rbegin = sqrt(this_g4hit->get_x(0)*this_g4hit->get_x(0) + this_g4hit->get_y(0)*this_g4hit->get_y(0));
		    float rend = sqrt(this_g4hit->get_x(1)*this_g4hit->get_x(1) + this_g4hit->get_y(1)*this_g4hit->get_y(1));
		    //cout << " Eval: g4hit " << this_g4hit->get_hit_id() <<  " rbegin " << rbegin << " rend " << rend << endl;

		    // make sure the entry point is at lower radius
		    float xl[2];
		    float yl[2];
		    float zl[2];

		    if(rbegin < rend)
		      {
			xl[0] = this_g4hit->get_x(0);
			yl[0] = this_g4hit->get_y(0);
			zl[0] = this_g4hit->get_z(0); 
			xl[1] = this_g4hit->get_x(1);
			yl[1] = this_g4hit->get_y(1);
			zl[1] = this_g4hit->get_z(1); 
		      }
		    else
		      {
			xl[0] = this_g4hit->get_x(1);
			yl[0] = this_g4hit->get_y(1);
			zl[0] = this_g4hit->get_z(1); 
			xl[1] = this_g4hit->get_x(0);
			yl[1] = this_g4hit->get_y(0);
			zl[1] = this_g4hit->get_z(0); 
			swap(rbegin,rend);
			//cout << "swapped in and out " << endl;
		      }

		    // check that the g4hit is not completely outside the cluster layer. Just skip this g4hit if it is
		    // this can happen because an electron moves across a layer boundary during drift and readout
		    // so the g4hit is recorded in the cell as contributing to that layer, even though it was outside the boundaries
		    if(rbegin < rbin && rend < rbin)
		      continue;
		    if(rbegin > rbout && rend > rbout)
		      continue;

		    float xin = xl[0];
		    float yin = yl[0];
		    float zin = zl[0];
		    float xout = xl[1];
		    float yout = yl[1];
		    float zout = zl[1];
		    
		    float t = NAN;
		    
		    if(rbegin < rbin)
		      {
			// line segment begins before boundary, find where it crosses
			t = line_circle_intersection(xl, yl, zl, rbin);
			if(t > 0)
			  {
			    xin = xl[0] + t * (xl[1]-xl[0]);
			    yin = yl[0] + t * (yl[1]-yl[0]);
			    zin = zl[0] + t * (zl[1]-zl[0]);
			  }
		      }
	    
		    if(rend > rbout)
		      {
			// line segment ends after boundary, find where it crosses
			t = line_circle_intersection(xl, yl, zl, rbout);
			if(t > 0)
			  {
			    xout = xl[0] + t * (xl[1]-xl[0]);
			    yout = yl[0] + t * (yl[1]-yl[0]);
			    zout = zl[0] + t * (zl[1]-zl[0]);
			  }
		      }
		    
		    // we want only the fraction of edep inside the layer	    
		    gx +=  (xin+xout) * 0.5 * this_g4hit->get_edep() * (xout-xin) / (xl[1]-xl[0]);
		    gy +=  (yin+yout) * 0.5 * this_g4hit->get_edep() * (yout-yin) / (yl[1]-yl[0]);
		    gz +=  (zin+zout) * 0.5 * this_g4hit->get_edep() * (zout-zin) / (zl[1]-zl[0]);
		    gt  += this_g4hit->get_avg_t() * this_g4hit->get_edep() * (zout-zin) / (zl[1]-zl[0]);
		    gwt +=  this_g4hit->get_edep() * (zout-zin) / (zl[1]-zl[0]);
		  }   // loop over this_g4hit
		gx /= gwt;  
		gy /= gwt; 
		gz /= gwt; 
		gt /= gwt;
		//cout << " weighted means: gx " << gx << " gy " << gy << " gz " << gz << endl;
	      }  // if TPC
	    else
	      {
		// not TPC, one g4hit per cluster
		gx = g4hit->get_avg_x();
		gy = g4hit->get_avg_y();
		gz = g4hit->get_avg_z();
	      }  // not TPC
	    
	    g4hitID  = g4hit->get_hit_id();
	    TVector3 gpos(gx,gy,gz);
	    gr = gpos.Perp();
	    gphi = gpos.Phi();
	    geta = gpos.Eta();
	
	    if (g4particle) {
	      
	      gtrackID = g4particle->get_track_id();
	      gflavor  = g4particle->get_pid();
	      gpx      = g4particle->get_px();
	      gpy      = g4particle->get_py();
	      gpz      = g4particle->get_pz();
	      
	      PHG4VtxPoint* vtx = trutheval->get_vertex(g4particle);
	      if (vtx) {
		gvx      = vtx->get_x();
		gvy      = vtx->get_y();
		gvz      = vtx->get_z();
		gvt      = vtx->get_t();
	      }
	      
	      PHG4Hit* outerhit = nullptr;
	      if(_do_eval_light == false)
		outerhit = trutheval->get_outermost_truth_hit(g4particle);	
	      if (outerhit) {
		gfpx     = outerhit->get_px(1);
		gfpy     = outerhit->get_py(1);
		gfpz     = outerhit->get_pz(1);
		gfx      = outerhit->get_x(1);
		gfy      = outerhit->get_y(1);
		gfz      = outerhit->get_z(1);
	      }
	      
	      gembed   = trutheval->get_embed(g4particle);
	      gprimary = trutheval->is_primary(g4particle);
	      
	    }      //   if (g4particle){
	  } //  if (g4hit) {
	
	if (g4particle){
	  efromtruth = clustereval->get_energy_contribution(cluster,g4particle);
	}
	
	float nparticles = clustereval->all_truth_particles(cluster).size();
	float cluster_data[50] = {(float) _ievent,
				  hitID,
				  x,
				  y,
				  z,
				  r,
				  phi,
				  eta,
				  ex,
				  ey,
				  ez,
				  ephi,
				  e,
				  adc,
				  layer,
				  size,
				  phisize,
				  zsize,
				  trackID,
				  g4hitID,
				  gx,
				  gy,
				  gz,
				  gr,
				  gphi,
				  geta,
				  gt,
				  gtrackID,
				  gflavor,
				  gpx,
				  gpy,
				  gpz,
				  gvx,
				  gvy,
				  gvz,
				  gvt,
				  gfpx,
				  gfpy,
				  gfpz,
				  gfx,
				  gfy,
				  gfz,
				  gembed,
				  gprimary,
				  efromtruth,
				  nparticles,
				  nhit_tpc_all,
				  nhit_tpc_in,
				  nhit_tpc_mid,
				  nhit_tpc_out
	};

	_ntp_cluster->Fill(cluster_data);
      }		  
    }
    
  } else if (_ntp_cluster && _scan_for_embedded) {

    if (Verbosity() > 0) cout << "Filling ntp_cluster (embedded only) " << endl;

    // if only scanning embedded signals, loop over all the tracks from
    // embedded particles and report all of their clusters, including those
    // from other sources (noise hits on the embedded track)
    
    // need things off of the DST...
    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode,_trackmapname.c_str());
    SvtxClusterMap* clustermap = findNode::getClass<SvtxClusterMap>(topNode,"SvtxClusterMap");
    if (trackmap) {
      
      for (SvtxTrackMap::Iter iter = trackmap->begin();
	   iter != trackmap->end();
	   ++iter) {
	
	SvtxTrack* track = iter->second;
	PHG4Particle* truth = trackeval->max_truth_particle_by_nclusters(track);
	if (truth) {	  
	  if (trutheval->get_embed(truth) <= 0) continue;
	}
	
	for (SvtxTrack::ConstClusterIter iter = track->begin_clusters();
	     iter != track->end_clusters();
	     ++iter) {

	  unsigned int cluster_id = *iter;
	  SvtxCluster* cluster = clustermap->get(cluster_id);

	  PHG4Hit *g4hit           = clustereval->max_truth_hit_by_energy(cluster); 
	  PHG4Particle *g4particle = trutheval->get_particle(g4hit);
    
	  float hitID    = cluster->get_id();
	  float x        = cluster->get_x();
	  float y        = cluster->get_y();
	  float z        = cluster->get_z();
	  TVector3 pos(x,y,z);
	  float r = pos.Perp();
	  float phi = pos.Phi();
	  float eta = pos.Eta();
	  float ex       = sqrt(cluster->get_error(0,0));
	  float ey       = sqrt(cluster->get_error(1,1));
	  float ez       = cluster->get_z_error();

	  float ephi     = cluster->get_rphi_error();
	  
	  float e        = cluster->get_e();
	  float adc      = cluster->get_adc();
	  float layer    = cluster->get_layer();
	  float size     = cluster->size_hits();
	  float phisize  = cluster->get_phi_size();
	  float zsize    = cluster->get_z_size();

	  float trackID  = NAN;
	  if (track) trackID = track->get_id();

	  float g4hitID  = NAN;
	  float gx       = NAN;
	  float gy       = NAN;
	  float gz       = NAN;
	  float gt       = NAN;
	  float gtrackID = NAN;
	  float gflavor  = NAN;
	  float gpx      = NAN;
	  float gpy      = NAN;
	  float gpz      = NAN;
	  float gvx      = NAN;
	  float gvy      = NAN;
	  float gvz      = NAN;
	  float gfpx     = NAN;
	  float gfpy     = NAN;
	  float gfpz     = NAN;
	  float gfx      = NAN;
	  float gfy      = NAN;
	  float gfz      = NAN;
	  float gembed   = NAN;
	  float gprimary = NAN;
    
	  float efromtruth = NAN;
      
	  if (g4hit) {
	    g4hitID  = g4hit->get_hit_id();
	    gx       = g4hit->get_avg_x();
	    gy       = g4hit->get_avg_y();
	    gz       = g4hit->get_avg_z();
	    gt       = g4hit->get_avg_t();

	    if (g4particle) {
	    
	      gtrackID = g4particle->get_track_id();
	      gflavor  = g4particle->get_pid();
	      gpx      = g4particle->get_px();
	      gpy      = g4particle->get_py();
	      gpz      = g4particle->get_pz();

	      PHG4VtxPoint* vtx = trutheval->get_vertex(g4particle);
	      if (vtx) {
		gvx      = vtx->get_x();
		gvy      = vtx->get_y();
		gvz      = vtx->get_z();
	      }
	      PHG4Hit* outerhit = nullptr;
	      if(_do_eval_light == false)
		outerhit = trutheval->get_outermost_truth_hit(g4particle);	
	      if (outerhit) {
		gfpx     = outerhit->get_px(1);
		gfpy     = outerhit->get_py(1);
		gfpz     = outerhit->get_pz(1);
		gfx      = outerhit->get_x(1);
		gfy      = outerhit->get_y(1);
		gfz      = outerhit->get_z(1);
	      }
	    
	      gembed   = trutheval->get_embed(g4particle);
	      gprimary = trutheval->is_primary(g4particle);
	    }      //   if (g4particle){
	  } //  if (g4hit) {

	  if (g4particle){
	    efromtruth = clustereval->get_energy_contribution(cluster,g4particle);
	  }

	  float nparticles = clustereval->all_truth_particles(cluster).size();

	  float cluster_data[46] = {(float) _ievent,
				    hitID,
				    x,
				    y,
				    z,
				    r,
				    phi,
				    eta,
				    ex,
				    ey,
				    ez,
				    ephi,
				    e,
				    adc,
				    layer,
				    size,
				    phisize,
				    zsize,
				    trackID,
				    g4hitID,
				    gx,
				    gy,
				    gz,
				    gt,
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
				    gembed,
				    gprimary,
				    efromtruth,
				    nparticles,
				    nhit_tpc_all,
				    nhit_tpc_in,
				    nhit_tpc_mid,
				    nhit_tpc_out
	  };

	  _ntp_cluster->Fill(cluster_data);
	}		  
      }
    }
  }
  if(Verbosity() >= 1){
    _timer->stop();
    cout << "cluster time:                "<<_timer->get_accumulated_time()/1000. << " sec" <<endl;
  }
  //------------------------
  // fill the Gtrack NTuple
  //------------------------

  // need things off of the DST...

  //cout << "check for ntp_gtrack" << endl;

  if (_ntp_gtrack) {
    if (Verbosity() > 0){ cout << "Filling ntp_gtrack " << endl;_timer->restart();}

    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");   
    SvtxClusterMap* clustermap = findNode::getClass<SvtxClusterMap>(topNode,"SvtxClusterMap");
    if (truthinfo) {

      PHG4TruthInfoContainer::ConstRange range = truthinfo->GetPrimaryParticleRange();
      for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
	   iter != range.second; 
	 ++iter) {
	
	PHG4Particle* g4particle = iter->second;

	if (_scan_for_embedded) {
	  if (trutheval->get_embed(g4particle) <= 0) continue;
	}
	
	float gtrackID = g4particle->get_track_id();
	float gflavor  = g4particle->get_pid();
      
	//std::set<PHG4Hit*> g4hits = trutheval->all_truth_hits(g4particle);
	std::set<SvtxCluster*> g4clusters = clustereval->all_clusters_from(g4particle);
    
	float ng4hits = g4clusters.size();
	unsigned int ngmaps  = 0;
	unsigned int ngintt  = 0;
	unsigned int ngtpc   = 0;
	unsigned int nglmaps = 0;
	unsigned int nglintt = 0;
	unsigned int ngltpc  = 0;

	int lmaps[_nlayers_maps+1]; 
	if(_nlayers_maps>0) for(unsigned int i = 0;i<_nlayers_maps;i++) lmaps[i] = 0;

	int lintt[_nlayers_intt+1];
	if(_nlayers_intt>0) for(unsigned int i = 0;i<_nlayers_intt;i++) lintt[i] = 0; 

	int ltpc[_nlayers_tpc+1];
	if(_nlayers_tpc>0) for(unsigned int i = 0;i<_nlayers_tpc;i++) ltpc[i] = 0; 

	for(const SvtxCluster* g4cluster : g4clusters){
	  unsigned int layer = g4cluster->get_layer();
	  //cout<<__LINE__<<": " << _ievent <<": " <<gtrackID << ": " << layer <<": " <<g4cluster->get_id() <<endl;
	  if(_nlayers_maps>0&&layer<_nlayers_maps) {
		  lmaps[layer] = 1;
		  ngmaps++;
	  }

	  if(_nlayers_intt>0&&layer>=_nlayers_maps&&layer<_nlayers_maps+_nlayers_intt){ 
	    lintt[layer-_nlayers_maps] = 1;
	    ngintt++;
	  }

	  if(_nlayers_tpc>0&&layer>=_nlayers_maps+_nlayers_intt && layer<_nlayers_maps+_nlayers_intt+_nlayers_tpc){ 
	    ltpc[layer-(_nlayers_maps+_nlayers_intt)] = 1;
	    ngtpc++;
	  }
	}
	if(_nlayers_maps>0) for(unsigned int i = 0;i<_nlayers_maps;i++) nglmaps+=lmaps[i];
	if(_nlayers_intt>0) for(unsigned int i = 0;i<_nlayers_intt;i++) nglintt+=lintt[i];
	if(_nlayers_tpc>0)  for(unsigned int i = 0;i<_nlayers_tpc;i++)  ngltpc+=ltpc[i];

	float gpx      = g4particle->get_px();
	float gpy      = g4particle->get_py();
	float gpz      = g4particle->get_pz();
	float gpt      = NAN;
	float geta     = NAN;
	float gphi     = NAN;
	if(gpx!=0&&gpy!=0){
	  TVector3 gv(gpx,gpy,gpz);
	  gpt  = gv.Pt();
	  geta = gv.Eta();
	  gphi = gv.Phi();
	}
	PHG4VtxPoint* vtx = trutheval->get_vertex(g4particle);	
	float gvx      = vtx->get_x();
	float gvy      = vtx->get_y();
	float gvz      = vtx->get_z();
	float gvt      = vtx->get_t();

	float gfpx      = 0.;
	float gfpy      = 0.;
	float gfpz      = 0.;
	float gfx       = 0.;
	float gfy       = 0.;
	float gfz       = 0.;
    
	PHG4Hit* outerhit = nullptr;
	if(_do_eval_light == false)
	  outerhit = trutheval->get_outermost_truth_hit(g4particle);	

	if (outerhit) {
	  gfpx      = outerhit->get_px(1);
	  gfpy      = outerhit->get_py(1);
	  gfpz      = outerhit->get_pz(1);
	  gfx       = outerhit->get_x(1);
	  gfy       = outerhit->get_y(1);
	  gfz       = outerhit->get_z(1);
	}
      
	float gembed   = trutheval->get_embed(g4particle);
	float gprimary = trutheval->is_primary(g4particle);

	float trackID       = NAN;
	float charge        = NAN;
	float quality       = NAN;
	float chisq         = NAN;
	float ndf           = NAN;
	float nhits         = NAN;
	float nmaps         =  0;
	float nintt         =  0;
	float ntpc          =  0;
	float nlintt        =  0;
	float nlmaps        =  0;
	float nltpc         =  0;
	unsigned int layers = 0x0;
	float dca2d         = NAN;
	float dca2dsigma    = NAN;
	float dca3dxy		 = NAN;
	float dca3dxysigma	 = NAN;
	float dca3dz		 = NAN;
	float dca3dzsigma	 = NAN;
	float px            = NAN;
	float py            = NAN;
	float pz            = NAN;
	float pt            = NAN;
	float eta           = NAN;
	float phi           = NAN;
	float pcax          = NAN;
	float pcay          = NAN;
	float pcaz          = NAN;

	float nfromtruth    = NAN;
	float nwrong        = NAN;
	float ntrumaps      = NAN;
	float ntruintt      = NAN;
	float ntrutpc       = NAN;
	float layersfromtruth = NAN;

	if(_do_track_match){
	  SvtxTrack* track = trackeval->best_track_from(g4particle);
	  
	  if (track) {
	    trackID   = track->get_id();     
	    charge    = track->get_charge();
	    quality   = track->get_quality();
	    chisq     = track->get_chisq();
	    ndf       = track->get_ndf();
	    nhits     = track->size_clusters();
	    int maps[_nlayers_maps];
	    int intt[_nlayers_intt];
	    int tpc[_nlayers_tpc];
	    
	    if(_nlayers_maps>0){
	      for(unsigned int i = 0;i<_nlayers_maps;i++) maps[i] = 0;
	    }
	    if(_nlayers_intt>0){
	      for(unsigned int i = 0;i<_nlayers_intt;i++) intt[i] = 0; 
	    }
	    if(_nlayers_tpc>0){
	      for(unsigned int i = 0;i<_nlayers_tpc;i++) tpc[i] = 0; 
	    }
	    
	    for (SvtxTrack::ConstClusterIter iter = track->begin_clusters();
		 iter != track->end_clusters();
		 ++iter) {
	      unsigned int cluster_id = *iter;
	      SvtxCluster* cluster = clustermap->get(cluster_id);
	      unsigned int layer = cluster->get_layer();
	      if(_nlayers_maps>0&&layer<_nlayers_maps){ maps[layer] = 1; nmaps++;}
	      if(_nlayers_intt>0&&layer>=_nlayers_maps&&layer<_nlayers_maps+_nlayers_intt){ intt[layer-_nlayers_maps] = 1;nintt++;}
	      if(_nlayers_tpc>0&&
		 layer>=(_nlayers_maps+_nlayers_intt)&&
		 layer<(_nlayers_maps+_nlayers_intt+_nlayers_tpc)){ 
		tpc[layer-(_nlayers_maps+_nlayers_intt)] = 1;ntpc++;
	      }
	    }
	    if(_nlayers_maps>0) for(unsigned int i = 0;i<_nlayers_maps;i++) nlmaps+=maps[i];
	    if(_nlayers_intt>0) for(unsigned int i = 0;i<_nlayers_intt;i++) nlintt+=intt[i];
	    if(_nlayers_tpc>0 ) for(unsigned int i = 0;i<_nlayers_tpc;i++)  nltpc+=tpc[i];
	    
	    layers = nlmaps+nlintt+nltpc;
	    
	    dca2d     = track->get_dca2d();
	    dca2dsigma = track->get_dca2d_error();
	    dca3dxy     = track->get_dca3d_xy();
	    dca3dxysigma = track->get_dca3d_xy_error();
	    dca3dz     = track->get_dca3d_z();
	    dca3dzsigma = track->get_dca3d_z_error();
	    px        = track->get_px();
	    py        = track->get_py();
	    pz        = track->get_pz();
	    TVector3 v(px,py,pz);
	    pt = v.Pt();
	    eta = v.Eta();
	    phi = v.Phi();
	    pcax      = track->get_x();
	    pcay      = track->get_y();
	    pcaz      = track->get_z();
	    
	    nfromtruth = trackeval->get_nclusters_contribution(track,g4particle);
	    nwrong     = trackeval->get_nwrongclusters_contribution(track,g4particle);
	    
	    if(_nlayers_maps==0){
	      ntrumaps   = 0;
	    }else{
	      ntrumaps   = trackeval->get_layer_range_contribution(track,g4particle,0,_nlayers_maps);
	    }
	    if(_nlayers_intt==0){
	      ntruintt   = 0;
	    }else{
	      ntruintt   = trackeval->get_layer_range_contribution(track,g4particle,_nlayers_maps,_nlayers_maps+_nlayers_intt);
	    }
	    ntrutpc    = trackeval->get_layer_range_contribution(track,g4particle,_nlayers_maps+_nlayers_intt,_nlayers_maps+_nlayers_intt+_nlayers_tpc);
	    
	    layersfromtruth = trackeval->get_nclusters_contribution_by_layer(track,g4particle);
	  }
	}
	float gtrack_data[66] = {(float) _ievent,
				 gtrackID,
				 gflavor,
				 ng4hits,
				 (float)ngmaps,
				 (float)ngintt,
				 (float)ngtpc,
				 (float)nglmaps,
				 (float)nglintt,
				 (float)ngltpc,
				 gpx,
				 gpy,
				 gpz,
				 gpt,
				 geta,
				 gphi,
				 gvx,
				 gvy,
				 gvz,
				 gvt,
				 gfpx,
				 gfpy,
				 gfpz,
				 gfx,
				 gfy,
				 gfz,
				 gembed,
				 gprimary,
				 trackID,
				 px,         
				 py,         
				 pz,  
				 pt,         
				 eta,         
				 phi,  
				 charge,     
				 quality,    
				 chisq,      
				 ndf,        
				 nhits,      
				 (float) layers,
				 nmaps,
				 nintt,
				 ntpc,
				 nlmaps,
				 nlintt,
				 nltpc,
				 dca2d,      
				 dca2dsigma,
				 dca3dxy,
				 dca3dxysigma,
				 dca3dz,
				 dca3dzsigma,
				 pcax,       
				 pcay,       
				 pcaz,
				 nfromtruth,
				 nwrong,
				 ntrumaps,
				 ntruintt,
				 ntrutpc,
				 layersfromtruth,
				 nhit_tpc_all,
				 nhit_tpc_in,
				 nhit_tpc_mid,
				 nhit_tpc_out
	};

	/*
	cout << " ievent " << _ievent
	     << " gtrackID " << gtrackID
	     << " gflavor " << gflavor
	     << " ng4hits " << ng4hits
	     << endl;
	*/

	_ntp_gtrack->Fill(gtrack_data);

      }	     
    }
    if(Verbosity() >= 1){
      _timer->stop();
      cout << "gtrack time:                "<<_timer->get_accumulated_time()/1000. << " sec" <<endl;
    }
  }
  
  //------------------------
  // fill the Track NTuple
  //------------------------



  if (_ntp_track) {
    if (Verbosity() > 0){ cout << "Filling ntp_track " << endl;_timer->restart();}

    // need things off of the DST...
    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode,_trackmapname.c_str());
    SvtxClusterMap* clustermap = findNode::getClass<SvtxClusterMap>(topNode,"SvtxClusterMap");
    if (trackmap) {

      for (SvtxTrackMap::Iter iter = trackmap->begin();
	   iter != trackmap->end();
	   ++iter) {
    
	SvtxTrack* track = iter->second;

	float trackID   = track->get_id();     
	float charge    = track->get_charge();
	float quality   = track->get_quality();
	float chisq     = track->get_chisq();
	float ndf       = track->get_ndf();
	float nhits     = track->size_clusters();
	unsigned int layers = 0x0;
	int maps[_nlayers_maps];
	int intt[_nlayers_intt];
	int tpc[_nlayers_tpc];
	if(_nlayers_maps>0){
	  for(unsigned int i = 0;i<_nlayers_maps;i++) maps[i] = 0;
	}
	if(_nlayers_intt>0){
	  for(unsigned int i = 0;i<_nlayers_intt;i++) intt[i] = 0; 
	}
	if(_nlayers_tpc>0){
	  for(unsigned int i = 0;i<_nlayers_tpc;i++) tpc[i] = 0; 
	}

	float nmaps  = 0;
	float nintt  = 0;
	float ntpc   = 0;
	float nlmaps = 0;
	float nlintt = 0;
	float nltpc  = 0;


	for (SvtxTrack::ConstClusterIter iter = track->begin_clusters();
	     iter != track->end_clusters();
	     ++iter) {
	  unsigned int cluster_id = *iter;
	  SvtxCluster* cluster = clustermap->get(cluster_id);
	  unsigned int layer = cluster->get_layer();

	  if(_nlayers_maps>0&&layer<_nlayers_maps){ maps[layer] = 1; nmaps++;}
	  if(_nlayers_intt>0&&layer>=_nlayers_maps&&layer<_nlayers_maps+_nlayers_intt){ intt[layer-_nlayers_maps] = 1;nintt++;}
	  if(_nlayers_tpc>0&&layer>=(_nlayers_maps+_nlayers_intt)&&layer<(_nlayers_maps+_nlayers_intt+_nlayers_tpc)){
	    tpc[layer - (_nlayers_maps+_nlayers_intt)] = 1;
	    ntpc++;
	  }
	}
	if(_nlayers_maps>0) for(unsigned int i = 0;i<_nlayers_maps;i++) nlmaps+=maps[i];
	if(_nlayers_intt>0) for(unsigned int i = 0;i<_nlayers_intt;i++) nlintt+=intt[i];
	if(_nlayers_tpc>0)  for(unsigned int i = 0;i<_nlayers_tpc;i++)  nltpc+=tpc[i];
	layers = nlmaps+nlintt+nltpc;
	float dca2d     = track->get_dca2d();
	float dca2dsigma = track->get_dca2d_error();
    float dca3dxy     = track->get_dca3d_xy();
    float dca3dxysigma = track->get_dca3d_xy_error();
    float dca3dz     = track->get_dca3d_z();
    float dca3dzsigma = track->get_dca3d_z_error();
	float px        = track->get_px();
	float py        = track->get_py();
	float pz        = track->get_pz();
	TVector3 v(px,py,pz);
	float pt        = v.Pt();
	float eta       = v.Eta();
	float phi       = v.Phi();
	float pcax      = track->get_x();
	float pcay      = track->get_y();
	float pcaz      = track->get_z();

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

	float gtrackID = NAN;
	float gflavor  = NAN;     
	float ng4hits  = NAN;
	unsigned int ngmaps   = 0;
	unsigned int ngintt   = 0;
	unsigned int ngtpc    = 0;
	unsigned int nglmaps  = 0;
	unsigned int nglintt  = 0;
	unsigned int ngltpc   = 0;
	float gpx      = NAN;
	float gpy      = NAN;
	float gpt      = NAN;
	float geta     = NAN;
	float gphi     = NAN;
	float gpz      = NAN;
	float gvx      = NAN;
	float gvy      = NAN;
	float gvz      = NAN;
	float gvt      = NAN;
	float gfpx     = NAN;
	float gfpy     = NAN;
	float gfpz     = NAN;
	float gfx      = NAN;
	float gfy      = NAN;
	float gfz      = NAN;
	float gembed   = NAN;
	float gprimary = NAN;

	float nfromtruth = NAN;
	float nwrong     = NAN;
	float ntrumaps   = NAN;
	float ntruintt   = NAN;
	float ntrutpc    = NAN;
	float layersfromtruth = NAN;
      
	if(_do_track_match){
	  PHG4Particle* g4particle = trackeval->max_truth_particle_by_nclusters(track);	
	  
	  if (g4particle) {
	    
	    if (_scan_for_embedded) {
	      if (trutheval->get_embed(g4particle) <= 0) continue;
	    }
	    
	    gtrackID = g4particle->get_track_id();
	    gflavor  = g4particle->get_pid();
	    
	    std::set<SvtxCluster*> g4clusters = clustereval->all_clusters_from(g4particle);
	    ng4hits = g4clusters.size();
	    gpx      = g4particle->get_px();
	    gpy      = g4particle->get_py();
	    gpz      = g4particle->get_pz();
	    
	    int lmaps[_nlayers_maps+1]; 
	    if(_nlayers_maps>0) for(unsigned int i = 0;i<_nlayers_maps;i++) lmaps[i] = 0;
	    
	    int lintt[_nlayers_intt+1];
	    if(_nlayers_intt>0) for(unsigned int i = 0;i<_nlayers_intt;i++) lintt[i] = 0; 
	    
	    int ltpc[_nlayers_tpc+1];
	    if(_nlayers_tpc>0) for(unsigned int i = 0;i<_nlayers_tpc;i++) ltpc[i] = 0; 
	    
	    for(const SvtxCluster* g4cluster : g4clusters){
	  	  unsigned int layer = g4cluster->get_layer();
	  	  if(_nlayers_maps>0&&layer<_nlayers_maps) {
	  		  lmaps[layer] = 1;
	  		  ngmaps++;
	  	  }

	  	  if(_nlayers_intt>0&&layer>=_nlayers_maps&&layer<_nlayers_maps+_nlayers_intt){
	  	    lintt[layer-_nlayers_maps] = 1;
	  	    ngintt++;
	  	  }

	  	  if(_nlayers_tpc>0&&layer>=_nlayers_maps+_nlayers_intt && layer<_nlayers_maps+_nlayers_intt+_nlayers_tpc){
	  	    ltpc[layer-(_nlayers_maps+_nlayers_intt)] = 1;
	  	    ngtpc++;
	  	  }
	    }
	    if(_nlayers_maps>0) for(unsigned int i = 0;i<_nlayers_maps;i++) nglmaps+=lmaps[i];
	    if(_nlayers_intt>0) for(unsigned int i = 0;i<_nlayers_intt;i++) nglintt+=lintt[i];
	    if(_nlayers_tpc>0)  for(unsigned int i = 0;i<_nlayers_tpc;i++)  ngltpc+=ltpc[i];
	    
	    TVector3 gv(gpx,gpy,gpz);
	    gpt      = gv.Pt();
	    geta     = gv.Eta();
	    gphi     = gv.Phi();
	    PHG4VtxPoint* vtx = trutheval->get_vertex(g4particle);	
	    gvx      = vtx->get_x();
	    gvy      = vtx->get_y();
	    gvz      = vtx->get_z();
	    gvt      = vtx->get_t();

	    PHG4Hit* outerhit = nullptr;
	    if(_do_eval_light == false)
	      outerhit = trutheval->get_outermost_truth_hit(g4particle);	
	    if (outerhit) {
	      gfpx     = outerhit->get_px(1);
	      gfpy     = outerhit->get_py(1);
	      gfpz     = outerhit->get_pz(1);
	      gfx      = outerhit->get_x(1);
	      gfy      = outerhit->get_y(1);
	      gfz      = outerhit->get_z(1);
	    }
	    gembed   = trutheval->get_embed(g4particle);
	    gprimary = trutheval->is_primary(g4particle);
	    
	    nfromtruth = trackeval->get_nclusters_contribution(track,g4particle);
	    nwrong     = trackeval->get_nwrongclusters_contribution(track,g4particle);
	    if(_nlayers_maps==0){
	      ntrumaps   = 0;
	    }else{
	      ntrumaps   = trackeval->get_layer_range_contribution(track,g4particle,0,_nlayers_maps);
	    }
	    if(_nlayers_intt==0){
	      ntruintt   = 0;
	    }else{
	      ntruintt   = trackeval->get_layer_range_contribution(track,g4particle,_nlayers_maps,_nlayers_maps+_nlayers_intt);
	    }
	    ntrutpc    = trackeval->get_layer_range_contribution(track,g4particle,_nlayers_maps+_nlayers_intt,_nlayers_maps+_nlayers_intt+_nlayers_tpc);
	    layersfromtruth = trackeval->get_nclusters_contribution_by_layer(track,g4particle);
	  }
	}
      
	float track_data[82] = {(float) _ievent,
				trackID, 
				px,        
				py,        
				pz,      
				pt,        
				eta,        
				phi, 
     				charge,  
				quality, 
				chisq,   
				ndf,     
				nhits,nmaps,nintt,ntpc,nlmaps,nlintt,nltpc,   
				(float) layers,
				dca2d,     
				dca2dsigma,
				dca3dxy,
				dca3dxysigma,
				dca3dz,
				dca3dzsigma,
				pcax,      
				pcay,      
				pcaz,      
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
				ng4hits,
				(float)ngmaps,
				(float)ngintt,
				(float)ngtpc,
				(float)nglmaps,
				(float)nglintt,
				(float)ngltpc,
				gpx,
				gpy,
				gpz,
				gpt,
				geta,
				gphi,
				gvx,
				gvy,
				gvz,
				gvt,
				gfpx,
				gfpy,
				gfpz,
				gfx,
				gfy,
				gfz,
				gembed,
				gprimary,
				nfromtruth,
				nwrong,
				ntrumaps,
				ntruintt,
				ntrutpc,
				layersfromtruth,
				nhit_tpc_all,
				nhit_tpc_in,
				nhit_tpc_mid,
				nhit_tpc_out
	};

	/*
	cout << "ievent " << _ievent
	     << " trackID " << trackID
	     << " nhits " << nhits
	     << " px " << px
	     << " py " << py
	     << " pz " << pz
	     << " gembed " << gembed
	     << " gprimary " << gprimary 
	     << endl;
	*/
	_ntp_track->Fill(track_data);
      }
    }
    if(Verbosity() >= 1){
      _timer->stop();
      cout << "track time:                "<<_timer->get_accumulated_time()/1000. << " sec" <<endl;
    }
  }

  //---------------------
  // fill the Gseed NTuple
  //---------------------

  if (_ntp_gseed ) {
    if (Verbosity() > 0) {cout << "Filling ntp_gseed " << endl;_timer->restart();}
    
    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
    
    float ntrk = 0;
    float gx = NAN;
    float gy = NAN;
    float gz = NAN;
    float gr = NAN;
    float geta = NAN;
    float gphi = NAN;
    float glayer = NAN;
    float gpx = NAN;
    float gpy = NAN;
    float gpz = NAN;
    float gtpt = NAN;
    float gtphi = NAN;
    float gteta = NAN;
    float gvx = NAN;
    float gvy = NAN;
    float gvz = NAN;
    float gembed = NAN;
    float gprimary = NAN;
    float gflav = NAN;
    float dphiprev = NAN;
    float detaprev = NAN;

    float xval[_nlayers_maps+_nlayers_intt+_nlayers_tpc];
    float yval[_nlayers_maps+_nlayers_intt+_nlayers_tpc];
    float zval[_nlayers_maps+_nlayers_intt+_nlayers_tpc];

    if (truthinfo) {
      PHG4TruthInfoContainer::ConstRange range = truthinfo->GetPrimaryParticleRange();
      for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
	   iter != range.second; 
	   ++iter) {
	ntrk++;
	PHG4Particle* g4particle = iter->second;
	for(unsigned int i = 0;i<_nlayers_maps+_nlayers_intt+_nlayers_tpc;i++){
	  xval[i] = 0;
	  yval[i] = 0;
	  zval[i] = 0;
	}
	std::set<PHG4Hit*> truth_hits = trutheval->all_truth_hits(g4particle);
	for (std::set<PHG4Hit*>::iterator iter = truth_hits.begin();
	     iter != truth_hits.end();
	     ++iter) {
	  PHG4Hit* g4hit = *iter;
	  int layer = g4hit->get_layer();
	  if (layer < 0)
	  {
	    cout << PHWHERE << " skipping negative detector id " << layer << endl;
	    continue;
	  }
	  xval[layer] = g4hit->get_avg_x();
	  yval[layer] = g4hit->get_avg_y();
	  zval[layer] = g4hit->get_avg_z();
	}
	
	for(unsigned int i = 0;i<_nlayers_maps+_nlayers_intt+_nlayers_tpc;i++){
	  gx = xval[i];
	  gy = yval[i];
	  gz = zval[i];
	  if(gx==0 && gy==0) continue;
	  
	  TVector3 vg4(gx,gy,gz);
	  glayer    = i;
	  gr = vg4.Perp();
	  geta = vg4.Eta();
	  gphi = vg4.Phi();
	  gpx       = g4particle->get_px();
	  gpy       = g4particle->get_py();
	  gpz       = g4particle->get_pz();
	  TVector3 vg4p(gpx,gpy,gpz);
	  
	  gtpt = vg4p.Perp();
	  gtphi = vg4p.Phi();
	  gteta = vg4p.Eta();
	  
	  PHG4VtxPoint* vtx = trutheval->get_vertex(g4particle);	
	  
	  if (vtx) {
	    gvx       = vtx->get_x();
	    gvy       = vtx->get_y();
	    gvz       = vtx->get_z();
	  }
	  
	  gembed    = trutheval->get_embed(g4particle);
	  gprimary  = trutheval->is_primary(g4particle);
	  gflav     = g4particle->get_pid();
	  if(i>=1){
	    if(xval[i-1]!=0&&yval[i-1]!=0){
	      TVector3 vg4prev(xval[i-1],yval[i-1],zval[i-1]);
	      dphiprev = vg4.DeltaPhi(vg4prev);
	      detaprev = geta - vg4prev.Eta();
	    }
	  }
	  
	
	  
	  float gseed_data[30] = {(float) _ievent,
				  ntrk,
				  gx,
				  gy,
				  gz,
				  gr,
				  geta,
				  gphi,
				  glayer,
				  gpx,
				  gpy,
				  gpz,
				  gtpt,
				  gtphi,
				  gteta,
				  gvx,
				  gvy,
				  gvz,
				  gembed,
				  gprimary,
				  gflav,
				  dphiprev,
				  detaprev,
				  nhit_tpc_all,
				  nhit_tpc_in,
				  nhit_tpc_mid,
				  nhit_tpc_out
	  };

	  
	  _ntp_gseed->Fill(gseed_data);
	}
      }
    }

    if(Verbosity() >= 1){
      _timer->stop();
      cout << "g4hit time:                "<<_timer->get_accumulated_time()/1000. << " sec" <<endl;
    }
  }
  
  return;
}

float SvtxEvaluator::line_circle_intersection(float x[], float y[], float z[], float radius) 
{
  // parameterize the line in terms of t (distance along the line segment, from 0-1) as
  // x = x0 + t * (x1-x0); y=y0 + t * (y1-y0); z = z0 + t * (z1-z0)
  // parameterize the cylinder (centered at x,y = 0,0) as  x^2 + y^2 = radius^2,   then
  // (x0 + t*(x1-z0))^2 + (y0+t*(y1-y0))^2 = radius^2
   // (x0^2 + y0^2 - radius^2) + (2x0*(x1-x0) + 2y0*(y1-y0))*t +  ((x1-x0)^2 + (y1-y0)^2)*t^2 = 0 = C + B*t + A*t^2
  // quadratic with:  A = (x1-x0)^2+(y1-y0)^2 ;  B = 2x0*(x1-x0) + 2y0*(y1-y0);  C = x0^2 + y0^2 - radius^2
  // solution: t = (-B +/- sqrt(B^2 - 4*A*C)) / (2*A) 
  
  float A = (x[1]-x[0])*(x[1]-x[0]) + (y[1]-y[0])*(y[1]-y[0]);
  float B = 2.0*x[0]*(x[1]-x[0]) + 2.0*y[0]*(y[1]-y[0]); 
  float C = x[0]*x[0] + y[0]*y[0] - radius*radius;
  float tup = (-B + sqrt(B*B - 4.0*A*C)) / (2.0*A); 
  float tdn = (-B - sqrt(B*B - 4.0*A*C)) / (2.0*A) ;

  // The limits are 0 and 1, but we allow a little for floating point precision
  float t;
  if(tdn >= -0.0e-4 && tdn <= 1.0004)
    t = tdn;
  else if(tup >= -0.0e-4 && tup <= 1.0004)
    t = tup;
  else
    {
      cout << PHWHERE << "   **** Oops! No valid solution for tup or tdn, tdn = " << tdn << " tup = " << tup << endl;
      cout << "   radius " << radius << " rbegin " << sqrt(x[0]*x[0]+y[0]*y[0]) << " rend " << sqrt(x[1]*x[1]+y[1]*y[1]) << endl;
      cout  << "   x0 " << x[0] << " x1 " << x[1] << endl;
      cout  << "   y0 " << y[0] << " y1 " << y[1] << endl;
      cout  << "   z0 " << z[0] << " z1 " << z[1] << endl;
      cout  << "   A " << A << " B " << B << " C " << C << endl;
      
      t = -1;
    }

  return t;


}
 
