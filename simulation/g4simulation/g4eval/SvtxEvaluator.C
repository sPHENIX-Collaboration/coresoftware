
#include "SvtxEvaluator.h"

#include "SvtxEvalStack.h"

#include <phool/PHCompositeNode.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

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

#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4detectors/PHG4CylinderCell.h>
#include <g4detectors/PHG4CylinderCellDefs.h>

#include <TNtuple.h>

#include <iostream>
#include <set>
#include <cmath>
#include <cassert>

using namespace std;

SvtxEvaluator::SvtxEvaluator(const string &name, const string &filename) :
  SubsysReco("SvtxEvaluator"),
  _ievent(0),
  _svtxevalstack(NULL),
  _strict(false),
  _errors(0),
  _do_vertex_eval(true),
  _do_gpoint_eval(true),
  _do_g4hit_eval(true),
  _do_hit_eval(true),
  _do_cluster_eval(true),
  _do_gtrack_eval(true),
  _do_track_eval(true),
  _ntp_vertex(NULL),
  _ntp_gpoint(NULL),
  _ntp_g4hit(NULL),
  _ntp_hit(NULL),
  _ntp_cluster(NULL),
  _ntp_gtrack(NULL),
  _ntp_track(NULL),
  _filename(filename),
  _tfile(NULL) {
  verbosity = 0;
}

int SvtxEvaluator::Init(PHCompositeNode *topNode) {
  
  _ievent = 0;

  _tfile = new TFile(_filename.c_str(), "RECREATE");

  if (_do_vertex_eval) _ntp_vertex = new TNtuple("ntp_vertex","vertex => max truth",
						 "event:vx:vy:vz:ntracks:"
						 "gvx:gvy:gvz:gntracks:"
						 "nfromtruth");

  if (_do_gpoint_eval) _ntp_gpoint = new TNtuple("ntp_gpoint","g4point => best vertex",
						 "event:gvx:gvy:gvz:gntracks:"
						 "vx:vy:vz:ntracks:"
						 "nfromtruth");
  
  if (_do_g4hit_eval) _ntp_g4hit = new TNtuple("ntp_g4hit","g4hit => best svtxcluster",
					       "event:g4hitID:gx:gy:gz:gedep:"
					       "glayer:gtrackID:gflavor:"
					       "gpx:gpy:gpz:gvx:gvy:gvz:"
					       "gfpx:gfpy:gfpz:gfx:gfy:gfz:"
					       "gembed:gprimary:nclusters:"
					       "clusID:x:y:z:e:adc:layer:size:"
					       "phisize:zsize:efromtruth");

  if (_do_hit_eval) _ntp_hit = new TNtuple("ntp_hit","svtxhit => max truth",
					   "event:hitID:e:adc:layer:"
					   "cellID:ecell:"
					   "g4hitID:gedep:gx:gy:gz:"
					   "gtrackID:gflavor:"
					   "gpx:gpy:gpz:gvx:gvy:gvz:"
					   "gfpx:gfpy:gfpz:gfx:gfy:gfz:"
					   "gembed:gprimary:efromtruth");

  if (_do_cluster_eval) _ntp_cluster = new TNtuple("ntp_cluster","svtxcluster => max truth",
						   "event:hitID:x:y:z:"
						   "e:adc:layer:size:phisize:"
						   "zsize:g4hitID:gx:"
						   "gy:gz:gtrackID:gflavor:"
						   "gpx:gpy:gpz:gvx:gvy:gvz:"
						   "gfpx:gfpy:gfpz:gfx:gfy:gfz:"
						   "gembed:gprimary:nhits:efromtruth");

  if (_do_gtrack_eval) _ntp_gtrack  = new TNtuple("ntp_gtrack","g4particle => best svtxtrack",
						  "event:gtrackID:gflavor:gnhits:"
						  "gpx:gpy:gpz:"
						  "gvx:gvy:gvz:"
						  "gfpx:gfpy:gfpz:gfx:gfy:gfz:"
						  "gembed:gprimary:"
						  "trackID:px:py:pz:charge:quality:chisq:ndf:nhits:layers:"
						  "dca2d:dca2dsigma:pcax:pcay:pcaz:nfromtruth");
  
  if (_do_track_eval) _ntp_track = new TNtuple("ntp_track","svtxtrack => max truth",
					       "event:trackID:px:py:pz:charge:"
					       "quality:chisq:ndf:nhits:layers:"
					       "dca2d:dca2dsigma:pcax:pcay:pcaz:"
					       "presdphi:presdeta:prese3x3:prese:"   
					       "cemcdphi:cemcdeta:cemce3x3:cemce:"
					       "hcalindphi:hcalindeta:hcaline3x3:hcaline:"
					       "hcaloutdphi:hcaloutdeta:hcaloute3x3:hcaloute:"
					       "gtrackID:gflavor:gnhits:"
					       "gpx:gpy:gpz:"
					       "gvx:gvy:gvz:"
					       "gfpx:gfpy:gfpz:gfx:gfy:gfz:"
					       "gembed:gprimary:nfromtruth");
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int SvtxEvaluator::InitRun(PHCompositeNode *topNode) {
  return Fun4AllReturnCodes::EVENT_OK;
}
  
int SvtxEvaluator::process_event(PHCompositeNode *topNode) {
  
  if ((verbosity > 0)&&(_ievent%100==0)) {
    cout << "SvtxEvaluator::process_event - Event = " << _ievent << endl;
  }

  if (!_svtxevalstack) {
    _svtxevalstack = new SvtxEvalStack(topNode);
    _svtxevalstack->set_strict(_strict);
    _svtxevalstack->set_verbosity(verbosity+1);
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

  _tfile->Close();

  delete _tfile;

  if (verbosity >  0) {
    cout << "========================= SvtxEvaluator::End() ============================" << endl;
    cout << " " << _ievent << " events of output written to: " << _filename << endl;
    cout << "===========================================================================" << endl;
  }

  _errors += _svtxevalstack->get_errors();
  
  if (verbosity > -1) {
    if ((_errors > 0)||(verbosity > 0)) {
      cout << "SvtxEvaluator::End() - Error Count: " << _errors << endl;
    }
  }
  
  delete _svtxevalstack;
  
  return Fun4AllReturnCodes::EVENT_OK;
}

void SvtxEvaluator::printInputInfo(PHCompositeNode *topNode) {
  
  if (verbosity > 1) cout << "SvtxEvaluator::printInputInfo() entered" << endl;

  if (verbosity > 3) {
    
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
    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
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
  
  if (verbosity > 1) cout << "SvtxEvaluator::printOutputInfo() entered" << endl;

  //==========================================
  // print out some useful stuff for debugging
  //==========================================

  if (verbosity > 0) {
    
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

    for (unsigned int ilayer = 0; ilayer < 100; ++ilayer) {
      cout << "layer " << ilayer << ": nG4hits = " << ng4hits[ilayer]
	   << " => nHits = " << nhits[ilayer]
	   << " => nClusters = " << nclusters[ilayer] << endl;
    }

    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
    
    cout << "nGtracks = " << std::distance(truthinfo->GetPrimaryParticleRange().first,
					  truthinfo->GetPrimaryParticleRange().second);
    cout << " => nTracks = ";
    if (trackmap) cout << trackmap->size() << endl;
    else cout << 0 << endl;

    // cluster wise information
    if (verbosity > 1) {
 
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

  } // if verbosity

  return;
}

void SvtxEvaluator::fillOutputNtuples(PHCompositeNode *topNode) {

  if (verbosity > 1) cout << "SvtxEvaluator::fillOutputNtuples() entered" << endl;

  SvtxVertexEval*   vertexeval = _svtxevalstack->get_vertex_eval();
  SvtxTrackEval*     trackeval = _svtxevalstack->get_track_eval();
  SvtxClusterEval* clustereval = _svtxevalstack->get_cluster_eval();
  SvtxHitEval*         hiteval = _svtxevalstack->get_hit_eval();
  SvtxTruthEval*     trutheval = _svtxevalstack->get_truth_eval();
  
  //-----------------------
  // fill the Vertex NTuple
  //-----------------------

  if (_ntp_vertex) {
    SvtxVertexMap* vertexmap = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMap");
    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
    if (vertexmap && truthinfo) {
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
	float gntracks   = truthinfo->GetNumPrimaryVertexParticles();
	float nfromtruth = NAN;
	
	if (point) {
	  gvx        = point->get_x();
	  gvy        = point->get_y();
	  gvz        = point->get_z();
	  gntracks   = truthinfo->GetNumPrimaryVertexParticles();
	  nfromtruth = vertexeval->get_ntracks_contribution(vertex,point);
	}
	  
	float vertex_data[10] = {(float) _ievent,
				 vx,
				 vy,
				 vz,
				 ntracks,
				 gvx,
				 gvy,
				 gvz,
				 gntracks,
				 nfromtruth
	};

	_ntp_vertex->Fill(vertex_data);      
      }
    }
  }
  
  //-----------------------
  // fill the gpoint NTuple
  //-----------------------

  if (_ntp_gpoint) {
    SvtxVertexMap* vertexmap = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMap");
    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
    if (vertexmap && truthinfo) {

      PHG4VtxPoint* point =  truthinfo->GetPrimaryVtx(truthinfo->GetPrimaryVertexIndex());
      SvtxVertex* vertex = vertexeval->best_vertex_from(point);
    
      float gvx        = point->get_x();
      float gvy        = point->get_y();
      float gvz        = point->get_z();
      float gntracks   = truthinfo->GetNumPrimaryVertexParticles();
      float vx         = NAN;
      float vy         = NAN;
      float vz         = NAN;
      float ntracks    = NAN;
      float nfromtruth = NAN;

      if (vertex) {
	vx         = vertex->get_x();
	vy         = vertex->get_y();
	vz         = vertex->get_z();
	ntracks    = vertex->size_tracks();
	nfromtruth = vertexeval->get_ntracks_contribution(vertex,point);
      }
	
      float gpoint_data[10] = {(float) _ievent,
			       gvx,
			       gvy,
			       gvz,
			       gntracks,
			       vx,
			       vy,
			       vz,
			       ntracks,
			       nfromtruth
      };

      _ntp_gpoint->Fill(gpoint_data);      
    }
  }
  
  //---------------------
  // fill the G4hit NTuple
  //---------------------

  if (_ntp_g4hit) {
    std::set<PHG4Hit*> g4hits = trutheval->all_truth_hits();
    for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
	 iter != g4hits.end();
	 ++iter) {
            
      PHG4Hit *g4hit = *iter;
      PHG4Particle *g4particle = trutheval->get_particle(g4hit);
      
      float g4hitID   = g4hit->get_hit_id();
      float gx        = 0.5*(g4hit->get_x(0)+g4hit->get_x(1));
      float gy        = 0.5*(g4hit->get_y(0)+g4hit->get_y(1));
      float gz        = 0.5*(g4hit->get_z(0)+g4hit->get_z(1));
      float gedep     = g4hit->get_edep();
      float glayer    = g4hit->get_layer();
  
      float gtrackID  = g4hit->get_trkid();


      float gflavor   = NAN;
      float gpx       = NAN;
      float gpy       = NAN;
      float gpz       = NAN;

      float gvx       = NAN;
      float gvy       = NAN;
      float gvz       = NAN;

      float gembed  = NAN;
      float gprimary = NAN;

      float gfpx      = 0.;
      float gfpy      = 0.;
      float gfpz      = 0.;
      float gfx       = 0.;
      float gfy       = 0.;
      float gfz       = 0.;

      if (g4particle) {

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
    
	PHG4Hit* outerhit = trutheval->get_outermost_truth_hit(g4particle);	
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
      SvtxCluster* cluster = clustereval->best_cluster_from(g4hit);

      float clusID     = NAN;
      float x          = NAN;
      float y          = NAN;
      float z          = NAN;
      float e          = NAN;
      float adc        = NAN;
      float layer      = NAN;
      float size       = NAN;
      float phisize    = NAN;
      float zsize      = NAN;
      float efromtruth = NAN;

      if (cluster) {
	clusID     = cluster->get_id();
	x          = cluster->get_x();
	y          = cluster->get_y();
	z          = cluster->get_z();
	e          = cluster->get_e();
	adc        = cluster->get_adc();
	layer      = cluster->get_layer();
	size       = cluster->size_hits();
	phisize    = cluster->get_phi_size();
	zsize      = cluster->get_z_size();
	if (g4particle) {
	  efromtruth = clustereval->get_energy_contribution(cluster,g4particle);
	}
      }

      float g4hit_data[35] = {(float) _ievent,
			      g4hitID,
			      gx,
			      gy,
			      gz,
			      gedep,
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
			      e,      
			      adc,    
			      layer,  
			      size,   
			      phisize,
			      zsize,  
			      efromtruth
      };

      _ntp_g4hit->Fill(g4hit_data);
    }
  }
  
  //--------------------
  // fill the Hit NTuple
  //--------------------

  if (_ntp_hit) {
    // need things off of the DST...
    SvtxHitMap* hitmap = findNode::getClass<SvtxHitMap>(topNode,"SvtxHitMap");
    if (hitmap) {

      for (SvtxHitMap::Iter iter = hitmap->begin();
	   iter != hitmap->end();
	   ++iter) {

	SvtxHit* hit             = iter->second;
	PHG4Hit* g4hit           = hiteval->max_truth_hit_by_energy(hit);
	PHG4CylinderCell* g4cell = hiteval->get_cell(hit);
	PHG4Particle* g4particle = trutheval->get_particle(g4hit);

	float event  = _ievent;
	float hitID  = hit->get_id();
	float e      = hit->get_e();
	float adc    = hit->get_adc();
	float layer  = hit->get_layer();
	float cellID = hit->get_cellid();
	float ecell  = g4cell->get_edep();

	float g4hitID  = NAN;
	float gedep    = NAN;
	float gx       = NAN;
	float gy       = NAN;
	float gz       = NAN;
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
	  gedep    = g4hit->get_edep();
	  gx       = 0.5*(g4hit->get_x(0)+g4hit->get_x(1));
	  gy       = 0.5*(g4hit->get_y(0)+g4hit->get_y(1));
	  gz       = 0.5*(g4hit->get_z(0)+g4hit->get_z(1));
	  gx       = g4hit->get_x(0);
	  gy       = g4hit->get_y(0);
	  gz       = g4hit->get_z(0);

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

	    PHG4Hit* outerhit = trutheval->get_outermost_truth_hit(g4particle);	
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

	float hit_data[32] = {
	  event,
	  hitID,
	  e,
	  adc,
	  layer,
	  cellID,
	  ecell,
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
	  gembed,
	  gprimary,
	  efromtruth
	};
	  
	_ntp_hit->Fill(hit_data);     
      }
    }
  }
  
  //------------------------
  // fill the Cluster NTuple
  //------------------------

  if (_ntp_cluster) {
    // need things off of the DST...
    SvtxClusterMap* clustermap = findNode::getClass<SvtxClusterMap>(topNode,"SvtxClusterMap");
    if (clustermap) {

      for (SvtxClusterMap::Iter iter = clustermap->begin();
	   iter != clustermap->end();
	   ++iter) {
    
	SvtxCluster* cluster     = iter->second;   
	PHG4Hit *g4hit           = clustereval->max_truth_hit_by_energy(cluster); 
	PHG4Particle *g4particle = trutheval->get_particle(g4hit);
    
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

	float g4hitID  = NAN;
	float gx       = NAN;
	float gy       = NAN;
	float gz       = NAN;
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
    
	float nhits      = NAN;
	float efromtruth = NAN;
      
	if (g4hit) {
	  g4hitID  = g4hit->get_hit_id();
	  gx       = 0.5*(g4hit->get_x(0)+g4hit->get_x(1));
	  gy       = 0.5*(g4hit->get_y(0)+g4hit->get_y(1));
	  gz       = 0.5*(g4hit->get_z(0)+g4hit->get_z(1));
	  gx       = g4hit->get_x(0);
	  gy       = g4hit->get_y(0);
	  gz       = g4hit->get_z(0);

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

	    PHG4Hit* outerhit = trutheval->get_outermost_truth_hit(g4particle);	
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

	float cluster_data[33] = {(float) _ievent,
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
				  gembed,
				  gprimary,
				  nhits,
				  efromtruth};

	_ntp_cluster->Fill(cluster_data);
      }		  
    }
  }
  
  //------------------------
  // fill the Gtrack NTuple
  //------------------------

  // need things off of the DST...

  if (_ntp_gtrack) {
    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");   
    SvtxClusterMap* clustermap = findNode::getClass<SvtxClusterMap>(topNode,"SvtxClusterMap");
    if (truthinfo) {

      PHG4TruthInfoContainer::ConstRange range = truthinfo->GetPrimaryParticleRange();
      for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
	   iter != range.second; 
	 ++iter) {
	
	PHG4Particle* g4particle = iter->second;
    
	float gtrackID = g4particle->get_track_id();
	float gflavor  = g4particle->get_pid();
      
	std::set<PHG4Hit*> g4hits = trutheval->all_truth_hits(g4particle);     
	float ng4hits = g4hits.size();  
	float gpx      = g4particle->get_px();
	float gpy      = g4particle->get_py();
	float gpz      = g4particle->get_pz();

	PHG4VtxPoint* vtx = trutheval->get_vertex(g4particle);	
	float gvx      = vtx->get_x();
	float gvy      = vtx->get_y();
	float gvz      = vtx->get_z();

	float gfpx      = 0.;
	float gfpy      = 0.;
	float gfpz      = 0.;
	float gfx       = 0.;
	float gfy       = 0.;
	float gfz       = 0.;
    
	PHG4Hit* outerhit = trutheval->get_outermost_truth_hit(g4particle);	
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

	SvtxTrack* track = trackeval->best_track_from(g4particle);

	float trackID       = NAN;
	float charge        = NAN;
	float quality       = NAN;
	float chisq         = NAN;
	float ndf           = NAN;
	float nhits         = NAN;
	unsigned int layers = 0x0;
	float dca2d         = NAN;
	float dca2dsigma    = NAN;
	float px            = NAN;
	float py            = NAN;
	float pz            = NAN;
	float pcax          = NAN;
	float pcay          = NAN;
	float pcaz          = NAN;

	float nfromtruth    = NAN;

	if (track) {
	  trackID   = track->get_id();     
	  charge    = track->get_charge();
	  quality   = track->get_quality();
	  chisq     = track->get_chisq();
	  ndf       = track->get_ndf();
	  nhits     = track->size_clusters();
	  
	  for (SvtxTrack::ConstClusterIter iter = track->begin_clusters();
	       iter != track->end_clusters();
	       ++iter) {
	    unsigned int cluster_id = *iter;
	    SvtxCluster* cluster = clustermap->get(cluster_id);
	    unsigned int layer = cluster->get_layer();
	    if (layer < 32) layers |= (0x1 << layer);
	  }

	  dca2d     = track->get_dca2d();
	  dca2dsigma = track->get_dca2d_error();
	  px        = track->get_px();
	  py        = track->get_py();
	  pz        = track->get_pz();
	  pcax      = track->get_x();
	  pcay      = track->get_y();
	  pcaz      = track->get_z();

	  nfromtruth = trackeval->get_nclusters_contribution(track,g4particle);
	}
      
	float gtrack_data[34] = {(float) _ievent,
				 gtrackID,
				 gflavor,
				 ng4hits,
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
				 trackID,
				 px,         
				 py,         
				 pz,  
				 charge,     
				 quality,    
				 chisq,      
				 ndf,        
				 nhits,      
				 (float) layers,     
				 dca2d,      
				 dca2dsigma, 			       
				 pcax,       
				 pcay,       
				 pcaz,
				 nfromtruth
	};

	_ntp_gtrack->Fill(gtrack_data);
      }	     
    }
  }
  
  //------------------------
  // fill the Track NTuple
  //------------------------

  if (_ntp_track) {
    // need things off of the DST...
    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
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
	for (SvtxTrack::ConstClusterIter iter = track->begin_clusters();
	     iter != track->end_clusters();
	     ++iter) {
	  unsigned int cluster_id = *iter;
	  SvtxCluster* cluster = clustermap->get(cluster_id);
	  unsigned int layer = cluster->get_layer();
	  if (layer < 32) layers |= (0x1 << layer);
	}
      
	float dca2d     = track->get_dca2d();
	float dca2dsigma = track->get_dca2d_error();
	float px        = track->get_px();
	float py        = track->get_py();
	float pz        = track->get_pz();
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

	float nfromtruth = NAN;
      
	PHG4Particle* g4particle = trackeval->max_truth_particle_by_nclusters(track);
      
	if (g4particle) {
	  gtrackID = g4particle->get_track_id();
	  gflavor  = g4particle->get_pid();
      
	  std::set<PHG4Hit*> g4hits = trutheval->all_truth_hits(g4particle);     
	  ng4hits = g4hits.size();  
	  gpx      = g4particle->get_px();
	  gpy      = g4particle->get_py();
	  gpz      = g4particle->get_pz();
	
	  PHG4VtxPoint* vtx = trutheval->get_vertex(g4particle);	
	  gvx      = vtx->get_x();
	  gvy      = vtx->get_y();
	  gvz      = vtx->get_z();

	  PHG4Hit* outerhit = trutheval->get_outermost_truth_hit(g4particle);	      
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
	}
      
	float track_data[50] = {(float) _ievent,
				trackID, 
				px,        
				py,        
				pz,      
				charge,  
				quality, 
				chisq,   
				ndf,     
				nhits,   
				(float) layers,
				dca2d,     
				dca2dsigma,      
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
				nfromtruth
	};
      
	_ntp_track->Fill(track_data);
      }
    }
  }
  
  return;
}
