
#include "SvtxEvaluator.h"

#include "SvtxVertexEval.h"
#include "SvtxTrackEval.h"
#include "SvtxClusterEval.h"
#include "SvtxHitEval.h"
#include "SvtxTruthEval.h"

#include <phool/PHCompositeNode.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/getClass.h>

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

using namespace std;

SvtxEvaluator::SvtxEvaluator(const string &name, const string &filename) :
  SubsysReco("SvtxEvaluator"),
  _ievent(0),
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

  //-----------------------------------
  // print what is coming into the code
  //-----------------------------------
  
  printInputInfo();
  
  //---------------------------
  // fill the Evaluator NTuples
  //---------------------------
  
  fillOutputNtuples(topNode);
  
  //--------------------------------------------------
  // Print out the ancestry information for this event
  //--------------------------------------------------
  
  printOutputInfo();
  
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

  if (verbosity >= 0) {
    cout << "========================= SvtxEvaluator::End() ============================" << endl;
    cout << " " << _ievent << " events of output written to: " << _filename << endl;
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void SvtxEvaluator::printInputInfo() {
  
  // if(verbosity > 1) cout << "SvtxEvaluator::fillClusterToG4hitMap() entered" << endl;

  // if(verbosity > 3)
  //   {
  //     // event information
  //     cout << endl;
  //     cout << PHWHERE << "   INPUT FOR EVENT " << _ievent << endl;

  //     cout << endl;
  //     cout << "---PHG4HITCONTAINER-------------" << endl;
  //     unsigned int ig4hit=0;
  //     for(HitMap::const_iterator g4hit_iter = _g4hitList.begin();
  // 	  g4hit_iter != _g4hitList.end();
  // 	  g4hit_iter++)
  // 	{
  // 	  PHG4Hit *g4hit = g4hit_iter->second;
  // 	  cout << ig4hit << " of " << _g4hitList.size();
  // 	  cout << ": PHG4Hit: " << g4hit << endl;
  // 	  ++ig4hit;
  // 	}

  //     cout << "---PHG4SVTXSIMPLECLUSTERLIST-------------" << endl;
  //     unsigned int icluster = 0;
  //     for (SvtxClusterMap::Iter iter = _clusterList->begin();
  // 	   iter != _clusterList->end();
  // 	   ++iter) {
  // 	SvtxCluster* cluster = &iter->second;
  // 	cout << icluster << " of " << _clusterList->size();	  
  // 	cout << ": SvtxCluster: " << cluster << endl;
  // 	++icluster;
  //     }

  //     if(_trackingWasRun)
  // 	{
  // 	  cout << "---SVXTRACKLIST-------------" << endl;
  // 	  for (SvtxTrackMap::Iter iter = _trackList->begin();
  // 	       iter != _trackList->end();
  // 	       ++iter) {
  // 	    SvtxTrack *track = &iter->second;
  // 	    cout << "SvtxTrack:" << endl;
  // 	    track->identify(cout);
  // 	    cout << endl;
  // 	  }
  // 	}

  //   }

  return;
}

void SvtxEvaluator::printOutputInfo() {
  
  // if(verbosity > 1) cout << "SvtxEvaluator::printLogInfo() entered" << endl;

  // //==========================================
  // // print out some useful stuff for debugging
  // //==========================================

  // if(verbosity > 0)
  //   {
  //     // event information
  //     cout << endl;
  //     cout << PHWHERE << "   NEW OUTPUT FOR EVENT " << _ievent << endl;
  //     cout << endl;

  //     PHG4VtxPoint *gvertex = _truth_info_container->GetPrimaryVtx( _truth_info_container->GetPrimaryVertexIndex() );
  //     float gvx = gvertex->get_x();
  //     float gvy = gvertex->get_y();
  //     float gvz = gvertex->get_z();

  //     float vx = NAN;
  //     float vy = NAN;
  //     float vz = NAN;
  //     if (_vertexList) {
  // 	if (!_vertexList->empty()) {
  // 	  SvtxVertex* vertex = &(_vertexList->begin()->second);
	
  // 	  vx = vertex->get_x();
  // 	  vy = vertex->get_y();
  // 	  vz = vertex->get_z();
  // 	}
  //     }

  //     cout << "===Vertex Reconstruction=======================" << endl;
  //     cout << "vtrue = (" << gvx << "," << gvy << "," << gvz << ") => vreco = (" << vx << "," << vy << "," << vz << ")" << endl;
  //     cout << endl;

  //     cout << "===Tracking Summary============================" << endl;
  //     unsigned int ng4hits[10] = {0};
  //     for(unsigned int ilayer=0; ilayer<_nlayers; ++ilayer)
  // 	{
  // 	  for(HitMap::const_iterator g4hit_iter = _g4hitList.begin();
  // 	      g4hit_iter != _g4hitList.end();
  // 	      g4hit_iter++)
  // 	    {
  // 	      PHG4Hit *g4hit = g4hit_iter->second;
  // 	      if (g4hit->get_layer() == ilayer) ++ng4hits[ilayer];
  // 	    }
  // 	}

  //     unsigned int nclusters[10] = {0};
  //     for (SvtxClusterMap::Iter iter = _clusterList->begin();
  // 	   iter != _clusterList->end();
  // 	   ++iter) {

  // 	SvtxCluster* cluster = &iter->second;
  // 	++nclusters[cluster->get_layer()];
  //     }

  //     for(unsigned int ilayer=0; ilayer<_nlayers; ++ilayer) {
  // 	cout << "layer " << ilayer << ": nG4hits = " << ng4hits[ilayer]
  // 	     << " => nCells = " << _nhits_per_layer[ilayer]
  // 	     << " => nClusters = " << nclusters[ilayer] << endl;
  //     }
    
  //     cout << "nGtracks = " << _gtrack_list.size();
  //     cout << " => nTracks = ";
  //     if(_trackingWasRun) 
  // 	{
  // 	  cout << _trackList->size() << endl;
  // 	}
  //     else
  // 	{
  // 	  cout << 0 << endl;
  // 	}

  //     // cluster wise information
  //     if(verbosity > 1)
  // 	{
 
  // 	  unsigned int ig4hit=0;
  // 	  for(HitMap::const_iterator g4hit_iter = _g4hitList.begin();
  // 	      g4hit_iter != _g4hitList.end();
  // 	      g4hit_iter++)
  // 	    {
  // 	      PHG4Hit *g4hit = g4hit_iter->second;

  // 	      cout << endl;
  // 	      cout << "===PHG4Hit===================================" << endl;
  // 	      cout << " PHG4Hit: " << g4hit;
  
  // 	      typedef multimap<PHG4Hit*,SvtxCluster*>::iterator mapiter2;
  // 	      typedef pair<mapiter2,mapiter2> maprange2;
  // 	      maprange2 therange2 = _allg4hits_cluster_mmap.equal_range( g4hit );
  // 	      for(mapiter2 theiter2=therange2.first; theiter2!=therange2.second; theiter2++) 
  // 		{
  // 		  SvtxCluster *cluster = theiter2->second;	  
  // 		  cout << "===Created-SvtxCluster================" << endl;      
  // 		  cout << "SvtxCluster: "; cluster->identify();
  // 		}

  // 	      ++ig4hit;
  // 	    } 
  // 	}
      
  //     for(unsigned igtrack = 0; igtrack < _gtrack_list.size(); igtrack++)
  // 	{
  // 	  SvxGtrack *gtrack = &_gtrack_list[igtrack];

  // 	  // don't print out the non-primary tracks
  // 	  if (!gtrack->get_primary()) continue;
	  
  // 	  // track-wise information
  // 	  cout << endl;

  // 	  cout << "=== Gtrack ===================================================" << endl;
  // 	  cout << " Gtrack id = " << gtrack->get_track_id() << endl;
  // 	  cout << " match = " << gtrack->get_standalone_match() << endl;
  // 	  cout << " best purity = " << gtrack->get_best_purity() << endl;
  // 	  cout << " best dp/p = " << gtrack->get_best_dpp() << endl;

  // 	  cout << " PHG4Particle: ";
  // 	  gtrack->get_particle()->identify(cout);
  // 	  cout << " ptrue = (";
  // 	  cout.width(5); cout << gtrack->get_px();
  // 	  cout << ",";
  // 	  cout.width(5); cout << gtrack->get_py();
  // 	  cout << ",";
  // 	  cout.width(5); cout << gtrack->get_pz();
  // 	  cout << ")" << endl;

  // 	  cout << " vtrue = (";
  // 	  cout.width(5); cout << gtrack->get_vx();
  // 	  cout << ",";
  // 	  cout.width(5); cout << gtrack->get_vy();
  // 	  cout << ",";
  // 	  cout.width(5); cout << gtrack->get_vz();
  // 	  cout << ")" << endl;
	  
  // 	  cout << " pt = " << sqrt(pow(gtrack->get_px(),2)+pow(gtrack->get_py(),2)) << endl;
  // 	  cout << " phi = " << atan2(gtrack->get_py(),gtrack->get_px()) << endl;
  // 	  cout << " eta = " << asinh(gtrack->get_pz()/sqrt(pow(gtrack->get_px(),2)+pow(gtrack->get_py(),2))) << endl;
	  
  // 	  cout << " chi^2 = " << gtrack->get_chisq() << endl;
  // 	  cout << " chi^2 w/ vertex = " << gtrack->get_chisqv() << endl;
  // 	  cout << " embed flag = " << gtrack->get_embed() << endl;
  // 	  cout << " primary flag = " << gtrack->get_primary() << endl;
  // 	  cout << " ---Associated-PHG4Hits-----------------------------------------" << endl;
	  
  // 	  for(unsigned int ig4hit = 0; ig4hit < gtrack->get_ng4hits(); ig4hit++)
  // 	    {
  // 	      PHG4Hit *g4hit = gtrack->get_g4hit(ig4hit);
	      
  // 	      float x = g4hit->get_x(0);
  // 	      float y = g4hit->get_y(0);
  // 	      float z = g4hit->get_z(0);
	      
  // 	      cout << " #" << g4hit->get_hit_id() << " xtrue = (";
  // 	      cout.width(5); cout << x;
  // 	      cout << ",";
  // 	      cout.width(5); cout << y;
  // 	      cout << ",";
  // 	      cout.width(5); cout << z;
  // 	      cout << ")";

  // 	      typedef multimap<PHG4Hit*,SvtxCluster*>::iterator mapiter2;
  // 	      typedef pair<mapiter2,mapiter2> maprange2;
  // 	      maprange2 therange2 = _allg4hits_cluster_mmap.equal_range( g4hit );
  // 	      for(mapiter2 theiter2=therange2.first; theiter2!=therange2.second; theiter2++) 
  // 		{
  // 		  SvtxCluster *cluster = theiter2->second;

  // 		  float x = cluster->get_x();
  // 		  float y = cluster->get_y();
  // 		  float z = cluster->get_z();
		 
  // 		  cout << " => #" << cluster->get_id(); 
  // 		  cout << " xreco = (";
  // 		  cout.width(5); cout << x;
  // 		  cout << ",";
  // 		  cout.width(5); cout << y;
  // 		  cout << ",";
  // 		  cout.width(5); cout << z;
  // 		  cout << ")";
  // 		}

  // 	      cout << endl;
  // 	    }

  // 	  if(_trackingWasRun)
  // 	    {
  // 	      typedef multimap<SvxGtrack*,SvtxTrack*>::iterator mapiter;
  // 	      typedef pair<mapiter,mapiter> maprange;
  // 	      maprange therange = _gtrack_track_mmap.equal_range( gtrack );
  // 	      for(mapiter theiter=therange.first; theiter!=therange.second; theiter++) 
  // 		{
  // 		  SvtxTrack *track = theiter->second;

  // 		  float px = track->get3Momentum(0);
  // 		  float py = track->get3Momentum(1);
  // 		  float pz = track->get3Momentum(2);

  // 		  cout << "===Created-SvtxTrack==========================================" << endl;
  // 		  cout << " SvtxTrack id = " << track->getTrackID() << endl;
  // 		  cout << " preco = (";
  // 		  cout.width(5); cout << px;
  // 		  cout << ",";
  // 		  cout.width(5); cout << py;
  // 		  cout << ",";
  // 		  cout.width(5); cout << pz;
  // 		  cout << ")" << endl;
  // 		  cout << " quality = " << track->getQuality() << endl;
  // 		  cout << " purity = " << _track_purity_map[track] << endl;

  // 		  cout << " ---Associated-SvtxClusters-to-PHG4Hits-------------------------" << endl;
	  
  // 		  // loop over the associated track clusters
  // 		  typedef multimap<SvtxTrack*,SvtxCluster*>::iterator mapiterS2C;
  // 		  typedef pair<mapiterS2C,mapiterS2C> maprangeS2C;
  // 		  maprangeS2C rangeS2C = _track_cluster_mmap.equal_range( track );
  // 		  for(mapiterS2C iterS2C = rangeS2C.first; iterS2C != rangeS2C.second; iterS2C++) 
  // 		    {
  // 		      SvtxCluster *cluster = iterS2C->second;
		  
  // 		      float x = cluster->get_x();
  // 		      float y = cluster->get_y();
  // 		      float z = cluster->get_z();
			  
  // 		      cout << " #" << cluster->get_id() << " xreco = (";
  // 		      cout.width(5); cout << x;
  // 		      cout << ",";
  // 		      cout.width(5); cout << y;
  // 		      cout << ",";
  // 		      cout.width(5); cout << z;
  // 		      cout << ") =>";
			  
  // 		      typedef multimap<SvtxCluster*,PHG4Hit*>::const_iterator mmapiter;
  // 		      pair<mmapiter,mmapiter> hitrange = _cluster_allg4hits_mmap.equal_range(cluster);
  // 		      for(mmapiter hititer = hitrange.first;
  // 			  hititer!=hitrange.second;
  // 			  hititer++) {
  // 			PHG4Hit *g4hit = hititer->second;
			
  // 			if (g4hit->get_trkid() == gtrack->get_track_id()) {
			  
  // 			  if (g4hit) {
  // 			    x = g4hit->get_x(0);
  // 			    y = g4hit->get_y(0);
  // 			    z = g4hit->get_z(0);
			    
  // 			    cout << " #" << g4hit->get_hit_id()
  // 				 << " xtrue = (";
  // 			    cout.width(5); cout << x;
  // 			    cout << ",";
  // 			    cout.width(5); cout << y;
  // 			    cout << ",";
  // 			    cout.width(5); cout << z;
  // 			    cout << ") => Gtrack id = " << g4hit->get_trkid();
  // 			    break; // print only the first match
  // 			  } else {
  // 			    cout << " noise hit";
  // 			  }
  // 			} 
  // 		      }
  // 		      cout << endl;
  // 		    }
  // 		}
  // 	    }
  // 	}
      
  //     cout << endl;

  //   } // if verbosity

  return;
}

void SvtxEvaluator::fillOutputNtuples(PHCompositeNode *topNode) {

  if (verbosity > 1) cout << "SvtxEvaluator::fillOutputNtuples() entered" << endl;

  SvtxTruthEval otrutheval(topNode);
  SvtxTruthEval* trutheval = &otrutheval;

  SvtxVertexEval overtexeval(topNode);
  SvtxVertexEval* vertexeval = &overtexeval;
  SvtxTrackEval* trackeval = vertexeval->get_track_eval();
  SvtxClusterEval* clustereval = vertexeval->get_cluster_eval();
  SvtxHitEval* hiteval = vertexeval->get_hit_eval();
  
  //-----------------------
  // fill the Vertex NTuple
  //-----------------------

  SvtxVertexMap* vertexmap = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMap");
  PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if (vertexmap && truthinfo) {
    for (SvtxVertexMap::Iter iter = vertexmap->begin();
	 iter != vertexmap->end();
	 ++iter) {
      SvtxVertex* vertex = &iter->second;
      PHG4VtxPoint* point = vertexeval->max_truth_point_by_ntracks(vertex);

      float vx         = vertex->get_x();
      float vy         = vertex->get_y();
      float vz         = vertex->get_z();
      float ntracks    = vertex->size_tracks();
      float gvx        = point->get_x();
      float gvy        = point->get_y();
      float gvz        = point->get_z();
      float gntracks   = truthinfo->GetNumPrimaryVertexParticles();
      float nfromtruth = vertexeval->get_ntracks_contribution(vertex,point);
      
      float event_data[10] = {_ievent,
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

      _ntp_vertex->Fill(event_data);      
    }
  }

  //-----------------------
  // fill the gpoint NTuple
  //-----------------------

  if (vertexmap && truthinfo) {
    PHG4VtxPoint* point = truthinfo->GetPrimaryVtx(truthinfo->GetVtxRange().first->second->get_id());
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
	
    float gpoint_data[10] = {_ievent,
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

  //---------------------
  // fill the G4hit NTuple
  //---------------------

  unsigned int i = 0;
  std::set<PHG4Hit*> g4hits = trutheval->all_truth_hits();
  for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
       iter != g4hits.end();
       ++iter) {

    if ((_ievent==0)&&(i==0)) cout << "SvtxEvaluator:: WARNING - g4hit eval limited to 100 entries" << endl;
    if (i > 100) break;
    ++i;
    
    PHG4Hit *g4hit = *iter;
    PHG4Particle *g4particle = trutheval->get_particle(g4hit);
    
    float g4hitID   = g4hit->get_hit_id();
    float gx        = 0.5*(g4hit->get_x(0)+g4hit->get_x(1));
    float gy        = 0.5*(g4hit->get_y(0)+g4hit->get_y(1));
    float gz        = 0.5*(g4hit->get_z(0)+g4hit->get_z(1));
    float gedep     = g4hit->get_edep();
    float glayer    = g4hit->get_layer();
  
    float gtrackID  = g4hit->get_trkid();
    float gflavor   = g4particle->get_pid();
    float gpx       = g4particle->get_px();
    float gpy       = g4particle->get_py();
    float gpz       = g4particle->get_pz();

    PHG4VtxPoint* vtx = trutheval->get_vertex(g4particle);	
    float gvx       = vtx->get_x();
    float gvy       = vtx->get_y();
    float gvz       = vtx->get_z();

    float gfpx      = NULL;
    float gfpy      = NULL;
    float gfpz      = NULL;
    float gfx       = NULL;
    float gfy       = NULL;
    float gfz       = NULL;
    
    PHG4Hit* outerhit = trutheval->get_outermost_truth_hit(g4particle);	
    if (outerhit) {
      gfpx      = outerhit->get_px(1);
      gfpy      = outerhit->get_py(1);
      gfpz      = outerhit->get_pz(1);
      gfx       = outerhit->get_x(1);
      gfy       = outerhit->get_y(1);
      gfz       = outerhit->get_z(1);
    }
    float gembed    = trutheval->get_embed(g4particle);
    float gprimary  = trutheval->is_primary(g4particle);

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
      efromtruth = clustereval->get_energy_contribution(cluster,g4particle);
    }

    float g4hit_data[35] = {_ievent,
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
  
  //--------------------
  // fill the Hit NTuple
  //--------------------

  // need things off of the DST...
  SvtxHitMap* hitmap = findNode::getClass<SvtxHitMap>(topNode,"SvtxHitMap");
  if (hitmap) {

    unsigned int i = 0;
    for (SvtxHitMap::Iter iter = hitmap->begin();
  	 iter != hitmap->end();
  	 ++iter) {

      if ((_ievent==0)&&(i==0)) cout << "SvtxEvaluator:: WARNING - hit eval limited to 100 entries" << endl;
      if (i > 100) break;
      ++i;
      
      SvtxHit* hit             = &iter->second;
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
      float glast    = NAN;
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
	gtrackID = g4particle->get_track_id();
	gflavor  = g4particle->get_pid();
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
	glast    = NAN;
	gembed   = trutheval->get_embed(g4particle);
	gprimary = trutheval->is_primary(g4particle);
      }      

      efromtruth = hiteval->get_energy_contribution(hit,g4particle);
      
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

  //------------------------
  // fill the Cluster NTuple
  //------------------------

  // need things off of the DST...
  SvtxClusterMap* clustermap = findNode::getClass<SvtxClusterMap>(topNode,"SvtxClusterMap");
  if (clustermap) {

    for (SvtxClusterMap::Iter iter = clustermap->begin();
	 iter != clustermap->end();
	 ++iter) {
    
      SvtxCluster* cluster     = &iter->second;   
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
      float glast    = NAN;
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
	gtrackID = g4particle->get_track_id();
	gflavor  = g4particle->get_pid();
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
	glast    = NAN;
	gembed   = trutheval->get_embed(g4particle);
	gprimary = trutheval->is_primary(g4particle);
      }      

      efromtruth = clustereval->get_energy_contribution(cluster,g4particle);

      float cluster_data[33] = {_ievent,
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
  
  //------------------------
  // fill the Gtrack NTuple
  //------------------------

  // need things off of the DST...

  if (truthinfo) {
    
    PHG4TruthInfoContainer::Map map = truthinfo->GetPrimaryMap();
    for (PHG4TruthInfoContainer::ConstIterator iter = map.begin(); 
	 iter != map.end(); 
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

      float gfpx      = NULL;
      float gfpy      = NULL;
      float gfpz      = NULL;
      float gfx       = NULL;
      float gfy       = NULL;
      float gfz       = NULL;
    
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
      float dca           = NAN;
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
	trackID   = track->getTrackID();     
	charge    = track->getCharge();
	quality   = track->getQuality();
	chisq     = track->getChisq();
	ndf       = track->getNDF();
	nhits     = track->getNhits();

	for (unsigned int i = 0; i < 32; ++i){ // only 32 bits available	
	  if (track->hasCluster(i)) {
	    layers |= (0x1 << i);
	  }
	}

	dca       = track->getDCA();
	dca2d     = track->getDCA2D();
	dca2dsigma = track->getDCA2Dsigma();
	px        = track->get3Momentum(0);
	py        = track->get3Momentum(1);
	pz        = track->get3Momentum(2);
	pcax      = track->d * sin(track->phi);
	pcay      = track->d * cos(track->phi);
	pcaz      = track->z0;

	nfromtruth = trackeval->get_nclusters_contribution(track,g4particle);
      }
      
      float gtrack_data[34] = {_ievent,
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
			       layers,     
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
  	
  //------------------------
  // fill the Track NTuple
  //------------------------

  // need things off of the DST...
  SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
  if (trackmap) {

    for (SvtxTrackMap::Iter iter = trackmap->begin();
	 iter != trackmap->end();
	 ++iter) {
    
      SvtxTrack* track         = &iter->second;

      float trackID   = track->getTrackID();     
      float charge    = track->getCharge();
      float quality   = track->getQuality();
      float chisq     = track->getChisq();
      float ndf       = track->getNDF();
      float nhits     = track->getNhits();

      unsigned int layers = 0x0;
      for (unsigned int i = 0; i < 32; ++i){ // only 32 bits available	
	if (track->hasCluster(i)) {
	  layers |= (0x1 << i);
	}
      }
      
      float dca2d     = track->getDCA2D();
      float dca2dsigma = track->getDCA2Dsigma();
      float px        = track->get3Momentum(0);
      float py        = track->get3Momentum(1);
      float pz        = track->get3Momentum(2);
      float pcax      = track->d * sin(track->phi);
      float pcay      = track->d * cos(track->phi);
      float pcaz      = track->z0;

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
      
      float track_data[50] = {_ievent,
			      trackID, 
			      px,        
			      py,        
			      pz,      
			      charge,  
			      quality, 
			      chisq,   
			      ndf,     
			      nhits,   
			      layers,
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
  
  return;
}
