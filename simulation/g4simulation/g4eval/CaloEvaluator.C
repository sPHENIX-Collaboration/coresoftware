
#include "CaloEvaluator.h"

#include "CaloEvalStack.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/getClass.h>
#include <fun4all/SubsysReco.h>
#include <phool/PHCompositeNode.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4Particle.h>
#include <g4hough/SvtxVertexMap.h>
#include <g4cemc/RawTowerContainer.h>
#include <g4cemc/RawTowerGeom.h>
#include <g4cemc/RawTower.h>
#include <g4cemc/RawClusterContainer.h>
#include <g4cemc/RawCluster.h>

#include <TNtuple.h>
#include <TFile.h>

#include <iostream>
//#include <vector>
//#include <algorithm>
#include <cmath>

using namespace std;

CaloEvaluator::CaloEvaluator(const string &name, const string &caloname, const string &filename) 
  : SubsysReco(name),
    _caloname(caloname),
    _ievent(0),
    _do_gpoint_eval(true),
    _do_gshower_eval(true),
    _do_tower_eval(true),
    _do_cluster_eval(true),
    _ntp_gpoint(NULL),
    _ntp_gshower(NULL),
    _ntp_tower(NULL),
    _ntp_cluster(NULL),
    _filename(filename),
    _tfile(NULL) {
  verbosity = 0;
}

int CaloEvaluator::Init(PHCompositeNode *topNode) {
  
  _ievent = 0;

  _tfile = new TFile(_filename.c_str(), "RECREATE");


  if (_do_gpoint_eval) _ntp_gpoint = new TNtuple("ntp_gpoint","primary vertex => best (first) vertex",
						 "event:gvx:gvy:gvz:"
						 "vx:vy:vz:");
  
  if (_do_gshower_eval) _ntp_gshower = new TNtuple("ntp_gshower","truth shower => best cluster",
						   "event:gparticleID:gflavor:gnhits:"
						   "geta:gphi:ge:gpt:gvx:gvy:gvz:gembed:gedep:gmrad:"
						   "clusterID:ntowers:eta:phi:e:efromtruth");
  
  if (_do_tower_eval) _ntp_tower = new TNtuple("ntp_tower","tower => max truth primary",
					       "event:towerID::ieta:iphi:eta:phi:e:"
					       "gparticleID:gflavor:gnhits:"
					       "geta:gphi:ge:gpt:gvx:gvy:gvz:"
					       "gembed:gedep:gmrad:"
					       "efromtruth");

  if (_do_cluster_eval) _ntp_cluster = new TNtuple("ntp_cluster","cluster => max truth primary",
						   "event:clusterID:ntowers:eta:phi:e:"
						   "gparticleID:gflavor:gnhits:"
						   "geta:gphi:ge:gpt:gvx:gvy:gvz:"
						   "gembed:gedep:gmrad:"
						   "efromtruth");

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloEvaluator::process_event(PHCompositeNode *topNode) {
  
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

int CaloEvaluator::End(PHCompositeNode *topNode) {
  
  _tfile->cd();

  if (_do_gpoint_eval)  _ntp_gpoint->Write();
  if (_do_gshower_eval) _ntp_gshower->Write();
  if (_do_tower_eval)   _ntp_tower->Write();
  if (_do_cluster_eval) _ntp_cluster->Write();
  
  _tfile->Close();

  delete _tfile;

  if (verbosity > 0) {
    cout << "========================= CaloEvaluator::End() ============================" << endl;
    cout << " " << _ievent << " events of output written to: " << _filename << endl;
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloEvaluator::printInputInfo(PHCompositeNode *topNode) {
  
  if (verbosity > 2) cout << "CaloEvaluator::printInputInfo() entered" << endl;

  // print out the truth container

  if (verbosity > 1) { 

    cout << endl;
    cout << PHWHERE << "   NEW INPUT FOR EVENT " << _ievent << endl;
    cout << endl;

    // need things off of the DST...
    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
    if (!truthinfo) {
      cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
      exit(-1);
    }
    
    cout << "PHG4TruthInfoContainer contents: " << endl; 

    PHG4TruthInfoContainer::Range truthrange = truthinfo->GetHitRange();
    for(PHG4TruthInfoContainer::Iterator truthiter = truthrange.first;
	truthiter != truthrange.second;
	++truthiter) {
      PHG4Particle *particle = truthiter->second;

      cout << truthiter->first << " => pid: " << particle->get_pid()
	   << " pt: " << sqrt(pow(particle->get_px(),2)+pow(particle->get_py(),2)) << endl;
    }
  }

  return;
}

void CaloEvaluator::printOutputInfo(PHCompositeNode *topNode) {
  
  if (verbosity > 2) cout << "CaloEvaluator::printOutputInfo() entered" << endl;

  //==========================================
  // print out some useful stuff for debugging
  //==========================================

  // if (verbosity > 1)
  //   {
  //     // event information
  //     cout << endl;
  //     cout << PHWHERE << "   NEW OUTPUT FOR EVENT " << _ievent << endl;
  //     cout << endl;

  //     PHG4VtxPoint *gvertex = truthinfo->GetPrimaryVtx( truthinfo->GetPrimaryVertexIndex() );
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

  //     cout << "vtrue = (" << gvx << "," << gvy << "," << gvz << ") => vreco = (" << vx << "," << vy << "," << vz << ")" << endl;
    
  //     float ngshowers = _gshower_list.size();
  //     float ng4hits  = _g4hitList->size();
  //     float ntowers  = _towerList->size();
  //     float nclusters = _clusterList->size();

  //     cout << "nGshowers = " << ngshowers << endl;
  //     cout << " => nGhits = " << ng4hits << endl;
  //     cout << " => nTowers = " << ntowers << endl;
  //     cout << " => nClusters = " << nclusters << endl;

  //     if(verbosity > 2)
  // 	{
  // 	  for(unsigned igshower = 0; igshower < _gshower_list.size(); igshower++)
  // 	    {
  // 	      CalGshower *gshower = &_gshower_list[igshower];
	      
  // 	      // track-wise information
  // 	      cout << endl;
      
  // 	      cout << "===CalGshower===================================================" << endl;
  // 	      cout << " CalGshower id = " << gshower->get_particle_id() << endl;
  // 	      cout << " flavor = " << gshower->get_flavor() << endl;
  // 	      cout << " ptrue = (";
  // 	      cout.width(5); cout << gshower->get_px();
  // 	      cout << ",";
  // 	      cout.width(5); cout << gshower->get_py();
  // 	      cout << ",";
  // 	      cout.width(5); cout << gshower->get_pz();
  // 	      cout << ")" << endl;
  // 	      cout << " vtrue = (";
  // 	      cout.width(5); cout << gshower->get_vx();
  // 	      cout << ",";
  // 	      cout.width(5); cout << gshower->get_vy();
  // 	      cout << ",";
  // 	      cout.width(5); cout << gshower->get_vz();
  // 	      cout << ")" << endl;
  // 	      cout << " ---Associated-PHG4Hits-----------------------------------------" << endl;
	      
  // 	      for(unsigned int ig4hit = 0; ig4hit < gshower->get_ng4hits(); ig4hit++)
  // 		{
  // 		  if((ig4hit > 5)&&(ig4hit < gshower->get_ng4hits() - 5)) continue;

  // 		  PHG4Hit *g4hit = gshower->get_g4hit(ig4hit);
		  
  // 		  float x = 0.5*(g4hit->get_x(1)+g4hit->get_x(0));
  // 		  float y = 0.5*(g4hit->get_y(1)+g4hit->get_y(0));
  // 		  float z = 0.5*(g4hit->get_z(1)+g4hit->get_z(0));
		  
  // 		  cout << " #" << ig4hit << " xtrue = (";
  // 		  cout.width(5); cout << x;
  // 		  cout << ",";
  // 		  cout.width(5); cout << y;
  // 		  cout << ",";
  // 		  cout.width(5); cout << z;
  // 		  cout << ")";
  // 		  cout << " e = " << g4hit->get_edep();
		  
  // 		  /*
  // 		    typedef multimap<PHG4Hit*,SvtxCluster*>::iterator mapiter2;
  // 		    typedef pair<mapiter2,mapiter2> maprange2;
  // 		    maprange2 therange2 = _g4hit_cluster_mmap.equal_range( g4hit );
  // 		    for(mapiter2 theiter2=therange2.first; theiter2!=therange2.second; theiter2++) 
  // 		    {
  // 		    SvtxCluster *cluster = theiter2->second;
		    
  // 		    float x = cluster->getHitPosition(0);
  // 		    float y = cluster->getHitPosition(1);
  // 		    float z = cluster->getHitPosition(2);
	    
  // 		    cout << " => #" << cluster->getClusterID() << " xreco = (";
  // 		    cout.width(5); cout << x;
  // 		    cout << ",";
  // 		    cout.width(5); cout << y;
  // 		    cout << ",";
  // 		    cout.width(5); cout << z;
  // 		    cout << ")";
  // 		    }
  // 		  */

  // 		  cout << endl;
  // 		}
  // 	    }      
  // 	}
  //   }

  return;
}

void CaloEvaluator::fillOutputNtuples(PHCompositeNode *topNode) {
  
  if (verbosity > 2) cout << "CaloEvaluator::fillOutputNtuples() entered" << endl;

  CaloEvalStack caloevalstack(topNode,_caloname); 
  CaloRawClusterEval* clustereval = caloevalstack.get_rawcluster_eval();
  CaloRawTowerEval*     towereval = caloevalstack.get_rawtower_eval();
  CaloTruthEval*        trutheval = caloevalstack.get_truth_eval();
  
  //----------------------
  // fill the Event NTuple
  //----------------------

  if (_do_gpoint_eval) {
    // need things off of the DST...
    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
    if (!truthinfo) {
      cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
      exit(-1);
    }

    // need things off of the DST...
    SvtxVertexMap* vertexmap = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMap");

    PHG4VtxPoint *gvertex = truthinfo->GetPrimaryVtx( truthinfo->GetPrimaryVertexIndex() );
    float gvx = gvertex->get_x();
    float gvy = gvertex->get_y();
    float gvz = gvertex->get_z();

    float vx = NAN;
    float vy = NAN;
    float vz = NAN;
    if (vertexmap) {
      if (!vertexmap->empty()) {
	SvtxVertex* vertex = &(vertexmap->begin()->second);
	
	vx = vertex->get_x();
	vy = vertex->get_y();
	vz = vertex->get_z();
      }
    }
	
    float gpoint_data[7] = {_ievent,
			    gvx,
			    gvy,
			    gvz,
			    vx,
			    vy,
			    vz
    };

    _ntp_gpoint->Fill(gpoint_data);   
  }
  
  //------------------------
  // fill the Gshower NTuple
  //------------------------
  
  if (_ntp_gshower) {

    if (verbosity > 1) cout << "CaloEvaluator::filling gshower ntuple..." << endl;
    
    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");   
    if (!truthinfo) {
      cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
      exit(-1);
    }
    
    PHG4TruthInfoContainer::Map map = truthinfo->GetPrimaryMap();
    for (PHG4TruthInfoContainer::ConstIterator iter = map.begin(); 
	 iter != map.end(); 
	 ++iter) {
      PHG4Particle* primary = iter->second;
      
      float gparticleID = primary->get_track_id();
      float gflavor     = primary->get_pid();
      
      std::set<PHG4Hit*> g4hits = trutheval->get_shower_from_primary(primary);     
      float gnhits   = g4hits.size();	
      float gpx      = primary->get_px();
      float gpy      = primary->get_py();
      float gpz      = primary->get_pz();
      float ge       = primary->get_e();

      float gpt = sqrt(gpx*gpx+gpy*gpy);
      float geta = NAN;
      if (gpt != 0.0) geta = asinh(gpz/gpt);
      float gphi = atan2(gpy,gpx);
      
      PHG4VtxPoint* vtx = trutheval->get_vertex(primary);	
      float gvx      = vtx->get_x();
      float gvy      = vtx->get_y();
      float gvz      = vtx->get_z();
      
      float gembed   = trutheval->get_embed(primary);
      float gedep    = trutheval->get_shower_energy_deposit(primary);
      float gmrad    = trutheval->get_shower_moliere_radius(primary);

      RawCluster* cluster = clustereval->best_cluster_from(primary);

      float clusterID = cluster->get_id();
      float ntowers   = cluster->getNTowers();
      float eta       = cluster->get_eta();
      float phi       = cluster->get_phi();
      float e         = cluster->get_energy();
	
      float efromtruth     = clustereval->get_energy_contribution(cluster, primary);
	
      float shower_data[20] = {_ievent,
			       gparticleID,
			       gflavor,
			       gnhits,
			       geta,
			       gphi,
			       ge,
			       gpt,
			       gvx,
			       gvy,
			       gvz,
			       gembed,
			       gedep,
			       gmrad,
			       clusterID,
			       ntowers,
			       eta,
			       phi,
			       e,
			       efromtruth
      };

      _ntp_gshower->Fill(shower_data);
    }
  }

  //----------------------
  // fill the Tower NTuple
  //----------------------
 
  if (_do_tower_eval) {

    if (verbosity > 1) cout << "CaloEvaluator::filling tower ntuple..." << endl;
    
    string towernode = "TOWER_" + _caloname;
    RawTowerContainer* towers = findNode::getClass<RawTowerContainer>(topNode,towernode.c_str());
    if (!towers) {
      cerr << PHWHERE << " ERROR: Can't find " << towernode << endl;
      exit(-1);
    }
    
    string towergeomnode = "TOWERGEOM_" + _caloname;
    RawTowerGeom* towergeom = findNode::getClass<RawTowerGeom>(topNode,towergeomnode.c_str());
    if (!towergeom) {
      cerr << PHWHERE << " ERROR: Can't find " << towergeomnode << endl;
      exit(-1);
    }
  
    RawTowerContainer::ConstRange begin_end = towers->getTowers();
    RawTowerContainer::ConstIterator rtiter;
    for (rtiter = begin_end.first; rtiter !=  begin_end.second; ++rtiter) {
      RawTower *tower = rtiter->second;

      float ieta    = tower->get_bineta();
      float iphi    = tower->get_binphi();
      float eta     = towergeom->get_etacenter(tower->get_bineta());
      float phi     = towergeom->get_phicenter(tower->get_binphi());
      float e       = tower->get_energy();

      PHG4Particle* primary = towereval->max_truth_primary_by_energy(tower);
    
      float gparticleID = primary->get_track_id();
      float gflavor     = primary->get_pid();
      
      std::set<PHG4Hit*> g4hits = trutheval->get_shower_from_primary(primary);     
      float gnhits   = g4hits.size();	
      float gpx      = primary->get_px();
      float gpy      = primary->get_py();
      float gpz      = primary->get_pz();
      float ge       = primary->get_e();

      float gpt = sqrt(gpx*gpx+gpy*gpy);
      float geta = NAN;
      if (gpt != 0.0) geta = asinh(gpz/gpt);
      float gphi = atan2(gpy,gpx);
      
      PHG4VtxPoint* vtx = trutheval->get_vertex(primary);	
      float gvx      = vtx->get_x();
      float gvy      = vtx->get_y();
      float gvz      = vtx->get_z();
      
      float gembed   = trutheval->get_embed(primary);
      float gedep    = trutheval->get_shower_energy_deposit(primary);
      float gmrad    = trutheval->get_shower_moliere_radius(primary);

      float efromtruth = towereval->get_energy_contribution(tower,primary);

      float tower_data[20] = {_ievent,
			      ieta,
			      iphi,
			      eta,
			      phi,
			      e,
			      gparticleID,
			      gflavor,
			      gnhits,
			      geta,
			      gphi,
			      ge,
			      gpt,
			      gvx,
			      gvy,
			      gvz,
			      gembed,
			      gedep,
			      gmrad,
			      efromtruth
      };
      
      _ntp_tower->Fill(tower_data);
    }
  }

  //------------------------
  // fill the Cluster NTuple
  //------------------------

  if (_do_cluster_eval) {
    if (verbosity > 1) cout << "CaloEvaluator::filling gcluster ntuple..." << endl;

    string clusternode = "CLUSTER_" + _caloname;
    RawClusterContainer* clusters = findNode::getClass<RawClusterContainer>(topNode,clusternode.c_str());
    if (!clusters) {
      cerr << PHWHERE << " ERROR: Can't find " << clusternode << endl;
      exit(-1);
    }
  
    // for every cluster
    for (unsigned int icluster = 0; icluster < clusters->size(); icluster++) {
      RawCluster *cluster = clusters->getCluster(icluster);

      float clusterID = cluster->get_id();
      float ntowers   = cluster->getNTowers();
      float eta       = cluster->get_eta();
      float phi       = cluster->get_phi();
      float e         = cluster->get_energy();
      
      PHG4Particle* primary = clustereval->max_truth_primary_by_energy(cluster);
    
      float gparticleID = primary->get_track_id();
      float gflavor     = primary->get_pid();
      
      std::set<PHG4Hit*> g4hits = trutheval->get_shower_from_primary(primary);     
      float gnhits   = g4hits.size();	
      float gpx      = primary->get_px();
      float gpy      = primary->get_py();
      float gpz      = primary->get_pz();
      float ge       = primary->get_e();

      float gpt = sqrt(gpx*gpx+gpy*gpy);
      float geta = NAN;
      if (gpt != 0.0) geta = asinh(gpz/gpt);
      float gphi = atan2(gpy,gpx);
      
      PHG4VtxPoint* vtx = trutheval->get_vertex(primary);	
      float gvx      = vtx->get_x();
      float gvy      = vtx->get_y();
      float gvz      = vtx->get_z();
      
      float gembed   = trutheval->get_embed(primary);
      float gedep    = trutheval->get_shower_energy_deposit(primary);
      float gmrad    = trutheval->get_shower_moliere_radius(primary);

      float efromtruth = clustereval->get_energy_contribution(cluster,primary);
      
      float cluster_data[20] = {_ievent,
				clusterID,
				ntowers,
				eta,
				phi,
				e,
				gparticleID,
				gflavor,
				gnhits,
				geta,
				gphi,
				ge,
				gpt,
				gvx,
				gvy,
				gvz,
				gembed,
				gedep,
				gmrad,
				efromtruth
      };

    _ntp_cluster->Fill(cluster_data);
    }
  }

  return;
}
