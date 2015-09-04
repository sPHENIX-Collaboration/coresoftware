
#include "JetEvaluator.h"

#include "JetEvalStack.h"
#include "JetRecoEval.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/getClass.h>
#include <fun4all/SubsysReco.h>
#include <phool/PHCompositeNode.h>
#include <g4jets/JetMap.h>
#include <g4jets/Jet.h>

#include <TNtuple.h>
#include <TFile.h>

#include <iostream>
#include <cmath>

using namespace std;

JetEvaluator::JetEvaluator(const string &name,
			   const string &recojetname,
			   const string &truthjetname,
			   const string &filename) 
  : SubsysReco(name),
    _recojetname(recojetname),
    _truthjetname(truthjetname),
    _ievent(0),
    _do_recojet_eval(true),
    _do_truthjet_eval(true),
    _ntp_recojet(NULL),
    _ntp_truthjet(NULL),
    _filename(filename),
    _tfile(NULL) {
  verbosity = 0;
}

int JetEvaluator::Init(PHCompositeNode *topNode) {
  
  _ievent = 0;

  _tfile = new TFile(_filename.c_str(), "RECREATE");
 
  if (_do_recojet_eval) _ntp_recojet = new TNtuple("ntp_recojet","reco jet => max truth jet",
						   "event:id:ncomp:eta:phi:e:pt:"
						   "gid:gncomp:geta:gphi:ge:gpt:"
						   "efromtruth");

  if (_do_truthjet_eval) _ntp_truthjet = new TNtuple("ntp_truthjet","truth jet => best reco jet",
						     "event:gid:gncomp:geta:gphi:ge:gpt:"
						     "id:ncomp:eta:phi:e:pt:"
						     "efromtruth");
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int JetEvaluator::process_event(PHCompositeNode *topNode) {
  
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

int JetEvaluator::End(PHCompositeNode *topNode) {
  
  _tfile->cd();

  if (_do_recojet_eval) _ntp_recojet->Write();
  if (_do_truthjet_eval) _ntp_truthjet->Write();
  
  _tfile->Close();

  delete _tfile;

  if (verbosity > 0) {
    cout << "========================== JetEvaluator::End() ============================" << endl;
    cout << " " << _ievent << " events of output written to: " << _filename << endl;
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void JetEvaluator::printInputInfo(PHCompositeNode *topNode) {
  
  // if (verbosity > 2) cout << "JetEvaluator::printInputInfo() entered" << endl;

  // // print out the truth container

  // if (verbosity > 1) { 

  //   cout << endl;
  //   cout << PHWHERE << "   NEW INPUT FOR EVENT " << _ievent << endl;
  //   cout << endl;

  //   // need things off of the DST...
  //   PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  //   if (!truthinfo) {
  //     cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
  //     exit(-1);
  //   }
    
  //   cout << "PHG4TruthInfoContainer contents: " << endl; 

  //   PHG4TruthInfoContainer::Range truthrange = truthinfo->GetHitRange();
  //   for(PHG4TruthInfoContainer::Iterator truthiter = truthrange.first;
  // 	truthiter != truthrange.second;
  // 	++truthiter) {
  //     PHG4Particle *particle = truthiter->second;

  //     cout << truthiter->first << " => pid: " << particle->get_pid()
  // 	   << " pt: " << sqrt(pow(particle->get_px(),2)+pow(particle->get_py(),2)) << endl;
  //   }
  // }

  return;
}

void JetEvaluator::printOutputInfo(PHCompositeNode *topNode) {
  
  // if (verbosity > 2) cout << "JetEvaluator::printOutputInfo() entered" << endl;

  // JetEvalStack jetevalstack(topNode,_jetname); 
  // JetRawClusterEval* clustereval = jetevalstack.get_rawcluster_eval();
  // JetTruthEval*        trutheval = jetevalstack.get_truth_eval();
  
  // //==========================================
  // // print out some useful stuff for debugging
  // //==========================================

  // if (verbosity > 1) {
    
  //   // event information
  //   cout << endl;
  //   cout << PHWHERE << "   NEW OUTPUT FOR EVENT " << _ievent << endl;
  //   cout << endl;

  //   // need things off of the DST...
  //   PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  //   if (!truthinfo) {
  //     cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
  //     exit(-1);
  //   }

  //   // need things off of the DST...
  //   SvtxVertexMap* vertexmap = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMap");
      
  //   PHG4VtxPoint *gvertex = truthinfo->GetPrimaryVtx( truthinfo->GetPrimaryVertexIndex() );
  //   float gvx = gvertex->get_x();
  //   float gvy = gvertex->get_y();
  //   float gvz = gvertex->get_z();

  //   float vx = NAN;
  //   float vy = NAN;
  //   float vz = NAN;
  //   if (vertexmap) {
  //     if (!vertexmap->empty()) {
  // 	SvtxVertex* vertex = &(vertexmap->begin()->second);
	
  // 	vx = vertex->get_x();
  // 	vy = vertex->get_y();
  // 	vz = vertex->get_z();
  //     }
  //   }

  //   cout << "vtrue = (" << gvx << "," << gvy << "," << gvz << ") => vreco = (" << vx << "," << vy << "," << vz << ")" << endl;

  //   PHG4TruthInfoContainer::Map map = truthinfo->GetPrimaryMap();
  //   for (PHG4TruthInfoContainer::ConstIterator iter = map.begin(); 
  // 	 iter != map.end(); 
  // 	 ++iter) {
  //     PHG4Particle* primary = iter->second;
      
  //     cout << endl;
      
  //     cout << "===Primary PHG4Particle=========================================" << endl;
  //     cout << " particle id = " << primary->get_track_id() << endl;
  //     cout << " flavor = " << primary->get_pid() << endl;
  //     cout << " (px,py,pz,e) = (";

  //     float gpx = primary->get_px();
  //     float gpy = primary->get_py();
  //     float gpz = primary->get_pz();
  //     float ge = primary->get_e();
	
  //     cout.width(5); cout << gpx;
  //     cout << ",";
  //     cout.width(5); cout << gpy;
  //     cout << ",";
  //     cout.width(5); cout << gpz;
  //     cout << ",";
  //     cout.width(5); cout << ge;
  //     cout << ")" << endl;

  //     float gpt = sqrt(gpx*gpx+gpy*gpy);
  //     float geta = NAN;
  //     if (gpt != 0.0) geta = asinh(gpz/gpt);
  //     float gphi = atan2(gpy,gpx);
      
  //     cout << "(eta,phi,e,pt) = (";
  //     cout.width(5); cout << geta;
  //     cout << ",";
  //     cout.width(5); cout << gphi;
  //     cout << ",";
  //     cout.width(5); cout << ge;
  //     cout << ",";
  //     cout.width(5); cout << gpt;
  //     cout << ")" << endl;

  //     PHG4VtxPoint* vtx = trutheval->get_vertex(primary);	
  //     float gvx      = vtx->get_x();
  //     float gvy      = vtx->get_y();
  //     float gvz      = vtx->get_z();
      
  //     cout << " vtrue = (";
  //     cout.width(5); cout << gvx;
  //     cout << ",";
  //     cout.width(5); cout << gvy;
  //     cout << ",";
  //     cout.width(5); cout << gvz;
  //     cout << ")" << endl;

  //     cout << " embed = " << trutheval->get_embed(primary) << endl;
  //     cout << " edep = " << trutheval->get_shower_energy_deposit(primary) << endl;
  //     cout << " mrad = " << trutheval->get_shower_moliere_radius(primary) << endl;

  //     std::set<RawCluster*> clusters = clustereval->all_clusters_from(primary);
  //     for (std::set<RawCluster*>::iterator clusiter = clusters.begin();
  // 	   clusiter != clusters.end();
  // 	   ++clusiter) {
  // 	RawCluster* cluster = (*clusiter);
	   
  // 	float ntowers   = cluster->getNTowers();
  // 	float eta       = cluster->get_eta();
  // 	float phi       = cluster->get_phi();
  // 	float e         = cluster->get_energy();
	
  // 	float efromtruth     = clustereval->get_energy_contribution(cluster, primary);
	
  // 	cout << " => #" << cluster->get_id() << " (eta,phi,e) = (";
  // 	cout.width(5); cout << eta;
  // 	cout << ",";
  // 	cout.width(5); cout << phi;
  // 	cout << ",";
  // 	cout.width(5); cout << e;
  // 	cout << "), ntowers = "<< ntowers <<", efromtruth = " << efromtruth << endl;
  //     }
  //   }
  //   cout << endl;
  // }

  return;
}

void JetEvaluator::fillOutputNtuples(PHCompositeNode *topNode) {
  
  if (verbosity > 2) cout << "JetEvaluator::fillOutputNtuples() entered" << endl;

  JetEvalStack jetevalstack(topNode,_recojetname,_truthjetname); 
  JetRecoEval*   recoeval = jetevalstack.get_reco_eval();
  //JetTruthEval* trutheval = jetevalstack.get_truth_eval();
 
  //-------------------------
  // fill the reco jet ntuple
  //-------------------------

  if (_do_recojet_eval) {
    if (verbosity > 1) cout << "JetEvaluator::filling recojet ntuple..." << endl;

    JetMap* recojets = findNode::getClass<JetMap>(topNode,_recojetname.c_str());
    if (!recojets) {
      cerr << PHWHERE << " ERROR: Can't find " << _recojetname << endl;
      exit(-1);
    }
  
    // for every recojet
    for (JetMap::Iter iter = recojets->begin();
	 iter != recojets->end();
	 ++iter) {
      Jet* recojet = iter->second;
      Jet* truthjet = recoeval->max_truth_jet_by_energy(recojet);

      float id    = recojet->get_id();
      float ncomp = recojet->size_comp();
      float eta   = recojet->get_eta();
      float phi   = recojet->get_phi();
      float e     = recojet->get_e();
      float pt    = recojet->get_pt();

      float gid        = NAN;
      float gncomp     = NAN;
      float geta       = NAN;
      float gphi       = NAN;
      float ge         = NAN;
      float gpt        = NAN;
      float efromtruth = NAN;

      if (truthjet) {
	gid    = truthjet->get_id();
	gncomp = truthjet->size_comp();
	geta   = truthjet->get_eta();
	gphi   = truthjet->get_phi();
	ge     = truthjet->get_e();
	gpt    = truthjet->get_pt();
	efromtruth = recoeval->get_energy_contribution(recojet,truthjet);
      }
      
      float recojet_data[14] = {_ievent,
				id,
				ncomp,
				eta,
				phi,
				e,
				pt,
				gid,
				gncomp,
				geta,
				gphi,
				ge,
				gpt,
				efromtruth
      };

      _ntp_recojet->Fill(recojet_data);
    }
  }

  //-------------------------
  // fill the truth jet ntuple
  //-------------------------

  if (_do_truthjet_eval) {
    if (verbosity > 1) cout << "JetEvaluator::filling truthjet ntuple..." << endl;

    JetMap* truthjets = findNode::getClass<JetMap>(topNode,_truthjetname.c_str());
    if (!truthjets) {
      cerr << PHWHERE << " ERROR: Can't find " << _truthjetname << endl;
      exit(-1);
    }
  
    // for every truthjet
    for (JetMap::Iter iter = truthjets->begin();
	 iter != truthjets->end();
	 ++iter) {
      Jet* truthjet = iter->second;
      Jet* recojet = recoeval->best_jet_from(truthjet);

      float gid    = truthjet->get_id();
      float gncomp = truthjet->size_comp();
      float geta   = truthjet->get_eta();
      float gphi   = truthjet->get_phi();
      float ge     = truthjet->get_e();
      float gpt    = truthjet->get_pt();

      float id         = NAN;
      float ncomp      = NAN;
      float eta        = NAN;
      float phi        = NAN;
      float e          = NAN;
      float pt         = NAN;
      float efromtruth = NAN;

      if (recojet) {
	id         = recojet->get_id();
	ncomp      = recojet->size_comp();
	eta        = recojet->get_eta();
	phi        = recojet->get_phi();
	e          = recojet->get_e();
	pt         = recojet->get_pt();
	efromtruth = recoeval->get_energy_contribution(recojet,truthjet);
      }
      
      float truthjet_data[14] = {_ievent,
				 gid,
				 gncomp,
				 geta,
				 gphi,
				 ge,
				 gpt,
				 id,
				 ncomp,
				 eta,
				 phi,
				 e,
				 pt,
				 efromtruth
      };

      _ntp_truthjet->Fill(truthjet_data);
    }
  }

  return;
}
