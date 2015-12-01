
#include "PHG4ShowerReco.h"

#include "CaloTruthEval.h"
#include "CaloRawTowerEval.h"

#include <g4main/PHG4ShowerMap.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4Shower_v1.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4cemc/RawTowerContainer.h>
#include <g4cemc/RawTower.h>

#include <phool/PHCompositeNode.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <Eigen/Dense>

#include <iostream>
#include <set>
#include <cmath>
#include <cstdlib>

using namespace std;

PHG4ShowerReco::PHG4ShowerReco(const string &name) :
  SubsysReco("PHG4ShowerReco"),
  _ievent(0),
  _truth_info(),
  _volume_names(),
  _volume_truthevals(),
  _volume_g4hits(),
  _volume_simtowers(),
  _volume_rawtowers(),
  _volume_calibtowers(),
  _shower_map() {
  verbosity = 0;

  _volume_names.insert(make_pair(PHG4Shower::CEMC_ELECTRONICS,"CEMC_ELECTRONICS"));
  _volume_names.insert(make_pair(PHG4Shower::CEMC            ,"CEMC"));
  _volume_names.insert(make_pair(PHG4Shower::ABSORBER_CEMC   ,"ABSORBER_CEMC"));
  _volume_names.insert(make_pair(PHG4Shower::CEMC_SPT        ,"CEMC_SPT"));

  _volume_names.insert(make_pair(PHG4Shower::ABSORBER_HCALIN ,"ABSORBER_HCALIN"));
  _volume_names.insert(make_pair(PHG4Shower::HCALIN          ,"HCALIN"));
  _volume_names.insert(make_pair(PHG4Shower::HCALIN_SPT      ,"HCALIN_SPT"));

  _volume_names.insert(make_pair(PHG4Shower::MAGNET          ,"MAGNET"));

  _volume_names.insert(make_pair(PHG4Shower::ABSORBER_HCALOUT,"ABSORBER_HCALOUT"));
  _volume_names.insert(make_pair(PHG4Shower::HCALOUT         ,"HCALOUT"));
  _volume_names.insert(make_pair(PHG4Shower::HCALOUT_SPT     ,"HCALOUT_SPT"));

  _volume_names.insert(make_pair(PHG4Shower::BH_1            ,"BH_1"));

  /// \todo expand this list to include forward detectors
}

int PHG4ShowerReco::Init(PHCompositeNode *topNode) {

  _ievent = 0;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4ShowerReco::InitRun(PHCompositeNode *topNode) {

  int code = CreateNodes(topNode);
  if (code != Fun4AllReturnCodes::EVENT_OK) return code;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4ShowerReco::process_event(PHCompositeNode *topNode) {

  int code = GetNodes(topNode);
  if (code != Fun4AllReturnCodes::EVENT_OK) return code;

  if ((verbosity > 0)&&(_ievent%100==0)) {
    cout << "PHG4ShowerReco::process_event - Event = " << _ievent << endl;
  }

  // create or update the tower evaluators for each volume
  if (_volume_towerevals.empty()) {
   for (std::map<PHG4Shower::VOLUME,RawTowerContainer*>::iterator iter = _volume_calibtowers.begin();
	 iter != _volume_calibtowers.end();
	 ++iter) {
     PHG4Shower::VOLUME volid = iter->first;
     if (_volume_towerevals.find(volid) != _volume_towerevals.end()) continue;
     _volume_towerevals.insert(make_pair(volid, new CaloRawTowerEval(topNode,_volume_names[volid])));
   }
  } else {
    for (std::map<PHG4Shower::VOLUME,CaloRawTowerEval*>::iterator iter = _volume_towerevals.begin();
	 iter != _volume_towerevals.end();
	 ++iter) {
      CaloRawTowerEval* towereval = iter->second;
      towereval->next_event(topNode);
   }
  }

  // create or update the truth evaluators for each volume
  if (_volume_truthevals.empty()) {
    for (std::map<PHG4Shower::VOLUME,PHG4HitContainer*>::iterator iter = _volume_g4hits.begin();
	 iter != _volume_g4hits.end();
	 ++iter) {
      PHG4Shower::VOLUME volid = iter->first;
      _volume_truthevals.insert(make_pair(volid, new CaloTruthEval(topNode,_volume_names[volid])));
    }
  } else {
    for (std::map<PHG4Shower::VOLUME,CaloTruthEval*>::iterator iter = _volume_truthevals.begin();
	 iter != _volume_truthevals.end();
	 ++iter) {
      CaloTruthEval* eval = iter->second;
      eval->next_event(topNode);
    }
  }

  // --- create the g4shower objects -------------------------------------------

  // which primary id will be stored on the g4shower?
  // the copy in the Map or the PrimaryMap?

  // typically we have been using
  // the id in the Map and confirming with the PrimaryMap since this is how
  // the primary ids on the secondary particles is tracked.

  // it is also faster to hop from the entry in the Map to the PrimaryMap
  // than in the reverse direction, also supporting recording the Map index

  // ~310 kilobytes / central hijing event right now
  // but far too slow due to multiple calls to the evaluator and
  // no sharing of caches between objects (determined via callgrind)

  // loop over all truth primary particles
  const PHG4TruthInfoContainer::Map &map = _truth_info->GetMap();
  for (PHG4TruthInfoContainer::ConstIterator primary_iter = map.begin();
       primary_iter != map.end();
       ++primary_iter) {
    PHG4Particle* primary = primary_iter->second;

    if (!_volume_truthevals.begin()->second->is_primary(primary)) continue;

    // does the output already contain a shower for this particle?
    // if so we don't need to create a new one
    bool exists = false;
    for (PHG4ShowerMap::Iter shower_iter = _shower_map->begin();
      	 shower_iter != _shower_map->end();
      	 ++shower_iter) {
       PHG4Shower *shower = shower_iter->second;
       PHG4Particle *candidate = _truth_info->GetParticle(shower->get_primary_id());
       if (_volume_truthevals.begin()->second->are_same_particle(primary,candidate)) {
	 exists = true;
	 break;
       }
    }
    if (exists) continue;

    PHG4Shower_v1 shower;
    shower.set_primary_id(primary->get_track_id());

    // Data structures to hold weighted
    std::vector<std::vector<float> > points;
    std::vector<float> weights;
    float sumw = 0.0;
    float sumw2 = 0.0;

    // loop over all volumes with evals
    for (std::map<PHG4Shower::VOLUME,CaloTruthEval*>::iterator volume_iter = _volume_truthevals.begin();
	 volume_iter != _volume_truthevals.end();
	 ++volume_iter) {
      PHG4Shower::VOLUME volid = volume_iter->first;
      CaloTruthEval* eval = volume_iter->second;

      float edep = 0.0;
      float eion = 0.0;
      float light_yield = 0.0;

      // get the g4hits from this particle in this volume

      std::set<PHG4Hit*> g4hits = eval->get_shower_from_primary(primary);
      for (std::set<PHG4Hit*>::iterator g4hit_iter = g4hits.begin();
	   g4hit_iter != g4hits.end();
	   ++g4hit_iter) {
	PHG4Hit* g4hit = *g4hit_iter;

	if (!isnan(g4hit->get_x(0)) &&
	    !isnan(g4hit->get_y(0)) &&
	    !isnan(g4hit->get_z(0))) {

	  std::vector<float> entry(3);
	  entry[0] = g4hit->get_x(0);
	  entry[1] = g4hit->get_y(0);
	  entry[2] = g4hit->get_z(0);

	  points.push_back(entry);
	  float w = g4hit->get_edep();
	  weights.push_back(w);
	  sumw += w;
	  sumw2 += w*w;
	}

	if (!isnan(g4hit->get_x(1)) &&
	    !isnan(g4hit->get_y(1)) &&
	    !isnan(g4hit->get_z(1))) {

	  std::vector<float> entry(3);
	  entry[0] = g4hit->get_x(1);
	  entry[1] = g4hit->get_y(1);
	  entry[2] = g4hit->get_z(1);

	  points.push_back(entry);
	  float w = g4hit->get_edep();
	  weights.push_back(w);
	  sumw += w;
	  sumw2 += w*w;
	}

	if (!isnan(g4hit->get_edep()))               edep += g4hit->get_edep();
	if (!isnan(g4hit->get_eion()))               eion += g4hit->get_eion();
	if (!isnan(g4hit->get_light_yield())) light_yield += g4hit->get_light_yield();
      } // g4hit loop

      if (edep != 0.0)        shower.set_edep(volid,edep);
      if (eion != 0.0)        shower.set_eion(volid,eion);
      if (light_yield != 0.0) shower.set_light_yield(volid,light_yield);
    } // volume loop

    if (points.empty()) continue;

    // fill Eigen matrices to compute wPCA
    // resizing these non-destructively is expensive
    // so I fill vectors and then copy
    Eigen::Matrix<double, Eigen::Dynamic, 3> X(points.size(),3);
    Eigen::Matrix<double, Eigen::Dynamic, 1> W(weights.size(),1);

    for (unsigned int i=0; i<points.size(); ++i) {
      for (unsigned int j=0; j<3; ++j)  {
	X(i,j) = points[i][j];
      }
      W(i,0) = weights[i];
    }

    // mean value of shower
    double prefactor = 1.0 / sumw;
    Eigen::Matrix<double, 1, 3> mean = prefactor * W.transpose() * X;

    // compute residual relative to the mean
    for (unsigned int i=0; i<points.size(); ++i) {
      for (unsigned int j=0; j<3; ++j) X(i,j) = points[i][j] - mean(0,j);
    }

    // weighted covariance matrix
    prefactor = sumw / (pow(sumw,2) - sumw2); // effectivelly 1/(N-1) when w_i = 1.0
    Eigen::Matrix<double, 3, 3> covar = prefactor * (X.transpose() * W.asDiagonal() * X);

    shower.set_x(mean(0,0));
    shower.set_y(mean(0,1));
    shower.set_z(mean(0,2));

    for (unsigned int i = 0; i < 3; ++i) {
      for (unsigned int j = 0; j <= i; ++j) {
	shower.set_covar(i,j,covar(i,j));
      }
    }

    PHG4Shower* ptr = _shower_map->insert(&shower);
    if (!ptr->isValid()) {
      static bool first = true;
      if (first) {
	cout << PHWHERE << "ERROR: Invalid PHG4Showers are being produced" << endl;
	ptr->identify();
	first = false;
      }
    }
  } // primary particle loop

  // loop over all showers and create a map to trace quickly between primary id and shower id
  std::map<int, unsigned int> _primaryid_showerid_map;
  for (PHG4ShowerMap::Iter lookup_iter = _shower_map->begin();
       lookup_iter != _shower_map->end();
       ++lookup_iter) {
    PHG4Shower *shower = lookup_iter->second;
    unsigned int showerid = shower->get_id();
    int primaryid = shower->get_primary_id();
    _primaryid_showerid_map.insert(make_pair(primaryid,showerid));
  }

  // --- update rawtower ancestry ----------------------------------------------

  std::set<std::map<PHG4Shower::VOLUME,RawTowerContainer*>*> volume_towers;
  volume_towers.insert(&_volume_simtowers);
  volume_towers.insert(&_volume_rawtowers);
  volume_towers.insert(&_volume_calibtowers);

  // loop over the different kinds of towers
  for (std::set<std::map<PHG4Shower::VOLUME,RawTowerContainer*>*>::iterator kind_iter = volume_towers.begin();
       kind_iter != volume_towers.end();
       ++kind_iter) {
    std::map<PHG4Shower::VOLUME,RawTowerContainer*> *volume_tower = *kind_iter;

    for (std::map<PHG4Shower::VOLUME,RawTowerContainer*>::iterator volume_iter = volume_tower->begin();
	 volume_iter != volume_tower->end();
	 ++volume_iter) {
      PHG4Shower::VOLUME volid = volume_iter->first;
      RawTowerContainer* towers = volume_iter->second;

      // loop over all towers...
      for (RawTowerContainer::Iterator tower_iter = towers->getTowers().first;
	   tower_iter != towers->getTowers().second;
	   ++tower_iter) {
	RawTower* tower = tower_iter->second;

	// get all primaries that contribute to tower
	std::set<PHG4Particle*> primaries = _volume_towerevals[volid]->all_truth_primaries(tower);

	// loop over primaries
	for (std::set<PHG4Particle*>::iterator primary_iter = primaries.begin();
	     primary_iter != primaries.end();
	     ++primary_iter) {
	  PHG4Particle* primary = *primary_iter;
	  unsigned int showerid = _primaryid_showerid_map[primary->get_track_id()];

	  bool exists = false;
	  for (RawTower::ShowerConstIterator shower_iter = tower->get_g4showers().first;
	       shower_iter != tower->get_g4showers().second;
	       ++shower_iter) {
	    if (showerid == shower_iter->first) {
	      exists = true;
	      break;
	    }
	  }

	  if (exists) continue;

	  float edep = _volume_towerevals[volid]->get_energy_contribution(tower,primary);

	  // insert this ancestry onto the tower
	  tower->add_eshower(showerid, edep);
	}
      }
    }
  }

  // remove all g4particles and vertexes that don't trace to an uncompressed volume
  // loop over all g4 secondary particles in Map
  // did this particle create any g4hit in the other containers
  
  ++_ievent;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4ShowerReco::End(PHCompositeNode *topNode) {

  for (std::map<PHG4Shower::VOLUME,CaloRawTowerEval*>::iterator iter = _volume_towerevals.begin();
       iter != _volume_towerevals.end();
       ++iter) {
    CaloRawTowerEval* eval = iter->second;
    if (eval) delete eval;
  }
  _volume_towerevals.clear();

  for (std::map<PHG4Shower::VOLUME,CaloTruthEval*>::iterator iter = _volume_truthevals.begin();
       iter != _volume_truthevals.end();
       ++iter) {
    CaloTruthEval* eval = iter->second;
    if (eval) delete eval;
  }
  _volume_truthevals.clear();

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4ShowerReco::CreateNodes(PHCompositeNode *topNode) {

  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode
    = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode","DST"));
  if (!dstNode) {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Create the Shower node if required
  PHG4ShowerMap *showers = findNode::getClass<PHG4ShowerMap>(topNode,"PHG4ShowerMap");
  if (!showers) {
    showers = new PHG4ShowerMap();
    PHIODataNode<PHObject> *PHG4ShowerMapNode =
      new PHIODataNode<PHObject>(showers, "PHG4ShowerMap", "PHObject");
    dstNode->addNode(PHG4ShowerMapNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4ShowerReco::GetNodes(PHCompositeNode *topNode) {

  for (std::map<PHG4Shower::VOLUME,std::string>::iterator iter = _volume_names.begin();
       iter != _volume_names.end();
       ++iter) {
    PHG4Shower::VOLUME volid = iter->first;
    std::string volname = iter->second;

    std::string nodename = "G4HIT_" + volname;
    PHG4HitContainer* g4hits = findNode::getClass<PHG4HitContainer>(topNode,nodename.c_str());
    if (g4hits) {
      _volume_g4hits.insert(make_pair(volid,g4hits));
    }

    nodename = "TOWER_SIM_" + volname;
    RawTowerContainer* sim_towers = findNode::getClass<RawTowerContainer>(topNode,nodename.c_str());
    if (sim_towers) {
      _volume_simtowers.insert(make_pair(volid,sim_towers));
    }

    nodename = "TOWER_RAW_" + volname;
    RawTowerContainer* raw_towers = findNode::getClass<RawTowerContainer>(topNode,nodename.c_str());
    if (raw_towers) {
      _volume_rawtowers.insert(make_pair(volid,raw_towers));
    }

    nodename = "TOWER_CALIB_" + volname;
    RawTowerContainer* calib_towers = findNode::getClass<RawTowerContainer>(topNode,nodename.c_str());
    if (calib_towers) {
      _volume_calibtowers.insert(make_pair(volid,calib_towers));
    }
  }

  _truth_info = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if (!_truth_info) {
    cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
    exit(-1);
  }

  _shower_map = findNode::getClass<PHG4ShowerMap>(topNode,"PHG4ShowerMap");
  if (!_shower_map) {
    cerr << PHWHERE << " ERROR: Can't find PHG4SHOWERMAP" << endl;
    exit(-1);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
