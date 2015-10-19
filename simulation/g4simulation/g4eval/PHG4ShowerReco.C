
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
#include <fun4all/getClass.h>

#include <TPrincipal.h>

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
  _volume_towers(),
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
   for (std::map<PHG4Shower::VOLUME,RawTowerContainer*>::iterator iter = _volume_towers.begin();
	 iter != _volume_towers.end();
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
  
  // loop over all truth primary particles
  const PHG4TruthInfoContainer::Map &map = _truth_info->GetPrimaryMap();
  for (PHG4TruthInfoContainer::ConstIterator iter = map.begin(); 
       iter != map.end(); 
       ++iter) {
    PHG4Particle* primary = iter->second;

    // does the output already contain a shower for this particle?
    // if so we don't need to create a new one
    bool exists = false;
    for (PHG4ShowerMap::Iter jter = _shower_map->begin();
	 jter != _shower_map->end();
	 ++jter) {
      PHG4Shower *shower = jter->second;
      if (shower->get_primary_id() == primary->get_track_id()) {
	exists = true;
	break;
      }
    }
    if (exists) continue;

    PHG4Shower_v1 shower;
    shower.set_primary_id(primary->get_track_id());    

    TPrincipal pca(3); // principal component analysis object

    // loop over all volumes with evals
    cout << endl;
    for (std::map<PHG4Shower::VOLUME,CaloTruthEval*>::iterator iter = _volume_truthevals.begin();
	 iter != _volume_truthevals.end();
	 ++iter) {
      PHG4Shower::VOLUME volid = iter->first;
      CaloTruthEval* eval = iter->second;

      float edep = 0.0;
      float eion = 0.0;
      float light_yield = 0.0;
      
      // get the g4hits from this particle in this volume
      
      std::set<PHG4Hit*> g4hits = eval->get_shower_from_primary(primary);     
      for (std::set<PHG4Hit*>::iterator jter = g4hits.begin();
	   jter != g4hits.end();
	   ++jter) {
	PHG4Hit* g4hit = *jter;

	Double_t data0[3] = {g4hit->get_x(0),
			     g4hit->get_y(0),
			     g4hit->get_z(0)};
	Double_t* pdata0 = &data0[0];
	
	if (!isnan(data0[0]) && !isnan(data0[1]) && !isnan(data0[2])) {
	  pca.AddRow(pdata0);
	}

	Double_t data1[3] = {g4hit->get_x(1),
			     g4hit->get_y(1),
			     g4hit->get_z(1)};
	Double_t* pdata1 = &data1[0];
	
	if (!isnan(data1[0]) && !isnan(data1[1]) && !isnan(data1[2])) {
	  pca.AddRow(pdata1);
	}
	
	if (!isnan(g4hit->get_edep()))               edep += g4hit->get_edep();
	if (!isnan(g4hit->get_eion()))               eion += g4hit->get_eion();
	if (!isnan(g4hit->get_light_yield())) light_yield += g4hit->get_light_yield();	
      } // g4hit loop

      shower.set_edep(volid,edep);
      shower.set_eion(volid,eion);
      shower.set_light_yield(volid,light_yield);     

    } // volume loop

    // fill shower with position and covariance information
    const TVectorD* MEAN  = pca.GetMeanValues();
    const TMatrixD* COVAR = pca.GetCovarianceMatrix();
    
    shower.set_x((*MEAN)[0]);
    shower.set_y((*MEAN)[1]);
    shower.set_z((*MEAN)[2]);

    for (unsigned int i = 0; i < 3; ++i) {
      for (unsigned int j = 0; j <= i; ++j) {
	shower.set_covar(i,j,(*COVAR)[i][j]);
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
  for (PHG4ShowerMap::Iter iter = _shower_map->begin();
       iter != _shower_map->end();
       ++iter) {
    PHG4Shower *shower = iter->second;
    unsigned int showerid = shower->get_id();
    int primaryid = shower->get_primary_id();
    _primaryid_showerid_map.insert(make_pair(primaryid,showerid));
  }
  
  // --- update rawtower ancestry ----------------------------------------------
  for (std::map<PHG4Shower::VOLUME,RawTowerContainer*>::iterator iter = _volume_towers.begin();
       iter != _volume_towers.end();
       ++iter) {
    PHG4Shower::VOLUME volid = iter->first;
    RawTowerContainer* towers = iter->second;

    // loop over all towers...
    for (RawTowerContainer::Iterator iter = towers->getTowers().first;
	 iter != towers->getTowers().second;
	 ++iter) {
      RawTower* tower = iter->second;

      // get all primaries that contribute to tower
      std::set<PHG4Particle*> primaries = _volume_towerevals[volid]->all_truth_primaries(tower);
           
      // loop over primaries
      for (std::set<PHG4Particle*>::iterator jter = primaries.begin();
	   jter != primaries.end();
	   ++jter) {
	PHG4Particle* primary = *jter;
	unsigned int showerid = _primaryid_showerid_map[primary->get_track_id()];

	float edep = _volume_towerevals[volid]->get_energy_contribution(tower,primary);
	
	// insert this ancestry onto the tower
	tower->add_eshower(showerid, edep);	
      }
    }
  }
      
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
    RawTowerContainer* calib_towers = findNode::getClass<RawTowerContainer>(topNode,nodename.c_str());
    if (calib_towers) {
      _volume_towers.insert(make_pair(volid,calib_towers));
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
