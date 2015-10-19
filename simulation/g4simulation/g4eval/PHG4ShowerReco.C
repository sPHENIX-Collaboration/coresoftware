
#include "PHG4SHowerReco.h"

#include "CaloTruthEval.h"

#include <g4main/PHG4ShowerMap.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <cemc/RawTowerContainer.h>
#include <cemc/RawTower.h>

#include <phool/PHCompositeNode.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/getClass.h>

#include <iostream>
#include <set>
#include <cmath>

using namespace std;

PHG4ShowerReco::PHG4ShowerReco(const string &name, const string &filename) :
  SubsysReco("PHG4ShowerReco"),
  _ievent(0),
  _volume_names(),
  _volume_evals() {
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

  // create or update the evaluators
  if (_volume_evals.empty()) {
    for (std::map<PHG4Shower::VOLUME,std::string>::iterator iter = _volume_names.begin();
	 iter != _volume_names.end();
	 ++iter) {
      std::string caloname = iter->second;
      _volume_evals.insert(make_pair(caloname,new CaloEvalStack(topNode,caloname)));
    }
  } else {
    for (std::map<PHG4Shower::VOLUME,CaloTruthEval*>::iterator iter = _volume_evals.begin();
	 iter != _volume_evals.end();
	 ++iter) {
      CaloTruthEval* eval = iter->second;
      eval->next_event(topNode);
    }
  }

  // --- create the g4shower objects -------------------------------------------
  
  // loop over all truth primary particles
  // loop over all volumes
  // get the g4hits from this particle in this volume

  // create a new g4shower from the g4hits
  // add the g4shower to the node tree

  // --- update rawtower ancestry ----------------------------------------------
  
  ++_ievent;
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4ShowerReco::End(PHCompositeNode *topNode) {

  for (std::map<PHG4Shower::VOLUME,CaloTruthEval*>::iterator iter = _volume_evals.begin();
       iter != _volume_evals.end();
       ++iter) {
    CaloEvalStack* eval = iter->second;
    if (eval) delete eval;
  }
  _volume_evals.clear();

  return;
}

int PHG4ShowerReco::CreateNodes(PHCompositeNode *topNode) {

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
    RawTowerContainer* towers = findNode::getClass<RawTowerContainer>(topNode,nodename.c_str());
    if (towers) {
      _volume_towers.insert(make_pair(volid,towers));
    }    
  }

  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if (!_truthinfo) {
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
