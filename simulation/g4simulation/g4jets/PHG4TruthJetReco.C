#include "PHG4TruthJetReco.h"

// PHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <fun4all/getClass.h>

// PHENIX Geant4 includes
#include <JetMap.h>
#include <Jet.h>

// standard includes
#include <iostream>

using namespace std;

PHG4TruthJetReco::PHG4TruthJetReco(const string &name)
  : SubsysReco(name) {
  verbosity = 0;
}

int PHG4TruthJetReco::Init(PHCompositeNode *topNode) {
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TruthJetReco::InitRun(PHCompositeNode *topNode) {
  
  if (verbosity >= 0) {
    cout << "===================== PHG4TruthJetReco::InitRun() =========================" << endl;
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TruthJetReco::process_event(PHCompositeNode *topNode) {
  
  if (verbosity > 0) cout << "PHG4TruthJetReco::process_event -- entered" << endl;

  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------
  /*
  // Pull the reconstructed track information off the node tree...
  _g4tracks = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if(!_g4tracks) 
    {
      cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap." << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  if (verbosity > 1) {
    _g4tracks->identify();
    for (SvtxTrackMap::Iter iter = _g4tracks->begin();
	 iter != _g4tracks->end();
	 ++iter) {
      SvtxTrack *track = &iter->second;
      track->identify();
    }
  }
  */
  
  if (verbosity > 0) cout << "PHG4TruthJetReco::process_event -- exited" << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TruthJetReco::End(PHCompositeNode *topNode) {
  return Fun4AllReturnCodes::EVENT_OK;
}

