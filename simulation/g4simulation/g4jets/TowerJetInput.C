
#include "TowerJetInput.h"

#include "JetInput.h"
#include "Jet.h"
#include "JetV1.h"

#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <fun4all/getClass.h>

// PHENIX Geant4 includes
//#include <g4cemc/SvtxTrackMap.h>

// standard includes
#include <iostream>
#include <vector>
#include <cstdlib>

using namespace std;

TowerJetInput::TowerJetInput(Jet::SRC input)
  : _verbosity(0),
    _input(input) {
}

std::vector<Jet*> TowerJetInput::get_input(PHCompositeNode *topNode) {
  
  if (_verbosity > 0) cout << "TowerJetInput::process_event -- entered" << endl;

  // // Pull the reconstructed tower information off the node tree...
  // SvtxTrackMap *trackmap = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
  // if (!trackmap) {
  //   cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap" << endl;
  //   return std::vector<Jet*>();
  // }

  std::vector<Jet*> pseudojets;
  // for (SvtxTrackMap::ConstIter iter = trackmap->begin(); 
  //      iter != trackmap->end(); 
  //      ++iter) {
  //   const SvtxTrack *track = &iter->second;

  //   Jet *jet = new JetV1();
  //   jet->set_px(track->get3Momentum(0));
  //   jet->set_py(track->get3Momentum(1));
  //   jet->set_pz(track->get3Momentum(2));
  //   jet->set_e(track->getMomentum());
  //   jet->insert_comp(Jet::TRACK,track->getTrackID());
  //   pseudojets.push_back(jet);
  // }

  if (_verbosity > 0) cout << "TowerJetInput::process_event -- exited" << endl;

  return pseudojets;
}
