
#include "TrackJetInput.h"

#include "JetInput.h"
#include "Jet.h"
#include "Jetv1.h"

#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

// PHENIX Geant4 includes
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack.h>

// standard includes
#include <iostream>
#include <vector>
#include <cstdlib>

using namespace std;

TrackJetInput::TrackJetInput(Jet::SRC input)
  : _verbosity(0),
    _input(input) {
}

void TrackJetInput::identify(std::ostream& os) {
  os << "   TrackJetInput: SvtxTrackMap to Jet::TRACK" << endl;
}

std::vector<Jet*> TrackJetInput::get_input(PHCompositeNode *topNode) {
  
  if (_verbosity > 0) cout << "TrackJetInput::process_event -- entered" << endl;

  // Pull the reconstructed track information off the node tree...
  SvtxTrackMap *trackmap = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
  if (!trackmap) {
    return std::vector<Jet*>();
  }

  std::vector<Jet*> pseudojets;
  for (SvtxTrackMap::ConstIter iter = trackmap->begin(); 
       iter != trackmap->end(); 
       ++iter) {
    const SvtxTrack *track = iter->second;

    Jet *jet = new Jetv1();
    jet->set_px(track->get_px());
    jet->set_py(track->get_py());
    jet->set_pz(track->get_pz());
    jet->set_e(track->get_p());
    jet->insert_comp(Jet::TRACK,track->get_id());
    pseudojets.push_back(jet);
  }

  if (_verbosity > 0) cout << "TrackJetInput::process_event -- exited" << endl;

  return pseudojets;
}
