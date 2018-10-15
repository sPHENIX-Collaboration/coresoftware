#include "PHG4SvtxMomentumRecal.h"
#include "SvtxTrackMap.h"
#include "SvtxTrack.h"

// PHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

// standard includes
#include <iostream>
#include <vector>

using namespace std;

PHG4SvtxMomentumRecal::PHG4SvtxMomentumRecal(const string &name,
					     TF1* corr) :
  SubsysReco(name),
  _corr(corr)
{}

int PHG4SvtxMomentumRecal::Init(PHCompositeNode *topNode) 
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4SvtxMomentumRecal::InitRun(PHCompositeNode *topNode) 
{
  if (Verbosity() > 0) {
    cout << "================== PHG4SvtxMomentumRecal::InitRun() =====================" << endl;
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4SvtxMomentumRecal::process_event(PHCompositeNode *topNode)
{
  if(Verbosity() > 1) cout << "PHG4SvtxMomentumRecal::process_event -- entered" << endl;

  if (!_corr) return Fun4AllReturnCodes::EVENT_OK;
  
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  // Pull the reconstructed track information off the node tree...
  SvtxTrackMap* _g4tracks = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!_g4tracks) {
    cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
    
  // loop over all tracks
  for (SvtxTrackMap::Iter iter = _g4tracks->begin();
       iter != _g4tracks->end();
       ++iter) {
    SvtxTrack *track = iter->second;
    
    double rescale = 1.0;

    double pt = track->get_pt();
    
    double xmin = 0.0;
    double xmax = 0.0;
    _corr->GetRange(xmin,xmax);
      
    if ((pt > xmin)&&(pt < xmax)) {
      rescale = _corr->Eval(pt);
    }

    track->set_px( track->get_px() * rescale );
    track->set_py( track->get_py() * rescale );
    track->set_pz( track->get_pz() * rescale );
  } // end track loop

  if (Verbosity() > 1) cout << "PHG4SvtxMomentumRecal::process_event -- exited" << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4SvtxMomentumRecal::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
