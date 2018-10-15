#include "PHG4TrackGhostRejection.h"

#include "SvtxVertexMap.h"
#include "SvtxVertex.h"
#include "SvtxTrackMap.h"
#include "SvtxTrack.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <iostream>
#include <algorithm>

using namespace std;

bool hit_sort(const unsigned int i, const unsigned int j) { return (i < j);}

PHG4TrackGhostRejection::PHG4TrackGhostRejection(const int nlayers, const string &name) :
  SubsysReco(name),
  _g4tracks(NULL),
  _nlayers(nlayers),
  _max_shared_hits(_nlayers)
{
  _layer_enabled.assign(_nlayers,true);
  _overlapping.clear();
  _candidates.clear();
}

int PHG4TrackGhostRejection::Init(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TrackGhostRejection::InitRun(PHCompositeNode *topNode)
{
  if (Verbosity() > 0) {
    cout << "================== PHG4TrackGhostRejection::InitRun() =====================" << endl;
    cout << " Maximum allowed shared hits: " << _max_shared_hits << endl;
    for (unsigned int i=0;i<_layer_enabled.size();++i) {
      cout << " Enabled for hits in layer #" << i << ": " << boolalpha << _layer_enabled[i] << noboolalpha << endl;
    }
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TrackGhostRejection::process_event(PHCompositeNode *topNode)
{
  if(Verbosity() > 0) cout << "PHG4TrackGhostRejection::process_event -- entered" << endl;

  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  // Pull the reconstructed track information off the node tree...
  _g4tracks = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if(!_g4tracks) 
    {
      cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap." << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  if (Verbosity() > 1) {
    _g4tracks->identify();
    for (SvtxTrackMap::Iter iter = _g4tracks->begin();
	 iter != _g4tracks->end();
	 ++iter) {
      SvtxTrack *track = iter->second;
      track->identify();
    }
  }

  //----------------------------
  // Sort the hits on each track
  //----------------------------

  _candidates.clear();
  
  for (SvtxTrackMap::Iter iter = _g4tracks->begin();
       iter != _g4tracks->end();
       ++iter) {

    SvtxTrack* track = iter->second;
  
    PHG4TrackCandidate combo;

    combo.trackid = track->get_id();
    combo.nhits = track->size_clusters();

    for (SvtxTrack::ConstClusterIter iter = track->begin_clusters();
	 iter != track->end_clusters();
	 ++iter) {
      unsigned int cluster_id = *iter;
      combo.hitids.push_back(cluster_id);
    }
      
    if (track->get_ndf() != 0) {
      combo.chisq = track->get_chisq()/track->get_ndf();
    }

    combo.keep = true;

    // sort the hits by index
    stable_sort(combo.hitids.begin(),combo.hitids.end(),hit_sort);

    _candidates.push_back(combo);
  }

  //---------------------
  // Fill the overlap map
  //---------------------

  _overlapping.clear();
  for(unsigned int i = 0; i < _candidates.size(); i++)
    {
      for(unsigned int j = i+1; j < _candidates.size(); j++)
  	{ 
	  // determine the maximum length of the anticipated storage
	  unsigned maxhits = _candidates[i].hitids.size();
	  if(_candidates[j].hitids.size() > maxhits)
	    {
	      maxhits = _candidates[j].hitids.size();
	    }

	  // create the difference storage
	  std::vector<unsigned int> diff;
	  diff.assign(maxhits,0);

	  // run the difference algorithm
	  std::vector<unsigned int>::iterator it = diff.begin();
	  it = std::set_difference(_candidates[i].hitids.begin(),_candidates[i].hitids.end(),
				   _candidates[j].hitids.begin(),_candidates[j].hitids.end(),
				   diff.begin());

	  // calculate the overlap
	  unsigned int overlap = maxhits - int(it - diff.begin());

	  // insert an overlapping pair into the map
	  if(overlap > _max_shared_hits) _overlapping.insert(std::make_pair(i,j));	  
   	}
    }

  //----------------------
  // Flag the ghost tracks
  //----------------------

  std::multimap<unsigned,unsigned int>::iterator iter;
  for (iter = _overlapping.begin(); iter != _overlapping.end(); iter++) {

    unsigned int key = iter->first;
    unsigned int value = iter->second;

    if (_candidates[key].nhits > _candidates[value].nhits) {
      // prefer longer track
      _candidates[value].keep = false;
    } else if (_candidates[key].nhits < _candidates[value].nhits) {
      // prefer longer track
      _candidates[key].nhits = false;
    } else {
      // choose between equal length tracks by chisq/dof
      if (_candidates[key].chisq < _candidates[value].chisq) {
	_candidates[value].keep = false;
      } else {
	_candidates[key].keep = false;
      }
    }
  }

  //------------------------
  // Remove the ghost tracks
  //------------------------

  SvtxVertexMap* vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  
  int initial_size = _g4tracks->size();

  // loop over container and delete!
  for (unsigned int i = 0; i < _candidates.size(); i++) {
    if (!_candidates[i].keep) {
      // look for the track to delete
      if (_g4tracks->find(_candidates[i].trackid) != _g4tracks->end()) {
	_g4tracks->erase(_candidates[i].trackid);

	// also remove the track id from any vertex that contains this track
	if (!vertexmap) continue;
	for (SvtxVertexMap::Iter iter = vertexmap->begin();
	     iter != vertexmap->end();
	     ++iter) {
	  SvtxVertex* vertex = iter->second;
	  vertex->erase_track(_candidates[i].trackid);
	}
      }
    }
  }

  if (Verbosity() > 1) {
    _g4tracks->identify();
    for (SvtxTrackMap::Iter iter = _g4tracks->begin();
	 iter != _g4tracks->end();
	 ++iter) {
      SvtxTrack *track = iter->second;
      track->identify();
    }
  }

  if(Verbosity() > 0)
    cout << "PHG4TrackGhostRejection - rejected and removed " 
         << initial_size - _g4tracks->size()
         << " tracks" << endl;;
  
  if(Verbosity() > 0) cout << "PHG4TrackGhostRejection::process_event -- exited" << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TrackGhostRejection::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

