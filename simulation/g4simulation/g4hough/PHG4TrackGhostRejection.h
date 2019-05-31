#ifndef G4HOUGH_PHG4TRACKGHOSTREJECTION_H
#define G4HOUGH_PHG4TRACKGHOSTREJECTION_H

//===========================================================
/// \file PHG4TrackGhostRejection.h
/// \brief Quick N-shared hit rejection routine
/// \author Mike McCumber
//===========================================================

// PHENIX includes
#include <fun4all/SubsysReco.h>

// standard includes
#include <map>
#include <string>
#include <vector>

// forward declarations
class PHCompositeNode;
class SvtxTrackMap;

class PHG4TrackCandidate
{
  
 public:

  int trackid;
  short nhits;
  std::vector<unsigned int> hitids;
  float chisq;
  bool keep;

};

/// \class PHG4TrackGhostRejection
///
/// \brief Quick N-shared hit rejection routine
///
/// This module runs after the pattern recognition to remove
/// track candidates with a user defined overlap. The steps are:
/// (1) Sort the hits on each track by index
/// (2) Fill a multi-map between a track and all of the overlapping tracks
/// (3) Keep only the best track from the set of overlapping tracks
///
class PHG4TrackGhostRejection : public SubsysReco
{

 public:
 
  PHG4TrackGhostRejection(const int nlayers, const std::string &name = "PHG4TrackGhostRejection");
  virtual ~PHG4TrackGhostRejection() {}
		
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  
  void set_max_shared_hits(const unsigned int nhits) { _max_shared_hits = nhits; }
  unsigned int get_max_shared_hits() { return _max_shared_hits; }

  void set_layer_enabled(const int layer, const bool enabled) {_layer_enabled[layer] = enabled;}
  bool get_layer_enabled(const int layer) const {return _layer_enabled[layer];}

 private:

  SvtxTrackMap *_g4tracks;
  std::vector< PHG4TrackCandidate > _candidates;

  unsigned int _nlayers;
  unsigned int _max_shared_hits;
  std::vector<bool> _layer_enabled;

  std::multimap< unsigned int, unsigned int > _overlapping;

};

#endif // __PHG4TRACKGHOSTREJECTION_H__
