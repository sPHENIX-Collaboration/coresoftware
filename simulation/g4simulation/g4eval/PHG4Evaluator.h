#ifndef __PHG4EVALUATOR_H__
#define __PHG4EVALUATOR_H__

//===============================================
/// \file PHG4Evaluator.h
/// \brief Compares reconstructed tracks to truth particles
/// \author Michael P. McCumber (original SVX version)
/// \author Matt Wysocki (translated to sPHENIX simulations)
//===============================================

#include "EvalLinks.h"

// PHENIX includes
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHTimeServer.h>

// PHG4 includes
#include <g4hough/SvtxVertexMap.h>
#include <g4hough/SvtxTrackMap.h>
#include <g4hough/SvtxTrack.h>
#include <g4hough/SvtxClusterMap.h>
#include <g4hough/SvtxCluster.h>
#include <g4hough/SvtxHitMap.h>
#include <g4main/PHG4HitContainer.h>
#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitDefs.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>

// ROOT includes
#include <TNtuple.h>
#include <TFile.h>
#include <TH2D.h>

// standard includes
#include <vector>
#include <map>
#include <float.h>

// forward declarations
class PHCompositeNode;
class PHG4TruthInfoContainer;

/// \class SvxGTrack
///
/// \brief A quick container class for truth particles
///
class SvxGtrack : public PHObject
{
 public:

  SvxGtrack() 
  {
    _is_last = false;
    _track_id = 0;
    _flavor = 0;
    _px = 0.0;
    _py = 0.0;
    _pz = 0.0;
    _vx = 0.0;
    _vy = 0.0;
    _vz = 0.0;
    _fpx = 0.0;
    _fpy = 0.0;
    _fpz = 0.0;
    _fx = 0.0;
    _fy = 0.0;
    _fz = 0.0;
    _standalone_match = 0.0;
    _best_purity = 0.0;
    _best_dpp = FLT_MAX;
    _g4hits.clear();
    _chisq = 0.0;
    _chisqv = 0.0;
    _embed = false;
    _primary = false;
  }
  virtual ~SvxGtrack() {}

  bool get_is_last() {return _is_last;}
  void set_is_last(bool last) {_is_last = last;}

  int get_track_id() {return _track_id;}
  void set_track_id(int track_id) {_track_id = track_id;}

  int get_flavor() {return _flavor;}
  void set_flavor(int flavor) {_flavor = flavor;}

  float get_px() {return _px;}
  void set_px(float px) {_px = px;}

  float get_py() {return _py;}
  void set_py(float py) {_py = py;}

  float get_pz() {return _pz;}
  void set_pz(float pz) {_pz = pz;}

  float get_vx() {return _vx;}
  void set_vx(float vx) {_vx = vx;}

  float get_vy() {return _vy;}
  void set_vy(float vy) {_vy = vy;}

  float get_vz() {return _vz;}
  void set_vz(float vz) {_vz = vz;}

  float get_fpx() {return _fpx;}
  void set_fpx(float fpx) {_fpx = fpx;}

  float get_fpy() {return _fpy;}
  void set_fpy(float fpy) {_fpy = fpy;}

  float get_fpz() {return _fpz;}
  void set_fpz(float fpz) {_fpz = fpz;}

  float get_fx() {return _fx;}
  void set_fx(float fx) {_fx = fx;}

  float get_fy() {return _fy;}
  void set_fy(float fy) {_fy = fy;}

  float get_fz() {return _fz;}
  void set_fz(float fz) {_fz = fz;}

  float get_standalone_match() {return _standalone_match;}
  void set_standalone_match(float match) {_standalone_match = match;}

  float get_best_purity() {return _best_purity;}
  void set_best_purity(float match) {_best_purity = match;}

  float get_best_dpp() {return _best_dpp;}
  void set_best_dpp(float dpp) {_best_dpp = dpp;}

  size_t get_ng4hits() {return _g4hits.size();}
  void clear_g4hits() {_g4hits.clear();}
  PHG4Hit* get_g4hit(int ig4hit) {return _g4hits.at(ig4hit);}
  void add_g4hit(PHG4Hit* g4hit) {_g4hits.push_back(g4hit);}

  PHG4Particle* get_particle() {return _particle;}
  void set_particle(PHG4Particle* particle) {_particle = particle;}

  float get_chisq() {return _chisq;}
  void set_chisq(float chisq) {_chisq = chisq;}

  float get_embed() {return _embed;}
  void set_embed(float embed) {_embed = embed;}

  float get_primary() {return _primary;}
  void set_primary(float primary) {_primary = primary;}

  float get_chisqv() {return _chisqv;}
  void set_chisqv(float chisqv) {_chisqv = chisqv;}

  size_t get_ncluster_chisqs() {return _cluster_chisq.size();}
  void clear_cluster_chisqs() {_cluster_chisq.clear();}
  float get_cluster_chisq(int icluster) {return _cluster_chisq.at(icluster);}
  void add_cluster_chisq(float chisq) {_cluster_chisq.push_back(chisq);}
  
 protected:

  bool _is_last; ///< last particle tossed in the event record (likely an embedded particle)

  int _track_id; ///< truth particle id

  int _flavor; ///< truth particle geant flavor code

  float _px; ///< truth particle x-momentum at creation
  float _py; ///< truth particle y-momentum at creation
  float _pz; ///< truth particle z-momentum at creation

  float _vx; ///< truth particle x-position at creation
  float _vy; ///< truth particle y-position at creation
  float _vz; ///< truth particle z-position at creation

  float _fpx; ///< truth particle x-momentum at exit of final tracking layer
  float _fpy; ///< truth particle y-momentum at exit of final tracking layer
  float _fpz; ///< truth particle z-momentum at exit of final tracking layer

  float _fx; ///< truth particle x-position at exit of final tracking layer
  float _fy; ///< truth particle y-position at exit of final tracking layer
  float _fz; ///< truth particle z-position at exit of final tracking layer

  float _standalone_match; ///< standalone reconstruction matching code
  float _best_purity; ///< best purity of associated reconstructed tracks
  float _best_dpp; ///< best delta momentum / momentum

  std::vector<PHG4Hit*> _g4hits; ///< storage of associated PHG4Hits
  PHG4Particle* _particle;

  float _chisq;      ///< chisq of just the clusters
  float _chisqv;     ///< chisq including the vertex
  std::vector<float> _cluster_chisq;

  bool _embed;
  bool _primary;
};

/// \class PHG4Evaluator
///
/// \brief Compares reconstructed tracks to truth particles
///
/// Plan: This module will trace the reconstructed tracks back to
/// the greatest contributor Monte Carlo particle and then
/// test one against the other.
///
class PHG4Evaluator : public SubsysReco
{
  typedef PHG4HitContainer::Map HitMap;

 public:
 
  PHG4Evaluator(const std::string &name = "PHG4EVALUATOR",
                const std::string &filename = "g4eval.root");
  virtual ~PHG4Evaluator() {};
		
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  //void Verbosity(int v); // inherited from SubsysReco

  // inclusive lower limit
  void set_min_hits(int n) {_min_hits=n;}
  int get_min_hits() {return _min_hits;}

  // exlusive upper limit
  void set_max_hits(int n) {_max_hits=n;}
  int get_max_hits() {return _max_hits;}

 private:

  enum LayerType {CylinderLayer, LadderLayer};
  
  unsigned int _ievent;
  PHTimeServer::timer _timer;   ///< Timer
  std::vector<PHTimeServer::timer> _internal_timer;

  //----------------------------------
  // evaluator output ntuples

  TNtuple *_ntp_event;
  TNtuple *_ntp_g4hit;
  TNtuple *_ntp_cell;
  TNtuple *_ntp_cluster;
  TNtuple *_ntp_gtrack;
  TNtuple *_ntp_track;

  // track query settings
  bool _trackingWasRun;

  unsigned int _nlayers;
  std::vector<unsigned int> _nhits_per_layer;
  std::vector<unsigned int> _nchannels_per_layer;
  unsigned int _min_hits;
  unsigned int _max_hits;

  // evaluator output file
  std::string _filename;
  TFile *_tfile;

  //----------------------------
  // internal evaluator storage

  // pointers to the object lists makes them availabe for all
  // the process_event subroutines
  SvtxTrackMap *_trackList;
  SvtxVertexMap *_vertexList;
  SvtxClusterMap *_clusterList;
  SvtxHitMap *_hitList;
  PHG4CylinderCellContainer *_cellList;
  PHG4CylinderCellGeomContainer *_cellGeos;
  // typedefed from PHG4HitContainer::Map
  HitMap _g4hitList;
  PHG4TruthInfoContainer* _truth_info_container;
  
 

  std::map <unsigned int, LayerType> _layer_type_map; // layer # => ladder or cylinder layer type
  PHG4CylinderGeomContainer *_ladderGeos;             // tilted ladder layer geos
  PHG4CylinderCellContainer *_ladderCellList;         // tilted ladder layer geos
  
  std::vector<SvxGtrack> _gtrack_list;              ///< holds a list of truth particles (convenience storage)
  std::map <PHG4Hit*,SvxGtrack*> _g4hit_gtrack_map; ///< reverse g4hit->gtrack lookup

  std::multimap <SvtxCluster*, PHG4Hit*> _cluster_g4hit_mmap;   ///< holds a map between the SvtxCluster and its principle truth contributor
  std::multimap <PHG4CylinderCell*, PHG4Hit*> _cell_g4hit_mmap;           ///< holds a map between the PHG4CylinderCell and its principle truth contributor
  std::map <SvtxTrack*, SvxGtrack*> _track_gtrack_map;                    ///< holds a map between the SvtxTrack and its principle truth contributor

  std::multimap <PHG4Hit*, SvtxCluster*> _g4hit_cluster_mmap;   ///< reverse lookup
  std::multimap <PHG4Hit*, PHG4CylinderCell*> _g4hit_cell_mmap;           ///< reverse lookup
  std::multimap <SvxGtrack*, SvtxTrack*> _gtrack_track_mmap;              ///< reverse lookup

  std::multimap <SvtxCluster*, PHG4Hit*> _cluster_allg4hits_mmap; ///< convenience storage for cluster to rawhit lookup
  std::multimap <PHG4Hit*, SvtxCluster*> _allg4hits_cluster_mmap; ///< reverse lookup
  std::multimap <SvtxTrack*, SvtxCluster*> _track_cluster_mmap;   ///< convenience storage for track to cluster lookup
	
  std::map <SvxGtrack*, SvtxTrack*> _gtrack_track_match;      ///< holds a map between the SvxGtrack and the properly reconstructed SvtxTrack if any
  std::map <SvtxTrack*, SvxGtrack*> _track_gtrack_match;      ///< holds a map between a properly reconstructed SvtxTrack and the source particle

  std::map <SvtxTrack*, unsigned int> _track_purity_map;


  std::multimap <int, PHG4HitDefs::keytype> _particleid_g4hitid_mmap; ///< forward look up between a truth particle id and a g4hit id
  EvalLinks* _cluster_g4hit_svtx_links;
  EvalLinks* _cluster_g4hit_silicon_tracker_links;
  EvalLinks* _track_particle_links;
  
  //-----------------------
  // evaluator subroutines

  // fill internal structure subroutines
  int  fillLayerTypeMap(PHCompositeNode* topNode); ///< determines which layer numbers correspond to what type of geometry info
  int  fillGhitList(PHCompositeNode* topNode); ///< fill a list to all the g4hits
  void fillGtrackObjects();         ///< creates the internal storage list of truth particles
  void fillCellToG4HitMap();        ///< creates the map between cells and g4hits pointers
  void fillClusterToG4HitMap();     ///< creates the map between clusters and g4hits pointers
  void fillTrackToGtrackMap();      ///< creates the map between tracks and gtracks pointers
  void fillTrackPurityMap();

  // fill DST output
  int fillClusterToG4HitLinks(PHCompositeNode *topNode);     // cluster ancestry eval
  int fillTrackToG4TruthInfoLinks(PHCompositeNode *topNode); // tracking ancestry eval
  int fillG4TruthInfoToTrackLinks(PHCompositeNode *topNode); // particle decendent eval
  
  // output subroutines
  void fillOutputNtuples();         ///< dump the evaluator information into ntuple for external analysis
  void printInputInfo();            ///< print out the input object information (debugging upstream components)
  void printOutputInfo();           ///< print out the ancestry information for detailed diagnosis

  void fitGtrackProducedClusters();
};

#endif // __PHG4EVALUATOR_H__
