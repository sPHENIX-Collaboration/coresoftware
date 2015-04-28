#ifndef __PHG4CALEVALUATOR_H__
#define __PHG4CALEVALUATOR_H__

//===============================================
/// \file PHG4CalEvaluator.h
/// \brief Compares reconstructed tracks to truth particles
/// \author Michael P. McCumber (original SVX version)
/// \author Richard Hollis (translated to sPHENIX CAL simulations)
//===============================================

// PHENIX includes
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>

// ROOT includes
#include <TNtuple.h>

// standard includes
#include <vector>
#include <map>

// forward declarations
class PHCompositeNode;
class PHG4TruthInfoContainer;
class SvtxTrackMap;
class SvtxVertexMap;
class PHG4CylinderCellContainer;
class PHG4Hit;
class PHG4HitContainer;
class RawTower;
class RawTowerContainer;
class RawTowerGeom;
class RawCluster;
class RawClusterContainer;
class TFile;

/// \class CalGshower
///
/// \brief A quick container class for truth particles
///        provides a mapping between an original particle (e.g. gamma)
///        and the end hits in the detector.
class CalGshower
{
 public:

  CalGshower() 
    {
      _embed = 0;
      _particle_id = 0;
      _flavor = 0;
      _px = 0.0;
      _py = 0.0;
      _pz = 0.0;
      _e  = 0.0;
      _vx = 0.0;
      _vy = 0.0;
      _vz = 0.0;
      _radius = 0.0;
      _g4hits.clear();
      _eabs = 0.0;
    }
  virtual ~CalGshower() {}

  int get_embed() const {return _embed;}
  void set_embed(const int embed) {_embed = embed;}

  int get_particle_id() const {return _particle_id;}
  void set_particle_id(const int particle_id) {_particle_id = particle_id;}

  int get_flavor() const {return _flavor;}
  void set_flavor(const int flavor) {_flavor = flavor;}

  float get_px() const {return _px;}
  void set_px(const float px) {_px = px;}

  float get_py() const {return _py;}
  void set_py(const float py) {_py = py;}

  float get_pz() const {return _pz;}
  void set_pz(const float pz) {_pz = pz;}

  float get_e() const {return _e;}
  void set_e(const float e) {_e = e;}

  float get_vx() const {return _vx;}
  void set_vx(const float vx) {_vx = vx;}

  float get_vy() const {return _vy;}
  void set_vy(const float vy) {_vy = vy;}

  float get_vz() const {return _vz;}
  void set_vz(const float vz) {_vz = vz;}

  float get_moliere_radius() const {return _radius;}
  void set_moliere_radius(const float r) {_radius = r;}

  size_t get_ng4hits() const {return _g4hits.size();}
  void clear_g4hits() {_g4hits.clear();}
  PHG4Hit* get_g4hit(const int ig4hit) {return _g4hits.at(ig4hit);}
  void add_g4hit(PHG4Hit* g4hit) {_g4hits.push_back(g4hit);}

  double get_edep() const;
  
  double get_eabs() const {return _eabs;}
  void set_eabs(float e) {_eabs = e;}
  
 protected:

  int _embed; ///< was the particle in the embedded map storage

  int _particle_id; ///< truth particle id which created the hit
  int _flavor; ///< truth particle geant flavor code which created the hit

  float _px; ///< truth particle x-momentum at creation
  float _py; ///< truth particle y-momentum at creation
  float _pz; ///< truth particle z-momentum at creation

  float _e; ///< truth particle energy at creation

  float _vx; ///< truth particle x-position at creation
  float _vy; ///< truth particle y-position at creation
  float _vz; ///< truth particle z-position at creation

  float _radius; ///< moliere radius for this shower

  std::vector<PHG4Hit*> _g4hits; ///< storage of associated PHG4Hits

  float _eabs; ///< energy in the absorber if available
};

/// \class PHG4CalEvaluator
///
/// \brief Compares reconstructed showers to truth particles
///
/// Plan: This module will trace the reconstructed tracks back to
/// the greatest contributor Monte Carlo particle and then
/// test one against the other.
///
class PHG4CalEvaluator : public SubsysReco
{

 public:
 
  PHG4CalEvaluator(const std::string &name = "PHG4CalEVALUATOR",
		    const std::string &filename = "g4eval_cal.root");
  virtual ~PHG4CalEvaluator() {};
		
  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void Verbosity(int verb){verbosity = verb;}

  void set_trace_truth(bool b) {_trace_truth = false;}

  void Detector(const std::string &name) {detector = name;}
	
 protected:
  
  unsigned long _ievent;

  bool _trace_truth;

  //----------------------------------
  // evaluator output ntuples

  TNtuple *_ntp_event;
  TNtuple *_ntp_gshower;
  TNtuple *_ntp_tower;
  TNtuple *_ntp_cluster;

  // evaluator output file
  std::string _filename;
  TFile *_tfile;

  //----------------------------
  // internal evaluator storage

  // pointers to the object lists makes them availabe for all
  // the process_event subroutines
  SvtxTrackMap *_trackList;
  SvtxVertexMap *_vertexList;
  PHG4TruthInfoContainer *_truth_info_container;
  PHG4HitContainer *_g4hitList;
  PHG4HitContainer *_g4hitAbsorberList;
  PHG4CylinderCellContainer *_cyl_cel_container;
  RawTowerContainer *_towerList;
  RawClusterContainer *_clusterList;
  RawTowerGeom *towergeom;

  std::string detector;
  std::vector<CalGshower> _gshower_list;              ///< holds a list of truth particles (convenience storage)
  std::map <PHG4Hit*,CalGshower*> _g4hit_gshower_map; ///< reverse g4hit->gshower lookup
  std::map <int,CalGshower*> _g4hitid_gshower_map;    ///< reverse g4hit->gshower lookup (uses id)
  
  std::map <RawTower*, CalGshower*> _tower_gshower_map;     ///< holds the greatest truth contributor map to a particular tower
  std::map <RawTower*, float> _tower_epurity_map;            ///< holds the greatest truth contributor map to a particular tower
  std::map <RawCluster*, CalGshower*> _cluster_gshower_map; ///< holds the greatest truth contributor map to a particular cluster
  std::map <CalGshower*, RawCluster*> _gshower_cluster_map; ///< holds the reverse lookup of the leading cluster produced by a gshower
  std::map <RawCluster*, float> _cluster_epurity_map;        ///< holds the greatest truth contributor map to a particular cluster

  std::multimap <RawTower*, PHG4Hit*> _tower_ghit_mmap; ///< holds anestry to skip nested loops second time around

  //-----------------------
  // evaluator subroutines

  // fill internal structure subroutines
  void fillGshowerObjects();      ///< creates the internal storage list of truth particles
  void fillTowerToGshowerMap();   ///< fills the map of the tower's greatest truth contributor
  void fillClusterToGshowerMap(); ///< fills the map of the cluster's greatest truth contributor

  // fill output subroutines
  void fillOutputNtuples();         ///< dump the evaluator information into ntuple for external analysis

  // diagnositics dump
  void printInputInfo();            ///< print out the input object information (debugging upstream components)
  void printOutputInfo();           ///< print out the ancestry information for detailed diagnosis
};

#endif // __PHG4CALEVALUATOR_H__
