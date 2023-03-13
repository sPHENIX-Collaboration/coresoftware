#ifndef G4EVAL_SVTXTRUTHEVAL_H
#define G4EVAL_SVTXTRUTHEVAL_H

#include "BaseTruthEval.h"

#include <trackbase/TrkrDefs.h>

class PHCompositeNode;

class PHG4Hit;
class PHG4HitContainer;
class PHG4Particle;
class PHG4TruthInfoContainer;
class PHG4CylinderGeomContainer;
class PHG4TpcCylinderGeomContainer;
class PHG4VtxPoint;
class TrkrCluster;
class ActsGeometry;

#include <map>
#include <set>
#include <vector>
#include <memory>

class SvtxTruthEval
{
 public:
  SvtxTruthEval(PHCompositeNode* topNode);
  virtual ~SvtxTruthEval();

  void next_event(PHCompositeNode* topNode);
  void do_caching(bool do_cache) { _do_cache = do_cache; }
  void set_strict(bool strict)
  {
    _strict = strict;
    _basetrutheval.set_strict(strict);
  }
  void set_verbosity(int verbosity)
  {
    _verbosity = verbosity;
    _basetrutheval.set_verbosity(verbosity);
  }

  std::set<PHG4Hit*> all_truth_hits();
  std::set<PHG4Hit*> all_truth_hits(PHG4Particle* particle);
  PHG4Particle* get_particle(PHG4Hit* g4hit);
  int get_embed(PHG4Particle* particle);
  PHG4VtxPoint* get_vertex(PHG4Particle* particle);
  bool is_primary(PHG4Particle* particle);
  PHG4Particle* get_primary_particle(PHG4Hit* g4hit);
  PHG4Particle* get_primary_particle(PHG4Particle* particle);
  PHG4Particle* get_particle(const int trackid); 

  std::map<TrkrDefs::cluskey, std::shared_ptr<TrkrCluster> > all_truth_clusters(PHG4Particle* particle);

  bool is_g4hit_from_particle(PHG4Hit* g4hit, PHG4Particle* particle);
  bool are_same_particle(PHG4Particle* p1, PHG4Particle* p2);
  bool are_same_vertex(PHG4VtxPoint* vtx1, PHG4VtxPoint* vtx2);

  PHG4Hit* get_innermost_truth_hit(PHG4Particle* particle);
  PHG4Hit* get_outermost_truth_hit(PHG4Particle* particle);

  unsigned int get_errors() { return _errors + _basetrutheval.get_errors(); }

  std::set<PHG4Hit*> get_truth_hits_from_truth_cluster(TrkrDefs::cluskey ckey);
  void  FillTruthHitsFromParticleCache();
 private:
  void get_node_pointers(PHCompositeNode* topNode);
  bool has_node_pointers();

  void LayerClusterG4Hits(std::set<PHG4Hit*> truth_hits, std::vector<PHG4Hit*> &contributing_hits, std::vector<double> &contributing_hits_energy, std::vector<std::vector<double>> &contributing_hits_entry, std::vector<std::vector<double>> &contributing_hits_exit, float layer, float &x, float &y, float &z,  float &t, float &e);

  float line_circle_intersection(float x[], float y[], float z[], float radius);

  void G4ClusterSize(TrkrDefs::cluskey ckey, unsigned int layer, std::vector<std::vector<double>> contributing_hits_entry,std::vector<std::vector<double>> contributing_hits_exit, float &g4phisize, float &g4zsize);

  unsigned int getAdcValue(double gedep);

  BaseTruthEval _basetrutheval;

  PHG4TruthInfoContainer* _truthinfo{nullptr};
  PHG4HitContainer* _g4hits_svtx{nullptr};
  PHG4HitContainer* _g4hits_mms{nullptr};
  PHG4HitContainer* _g4hits_tracker{nullptr};
  PHG4HitContainer* _g4hits_maps{nullptr};

  PHG4TpcCylinderGeomContainer* _tpc_geom_container;
  PHG4CylinderGeomContainer *_intt_geom_container;
  PHG4CylinderGeomContainer* _mvtx_geom_container;
  PHG4CylinderGeomContainer* _mms_geom_container;
  ActsGeometry* _tgeometry{nullptr};

  bool _strict;
  int _verbosity;
  unsigned int _errors;
  unsigned long iclus;
  
  const unsigned int _nlayers_maps = 3;
  const unsigned int _nlayers_intt = 4;
  const unsigned int _nlayers_tpc = 48;
  const unsigned int _nlayers_mms = 2;

  std::multimap<TrkrDefs::cluskey, PHG4Hit*> _truth_cluster_truth_hit_map;

  bool _do_cache;
  std::set<PHG4Hit*> _cache_all_truth_hits;
  std::map<PHG4Particle*, std::set<PHG4Hit*> > _cache_all_truth_hits_g4particle;
  std::map<PHG4Particle*, std::map<TrkrDefs::cluskey, std::shared_ptr<TrkrCluster> > > _cache_all_truth_clusters_g4particle;
  std::map<PHG4Particle*, PHG4Hit*> _cache_get_innermost_truth_hit;
  std::map<PHG4Particle*, PHG4Hit*> _cache_get_outermost_truth_hit;
  std::map<PHG4Hit*, PHG4Particle*> _cache_get_primary_particle_g4hit;

};

#endif  // G4EVAL_SVTXTRUTHEVAL_H
