#ifndef INTT_INTTCLUSTERIZER_H
#define INTT_INTTCLUSTERIZER_H

#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrCluster.h>

#include <climits>
#include <map>
#include <string>
#include <utility>

class PHCompositeNode;
class TrkrHitSetContainer;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrClusterCrossingAssoc;
class TrkrHit;
class RawHit;
class RawHitSet;
class RawHitSetContainer;

class InttClusterizer : public SubsysReco
{
 public:
  InttClusterizer(const std::string &name = "InttClusterizer",
                  unsigned int min_layer = 0, unsigned int max_layer = UINT_MAX);
  ~InttClusterizer() override {}

  //! module initialization
  int Init(PHCompositeNode */*topNode*/) override { return 0; }

  //! run initialization
  int InitRun(PHCompositeNode *topNode) override;

  //! event processing
  int process_event(PHCompositeNode *topNode) override;

  //! end of process
  int End(PHCompositeNode */*topNode*/) override { return 0; }

  //! set an energy requirement relative to the thickness MIP expectation
  void set_threshold(const float fraction_of_mip)
  {
    _fraction_of_mip = fraction_of_mip;
  }
  float get_threshold_by_layer(const int layer) const
  {
    if (_thresholds_by_layer.find(layer) == _thresholds_by_layer.end()) return 0.0;
    return _thresholds_by_layer.find(layer)->second;
  }

  //! option to turn off z-dimension clustering
  void set_z_clustering(const int layer, const bool make_z_clustering)
  {
    _make_z_clustering.insert(std::make_pair(layer, make_z_clustering));
  }
  bool get_z_clustering(const int layer) const
  {
    if (_make_z_clustering.find(layer) == _make_z_clustering.end()) return true;
    return _make_z_clustering.find(layer)->second;
  }

  //! option to turn on/off energy weighted clustering
  void set_energy_weighting(const int layer, const bool make_e_weights)
  {
    _make_e_weights.insert(std::make_pair(layer, make_e_weights));
  }
  bool get_energy_weighting(const int layer) const
  {
    if (_make_e_weights.find(layer) == _make_e_weights.end()) return false;
    return _make_e_weights.find(layer)->second;
  }
  void set_cluster_version(int value) { m_cluster_version = value; }
  void set_do_hit_association(bool do_assoc){do_hit_assoc = do_assoc;}
  void set_read_raw(bool read_raw){ do_read_raw = read_raw;}

 private:
  bool ladder_are_adjacent(const std::pair<TrkrDefs::hitkey, TrkrHit*> &lhs, const std::pair<TrkrDefs::hitkey, TrkrHit*> &rhs, const int layer);
  bool ladder_are_adjacent(RawHit* lhs,  RawHit* rhs, const int layer);

  void CalculateLadderThresholds(PHCompositeNode *topNode);
  void ClusterLadderCells(PHCompositeNode *topNode);
  void ClusterLadderCellsRaw(PHCompositeNode *topNode);
  void PrintClusters(PHCompositeNode *topNode);

  // node tree storage pointers
  TrkrHitSetContainer *m_hits;
  RawHitSetContainer *m_rawhits;
  TrkrClusterContainer *m_clusterlist; 
  TrkrClusterHitAssoc *m_clusterhitassoc;
  TrkrClusterCrossingAssoc *m_clustercrossingassoc{nullptr};

  // settings
  float _fraction_of_mip;
  std::map<int, float> _thresholds_by_layer;  // layer->threshold
  std::map<int, bool> _make_z_clustering;     // layer->z_clustering_option
  std::map<int, bool> _make_e_weights;        // layer->energy_weighting_option
  bool do_hit_assoc = true;
  bool do_read_raw = false;
  int m_cluster_version = 4;
};

#endif
