#ifndef OTRACK_OUTERTRACKERCLUSTERIZER_H
#define OTRACK_OUTERTRACKERCLUSTERIZER_H

#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrDefs.h>

#include <phool/PHTimeServer.h>
#include <limits>
#include <map>
#include <string>
#include <utility>

class PHCompositeNode;
class TrkrHitSetContainer;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHit;

class OuterTrackerClusterizer : public SubsysReco
{
 public:
  OuterTrackerClusterizer(const std::string &name = "OuterTrackerClusterizer",
                  unsigned int min_layer = 0, unsigned int max_layer = UINT_MAX);
  virtual ~OuterTrackerClusterizer() {}

  //! module initialization
  int Init(PHCompositeNode *topNode) { return 0; }

  //! run initialization
  int InitRun(PHCompositeNode *topNode);

  //! event processing
  int process_event(PHCompositeNode *topNode);

  //! end of process
  int End(PHCompositeNode *topNode) { return 0; }

  /*
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
  */

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

  /*
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
  */

 private:
  bool ladder_are_adjacent(const std::pair<TrkrDefs::hitkey, TrkrHit*> &lhs, const std::pair<TrkrDefs::hitkey, TrkrHit*> &rhs, const int layer);

  //void CalculateLadderThresholds(PHCompositeNode *topNode);
  void ClusterLadderCells(PHCompositeNode *topNode);

  void PrintClusters(PHCompositeNode *topNode);

  // node tree storage pointers
  TrkrHitSetContainer *m_hits;
  TrkrClusterContainer *m_clusterlist; 
  TrkrClusterHitAssoc *m_clusterhitassoc;


  // settings
  float _fraction_of_mip;
  std::map<int, float> _thresholds_by_layer;  // layer->threshold
  std::map<int, bool> _make_z_clustering;     // layer->z_clustering_option
  std::map<int, bool> _make_e_weights;        // layer->energy_weighting_option

  PHTimeServer::timer _timer;
};

#endif
