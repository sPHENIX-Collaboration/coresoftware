#ifndef G4INTT_TRUTHCLUSTERIZER__H
#define G4INTT_TRUTHCLUSTERIZER__H

#include <g4tracking/TruthClusterizerBase.h>
#include <vector>

class PHCompositeNode;
class TrkrHit;

class  PHG4InttTruthClusterizer : public TruthClusterizerBase {
  public:

  PHG4InttTruthClusterizer ( );

  //void init_run(PHCompositeNode*& _topNode);
  
  // regular functions (in common iwht PHG4MvtxTruthClusterizer and PHG4TpcTruthClusterizer)
  void init_run(PHCompositeNode*& _topnode, int _verbosity);
  int  clusterize_hits(TrkrClusterContainer*);
  void check_g4hit(PHG4Hit*);
  void end_of_event();
  ~PHG4InttTruthClusterizer(){};

  // all data below is implemented for the clusterize_hits() function
  // It does the work that the following modules do for SvtxTrack's:
  // (1) g4intt/PHG4InttDigitizer -- use prefix _D__
  // (2) intt/InttClusterizer     -- use prefix _C__
  
  // ---------------------------------------
  // (1) g4intt/PHG4InttDigitizer _D__
  //     note that no noise is added during the digitization
  // ---------------------------------------
  public:
  void set_adc_scale(const int &layer, const std::vector<double> &userrange);
  void _D__DigitizeLadderCells(PHCompositeNode *topNode);
  private:
  int  _D__InitRun(PHCompositeNode* topNode);
  void _D__CalculateLadderCellADCScale(PHCompositeNode* topNode);

  // settings
  std::map<int, unsigned int> _max_adc;
  std::map<int, std::vector<std::pair<double, double> > > _max_fphx_adc;
  std::map<int, float> _energy_scale;
  std::string detector = "INTT";
  const unsigned int nadcbins = 8;
  // unsigned int m_nCells = 0;     // not really used by PHG4InttDigitizer
  // unsigned int m_nDeadCells = 0;     // not really used by PHG4InttDigitizer
   
  // ---------------------------------------
  // (2) intt/InttClusterizer _C__
  // ---------------------------------------
  public:
  void _C__ClusterLadderCells(PHCompositeNode* topNode, TrkrClusterContainer*);
  void PrintClusters(PHCompositeNode* topNode);
  void set_cluster_version(int value) { m_cluster_version = value; }

  private:
  bool _C__ladder_are_adjacent(const std::pair<TrkrDefs::hitkey, TrkrHit*> &lhs, const std::pair<TrkrDefs::hitkey, TrkrHit*> &rhs, const int layer);
  int m_cluster_version = 4;

  void set_z_clustering(const int layer, const bool make_z_clustering)
  { _make_z_clustering.insert(std::make_pair(layer, make_z_clustering)); }
  bool get_z_clustering(const int layer) const
  {
    if (_make_z_clustering.find(layer) == _make_z_clustering.end()) return true;
    return _make_z_clustering.find(layer)->second;
  }

  std::map<int, float> _thresholds_by_layer;  // layer->threshold
  std::map<int, bool> _make_z_clustering;     // layer->z_clustering_option
  std::map<int, bool> _make_e_weights;        // layer->energy_weighting_option

};

#endif

