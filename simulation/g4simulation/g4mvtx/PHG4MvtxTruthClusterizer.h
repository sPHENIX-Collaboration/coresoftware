#ifndef G4MVTX_TRUTHCLUSTERIZER__H
#define G4MVTX_TRUTHCLUSTERIZER__H

#include <g4tracking/TruthClusterizerBase.h>

class PHCompositeNode;
class TrkrHit;

class  PHG4MvtxTruthClusterizer : public TruthClusterizerBase {
  public:

  PHG4MvtxTruthClusterizer ();

  //void init_run(PHCompositeNode*& _topNode);
  
  // regular functions (in common iwht PHG4InttTruthClusterizer and PHG4TpcTruthClusterizer)
  void init_run ( PHCompositeNode*& _topNode, int _verbosity );
  int  clusterize_hits(TrkrClusterContainer*);
  void check_g4hit(PHG4Hit*);
  void end_of_event();
  ~PHG4MvtxTruthClusterizer(){};

  // all data below is implemented for the clusterize_hits() function
  // It does the work that the following modules do for SvtxTrack's:
  // (1) g4mvtx/PHG4MvtxDigitizer -- use prefix _D__
  // (2) mvtx/MvtxHitPruner       -- use prefix _P__
  // (3) mvtx/MvtxClusterizer     -- use prefix _C__
  
  // ---------------------------------------
  // (1) g4mvtx/PHG4MvtxDigitizer _D__
  // ---------------------------------------
  private:
  std::map<int, unsigned short> _max_adc;
  std::map<int, float> _energy_scale;
  float _energy_threshold;
  public:
  // Need to update the G4INTT macro accordingly with the set_adc_scale values
  void set_adc_scale(const int layer, const unsigned short max_adc, const float energy_per_adc)
  {
    _max_adc.insert(std::make_pair(layer, max_adc));
    _energy_scale.insert(std::make_pair(layer, energy_per_adc));
  }
  void set_energy_threshold(const float threshold) { _energy_threshold = threshold; }
  float _D__get_energy_threshold() { return _energy_threshold; }

  int  _D__InitRun                         (PHCompositeNode *topNode);
  void _D__CalculateMvtxLadderCellADCScale (PHCompositeNode *topNode); // called by _D_InitRun
  void _D__DigitizeMvtxLadderCells         ();
 
  // ---------------------------------------
  // (2) mvtx/MvtxHitPruner _P__
  // ---------------------------------------
  public:
  int _P__MvtxHitPruner();

  // ---------------------------------------
  // (3) mvtx/MvtxClusterizer _C__
  // ---------------------------------------
  private:
  bool m_makeZClustering { true };  // z_clustering_option
  int  m_cluster_version { 4 };
  bool _C__are_adjacent(const std::pair<TrkrDefs::hitkey, TrkrHit*> &lhs, const std::pair<TrkrDefs::hitkey, TrkrHit*> &rhs);

  public:
  void set_cluster_version(int /*value*/) {
    std::cout << "None-implemented function in PHG4MvtxTruthClusterizer. " << std::endl
              << "Only TrkrClusterv4 currently implemented." << std::endl;
  };
  void _C__ClusterMvtx(TrkrClusterContainer* clusters);
  bool GetZClustering() const { return m_makeZClustering; };
  void SetZClustering(const bool make_z_clustering) { m_makeZClustering = make_z_clustering; };

};

#endif
