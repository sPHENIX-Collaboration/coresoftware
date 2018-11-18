#ifndef __PHG4SVTXCLUSTERIZER_H__
#define __PHG4SVTXCLUSTERIZER_H__

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>
#include <map>
#include <limits.h>

class SvtxHitMap;
class SvtxClusterMap;
class PHG4Cell;

class PHG4SvtxClusterizer : public SubsysReco {

public:

  PHG4SvtxClusterizer(const std::string &name = "PHG4SvtxClusterizer",
		      unsigned int min_layer = 0, unsigned int max_layer = UINT_MAX);
  virtual ~PHG4SvtxClusterizer(){}
  
  //! module initialization
  int Init(PHCompositeNode *topNode){return 0;}
  
  //! run initialization
  int InitRun(PHCompositeNode *topNode);
  
  //! event processing
  int process_event(PHCompositeNode *topNode);
  
  //! end of process
  int End(PHCompositeNode *topNode){return 0;}
  
  //! set an energy requirement relative to the thickness MIP expectation
  void set_threshold(const float fraction_of_mip) {
    _fraction_of_mip = fraction_of_mip;
  }
  float get_threshold_by_layer(const int layer) const {
    if (_thresholds_by_layer.find(layer) == _thresholds_by_layer.end()) return 0.0;
    return _thresholds_by_layer.find(layer)->second;
  }
  
  //! option to turn off z-dimension clustering
  void set_z_clustering(const int layer, const bool make_z_clustering) {
    _make_z_clustering.insert(std::make_pair(layer,make_z_clustering));
  }
  bool get_z_clustering(const int layer) const {
    if (_make_z_clustering.find(layer) == _make_z_clustering.end()) return true;
    return _make_z_clustering.find(layer)->second;
  }

  //! option to turn on/off energy weighted clustering
  void set_energy_weighting(const int layer, const bool make_e_weights) {
    _make_e_weights.insert(std::make_pair(layer,make_e_weights));
  }
  bool get_energy_weighting(const int layer) const {
    if (_make_e_weights.find(layer) == _make_e_weights.end()) return false;
    return _make_e_weights.find(layer)->second;
  }  

private:

  static bool lessthan(const PHG4Cell*, 
		       const PHG4Cell*);
  static bool ladder_lessthan(const PHG4Cell*, 
			      const PHG4Cell*);
  static bool mvtx_ladder_lessthan(const PHG4Cell*, 
			      const PHG4Cell*);
  bool are_adjacent(const PHG4Cell*, 
		    const PHG4Cell*, 
		    const int &);
  bool ladder_are_adjacent(const PHG4Cell*, 
			   const PHG4Cell*);
  bool mvtx_ladder_are_adjacent(const PHG4Cell*,
				const PHG4Cell*);

  void CalculateCylinderThresholds(PHCompositeNode *topNode);
  void CalculateLadderThresholds(PHCompositeNode *topNode);
  void CalculateMVTXLadderThresholds(PHCompositeNode *topNode);
  
  void ClusterCylinderCells(PHCompositeNode *topNode);
  void ClusterLadderCells(PHCompositeNode *topNode);
  void ClusterMVTXLadderCells(PHCompositeNode *topNode);

  void PrintClusters(PHCompositeNode *topNode);
  
  // node tree storage pointers
  SvtxHitMap* _hits;
  SvtxClusterMap* _clusterlist;

  // settings
  float _fraction_of_mip;
  std::map<int,float> _thresholds_by_layer; // layer->threshold
  std::map<int,bool> _make_z_clustering;    // layer->z_clustering_option
  std::map<int,bool> _make_e_weights;       // layer->energy_weighting_option

  unsigned int _min_layer;
  unsigned int _max_layer;
  
  PHTimeServer::timer _timer;
};

#endif
