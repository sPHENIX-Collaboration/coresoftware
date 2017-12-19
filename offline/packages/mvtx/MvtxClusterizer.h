#ifndef __PHG4SVTXCLUSTERIZER_H__
#define __PHG4SVTXCLUSTERIZER_H__

#include <tracker/TrackerDefs.h>
#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>
#include <map>
#include <limits.h>

class TrackerClusterContainer;
class TrackerHit;
class TrackerHitContainer;

class MvtxClusterizer : public SubsysReco {

public:

  MvtxClusterizer(const std::string &name = "MvtxClusterizer",
		      unsigned int min_layer = 0, unsigned int max_layer = UINT_MAX);
  virtual ~MvtxClusterizer(){}
  
  //! module initialization
  int Init(PHCompositeNode *topNode){return 0;}
  
  //! run initialization
  int InitRun(PHCompositeNode *topNode);
  
  //! event processing
  int process_event(PHCompositeNode *topNode);
  
  //! end of process
  int End(PHCompositeNode *topNode){return 0;}
  
  //! option to turn off z-dimension clustering
  void set_z_clustering(const int layer, const bool make_z_clustering) {
    _make_z_clustering.insert(std::make_pair(layer,make_z_clustering));
  }
  bool get_z_clustering(const int layer) const {
    if (_make_z_clustering.find(layer) == _make_z_clustering.end()) return true;
    return _make_z_clustering.find(layer)->second;
  }

private:

  static bool pixel_lessthan(const TrackerDefs::keytype, 
			      const TrackerDefs::keytype);
  bool pixel_are_adjacent(const TrackerDefs::keytype,
				const TrackerDefs::keytype);

  void ClusterMapsLadderHits(PHCompositeNode *topNode);

  void PrintClusters(PHCompositeNode *topNode);
  
  // node tree storage pointers
  TrackerHitContainer* _hits;
  TrackerClusterContainer* _clusters;

  // settings
  std::map<int,bool> _make_z_clustering;    // layer->z_clustering_option

  unsigned int _min_layer;
  unsigned int _max_layer;
  
  PHTimeServer::timer _timer;
};

#endif
