#ifndef CALORECO_RAWCLUSTERTOPO_H
#define CALORECO_RAWCLUSTERTOPO_H

//===========================================================
/// \file RawClusterBuilderTopo.h
/// \brief 3-D topoClustering across calorimeter layers
/// \author Dennis V. Perepelitsa
//===========================================================

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class PHCompositeNode;
class RawClusterContainer;

class RawClusterBuilderTopo : public SubsysReco
{
 public:
  RawClusterBuilderTopo(const std::string &name = "RawClusterBuilderTopo");
  virtual ~RawClusterBuilderTopo() {}

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void set_noise( float noise_0, float noise_1 ) {

    _noise_LAYER[ 0 ] = noise_0;
    _noise_LAYER[ 1 ] = noise_1;

  }

  void set_significance( float seed, float grow, float peri ) {

    _sigma_seed = seed;
    _sigma_grow = grow;
    _sigma_peri = peri;

  }

  void allow_corner_neighbor( bool allow ) { 

    _allow_corner_neighbor = allow;

  }

 private:
  void CreateNodes(PHCompositeNode *topNode);


  int get_ID( int ilayer, int ieta, int iphi ) {
    return ilayer * 24 * 64 + ieta * 64 + iphi;
  }
  
  int get_ilayer_from_ID( int ID ) {
    return ( (int) ( ID / ( 24 * 64 ) ) );
  }
  
  int get_ieta_from_ID( int ID ) {
    return ( (int) ( ( ID % ( 24 * 64 ) ) / ( 64 ) ) );
  }
  
  int get_iphi_from_ID( int ID ) {
    return ( (int) ( ID % 64 ) );
  }
    
  RawClusterContainer *_clusters;

  std::vector< std::vector< std::vector<float> > > _TOWERMAP_E_LAYER_ETA_PHI;
  std::vector< std::vector< std::vector<int> > > _TOWERMAP_STATUS_LAYER_ETA_PHI;

  int _HCAL_NETA;
  int _HCAL_NPHI;

  float _noise_LAYER[2];

  float _sigma_seed;
  float _sigma_grow;
  float _sigma_peri;

  bool _allow_corner_neighbor;

  std::string ClusterNodeName;
};

#endif
