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

  void set_noise( float noise_0 = 0.0025, float noise_1 = 0.006, float noise_2 = 0.03) {

    _noise_LAYER[ 0 ] = noise_0;
    _noise_LAYER[ 1 ] = noise_1;
    _noise_LAYER[ 2 ] = noise_2;

  }

  void set_significance( float seed, float grow, float peri ) {

    _sigma_seed = seed;
    _sigma_grow = grow;
    _sigma_peri = peri;

  }

  void allow_corner_neighbor( bool allow ) { 

    _allow_corner_neighbor = allow;

  }

  void set_enable_HCal( bool enable_HCal ) {

    _enable_HCal = enable_HCal;

  }

  void set_enable_EMCal( bool enable_EMCal ) {

    _enable_EMCal = enable_EMCal;

  }

 private:

  void CreateNodes(PHCompositeNode *topNode);

  std::vector< std::vector< std::vector<float> > > _TOWERMAP_E_LAYER_ETA_PHI;
  std::vector< std::vector< std::vector<int> > > _TOWERMAP_STATUS_LAYER_ETA_PHI;

  std::vector< std::vector<float> > _EMTOWERMAP_E_ETA_PHI;
  std::vector< std::vector<int> > _EMTOWERMAP_STATUS_ETA_PHI;

  // geometric constants to express IHCal<->EMCal overlap in eta
  static int RawClusterBuilderTopo_constants_EMCal_eta_start_given_IHCal[];

  static int RawClusterBuilderTopo_constants_EMCal_eta_end_given_IHCal[];

  static int RawClusterBuilderTopo_constants_IHCal_eta_given_EMCal[];

  // utility functions to express IHCal<->EMCal overlap in phi
  int get_first_matching_EMCal_phi_from_IHCal( int index_hcal_phi ) {
    return ( (68 + 4 * ( index_hcal_phi - 32) + 256 ) % 256 );
  }

  int get_matching_HCal_phi_from_EMCal( int index_emcal_phi ) { 
    return ( (32 + ( index_emcal_phi - 68 ) / 4 + 64 ) % 64 );
  }

  std::vector<int> get_adjacent_towers_by_ID( int ID );

  int get_ID( int ilayer, int ieta, int iphi ) {
    if ( ilayer < 2 ) return ilayer * 24 * 64 + ieta * 64 + iphi;
    else return 256 * 96 + ieta * 256 + iphi;
  }
  
  int get_ilayer_from_ID( int ID ) {
    if ( ID < 256 * 96 ) return ( (int) ( ID / ( 24 * 64 ) ) );
    else return 2;
  }
  
  int get_ieta_from_ID( int ID ) {
    if ( ID < 256 * 96 ) return ( (int) ( ( ID % ( 24 * 64 ) ) / ( 64 ) ) );
    else return ( (int) ( ( ID - 256 * 96 ) / 256 ) );
  }
  
  int get_iphi_from_ID( int ID ) {
    if ( ID < 256 * 96 ) return ( (int) ( ID % 64 ) );
    else return ( (int) ( ( ID - 256 * 96 ) % 256 ) );
  }

  int get_status_from_ID( int ID ) {
    if ( ID < 256 * 96 ) return _TOWERMAP_STATUS_LAYER_ETA_PHI[ get_ilayer_from_ID( ID ) ][ get_ieta_from_ID( ID ) ][ get_iphi_from_ID( ID ) ];
    else return _EMTOWERMAP_STATUS_ETA_PHI[ get_ieta_from_ID( ID ) ][ get_iphi_from_ID( ID ) ];
  }

  float get_E_from_ID( int ID ) {
    if ( ID < 256 * 96 ) return _TOWERMAP_E_LAYER_ETA_PHI[ get_ilayer_from_ID( ID ) ][ get_ieta_from_ID( ID ) ][ get_iphi_from_ID( ID ) ];
    else return _EMTOWERMAP_E_ETA_PHI[ get_ieta_from_ID( ID ) ][ get_iphi_from_ID( ID ) ];
  }

  void set_status_by_ID( int ID , int status ) {
    if ( ID < 256 * 96 ) _TOWERMAP_STATUS_LAYER_ETA_PHI[ get_ilayer_from_ID( ID ) ][ get_ieta_from_ID( ID ) ][ get_iphi_from_ID( ID ) ] = status;
    else _EMTOWERMAP_STATUS_ETA_PHI[ get_ieta_from_ID( ID ) ][ get_iphi_from_ID( ID ) ] = status;
  }
  
  RawClusterContainer *_clusters;

  int _EMCAL_NETA;
  int _EMCAL_NPHI;

  int _HCAL_NETA;
  int _HCAL_NPHI;

  float _noise_LAYER[3];

  float _sigma_seed;
  float _sigma_grow;
  float _sigma_peri;

  bool _allow_corner_neighbor;

  bool _enable_HCal;
  bool _enable_EMCal;

  std::string ClusterNodeName;
};

#endif
