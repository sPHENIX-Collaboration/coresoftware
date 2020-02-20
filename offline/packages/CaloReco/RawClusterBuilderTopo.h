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

int RawClusterBuilderTopo_constants_EMCal_eta_start_given_IHCal[24] = 
  {2,  6,  10, 14, 18, 22, 26, 30, 33, 37, 41, 44, 
   48, 52, 55, 59, 63, 66, 70, 74, 78, 82, 86, 90 };

int RawClusterBuilderTopo_constants_EMCal_eta_end_given_IHCal[24] = 
  {5,  9,  13, 17, 21, 25, 29, 32, 36, 40, 43, 47,
   51, 54, 58, 62, 65, 69, 73, 77, 81, 85, 89, 93 };

int RawClusterBuilderTopo_constants_IHCal_eta_given_EMCal[96] = {
  -1, -1,  0,  0,  0,  0,  1,  1,  1,  1,  2,  2, 
   2,  2,  3,  3,  3,  3,  4,  4,  4,  4,  5,  5, 
   5,  5,  6,  6,  6,  6,  7,  7,  7,  8,  8,  8, 
   8,  9,  9,  9,  9, 10, 10, 10, 11, 11, 11, 11, 
  12, 12, 12, 12, 13, 13, 13, 14, 14, 14, 14, 15, 
  15, 15, 15, 16, 16, 16, 17, 17, 17, 17, 18, 18, 
  18, 18, 19, 19, 19, 19, 20, 20, 20, 20, 21, 21, 
  21, 21, 22, 22, 22, 22, 23, 23, 23, 23, -1, -1 };

class RawClusterBuilderTopo : public SubsysReco
{
 public:
  RawClusterBuilderTopo(const std::string &name = "RawClusterBuilderTopo");
  virtual ~RawClusterBuilderTopo() {}

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void set_noise( float noise_0, float noise_1, float noise_2 ) {

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

  int get_first_matching_EMCal_phi_from_IHCal( int index_hcal_phi ) {
    return ( (68 + 4 * ( index_hcal_phi - 32) + 256 ) % 256 );
  }

  int get_matching_HCal_phi_from_EMCal( int index_emcal_phi ) { 
    return ( (32 + ( index_emcal_phi - 68 ) / 4 + 64 ) % 64 );
  }

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

  std::vector<int> get_adjacent_towers_by_ID( int ID ) {

    int this_layer = get_ilayer_from_ID( ID );
    int this_eta = get_ieta_from_ID( ID );
    int this_phi = get_iphi_from_ID( ID );
    
    std::vector<int> adjacent_towers;

    // for both IHCal and OHCal, add adjacent layers in the HCal
    if ( this_layer == 0 || this_layer == 1 ) {
      
      for (int delta_layer = 0; delta_layer <= 1; delta_layer++) {
	for (int delta_eta = -1; delta_eta <= 1; delta_eta++) {
	  for (int delta_phi = -1; delta_phi <= 1; delta_phi++) {
	    
	    if ( delta_layer == 0 && delta_eta == 0 && delta_phi == 0 ) continue; // this is the same tower
	    
	    int test_eta = this_eta + delta_eta;
	    if ( test_eta < 0 || test_eta >= _HCAL_NETA ) { continue; } // ignore if at the (eta) edge of calorimeter
	    
	    int test_layer = ( this_layer + delta_layer ) % 2; // wrap around in layer
	    int test_phi = ( this_phi + delta_phi + _HCAL_NPHI ) % _HCAL_NPHI; // wrap around in phi (add 64 to avoid -1)
	    
	    // disallow "corner" adjacency (diagonal in eta/phi plane and in different layer) if this option not enabled
	    if ( !_allow_corner_neighbor && delta_layer == 1 && abs( delta_phi ) == 1 && abs( delta_eta ) == 1 ) {
	      //if (Verbosity() > 10) std::cout << " --> --> --> corner growth not allowed " << std::endl;
	      continue;
	    }
	    
	    // add to list of adjacent towers
	    adjacent_towers.push_back( get_ID( test_layer, test_eta, test_phi ) );
	    
	  }
	}
      }
    }

    // for IHCal only, also add 4x4 group of EMCal towers
    if ( this_layer == 0 && _enable_EMCal ) {

      int EMCal_phi_start = get_first_matching_EMCal_phi_from_IHCal( this_phi );
      int EMCal_eta_start = RawClusterBuilderTopo_constants_EMCal_eta_start_given_IHCal[ this_eta ];
      int EMCal_eta_end = RawClusterBuilderTopo_constants_EMCal_eta_end_given_IHCal[ this_eta ];

      for (int new_eta = EMCal_eta_start; new_eta <= EMCal_eta_end; new_eta++) {
	for (int delta_phi = 0; delta_phi < 4; delta_phi++) {
	  int new_phi = ( EMCal_phi_start + delta_phi + 256 ) % 256;

	  int EMCal_tower = get_ID( 2, new_eta, new_phi );
	  //std::cout << " DVP : HCal tower with eta / phi = " << this_eta << " / " << this_phi << ", adding EMCal tower with eta / phi = " << new_eta << " / " << new_phi << std::endl;
	  adjacent_towers.push_back( EMCal_tower );	
	}
      }
    }

    // for EMCal, add adjacent EMCal towers and (usually) one IHCal tower
    if ( this_layer == 2 ) {

      // EMCal towers first
      for (int delta_eta = -1; delta_eta <= 1; delta_eta++) {
	for (int delta_phi = -1; delta_phi <= 1; delta_phi++) {
	  
	  if ( delta_eta == 0 && delta_phi == 0 ) continue; // this is the same tower
	  
	  int test_eta = this_eta + delta_eta;
	  if ( test_eta < 0 || test_eta >= _EMCAL_NETA ) { continue; } // ignore if at the (eta) edge of calorimeter
	  
	  int test_phi = ( this_phi + delta_phi + _EMCAL_NPHI ) % _EMCAL_NPHI; // wrap around in phi (add 256 to avoid -1)
	  	    
	  // add to list of adjacent towers
	  adjacent_towers.push_back( get_ID( this_layer, test_eta, test_phi ) );
	    
	}
      }

      // now add IHCal towers
      if ( _enable_HCal ) {
	int HCal_eta = RawClusterBuilderTopo_constants_IHCal_eta_given_EMCal[ this_eta ];
	int HCal_phi = get_matching_HCal_phi_from_EMCal( this_phi );
	
	if ( HCal_eta >= 0 ) {
	  int IHCal_tower = get_ID( 0, HCal_eta, HCal_phi );
	  //std::cout << " DVP : EMCal tower with eta / phi = " << this_eta << " / " << this_phi << ", adding IHCal tower with eta / phi = " << HCal_eta << " / " << HCal_phi << std::endl;
	  adjacent_towers.push_back( IHCal_tower );
	} else {
	  //std::cout << " DVP : EMCal tower with eta / phi = " << this_eta << " / " << this_phi << ", does not have matching IHCal due to large eta " << std::endl;
	}
      }

    }

    return adjacent_towers;

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
