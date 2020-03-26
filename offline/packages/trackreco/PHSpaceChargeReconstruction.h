#ifndef TRACKRECO_PHSPACECHARGERECONSTRUCTION_H
#define TRACKRECO_PHSPACECHARGERECONSTRUCTION_H

#include <fun4all/SubsysReco.h>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>

// forward declaration
class SvtxTrack;
class SvtxTrackMap;
class TrkrCluster;
class TrkrClusterContainer;

class PHSpaceChargeReconstruction: public SubsysReco
{
  public:

  /// constructor
  PHSpaceChargeReconstruction( const std::string& = "PHSPACECHARGERECONSTRUCTION" );

  ///@name configuration
  //@{

  /// set tpc layers
  void set_tpc_layers( unsigned int first_layer, unsigned int n_layers );

  /// set grid dimensions
  /**
  \param zbins the number of bins along z
  \param rbins the number of bins in the radial direction. It must be a divider of the number of tpc layers (by default 48), e.g. 1, 12, 24, 48
  \param phibins the number of bins in the azimuth direction
  */
  void set_grid_dimensions( int zbins, int rbins, int phibins );

  /// output file
  /**
  for now only TGraphs of space charge corrections vs r are stored
  ultimately one wants to save the data in a format that is suitable for use in the reconstruction
  */
  void set_outputfile( const std::string& filename );

  //@}

  /// global initialization
  virtual int Init(PHCompositeNode*);

  /// run initialization
  virtual int InitRun(PHCompositeNode*);

  /// event processing
  virtual int process_event(PHCompositeNode*);

  /// end of processing
  virtual int End(PHCompositeNode*);

  private:

  /// load nodes
  int load_nodes( PHCompositeNode* );

  /// process tracks
  void process_tracks();

  /// process track
  void process_track( SvtxTrack* );

  /// calculate distortions
  void calculate_distortions();

  /// get cell from z, r and phi index
  int get_cell( int iz, int ir, int iphi ) const;

  /// get relevant cell for a given cluster
  int get_cell( TrkrCluster* ) const;

  /// output file
  std::string m_outputfile = "PHSpaceChargeReconstruction.root";

  // tpc layers
  unsigned int m_firstlayer_tpc = 7;
  unsigned int m_nlayers_tpc = 48;

  ///@name grid size
  //@{
  int m_zbins = 50;
  int m_phibins = 72;
  int m_rbins = 48;
  int m_totalbins = m_zbins*m_phibins*m_rbins;
  //@}

  // shortcut for relevant eigen matrices
  static constexpr int m_ncoord = 3;
  using matrix_t = Eigen::Matrix<float, m_ncoord, m_ncoord >;
  using column_t = Eigen::Matrix<float, m_ncoord, 1 >;

  /// left hand side matrices for distortion inversions
  std::vector<matrix_t> m_lhs;

  /// right hand side matrices for distortion inversions
  std::vector<column_t> m_rhs;

  /// keep track of how many clusters are used per cell
  std::vector<int> m_cluster_count;

  ///@name nodes
  //@{
  SvtxTrackMap* m_track_map = nullptr;
  TrkrClusterContainer* m_cluster_map = nullptr;
  //@}

};

#endif
