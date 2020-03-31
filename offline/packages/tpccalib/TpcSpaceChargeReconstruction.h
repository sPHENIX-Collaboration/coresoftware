#ifndef TPCCALIB_TPCSPACECHARGERECONSTRUCTION_H
#define TPCCALIB_TPCSPACECHARGERECONSTRUCTION_H
/**
 * \file TpcSpaceChargeReconstruction.h
 * \brief performs space charge distortion reconstruction using tracks
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */
#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>

// forward declaration
class SvtxTrack;
class SvtxTrackMap;
class TrkrCluster;
class TrkrClusterContainer;

/**
 * \class TpcSpaceChargeReconstruction
 * \brief performs space charge distortion reconstruction using tracks
 * \detail To reconstruct the distortions dr0, drphi0 and dz0 in a given volume element, the following chisquare is minimized
 chisquare = sum_cluster (drphi - (drphi0 + dr0 tan alpha))**2/error**2 + sum_cluster ( dz - (dz0 + dr0 tan beta))**2/error**2
 with
 - drphi and dz the residuals (track - cluster) measured for a given cluster
 - alpha and beta the track angles at the cluster in the rphi,r plane and the z,r plane, respectively
 The chisquare being quadratic in dr0, drphi0 and dz0, it can be minimized analytically.
 This results in a linear equation lhs[i].[corrections] = rhs[i], and thus [corrections] = lhs[i]**(-1).rhs[i]
 The lhs and rhs matrices are filled in TpcSpaceChargeReconstruction::process_track
 The inversion is performed in TpcSpaceChargeReconstruction::calculate_distortions
 */

class TpcSpaceChargeReconstruction: public SubsysReco
{
  public:

  /// constructor
  TpcSpaceChargeReconstruction( const std::string& = "TPCSPACECHARGERECONSTRUCTION" );

  ///@name configuration
  //@{

  /// set tpc layers
  // TODO: use recoconst instead ?
  void set_tpc_layers( unsigned int first_layer, unsigned int n_layers );

  /// set grid dimensions
  // TODO: use recoconst instead ?
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
  void calculate_distortions( PHCompositeNode* );

  /// get cell from z, r and phi index
  int get_cell( int iz, int ir, int iphi ) const;

  /// get relevant cell for a given cluster
  int get_cell( TrkrDefs::cluskey, TrkrCluster* ) const;

  /// output file
  std::string m_outputfile = "TpcSpaceChargeReconstruction.root";

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
