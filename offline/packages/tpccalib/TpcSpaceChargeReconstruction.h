#ifndef TPCCALIB_TPCSPACECHARGERECONSTRUCTION_H
#define TPCCALIB_TPCSPACECHARGERECONSTRUCTION_H
/**
 * \file TpcSpaceChargeReconstruction.h
 * \brief performs space charge distortion reconstruction using tracks
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */
#include <fun4all/SubsysReco.h>
#include <phparameter/PHParameterInterface.h>

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

class TpcSpaceChargeReconstruction: public SubsysReco, public PHParameterInterface
{
  public:

  /// constructor
  TpcSpaceChargeReconstruction( const std::string& = "TPCSPACECHARGERECONSTRUCTION" );

  ///@name configuration
  //@{

  /// set whether to use only tracks with micromegas or not
  void set_use_micromegas( bool value )
  { m_use_micromegas = value; }

  /// set grid dimensions
  /**
  \param phibins the number of bins in the azimuth direction
  \param zbins the number of bins along z
  */
  void set_grid_dimensions( int phibins, int rbins, int zbins );

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

  //! parameters
  void SetDefaultParameters() override;

  private:

  /// load nodes
  int load_nodes( PHCompositeNode* );

  /// process tracks
  void process_tracks();

  /// returns true if track fulfills basic requirement for distortion calculations
  bool accept_track( SvtxTrack* ) const;

  /// process track
  void process_track( SvtxTrack* );

  /// calculate distortions
  void calculate_distortions( PHCompositeNode* );

  /// get cell from phi, r and z index
  int get_cell( int iphi, int ir, int iz ) const;

  /// get cell from phi, r and z values
  int get_cell( float phi, float r, float z ) const;

  /// get relevant cell for a given cluster
  int get_cell( TrkrCluster* ) const;

  /// output file
  std::string m_outputfile = "TpcSpaceChargeReconstruction.root";

  /// true if only tracks with micromegas must be used
  bool m_use_micromegas = true;

  ///@name grid dimensions
  //@{

  // phi range
  static constexpr float m_phimin = 0;
  static constexpr float m_phimax = 2.*M_PI;

  // TODO: could try to get the r and z range from TPC geometry
  // r range
  static constexpr float m_rmin = 20;
  static constexpr float m_rmax = 78;

  // z range
  static constexpr float m_zmin = -105.5;
  static constexpr float m_zmax = 105.5;

  //@}

  //!@name grid size
  //@{
  int m_phibins = 36;
  int m_rbins = 16;
  int m_zbins = 80;
  int m_totalbins = m_phibins*m_rbins*m_zbins;
  //@}

  //!@name selection parameters
  //@{
  // residual cuts in r, phi plane
  float m_max_talpha = 0.6;
  float m_max_drphi = 0.5;

  // residual cuts in r, z plane
  float m_max_tbeta = 1.5;
  float m_max_dz = 0.5;
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
