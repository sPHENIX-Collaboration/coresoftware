#ifndef TPCCALIB_TPCSPACECHARGERECONSTRUCTION_H
#define TPCCALIB_TPCSPACECHARGERECONSTRUCTION_H
/**
 * \file TpcSpaceChargeReconstruction.h
 * \brief performs space charge distortion reconstruction using tracks
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */
#include <fun4all/SubsysReco.h>
#include <phparameter/PHParameterInterface.h>

#include <TString.h>

#include <memory>
#include <vector>

// forward declaration
class SvtxTrack;
class SvtxTrackMap;
class TpcSpaceChargeMatrixContainer;
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
 The actual inversion is performed in TpcSpaceChargeMatrixInversion::calculate_distortions
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
   * this is the file where space charge matrix container is stored 
   */
  void set_outputfile( const std::string& filename )
  { m_outputfile = filename; }

  //@}

  /// global initialization
  int Init(PHCompositeNode*) override;

  /// run initialization
  int InitRun(PHCompositeNode*) override;

  /// event processing
  int process_event(PHCompositeNode*) override;

  /// end of processing
  int End(PHCompositeNode*) override;

  /// parameters
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

  /// get relevant cell for a given cluster
  int get_cell_index( TrkrCluster* ) const;

  /// output file
  std::string m_outputfile = "TpcSpaceChargeMatrices.root";

  /// true if only tracks with micromegas must be used
  bool m_use_micromegas = true;

  ///@name selection parameters
  //@{
  // residual cuts in r, phi plane
  float m_max_talpha = 0.6;
  float m_max_drphi = 0.5;

  // residual cuts in r, z plane
  float m_max_tbeta = 1.5;
  float m_max_dz = 0.5;
  //@}

  /// matrix container
  std::unique_ptr<TpcSpaceChargeMatrixContainer> m_matrix_container;
  
  ///@name counters
  //@{
  int m_total_tracks = 0;
  int m_accepted_tracks = 0;

  int m_total_clusters = 0;
  int m_accepted_clusters = 0;
  //@}
  
  ///@name nodes
  //@{
  SvtxTrackMap* m_track_map = nullptr;
  TrkrClusterContainer* m_cluster_map = nullptr;
  //@}

};

#endif
