#ifndef TPCCALIB_TPCDIRECTLASERRECONSTRUCTION_H
#define TPCCALIB_TPCDIRECTLASERRECONSTRUCTION_H

/**
 * \file TpcDirectLaserReconstruction.h
 * \brief performs the reconstruction of TPC direct laser tracks
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <fun4all/SubsysReco.h>
#include <phparameter/PHParameterInterface.h>
#include <trackbase_historic/ActsTransformations.h>

#include <memory>

class SvtxTrack;
class SvtxTrackMap;
class TpcSpaceChargeMatrixContainer;
class TrkrHitSetContainer;
class PHG4TpcCylinderGeomContainer;

class TFile;
class TH1;
class TH2;
class TH3;
class TH2;
class TVector3;

class TpcDirectLaserReconstruction: public SubsysReco, public PHParameterInterface
{

  public:

  /// constructor
  TpcDirectLaserReconstruction( const std::string& = "TpcDirectLaserReconstruction" );

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

  /// output file
  /**
   * this is the file where space charge matrix container is stored
   */
  void set_outputfile( const std::string& filename )
  { m_outputfile = filename; }

  /// set to true to store evaluation histograms and ntuples
  void set_savehistograms( bool value ) { m_savehistograms = value; }

  /// output file name for evaluation histograms
  void set_histogram_outputfile(const std::string &outputfile)
  {m_histogramfilename = outputfile;}

  /// set grid dimensions
  void set_grid_dimensions( int phibins, int rbins, int zbins );

  void set_max_zrange(float length)
  {m_max_zrange = length;}

  void set_max_dca(float length)
  {m_max_dca = length;}

  private:

  /// load nodes
  int load_nodes( PHCompositeNode* );

  /// create evaluation histograms
  void create_histograms();

  /// process tracks
  void process_tracks();

  /// process track
  void process_track( SvtxTrack* );

  /// get relevant cell for a given hit
  int get_cell_index( const TVector3& ) const;

  /// output file
  std::string m_outputfile = "TpcSpaceChargeMatrices.root";

  float m_max_zrange = 10.0; // cm

  ///@name selection parameters
  //@{
  // residual cuts in r, phi plane
  float m_max_dca = NAN;

  /// residual cuts in r, phi plane
  float m_max_drphi = NAN;

  /// residual cuts in r, z plane
  float m_max_dz = NAN;

  float m_pedestal = 74.4;  // pedestal for hit ASDC values

  //@}

  /// matrix container
  std::unique_ptr<TpcSpaceChargeMatrixContainer> m_matrix_container;

  ///@name counters
  //@{
  int m_total_hits = 0;
  int m_matched_hits = 0;
  int m_accepted_clusters = 0;
  //@}

  ///@name nodes
  //@{

  PHG4TpcCylinderGeomContainer *m_geom_container = nullptr;

  /// Acts geometry
  ActsGeometry *m_tGeometry = nullptr;

  /// acts transformation
  ActsTransformations m_transformer;

  /// tracks
  SvtxTrackMap* m_track_map = nullptr;
  
  TrkrHitSetContainer* m_hit_map = nullptr;

  //@}

  ///@name evaluation
  //@{
  bool m_savehistograms = false;
  std::string m_histogramfilename = "TpcDirectLaserReconstruction.root";
  std::unique_ptr<TFile> m_histogramfile = nullptr;

  /// dca vs layer number
  TH2* h_dca_layer = nullptr;

  /// delta rphi vs layer number
  TH2 *h_deltarphi_layer_south = nullptr;
  TH2 *h_deltarphi_layer_north = nullptr;

  /// delta z vs layer number
  TH2 *h_deltaz_layer = nullptr;

  TH2 *h_deltar_r = nullptr;

  /// number of entries per cell
  TH3 *h_entries = nullptr;

  TH2* h_xy = nullptr;
  TH2* h_xz = nullptr;
  TH2* h_xy_pca = nullptr;
  TH2* h_xz_pca = nullptr;
  TH2* h_dca_path = nullptr;
  TH2* h_zr = nullptr;
  TH2* h_zr_pca = nullptr;
  TH2* h_dz_z = nullptr;

  //@}

};

#endif
