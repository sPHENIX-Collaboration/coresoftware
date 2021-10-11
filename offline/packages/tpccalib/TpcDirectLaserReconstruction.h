#ifndef TPCCALIB_TPCDIRECTLASERRECONSTRUCTION_H
#define TPCCALIB_TPCDIRECTLASERRECONSTRUCTION_H

/**
 * \file TpcDirectLaserReconstruction.h
 * \brief performs the reconstruction of TPC direct laser tracks
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <fun4all/SubsysReco.h>
#include <phparameter/PHParameterInterface.h>

#include <memory>

class SvtxTrack;
class SvtxTrackMap;
class TpcSpaceChargeMatrixContainer;
class TrkrClusterContainer;
class TrkrHitSetContainer;

class TFile;
class TH1;
class TH2;
class TH3;

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

  private:

  /// load nodes
  int load_nodes( PHCompositeNode* );

  /// create evaluation histograms
  void create_histograms();

  /// process tracks
  void process_tracks();

  /// process track
  void process_track( SvtxTrack* );

  /// output file
  std::string m_outputfile = "TpcSpaceChargeMatrices.root";

  ///@name nodes
  //@{
  TrkrHitSetContainer* m_hitsetcontainer = nullptr;
  SvtxTrackMap* m_track_map = nullptr;
  TrkrClusterContainer* m_cluster_map = nullptr;
  //@}

  /// matrix container
  std::unique_ptr<TpcSpaceChargeMatrixContainer> m_matrix_container;

  ///@name selection parameters
  //@{
  // residual cuts in r, phi plane
  float m_max_dca = 1.5;
  //@}

  ///@name evaluation
  //@{
  bool m_savehistograms = false;
  std::string m_histogramfilename = "TpcDirectLaserReconstruction.root";
  std::unique_ptr<TFile> m_histogramfile = nullptr;

  /// dca vs layer number
  TH2* h_dca_layer = nullptr;

  /// delta rphi vs layer number
  TH2 *h_deltarphi_layer = nullptr;

  /// delta z vs layer number
  TH2 *h_deltaz_layer = nullptr;

  /// number of entries per cell
  TH3 *h_entries = nullptr;

  //@}

};

#endif
