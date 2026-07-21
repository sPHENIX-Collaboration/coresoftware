#ifndef QA_TRACKING_MICROMEGASDRIFTQA_H
#define QA_TRACKING_MICROMEGASDRIFTQA_H

/*
 * Bade Sayki June 10th, 2026 -- LANL
 * This QA module is created to monitor the calibration of the drift velocity in the TPC by fitting a helix to the clusters within a certain layer range, and projecting it to the TPOT z view module plane. The default layers in the TPC are set to be 39-55, which correspond to R3.
 * This is heavily inspired by and distilled from Dr. Hugo Pereira Da Costa's MicromegasTrackEvaluator_hp module. It is meant to be a more lightweight and specialized version.
 * If you have any questions, please feel free to message me on mattermost.
 * Claude Code tool was used to format and debug this module
 */

#include <fun4all/SubsysReco.h>
#include <tpc/TpcGlobalPositionWrapper.h>

#include <string>

class ActsGeometry;
class PHCompositeNode;
class PHG4CylinderGeomContainer;
class SvtxTrackMap;
class TrkrClusterContainer;
class TH1;
class TH2;

class MicromegasDriftQA : public SubsysReco
{
 public:
  explicit MicromegasDriftQA(const std::string& name = "MicromegasDriftQA");

  ~MicromegasDriftQA() override = default;

  //! run initialization: load nodes, create and register histograms
  int InitRun(PHCompositeNode* topNode) override;

  //! event processing: fill histograms
  int process_event(PHCompositeNode* topNode) override;

  //! end of processing: fit accumulated distributions, fill summary histogram
  int End(PHCompositeNode* topNode) override;

  //! track map name
  void set_trackmapname(const std::string& value) { m_trackmapname = value; }

  //! initial drift velocity (cm/ns); starting point for the fit. Use the drift velocity used at reconstruction.
  void set_drift_velocity(double value) { m_drift_velocity = value; }

  //! TPC layer range used for the helix fit (default: R3)
  void set_min_tpc_layer(unsigned int value) { m_min_tpc_layer = value; }
  void set_max_tpc_layer(unsigned int value) { m_max_tpc_layer = value; }

  //! reject track states near the tile edge (cm, local y)
  void set_y_local_cut(double value) { m_y_local_cut = value; }

  //! search window to match a Micromegas cluster to the prediction (cm)
  void set_z_search_window(double value) { m_z_search_win = value; }

  //! minimum entries per z slice required by FitSlicesY
  void set_min_slice_entries(int value) { m_min_slice_entries = value; }

 private:
  int load_nodes(PHCompositeNode* topNode);

  void createHistos();
  std::string getHistoPrefix() const;

  //!@name histograms (owned by the QA histogram manager)
  //@{

  //! z_track vs dz, one per z-view tile
  TH2* h_ztrk_dz[8]{nullptr};

  //! dz = z_track - z_cluster, all tiles
  TH1* h_dz{nullptr};

  //! matched track states per tile
  TH1* h_tile{nullptr};

  //! local y of the track state on the tile
  TH1* h_ylocal{nullptr};

  //! number of matched track states per event
  TH1* h_ntracks{nullptr};

  //! drift velocity fit summary, filled in End()
  TH1* h_driftSummary{nullptr};

  //@}

  //!@name nodes
  //@{
  ActsGeometry* m_tGeometry{nullptr};
  TpcGlobalPositionWrapper m_globalPositionWrapper;
  PHG4CylinderGeomContainer* m_micromegas_geomcontainer{nullptr};
  TrkrClusterContainer* m_cluster_map{nullptr};
  SvtxTrackMap* m_track_map{nullptr};
  //@}

  //! track map name
  std::string m_trackmapname{"SvtxTrackMap"};

  //! initial drift velocity (cm/ns)
  double m_drift_velocity{0.00745};

  //! TPC layer range used for the helix fit
  unsigned int m_min_tpc_layer{39};
  unsigned int m_max_tpc_layer{55};

  //! reject track states near the tile edge (cm)
  double m_y_local_cut{22.0};

  //! search window to match a Micromegas cluster to the prediction (cm)
  double m_z_search_win{3.0};

  //! minimum entries per z slice for FitSlicesY
  int m_min_slice_entries{10};
};

#endif  // QA_TRACKING_MICROMEGASDRIFTQA_H
