#ifndef G4EVAL_MicromegasDriftEvaluator_H
#define G4EVAL_MicromegasDriftEvaluator_H

/*
 * Bade Sayki June 10th, 2026 -- LANL
 * This module is created to calibrate the drift velocity in the TPC by fitting a helix to the clusters within a certain layer range, and projecting it to the TPOT z view module plane. The default layers in the TPC are set to be 39-55, which correspond to R3.
 * This is heavily inspired by and distilled from Dr. Hugo Pereira Da Costa's MicromegasTrackEvaluator_hp module. It is meant to be a more lightweight and specialized version.
 * It accumulates a TH3F(tile, z_track, dz) histogram during process_event, then in End() fits a piecewise function to suggest an updated drift velocity.
 * If you have any questions, please feel free to message me on mattermost. 
 * Claude Code tool was used to format and comment this module.
 */

#include <fun4all/SubsysReco.h>
#include <micromegas/MicromegasDefs.h>
#include <phool/PHObject.h>
#include <tpc/TpcGlobalPositionWrapper.h>
#include <trackbase/TrkrDefs.h>

#include <array>
#include <string>
#include <vector>

class ActsGeometry;
class PHG4CylinderGeomContainer;
class TH3F;
class TrkrCluster;
class TrkrClusterContainer;
class SvtxTrackMap;

class MicromegasDriftEvaluator : public SubsysReco
{
  public:

  explicit MicromegasDriftEvaluator( const std::string& name = "MicromegasDriftEvaluator" );

  int Init(PHCompositeNode*) override;
  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  struct TrackStateStruct
  {
    unsigned short _layer  = 0; 
    unsigned short _tile   = 0;
    double         _z      = 0;  
    double         _y_local = 0; 
  };

  struct ClusterStruct
  {
    unsigned short _layer = 0; 
    unsigned short _tile  = 0;
    double         _z     = 0;
  };

  struct TrackStruct
  {
    float _chisquare = 0;
    int   _ndf       = 0;

    unsigned int _nclusters_tpc        = 0;
    unsigned int _nclusters_mvtx       = 0;
    unsigned int _nclusters_intt       = 0;
    unsigned int _nclusters_micromegas = 0;

    TrackStateStruct _trk_state_z; 
    ClusterStruct    _found_cluster_z;

    using List = std::vector<TrackStruct>;
  };


  class Container : public PHObject
  {
    public:

    explicit Container() = default;
    Container(const Container&) = delete;
    Container& operator=(const Container&) = delete;

    void Reset() override { _tracks.clear(); }

    const TrackStruct::List& tracks() const { return _tracks; }
    void add_track(const TrackStruct& t)    { _tracks.push_back(t); }
    void clear_tracks()                     { _tracks.clear(); }

    private:

    TrackStruct::List _tracks;

    TrackStateStruct _unused_state;
    ClusterStruct    _unused_cluster;

    ClassDefOverride(Container, 1)
  };

  void set_trackmapname(const std::string& value) { m_trackmapname = value; }

  /// This function is specifically used to give the fitting function a starting point. Use the initial drift velocity you used when reconstructing.
  void set_drift_velocity(double value) { m_drift_velocity = value; }

  /// TPC layer range used for the helix fit. The default is R3, but this is an area with huge static distortions. It can easily be adjusted in the Fun4All macro with these functions.
  void set_min_tpc_layer(unsigned int value) { m_min_tpc_layer = value; }
  void set_max_tpc_layer(unsigned int value) { m_max_tpc_layer = value; }

  /// This one rejects track states near tile edge
  void set_y_local_cut(double value) { m_y_local_cut = value; }

  /// Search window to match a Micromegas cluster to the prediction
  void set_z_search_window(double value) { m_z_search_win = value; }

  /// Output filename for the QA plot. Make this a .png
  void set_plot_filename(const std::string& value) { m_plot_filename = value; }

  /// Output ROOT filename for histograms and fit results.
  void set_root_filename(const std::string& value) { m_root_filename = value; }

  /// If true (default), append -<runnumber>-<segment> to output filenames, following sPHENIX convention
  void set_add_run_segment(bool value) { m_add_run_segment = value; }

  /// Manually set the segment number used in output filenames (otherwise parsed from the input filename)
  void set_segment(int value) { m_segment = value; }

  private:

  int  load_nodes(PHCompositeNode*);
  std::string make_output_filename(const std::string&) const;
  void evaluate_tracks();

  Container*                 m_container             = nullptr;
  ActsGeometry*              m_tGeometry             = nullptr;
  TpcGlobalPositionWrapper   m_globalPositionWrapper;
  PHG4CylinderGeomContainer* m_micromegas_geomcontainer = nullptr;
  TrkrClusterContainer*      m_cluster_map           = nullptr;
  SvtxTrackMap*              m_track_map             = nullptr;

  std::string  m_trackmapname   = "SvtxTrackMap";

  //These are all adjustable in your F4A macro. You should probably put in a better m_plot_filename.
  double       m_drift_velocity = 0.00747; 
  unsigned int m_min_tpc_layer  = 39;
  unsigned int m_max_tpc_layer  = 55;
  double       m_y_local_cut    = 22.0;
  double       m_z_search_win   = 3.0;
  std::string  m_plot_filename  = "micromegas_drift_calib.png";
  std::string  m_root_filename  = "micromegas_drift_calib.root";
  bool         m_add_run_segment = true;
  int          m_segment        = -1;

  TH3F* m_hist3D = nullptr;
};

#endif