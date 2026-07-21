#pragma once

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>


#include <string>
#include <vector>

class PHCompositeNode;

class TFile;
class TTree;

class Tpc_ModuleTrackContainer;

class TrkrHitSetContainer;
class TrkrHitSet;

// ===================================================================
// Per-module thread data
// ===================================================================
struct InModuleThreadData
{
  InModuleThreadData();

  struct LayerHitSet
  {
    LayerHitSet();

    unsigned int layer;
    TrkrDefs::hitsetkey hitsetkey;
    TrkrHitSet* hitset;
  };

  struct RawHit
  {
    RawHit();

    unsigned int layer;
    TrkrDefs::hitsetkey hitsetkey;
    TrkrDefs::hitkey hitkey;

    unsigned short pad;
    unsigned short tbin;
    unsigned short adc;
  };

  struct Blob
  {
    Blob();

    unsigned int layer;

    double pad;
    double tbin;
    double adc;

    unsigned int nhits;
    int used;

    std::vector<unsigned int> raw_hit_indices;
  };

  struct Track
  {
    Track();

    unsigned int track_id;

    unsigned int first_layer;
    unsigned int last_layer;

    unsigned int nblobs;
    unsigned int nrawhits;

    // Internal temporary straight-line parameters used only for chain growth
    // and piece connection. They are not copied to Tpc_ModuleTrackContainer or
    // saved as track fit output.
    double pad_slope;
    double pad_intercept;
    double tbin_slope;
    double tbin_intercept;

    std::vector<unsigned int> blob_indices;
    std::vector<unsigned int> raw_hit_indices;
  };

  // Detector/module identity
  unsigned int region;
  unsigned int sector;
  int side;
  TrkrDefs::hitsetkey module_key;


  // General configuration
  double pedestal;
  int verbosity;

  // Noise rejection
  int noise_max_consecutive_timebins;
  int noise_keep_first_timebins;
  int noise_adc_tolerance;

  // Blob building
  int blob_dt;
  int blob_dp;

  // Initial chain growing
  int search_dt;
  int search_dp;
  unsigned int min_track_blobs;

  // Track-piece connection.
  unsigned int connect_max_layer_gap;

  double connect_dp;
  double connect_dt;

  double connect_dpad_slope;
  double connect_dtbin_slope;

  // ADC weighting for temporary internal fits
  double weight_power;
  double adc_weight_floor_frac;

  // Per-module containers
  std::vector<LayerHitSet> layer_hitsets;
  std::vector<RawHit> raw_hits;
  std::vector<Blob> blobs;
  std::vector<Track> tracks;
};

// ===================================================================
// Main Fun4All module
// ===================================================================
class Tpc_ModuleTrackReco : public SubsysReco
{
 public:
  explicit Tpc_ModuleTrackReco(const std::string& name = "Tpc_ModuleTrackReco",
                          const std::string& filename = "Tpc_ModuleTrackReco.root");

  virtual ~Tpc_ModuleTrackReco();

  int Init(PHCompositeNode*);
  int InitRun(PHCompositeNode*);
  int process_event(PHCompositeNode*);
  int End(PHCompositeNode*);

  void setMaxThreads(unsigned int n);

  void setPedestal(double p)
  {
    m_pedestal = p;
  }

  void setBlobWindow(int dt, int dp)
  {
    m_blob_dt = dt;
    m_blob_dp = dp;
  }

  void setSearchWindow(int dt, int dp)
  {
    m_search_dt = dt;
    m_search_dp = dp;
  }

  // Reject long same-pad tails before blob/track building.
  void setNoiseRejection(int max_consecutive_timebins = 10,
                         int keep_first_timebins = 3,
                         int adc_tolerance = 5)
  {
    m_noiseMaxConsecutiveTimebins = max_consecutive_timebins;
    m_noiseKeepFirstTimebins = keep_first_timebins;
    m_noiseAdcTolerance = adc_tolerance;
  }

  void setMinTrackBlobs(unsigned int n)
  {
    m_minTrackBlobs = n;
  }

  void setConnectMaxLayerGap(unsigned int n)
  {
    m_connectMaxLayerGap = n;
  }

  void setConnectWindow(double dt, double dp)
  {
    m_connect_dt = dt;
    m_connect_dp = dp;
  }

  void setConnectSlopeWindow(double dtbin_slope, double dpad_slope)
  {
    m_connect_dtbin_slope = dtbin_slope;
    m_connect_dpad_slope = dpad_slope;
  }

 private:
  int getNodes(PHCompositeNode*);
  void reset_tree_vars();
  int createNodes(PHCompositeNode*);

  std::string m_outputFileName;

  TFile* m_outputFile;
  TTree* m_tree;

  TrkrHitSetContainer* m_hits;

  Tpc_ModuleTrackContainer* m_tpcModuleTrackContainer;
  int m_event;
  unsigned int m_maxThreads;

  // General configuration
  double m_pedestal;

  // Noise rejection
  int m_noiseMaxConsecutiveTimebins;
  int m_noiseKeepFirstTimebins;
  int m_noiseAdcTolerance;

  // Blob building
  int m_blob_dt;
  int m_blob_dp;

  // Initial chain growing
  int m_search_dt;
  int m_search_dp;
  unsigned int m_minTrackBlobs;

  // Track-piece connection parameters
  unsigned int m_connectMaxLayerGap;

  double m_connect_dp;
  double m_connect_dt;

  double m_connect_dpad_slope;
  double m_connect_dtbin_slope;

  // Event number saved once per tree entry
  int m_tree_event;

  // One entry per found module-track. No fit branches are saved here.
  std::vector<unsigned int> m_tree_track_id;
  std::vector<unsigned int> m_tree_region;
  std::vector<unsigned int> m_tree_sector;
  std::vector<int> m_tree_side;

  std::vector<unsigned int> m_tree_nblobs;
  std::vector<unsigned int> m_tree_nrawhits;

  std::vector<unsigned int> m_tree_first_layer;
  std::vector<unsigned int> m_tree_last_layer;

  // Flat per-hit content for TTree reading.
  // Hits are identified by their TrkrHitSetContainer keys only;
  // no hit data is duplicated here.
  std::vector<unsigned int> m_tree_hit_event;
  std::vector<unsigned int> m_tree_hit_track_id;

  std::vector<unsigned int> m_tree_hit_region;
  std::vector<unsigned int> m_tree_hit_sector;
  std::vector<int> m_tree_hit_side;

  std::vector<unsigned int> m_tree_hit_layer;

  std::vector<unsigned long long> m_tree_hit_hitsetkey;
  std::vector<unsigned long long> m_tree_hit_hitkey;
};
