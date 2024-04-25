#ifndef G4MVTX_PHG4MVTXHITRECO_H
#define G4MVTX_PHG4MVTXHITRECO_H

#include <phparameter/PHParameterInterface.h>
#include <trackbase/TrkrDefs.h>

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_rng.h>

#include <map>
#include <memory>  // for unique_ptr
#include <string>
#include <vector>

class ClusHitsVerbosev1;
class PHCompositeNode;
class PHG4Hit;
class PHG4TruthInfoContainer;
class TrkrClusterContainer;
class TrkrHitSetContainer;
class TrkrTruthTrack;
class TrkrTruthTrackContainer;

////// Dead Pixels /////////////
class MvtxRawEvtHeader;

///////////////////////////////

class PHG4MvtxHitReco : public SubsysReco, public PHParameterInterface
{
 public:
  explicit PHG4MvtxHitReco(
      const std::string &name = "PHG4MvtxHitReco",
      const std::string &detector = "MVTX");

  ~PHG4MvtxHitReco() override;

  //! module initialization
  int InitRun(PHCompositeNode *topNode) override;

  //! event processing
  int process_event(PHCompositeNode *topNode) override;

  void Detector(const std::string &d) { m_detector = d; }

  //! TODO keep it for backward compatibily. remove after PR merged
  //In the future use the relevant set parameter function
  void set_timing_window(const int detid, const double tmin, const double tmax);

  //! parameters
  void SetDefaultParameters() override;

  void useRawEvtHeaderNodeName(const std::string& name) { m_MvtxRawEvtHeaderNodeName = name; }

 private:
  std::pair<double, double> generate_alpide_pulse(const double energy_deposited);

  double generate_strobe_zero_tm_start();

  int get_strobe_frame(double alpide_time, double strobe_zero_tm_start);

  TrkrDefs::hitsetkey zero_strobe_bits(TrkrDefs::hitsetkey hitsetkey);

  std::string m_detector;

  double m_tmin;
  double m_tmax;
  double m_strobe_width;
  double m_strobe_separation;
  //double crossing_period = 106.0;
  double m_extended_readout_time = 0.0;

  bool m_in_sphenix_srdo = false;

  class Deleter
  {
   public:
    //! delection operation
    void operator()(gsl_rng *rng) const { gsl_rng_free(rng); }
  };

  std::unique_ptr<gsl_rng, Deleter> m_rng;

  // needed for clustering truth tracks
  private:
  TrkrTruthTrackContainer* m_truthtracks     { nullptr }; // output truth tracks
  TrkrClusterContainer*    m_truthclusters   { nullptr }; // output clusters indexed to TrkrDefs::cluskeys in m_truthtracks
  PHG4TruthInfoContainer*  m_truthinfo       { nullptr };
  int                      m_trkid           { -1      };
  bool                     m_is_emb          { false   };
  TrkrTruthTrack*          m_current_track   { nullptr };
  const int                m_cluster_version { 4 };
  TrkrHitSetContainer*     m_truth_hits; // generate and delete a container for each truth track
  std::map<TrkrDefs::hitsetkey,unsigned int> m_hitsetkey_cnt {}; // counter for making ckeys form hitsetkeys

  std::string m_MvtxRawEvtHeaderNodeName = "MVTXRAWEVTHEADER";
  std::vector<std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey>> m_deadPixelMap;

  PHG4Hit* prior_g4hit { nullptr }; // used to check for jumps in g4hits for loopers;
  void addtruthhitset ( TrkrDefs::hitsetkey, TrkrDefs::hitkey, float neffelectrons );
  void truthcheck_g4hit       ( PHG4Hit*,        PHCompositeNode* topNode );
  void cluster_truthhits      ( PHCompositeNode* topNode          );
  void end_event_truthcluster ( PHCompositeNode* topNode          );

  double m_pixel_thresholdrat { 0.01 };
  float  max_g4hitstep        { 3.5  };
  bool  record_ClusHitsVerbose { false };

  public:
  void set_pixel_thresholdrat (double val) { m_pixel_thresholdrat = val; };

  void set_ClusHitsVerbose(bool set=true) { record_ClusHitsVerbose = set; };
  ClusHitsVerbosev1* mClusHitsVerbose { nullptr };
};

#endif
