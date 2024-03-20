// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef G4TPC_PHG4TPCELECTRONDRIFT_H
#define G4TPC_PHG4TPCELECTRONDRIFT_H

#include "TpcClusterBuilder.h"

#include <trackbase/ActsGeometry.h>

#include <g4main/PHG4HitContainer.h>

#include <phparameter/PHParameterInterface.h>

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_rng.h>

#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <memory>
#include <string>

class PHG4TpcPadPlane;
class PHG4TpcDistortion;
class PHCompositeNode;
class TH1;
class TH2;
class TNtuple;
class TFile;
class TrkrHitSetContainer;
class TrkrHitTruthAssoc;
class TrkrTruthTrackContainer;
class TrkrClusterContainer;
class TrkrTruthTrack;
class DistortedTrackContainer;
class TpcClusterBuilder;
class PHG4TpcCylinderGeomContainer;
class ClusHitsVerbose;

class PHG4TpcElectronDrift : public SubsysReco, public PHParameterInterface
{
 public:
  PHG4TpcElectronDrift(const std::string &name = "PHG4TpcElectronDrift");
  ~PHG4TpcElectronDrift() override{};
  int Init(PHCompositeNode *) override;
  int InitRun(PHCompositeNode *) override;
  int process_event(PHCompositeNode *) override;
  int End(PHCompositeNode *) override;

  void SetDefaultParameters() override;

  //! detector name
  void Detector(const std::string &d)
  {
    detector = d;
  }

  //! detector name
  std::string Detector() const
  {
    return detector;
  }

  //! random seed
  void set_seed(const unsigned int iseed);

  //! setup TPC distortion
  void setTpcDistortion(PHG4TpcDistortion *);

  //
  void set_flag_threshold_distortion(bool setflag, float setthreshold);

  //! setup readout plane
  void registerPadPlane(PHG4TpcPadPlane *padplane);

  // cluster the PHG4Tracks individually
  TpcClusterBuilder truth_clusterer{};
  void set_pixel_thresholdrat(double val) { truth_clusterer.set_pixel_thresholdrat(val); };
  void set_max_g4hitstep(float f) { max_g4hitstep = f; };
  void set_ClusHitsVerbose(bool set = true) { record_ClusHitsVerbose = set; };
  void set_zero_bfield_flag(bool flag) { zero_bfield = flag; };
  void set_zero_bfield_diffusion_factor(double f) { zero_bfield_diffusion_factor = f; };
  ClusHitsVerbosev1 *mClusHitsVerbose{nullptr};

 private:
  TrkrHitSetContainer *hitsetcontainer{nullptr};
  TrkrHitTruthAssoc *hittruthassoc{nullptr};
  TrkrTruthTrackContainer *truthtracks{nullptr};
  TrkrTruthTrack *truth_track{nullptr};
  TrkrClusterContainer *truthclustercontainer{nullptr};  // the TrkrClusterContainer for truth clusters
  ActsGeometry *m_tGeometry{nullptr};
  PHG4TpcCylinderGeomContainer *seggeo{nullptr};

  TNtuple *nt{nullptr};
  TNtuple *nthit{nullptr};
  TNtuple *ntfinalhit{nullptr};
  TNtuple *ntpad{nullptr};

  ///@name evaluation histograms
  //@{
  TH1 *dlong{nullptr};
  TH1 *dtrans{nullptr};
  TH1 *ratioElectronsRR{nullptr};
  TH2 *hitmapstart{nullptr};
  TH2 *hitmapend{nullptr};
  TH2 *hitmapstart_z{nullptr};
  TH2 *hitmapend_z{nullptr};
  TH2 *z_startmap{nullptr};
  TH2 *deltaphi{nullptr};
  TH2 *deltar{nullptr};
  TH2 *deltaphinodiff{nullptr};
  TH2 *deltaRphinodiff{nullptr};
  TH2 *deltaphivsRnodiff{nullptr};
  TH2 *deltaphinodist{nullptr};
  TH2 *deltarnodiff{nullptr};
  TH2 *deltarnodist{nullptr};
  TH2 *deltaz{nullptr};
  //@}

  int event_num{0};

  float max_g4hitstep{7.};
  float thresholdforreachesreadout{0.5};

  double diffusion_trans = std::numeric_limits<double>::signaling_NaN();
  double added_smear_sigma_trans = std::numeric_limits<double>::signaling_NaN();
  double diffusion_long = std::numeric_limits<double>::signaling_NaN();
  double added_smear_sigma_long = std::numeric_limits<double>::signaling_NaN();
  double drift_velocity = std::numeric_limits<double>::signaling_NaN();
  double tpc_length = std::numeric_limits<double>::signaling_NaN();
  double electrons_per_gev = std::numeric_limits<double>::signaling_NaN();
  double min_active_radius = std::numeric_limits<double>::signaling_NaN();
  double max_active_radius = std::numeric_limits<double>::signaling_NaN();
  double min_time = std::numeric_limits<double>::signaling_NaN();
  double max_time = std::numeric_limits<double>::signaling_NaN();
  double zero_bfield_diffusion_factor{3.5};  // at drift voltage of 400 V

  bool record_ClusHitsVerbose{false};
  bool do_ElectronDriftQAHistos{false};
  bool do_getReachReadout{false};
  bool zero_bfield{false};

  std::unique_ptr<TrkrHitSetContainer> temp_hitsetcontainer;
  std::unique_ptr<TrkrHitSetContainer> single_hitsetcontainer;
  std::unique_ptr<PHG4TpcPadPlane> padplane;
  std::unique_ptr<PHG4TpcDistortion> m_distortionMap;
  std::unique_ptr<TFile> m_outf;
  std::unique_ptr<TFile> EDrift_outf;

  std::string detector;
  std::string hitnodename;
  std::string seggeonodename;

  //! rng de-allocator
  class Deleter
  {
   public:
    //! deletion operator
    void operator()(gsl_rng *rng) const { gsl_rng_free(rng); }
  };
  std::unique_ptr<gsl_rng, Deleter> RandomGenerator;
};

#endif  // G4TPC_PHG4TPCELECTRONDRIFT_H
