// Tell emacs that this is a C++ source
// -*- C++ -*-.

#ifndef G4TPC_PHG4TPCELECTRONDRIFT_H
#define G4TPC_PHG4TPCELECTRONDRIFT_H


#include <phparameter/PHParameterInterface.h>

#include <g4main/PHG4HitContainer.h>

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_rng.h>

#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <trackbase/ActsGeometry.h>

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

class PHG4TpcElectronDrift : public SubsysReco, public PHParameterInterface
{
 public:
  PHG4TpcElectronDrift(const std::string &name = "PHG4TpcElectronDrift");
  ~PHG4TpcElectronDrift() override;
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

  //! setup readout plane
  void registerPadPlane(PHG4TpcPadPlane *padplane);

 private:
  //! map a given x,y,z coordinates to plane hits
  /* TpcClusterBuilder MapToPadPlane(const double x, const double y, const */
  /*     double z, const unsigned int side, PHG4HitContainer::ConstIterator hiter, */
  /*     TNtuple *ntpad, TNtuple *nthit); */

  TrkrHitSetContainer *hitsetcontainer = nullptr;
  TrkrHitTruthAssoc *hittruthassoc = nullptr;
  TrkrTruthTrackContainer *truthtracks = nullptr;
  TrkrTruthTrack *current_track = nullptr;
  TrkrClusterContainer *truthclustercontainer = nullptr; // the TrkrClusterContainer for truth clusters
  std::unique_ptr<TrkrHitSetContainer> temp_hitsetcontainer;
  std::unique_ptr<TrkrHitSetContainer> single_hitsetcontainer;
  std::unique_ptr<PHG4TpcPadPlane> padplane;

  std::unique_ptr<PHG4TpcDistortion> m_distortionMap;
  ActsGeometry* m_tGeometry;
  PHG4TpcCylinderGeomContainer *seggeo;

  int event_num = 0;
  bool do_ElectronDriftQAHistos = false;

  ///@name evaluation histograms
  //@{
  TH1 *dlong = nullptr;
  TH1 *dtrans = nullptr;
  TH2 *hitmapstart = nullptr;
  TH2 *hitmapend = nullptr;
  TH2 *z_startmap = nullptr;
  TH2 *deltaphi = nullptr;
  TH2 *deltar = nullptr;
  TH2 *deltaphinodiff = nullptr;
  TH2 *deltaRphinodiff = nullptr;
  TH2 *deltaphivsRnodiff = nullptr;
  TH2 *deltaphinodist = nullptr;
  TH2 *deltarnodiff = nullptr;
  TH2 *deltarnodist = nullptr;
  TH2 *deltaz = nullptr;
  //@}

  std::unique_ptr<TFile> m_outf;
  std::unique_ptr<TFile> EDrift_outf;

  TNtuple *nt = nullptr;
  TNtuple *nthit = nullptr;
  TNtuple *ntfinalhit = nullptr;
  TNtuple *ntpad = nullptr;

  std::string detector;
  std::string hitnodename;
  std::string seggeonodename;

  double diffusion_trans = NAN;
  double added_smear_sigma_trans = NAN;
  double diffusion_long = NAN;
  double added_smear_sigma_long = NAN;
  double drift_velocity = NAN;
  double tpc_length = NAN;
  double electrons_per_gev = NAN;
  double min_active_radius = NAN;
  double max_active_radius = NAN;
  double min_time = NAN;
  double max_time = NAN;


  /* std::array<TpcClusterBuilder,55> layer_clusterers; // Generate TrkrClusterv4's for TrkrTruthTracks */
  TpcClusterBuilder* truth_clusterer { nullptr };

  /* void buildTruthClusters(std::map<TrkrDefs::hitsetkey,unsigned int>&); */

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
