// Tell emacs that this is a C++ source
// -*- C++ -*-.

#ifndef G4TPC_PHG4TPCELECTRONDRIFT_H
#define G4TPC_PHG4TPCELECTRONDRIFT_H

#include <fun4all/SubsysReco.h>
#include <g4main/PHG4HitContainer.h>

#include <cmath>
#include <memory>
#include <phparameter/PHParameterInterface.h>

#include <gsl/gsl_rng.h>
#include <string>                              // for string

class PHG4TpcPadPlane;
class PHG4TpcDistortion;
class PHCompositeNode;
class TH1;
class TH2;
class TH3;
class TNtuple;
class TFile;
class TrkrHitSetContainer;
class TrkrHitTruthAssoc;
class DistortedTrackContainer;

class PHG4TpcElectronDrift : public SubsysReco, public PHParameterInterface
{
 public:
  PHG4TpcElectronDrift(const std::string &name = "PHG4TpcElectronDrift");
  virtual ~PHG4TpcElectronDrift() = default;
  virtual int Init(PHCompositeNode*);
  virtual int InitRun(PHCompositeNode*);
  virtual int process_event(PHCompositeNode*);
  virtual int End(PHCompositeNode*);

  void SetDefaultParameters();

  //! detector name
  void Detector(const std::string &d)
  { detector = d; }

  //! detector name
  std::string Detector() const
  { return detector; }

  //! random seed
  void set_seed(const unsigned int iseed);
  void set_time_ordered_distortions_on();
  void set_static_distortions_on();

  //! setup readout plane
  void registerPadPlane(PHG4TpcPadPlane *padplane);

 private:

  //! map a given x,y,z coordinates to plane hits
  void MapToPadPlane(const double x, const double y, const double z, PHG4HitContainer::ConstIterator hiter, TNtuple *ntpad, TNtuple *nthit);

  TrkrHitSetContainer *hitsetcontainer = nullptr;
  TrkrHitTruthAssoc *hittruthassoc = nullptr;
  std::unique_ptr<TrkrHitSetContainer> temp_hitsetcontainer;
  std::unique_ptr<PHG4TpcPadPlane> padplane;

  std::unique_ptr<PHG4TpcDistortion> m_distortionMap;
  int event_num = 0;
  bool do_static_distortion = false;
  bool do_time_ordered_distortion = false;
  bool do_ElectronDriftQAHistos = false;

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

  //! rng de-allocator
  class Deleter
  {
    public:
    //! deletion operator
    void operator() (gsl_rng* rng) const { gsl_rng_free(rng); }
  };
  std::unique_ptr<gsl_rng, Deleter> RandomGenerator;

};

#endif  // G4TPC_PHG4TPCELECTRONDRIFT_H
