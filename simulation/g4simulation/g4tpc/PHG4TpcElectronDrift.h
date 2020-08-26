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
class PHCompositeNode;
class TH1;
class TNtuple;
class TFile;
class TrkrHitSetContainer;
class TrkrHitTruthAssoc;

class PHG4TpcElectronDrift : public SubsysReco, public PHParameterInterface
{
 public:
  PHG4TpcElectronDrift(const std::string &name = "PHG4TpcElectronDrift");
  virtual ~PHG4TpcElectronDrift() = default;
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void SetDefaultParameters();

  void Detector(const std::string &d) { detector = d; }
  std::string Detector() const { return detector; }
  void set_seed(const unsigned int iseed);
  void MapToPadPlane(const double x, const double y, const double z, PHG4HitContainer::ConstIterator hiter, TNtuple *ntpad, TNtuple *nthit);
  void registerPadPlane(PHG4TpcPadPlane *padplane);

 private:
  TrkrHitSetContainer *hitsetcontainer = nullptr;
  TrkrHitTruthAssoc *hittruthassoc = nullptr;
  std::unique_ptr<TrkrHitSetContainer> temp_hitsetcontainer;
  std::unique_ptr<PHG4TpcPadPlane> padplane;
  TH1 *dlong = nullptr;
  TH1 *dtrans = nullptr;
  TFile *m_outf = nullptr;
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
