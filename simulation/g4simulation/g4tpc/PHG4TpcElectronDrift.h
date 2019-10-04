// Tell emacs that this is a C++ source
// -*- C++ -*-.

#ifndef G4TPC_PHG4TPCELECTRONDRIFT_H
#define G4TPC_PHG4TPCELECTRONDRIFT_H

#include <fun4all/SubsysReco.h>
#include <g4main/PHG4HitContainer.h>

#include <phparameter/PHParameterInterface.h>

// rootcint barfs with this header so we need to hide it
#if !defined(__CINT__) || defined(__CLING__)
#include <gsl/gsl_rng.h>
#endif

#include <string>                              // for string

class PHG4TpcPadPlane;
class PHCompositeNode;
class TH1;
class TNtuple;
class TrkrHitSetContainer;
class TrkrHitTruthAssoc;

class PHG4TpcElectronDrift : public SubsysReco, public PHParameterInterface
{
 public:
  PHG4TpcElectronDrift(const std::string &name = "PHG4TpcElectronDrift");
  virtual ~PHG4TpcElectronDrift();
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
  TrkrHitSetContainer *hitsetcontainer;
  TrkrHitSetContainer *temp_hitsetcontainer;
  TrkrHitTruthAssoc *hittruthassoc;
  PHG4TpcPadPlane *padplane;
  TH1 *dlong;
  TH1 *dtrans;
  TNtuple *nt;
  TNtuple *nthit;
  TNtuple *ntpad;
  std::string detector;
  std::string hitnodename;
  std::string seggeonodename;
  unsigned int seed;
  double diffusion_trans;
  double added_smear_sigma_trans;
  double diffusion_long;
  double added_smear_sigma_long;
  double drift_velocity;
  double tpc_length;
  double electrons_per_gev;
  double min_active_radius;
  double max_active_radius;
  double min_time;
  double max_time;

#if !defined(__CINT__) || defined(__CLING__)
  gsl_rng *RandomGenerator;
#endif
};

#endif  // G4TPC_PHG4TPCELECTRONDRIFT_H
