// Tell emacs that this is a C++ source
// -*- C++ -*-.

#ifndef G4TPC_PHG4TPCELECTRONDRIFT_H
#define G4TPC_PHG4TPCELECTRONDRIFT_H
#include "PHG4TpcDistortion.h"
#include <fun4all/SubsysReco.h>
#include <g4main/PHG4HitContainer.h>

#include <phparameter/PHParameterInterface.h>

#include <gsl/gsl_rng.h>

#include <string>                              // for string

class PHG4TpcPadPlane;
class PHG4TpcDistortion;
class PHCompositeNode;
class TH1;
class TH2;
class TH3F;
class TAxis;
class TGraph;
class TFile;
class TTree;
class TNtuple;
class TFile;
class TrkrHitSetContainer;
class TrkrHitTruthAssoc;
class Fun4AllHistoManager;

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
  void set_time_ordered_distortions_on();
  void set_static_distortions_on();
  void MapToPadPlane(const double x, const double y, const double z, PHG4HitContainer::ConstIterator hiter, TNtuple *ntpad, TNtuple *nthit);
  void registerPadPlane(PHG4TpcPadPlane *padplane);


 private:

  TrkrHitSetContainer *hitsetcontainer = nullptr;
  TrkrHitSetContainer *temp_hitsetcontainer = nullptr;
  TrkrHitTruthAssoc *hittruthassoc = nullptr;
  PHG4TpcPadPlane *padplane = nullptr;

  PHG4TpcDistortion* DistortionMap = nullptr;
  int event_num;
  bool do_static_distortion = false;
  bool do_time_ordered_distortion = false;
  bool do_ElectronDriftQAHistos;

  TH1 *dlong;
  TH1 *dtrans;
  TH2 *hitmapstart;
  TH2 *hitmapend;
  TH2 *z_startmap;
  TH2 *deltaphi;
  TH2 *deltar;
  TH2 *deltaphinodiff;
  TH2 *deltaRphinodiff;
  TH2 *deltaphivsRnodiff;
  TH2 *deltaphinodist;
  TH2 *deltarnodiff;
  TH2 *deltarnodist;
  TH2 *deltaz;
  TFile *m_outf;
  TFile *EDrift_outf;
  TNtuple *nt;
  TNtuple *nthit;
  TNtuple *ntfinalhit;
  TNtuple *ntpad;
  std::string detector;
  std::string hitnodename;
  std::string seggeonodename;
  unsigned int seed,print_layer;
  int nBinZ, nBinR,nBinP,e_num;
  double x_start,y_start,x_final,y_final;
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

  gsl_rng *RandomGenerator;
};

#endif  // G4TPC_PHG4TPCELECTRONDRIFT_H
