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
  int DistortionIntegral(double radstart,double phistart,double z_start,double* rad_final, double* phi_final);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  int event_num;
  TH3F *TimehDR;
  TH3F *TimehDP;
  TH3F *TimehDZ;
  TH3F *TimeInthDR;
  TH3F *TimeInthDP;
  TH3F *TimeInthDZ;
    
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
  TFile *DistFile;
  TFile *TimeDistFile;
  TFile *CMFile;
  TTree *TimeTree;
  TTree *CMTimeDists;
  TH3F *hDRint;
  TH3F *hDPint;
  TH3F *hDZint;
  TH3F *hDRdiff;
  TH3F *hDPdiff;
  TH3F *hDZdiff;
  TH3F *three_d_startmap;
  TGraph *Graph;
  TGraph *CM;
  TAxis *xaxis, *yaxis, *zaxis; 
  TH1 *dlong;
  TH1 *dtrans;
  TH1 *deltatime; 
  TH2 *hitmapstart;
  TH2 *hitmapend;
  TH2 *z_startmap;
  TH2 *deltaphi;
  TH2 *deltaphiint;
  TH2 *deltaphidiff;
  TH2 *deltaphidifference;
  TH2 *deltaphidifferencepercent;
  TH2 *deltar;
  TH2 *deltarint;
  TH2 *deltardiff;
  TH2 *deltardifference;
  TH2 *deltardifferencepercent;
  TH2 *deltaphinodiff;
  TH2 *deltaphinodist;
  TH2 *deltarnodiff;
  TH2 *deltarnodist;
  TH2 *deltaz;
  TFile *m_outf;
  TFile *EDrift_outf;
  TFile *CM_outf;
  TNtuple *nt;
  TNtuple *nthit;
  TNtuple *ntfinalhit;
  TNtuple *ntpad;
  std::string detector;
  std::string hitnodename;
  std::string seggeonodename;
  bool do_Centralmem,do_diff_SC_Distortion,do_Int_SC_Distortion;
  unsigned int seed,print_layer;
  int nBinZ, nBinR,nBinP,e_num;
  double x_start,y_start,x_final,y_final; 
  double Start_x;
  double Start_y;
  double Start_z;
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
