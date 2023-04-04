// Tell emacs that this is a C++ source
// -*- C++ -*-.

#ifndef G4TPC_PHG4TPCELECTRONDRIFT_H
#define G4TPC_PHG4TPCELECTRONDRIFT_H

#include <array>
#include <cmath>
#include <fstream>
#include <fun4all/SubsysReco.h>
#include <g4main/PHG4HitContainer.h>
#include <gsl/gsl_rng.h>
#include <memory>
#include <phparameter/PHParameterInterface.h>
#include <string>
#include <trackbase/ActsGeometry.h>

class PHG4TpcPadPlane;
class PHG4TpcDistortion;
class PHCompositeNode;
class TH1;
class TH2;
class TNtuple;
class TTree;
class TFile;
class TrkrHitSetContainer;
class TrkrHitTruthAssoc;
class TrkrTruthTrackContainer;
class TrkrClusterContainer;
class TrkrTruthTrack;
class DistortedTrackContainer;
class TpcClusterBuilder;
class PHG4TpcCylinderGeomContainer;

class genElecs{
 public:
  genElecs();
  ~genElecs();

  void set_xStart(double a_x){ m_xStart = a_x; };
  void set_yStart(double a_y){ m_yStart = a_y; };
  void set_zStart(double a_z){ m_zStart = a_z; };
  void set_tStart(double a_t){ m_yStart = a_t; };

  void set_xEnd(double a_x){ m_xEnd = a_x; };
  void set_yEnd(double a_y){ m_yEnd = a_y; };
  void set_zEnd(double a_z){ m_zEnd = a_z; };
  void set_tEnd(double a_t){ m_tEnd = a_t; };

  double get_xStart(){ return m_xStart; };
  double get_yStart(){ return m_yStart; };
  double get_zStart(){ return m_zStart; };
  double get_tStart(){ return m_tStart; };

  double get_xEnd(){ return m_xEnd; };
  double get_yEnd(){ return m_yEnd; };
  double get_zEnd(){ return m_zEnd; };
  double get_tEnd(){ return m_tEnd; };

  //void resetElec(){ m_xStart=0.0; m_yStart=0.0; m_zStart=0.0; m_tStart=0.0; m_xEnd=0.0; m_yEnd=0.0; m_zEnd=0.0; m_tEnd=0.0; };
 
 private:
  double m_xStart;
  double m_yStart;
  double m_zStart;
  double m_tStart;

  double m_xEnd;
  double m_yEnd;
  double m_zEnd;
  double m_tEnd;
  
};

class g4HitClass{
 public:
  g4HitClass();
  ~g4HitClass();

  void set_g4xStart(double a_x){ m_g4xStart = a_x; };
  void set_g4yStart(double a_y){ m_g4yStart = a_y; };
  void set_g4zStart(double a_z){ m_g4zStart = a_z; };
  void set_g4tStart(double a_t){ m_g4tStart = a_t; };

  void set_g4xEnd(double a_x){ m_g4xEnd = a_x; };
  void set_g4yEnd(double a_y){ m_g4yEnd = a_y; };
  void set_g4zEnd(double a_z){ m_g4zEnd = a_z; };
  void set_g4tEnd(double a_t){ m_g4tEnd = a_t; };

  void set_g4E(double a_E){ m_g4E = a_E; };

  void addElectron(genElecs *a_elec){ m_elecs.push_back(a_elec); };

  double get_g4xStart(){ return m_g4xStart; };
  double get_g4yStart(){ return m_g4yStart; };
  double get_g4zStart(){ return m_g4zStart; };
  double get_g4tStart(){ return m_g4tStart; };

  double get_g4xEnd(){ return m_g4xEnd; };
  double get_g4yEnd(){ return m_g4yEnd; };
  double get_g4zEnd(){ return m_g4zEnd; };
  double get_g4tEnd(){ return m_g4tEnd; };

  double get_g4E(){ return m_g4E; };
  int get_nElec(){ return (int)m_elecs.size(); };
   
  genElecs* getElectron(int elecIndex){ return m_elecs[elecIndex]; };

  void reset_g4_electrons(){ m_elecs.clear(); };

  //void reset_g4Hits(){ m_g4x=0.0; m_g4y=0.0; m_g4z=0.0; m_g4t=0.0; m_g4E=0.0; m_elecs.clear(); };


 private:

  double m_g4xStart;
  double m_g4yStart;
  double m_g4zStart;
  double m_g4tStart;

  double m_g4xEnd;
  double m_g4yEnd;
  double m_g4zEnd;
  double m_g4tEnd;

  double m_g4E;

  std::vector<genElecs*> m_elecs;

};


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

  std::ofstream f_out;

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


  TFile *m_outfile;
  TTree *m_hitTree;
  g4HitClass *m_g4Hits;


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
