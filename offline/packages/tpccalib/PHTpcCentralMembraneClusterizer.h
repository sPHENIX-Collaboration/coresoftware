// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHTPCCENTRALMEMBRANECLUSTERIZER_H
#define PHTPCCENTRALMEMBRANECLUSTERIZER_H

#include <string>

#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrClusterContainer.h>
#include <tpc/TpcDistortionCorrectionContainer.h>
#include <tpc/TpcDistortionCorrection.h>

#include <memory>

class PHCompositeNode;
class SvtxTrackMap;
class SvtxTrack;
class SvtxVertexMap;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class CMFlashClusterContainer;
class  PHG4TpcCylinderGeomContainer;

class TF1;
class TNtuple;
class TFile;
class TH1F;
class TH2F;

class PHTpcCentralMembraneClusterizer : public SubsysReco
{
 public:

  PHTpcCentralMembraneClusterizer(const std::string &name = "PHTpcCentralMembraneClusterizer");

  void set_histos_on(const bool val) {_histos = val;}

  /// output file name for evaluation histograms
  void set_histogram_outputfile(const std::string &outputfile) 
  {m_histogramfilename = outputfile;}
  
  void set_min_adc_value(const unsigned int val) {_min_adc_value = val;}
  void set_min_z_value(const double val) {_min_z_value = val;}
  void set_stripe_dr_values(const double dr1, const double dr2, const double dr3){ _cmclus_dr_inner = dr1; _cmclus_dr_mid = dr2; _cmclus_dr_outer = dr3;}

  void set_modulo_threshold( int val ) { m_moduloThreshold = val; }
  void set_metaCluster_threshold( int val ) { m_metaClusterThreshold = val; }

 //! run initialization
  int InitRun(PHCompositeNode *topNode);

  //! event processing
  int process_event(PHCompositeNode *topNode);

  //! end of process
  int End(PHCompositeNode * topNode);

 protected:
  
 private:

  int GetNodes(PHCompositeNode* topNode);

  std::string _track_map_name;

 /// tpc distortion correction utility class
  TpcDistortionCorrection _distortionCorrection;

  TrkrClusterContainer *_cluster_map{nullptr};
  CMFlashClusterContainer *_corrected_CMcluster_map{nullptr};
  PHG4TpcCylinderGeomContainer *_geom_container{nullptr};
  TpcDistortionCorrectionContainer* _dcc{nullptr};

  ///@name counters
  //@{
  int m_total_clusters = 0;
  int m_accepted_clusters = 0;
  int m_cm_clusters = 0;  
  int m_cm_clusters_size1 = 0;  
  int m_cm_clusters_size2 = 0;  
  //@}

  int m_moduloThreshold = 5;
  int m_metaClusterThreshold = 18;
  
  bool _histos = false;
  TH1F *henergy = nullptr;
  TH1F *hz = nullptr;
  TH1F *hz_pos = nullptr;
  TH1F *hz_neg = nullptr;
  TH2F *hxy = nullptr;
  TH1F *hDist = nullptr;
  TH2F *hDistRow = nullptr;
  TH1F *hDist2 = nullptr;
  TH2F *hDistRowAdj = nullptr;
  TH1F *hDist2Adj = nullptr;
  TH1F *hClustE[3] = {nullptr};

  TH2F *hrPhi_reco_petalModulo_pos = nullptr;
  TH2F *hrPhi_reco_petalModulo_neg = nullptr;
  
  TH1F *hphi_reco_pos[48] = {nullptr};
  TH1F *hphi_reco_neg[48] = {nullptr};

  TH1F *hphi_reco_pair_pos[47] = {nullptr};
  TH1F *hphi_reco_pair_neg[47] = {nullptr};

  int nPairAbove_pos[47] = {0};
  int nPairAbove_neg[47] = {0};

  double pairAboveContent_pos[47] = {0.0};
  double pairAboveContent_neg[47] = {0.0};

  std::string m_histogramfilename = "PHTpcCentralMembraneClusterizer.root";
  std::unique_ptr<TFile> m_histogramfile;

  unsigned int _min_adc_value = 0;
  double _min_z_value = 0.0;
  double _cmclus_dr_inner = 0.51;  //cm
  double _cmclus_dr_mid = 0.95;  //cm
  double _cmclus_dr_outer = 1.025;  //cm

};

#endif // PHTPCCENTRALMEMBRANECLUSTERIZER_H
