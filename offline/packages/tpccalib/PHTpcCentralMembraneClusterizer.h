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

  virtual ~PHTpcCentralMembraneClusterizer();

  void set_process(const int proc)  { _process = proc;  }
  void set_histos_on(const bool val) {_histos = val;}
  void set_min_adc_value(const unsigned int val) {_min_adc_value = val;}
  void set_min_z_value(const double val) {_min_z_value = val;}
  void set_stripe_dr_values(const double dr1, const double dr2, const double dr3){ _cmclus_dr_inner = dr1; _cmclus_dr_mid = dr2; _cmclus_dr_outer = dr3;}

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

  TH1F *henergy;
  TH1F *hz;
  TH2F *hxy;
  TH1F *hDist;
  TH2F *hDistRow;
  TH1F *hDist2;
  TH2F *hDistRowAdj;
  TH1F *hDist2Adj;
  TH1F *hClustE[3];
  
  int _process = 0;
  unsigned int _min_adc_value = 0;
  double _min_z_value = 1.0;
  double _cmclus_dr_inner = 0.51;  //cm
  double _cmclus_dr_mid = 0.95;  //cm
  double _cmclus_dr_outer = 1.025;  //cm

  bool _histos = false;

  TFile *fout;

};

#endif // PHTPCCENTRALMEMBRANECLUSTERIZER_H
