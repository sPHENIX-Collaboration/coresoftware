// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHTPCCENTRALMEMBRANECLUSTERIZER_H
#define PHTPCCENTRALMEMBRANECLUSTERIZER_H

#include <string>

#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrClusterContainer.h>

class PHCompositeNode;
class SvtxTrackMap;
class SvtxTrack;
class SvtxVertexMap;
class TrkrClusterHitAssoc;
class TrkrHitSetContainer;
class CMFlashClusterContainer;

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

  TrkrClusterContainer *_cluster_map{nullptr};
  CMFlashClusterContainer *_corrected_CMcluster_map{nullptr};
  TrkrClusterHitAssoc *_cluster_hit_map{nullptr};
  TrkrHitSetContainer *_hitset_map{nullptr};

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
  unsigned int _min_adc_value = 180;
  bool _histos = false;

  TFile *fout;

};

#endif // PHTPCCENTRALMEMBRANECLUSTERIZER_H
