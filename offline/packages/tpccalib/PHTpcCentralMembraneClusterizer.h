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

  void set_track_map_name(const std::string &map_name) { _track_map_name = map_name; }
  void set_process(const unsigned int proc)
  {
    process = proc;
  }

 //! run initialization
  int InitRun(PHCompositeNode *topNode);

  //! event processing
  int process_event(PHCompositeNode *topNode);

  int EndRun(PHCompositeNode *topNode);

  //! end of process
int End(PHCompositeNode *topNode);

 protected:
  
 private:

  int GetNodes(PHCompositeNode* topNode);

  std::string _track_map_name;

  SvtxTrackMap *_track_map{nullptr};
  SvtxTrack *_track{nullptr};
  SvtxVertexMap *_vertex_map{nullptr};
  TrkrClusterContainer *_cluster_map;
  CMFlashClusterContainer *_corrected_CMcluster_map{nullptr};
  TrkrClusterHitAssoc *_cluster_hit_map;
  TrkrHitSetContainer *_hitset_map;
  const bool clst_track = true;


  TH1F *henergy;
  TH1F *hz;
  TH2F *hxy;
  TH1F *hDist;
  TH2F *hDistRow;
  TH1F *hDist2;
  TH2F *hDistRowAdj;
  TH1F *hDist2Adj;
  TH1F *hClustE[3];
  
  unsigned int MAXPAD = 2400;
  unsigned int process = 0;



  TNtuple *ntp{nullptr};
  TFile *fout;

};

#endif // PHTPCCENTRALMEMBRANECLUSTERIZER_H
