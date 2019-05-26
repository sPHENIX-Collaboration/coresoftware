#ifndef G4EVAL_TRKREVALUATOR_H
#define G4EVAL_TRKREVALUATOR_H


#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class PHTimer;
class TFile;
class TNtuple;

class TrkrEvaluator : public SubsysReco
{
 public:
  TrkrEvaluator(const std::string &name = "TRKREVALUATOR",
                const std::string &filename = "g4eval.root",
                const std::string &trackmapname = "SvtxTrackMap",
                unsigned int nlayers_maps = 3,
                unsigned int nlayers_intt = 4,
                unsigned int nlayers_tpc = 48);
  virtual ~TrkrEvaluator() {}

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  
  void scan_for_embedded(bool b);
  
 private:
  unsigned int _ievent;

  //----------------------------------
  // evaluator output ntuples

  bool _do_vertex_eval;
  bool _do_cluster_eval;
  bool _do_gtrack_eval;
  bool _do_track_eval;
  bool _do_track_match;
  bool _scan_for_embedded;


  unsigned int _nlayers_maps;
  unsigned int _nlayers_intt;
  unsigned int _nlayers_tpc;

  TNtuple *_ntp_vertex;
  TNtuple *_ntp_cluster;
  TNtuple *_ntp_gtrack;
  TNtuple *_ntp_track;

  // evaluator output file
  std::string _filename;
  std::string _trackmapname;

  TFile *_tfile;

  PHTimer *_timer;

  float line_circle_intersection(float x[], float y[], float z[], float radius);

  // output subroutines
  void fillOutputNtuples(PHCompositeNode *topNode);  ///< dump the evaluator information into ntuple for external analysis
  void printInputInfo(PHCompositeNode *topNode);     ///< print out the input object information (debugging upstream components)
  void printOutputInfo(PHCompositeNode *topNode);    ///< print out the ancestry information for detailed diagnosis
};

#endif  // G4EVAL_SVTXEVALUATOR_H
