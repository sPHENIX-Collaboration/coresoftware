#ifndef G4EVAL_TRKREVALUATOR_H
#define G4EVAL_TRKREVALUATOR_H


#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>
#include <phool/PHTimer.h>
#include <string>

class PHCompositeNode;

class TFile;
class TNtuple;

class TrkrEvaluator : public SubsysReco
{
 public:
  TrkrEvaluator(const std::string &name = "SVTXEVALUATOR",
                const std::string &filename = "g4eval.root",
                unsigned int nlayers_maps = 3,
                unsigned int nlayers_intt = 8,
                unsigned int nlayers_tpc = 60);
  virtual ~TrkrEvaluator() {}

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void set_strict(bool b) { _strict = b; }

  void do_cluster_eval(bool b) { _do_cluster_eval = b; }

 private:
  unsigned int _ievent;

  // eval stack

  //----------------------------------
  // evaluator output ntuples

  bool _strict;
  unsigned int _errors;

  bool _do_cluster_eval;
  bool _do_track_eval;
  bool _do_gtrack_eval;
  bool _do_vertex_eval;

  bool _scan_for_embeded;
  bool _do_track_match;

  unsigned int _nlayers_maps;
  unsigned int _nlayers_intt;
  unsigned int _nlayers_tpc;

  TNtuple *_ntp_vertex;
  TNtuple *_ntp_gpoint;
  TNtuple *_ntp_g4hit;
  TNtuple *_ntp_hit;
  TNtuple *_ntp_cluster;
  TNtuple *_ntp_gtrack;
  TNtuple *_ntp_track;
  TNtuple *_ntp_gseed;

  // evaluator output file
  std::string _filename;
  //Track map name
  TFile *_tfile;

  PHTimer *_timer;

  float line_circle_intersection(float x[], float y[], float z[], float radius);

  // output subroutines
  void fillOutputNtuples(PHCompositeNode *topNode);  ///< dump the evaluator information into ntuple for external analysis
  void printInputInfo(PHCompositeNode *topNode);     ///< print out the input object information (debugging upstream components)
  void printOutputInfo(PHCompositeNode *topNode);    ///< print out the ancestry information for detailed diagnosis
};

#endif  // G4EVAL_SVTXEVALUATOR_H
