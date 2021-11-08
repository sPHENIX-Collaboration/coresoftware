#ifndef G4EVAL_MOMENTUMEVALUATOR_H
#define G4EVAL_MOMENTUMEVALUATOR_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class TNtuple;

class MomentumEvaluator : public SubsysReco
{
 public:
  MomentumEvaluator(const std::string &fname, float pt_s = 0.1, float pz_s = 0.2, unsigned int n_l = 62, unsigned int n_i = 2, unsigned int n_r = 50, float i_z = 10., float o_z = 80.);
  ~MomentumEvaluator() override;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 private:
  TNtuple *ntp_true;
  TNtuple *ntp_reco;
  float pt_search_scale;
  float pz_search_scale;
  unsigned int event_counter;
  std::string file_name;
  unsigned int n_inner_layers;
  unsigned int n_required_layers;
  float inner_z_length;
  float outer_z_length;
};

#endif
