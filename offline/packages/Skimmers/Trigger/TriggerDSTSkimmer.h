// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TRIGGERDSTSKIMMER_H
#define TRIGGERDSTSKIMMER_H

#include <fun4all/SubsysReco.h>

#include <vector>
#include <string>

class PHCompositeNode;

class TriggerDSTSkimmer : public SubsysReco
{
 public:

  TriggerDSTSkimmer(const std::string &name = "TriggerDSTSkimmer");

  ~TriggerDSTSkimmer() override = default;

  int process_event(PHCompositeNode *topNode) override;

  void SetTrigger(std::vector<int> &trigger_vector) {m_trigger_index = trigger_vector;}

  void set_accept_max(int max_events) 
  {
    use_max_accept = true;
    max_accept = max_events;
    return;
  }

 private:

  std::vector<int> m_trigger_index{10}; 
  int ievent{0};

  int accepted_events{0};
  int max_accept{0};
  bool use_max_accept{false};

};

#endif // JETDSTSKIMMER_H
