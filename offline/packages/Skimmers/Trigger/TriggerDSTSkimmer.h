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

 private:

  std::vector<int> m_trigger_index{10}; 
  int ievent{0};
};

#endif // JETDSTSKIMMER_H
