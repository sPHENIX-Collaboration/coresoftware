#ifndef TRIGGERVALID_TRIGGERVALID_H
#define TRIGGERVALID_TRIGGERVALID_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

// Forward declarations
class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TNtuple;
class TTree;
class TH2F;
class TH1F;
class TH1;
class TProfile2D;

class TriggerValid : public SubsysReco
{
 public:
  //! constructor
  TriggerValid(const std::string& name = "TriggerValid");

  //! destructor
  ~TriggerValid() override = default;

  //! full initialization
  int Init(PHCompositeNode*) override;

  //! event processing method
  int process_event(PHCompositeNode*) override;

  int process_towers(PHCompositeNode*);
  int process_primitives(PHCompositeNode*);
  int process_ll1out(PHCompositeNode*);

  void Trigger(const std::string& name) { trigger = name; }

  void set_debug(bool debug) { m_debug = debug; }

 private:
  int Getpeaktime(TH1* h);

  bool m_debug{0};
  std::string trigger;

  int _eventcounter{0};
  int _range{1};
};

#endif
