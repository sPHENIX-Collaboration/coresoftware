// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFARAWMODULES_EVENTCOMBINER_H
#define FFARAWMODULES_EVENTCOMBINER_H

#include <fun4all/SubsysReco.h>

#include <set>
#include <string>

class Event;
class Fun4AllInputManager;
class PHCompositeNode;

class EventCombiner : public SubsysReco
{
 public:
  EventCombiner(const std::string &name = "EventCombiner");

  ~EventCombiner() override {}

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int ResetEvent(PHCompositeNode *topNode) override;

  void AddPrdfInputNodeFromManager(const Fun4AllInputManager *in);
  void AddPrdfInputNodeName(const std::string &name);

 private:
  Event *m_Event = nullptr;
  int *m_OutArray = nullptr;

  std::string m_PrdfOutputNodeName = "PRDF";
  std::set<std::string> m_PrdfInputNodeNameSet;
};

#endif  // FFARAWMODULES_EVENTCOMBINER_H
