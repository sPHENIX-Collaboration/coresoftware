// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFARAWMODULES_EVENTNUMBERCHECK_H
#define FFARAWMODULES_EVENTNUMBERCHECK_H

#include <fun4all/SubsysReco.h>

#include <set>
#include <string>

class Event;
class Fun4AllInputManager;
class Packet;
class PHCompositeNode;

class EventNumberCheck : public SubsysReco
{
 public:
  EventNumberCheck(const std::string &name = "EventNumberCheck");

  ~EventNumberCheck() override {}

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

//  int ResetEvent(PHCompositeNode *topNode) override;

  void MyPrdfNode(const std::string &name) {m_MyPrdfNode = name;}

 private:
  void CheckFem(int nw);
  Packet *plist[10000]{0};
  int previous_event_clkdiff {0};
  std::set<int> m_EventSeen;
  std::string m_MyPrdfNode {"PRDF"};
};

#endif  // FFARAWMODULES_EVENTCOMBINER_H
