
#ifndef EVENTFLOWRECO_H
#define EVENTFLOWRECO_H

//===========================================================
/// \author Tanner
//===========================================================

#include <fun4all/SubsysReco.h>

#include "EventFlowInfo.h"

class PHCompositeNode;

class EventFlowReco : public SubsysReco {
public:
  
  EventFlowReco(const std::string &name = "EventFlowReco");
  ~EventFlowReco() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode * /*topNode*/) override;

  void SetFlowSrc(EventFlowInfo::EventFlowSrc src) { _src = src; }

private:

  int CreateNodes(PHCompositeNode *topNode);
  EventFlowInfo::EventFlowSrc _src = EventFlowInfo::VOID;
  
};

#endif // EVENTFLOWRECO_H
