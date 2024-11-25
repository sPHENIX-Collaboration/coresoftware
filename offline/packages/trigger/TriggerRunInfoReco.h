#ifndef TRIGGER_TRIGGERRUNINFORECO_H
#define TRIGGER_TRIGGERRUNINFORECO_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class TriggerRunInfo;

class TriggerRunInfoReco : public SubsysReco
{
 public:
  TriggerRunInfoReco(const std::string &name = "TriggerRunInfoReco");
  ~TriggerRunInfoReco() override = default;
  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;

  void UseEmulator(bool use) { m_useEmulator = use; }

  
 private:
  bool m_useEmulator{false};
  
  void SetTriggerEmulator(TriggerRunInfo *triggerRunInfo);

  int fetchTriggerPrescales(int runnumber, TriggerRunInfo *triggerRunInfo);
  int fetchTriggerScalers(int runnumber, TriggerRunInfo *triggerRunInfo);

};

#endif /* TRIGGER_TRIGGERRUNINFORECO_H */
