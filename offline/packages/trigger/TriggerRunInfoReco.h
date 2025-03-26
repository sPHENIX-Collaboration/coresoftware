#ifndef TRIGGER_TRIGGERRUNINFORECO_H
#define TRIGGER_TRIGGERRUNINFORECO_H

#include <fun4all/SubsysReco.h>

#include <string>

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

  static void SetTriggerEmulator(TriggerRunInfo *triggerRunInfo);

  static int fetchTriggerPrescales(int runnumber, TriggerRunInfo *triggerRunInfo);
  static int fetchTriggerScalers(int runnumber, TriggerRunInfo *triggerRunInfo);
};

#endif /* TRIGGER_TRIGGERRUNINFORECO_H */
