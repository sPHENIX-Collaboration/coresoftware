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

 private:
  int fetchTriggerPrescales(int runnumber, TriggerRunInfo *triggerRunInfo);
  int fetchTriggerScalers(int runnumber, TriggerRunInfo *triggerRunInfo);

  std::string m_dbConnectionStr;
};

#endif /* TRIGGER_TRIGGERRUNINFORECO_H */
