#ifndef CALOTRIGGER_TRIGGERANALYZER_H
#define CALOTRIGGER_TRIGGERANALYZER_H

#include <string>
#include <phool/PHCompositeNode.h>
#include <ffarawobjects/Gl1Packetv2.h>
#include <ffarawobjects/Gl1Packet.h>
#include "TriggerRunInfo.h"
#include "TriggerRunInfov1.h"

class TriggerAnalyzer
{
 public:
  TriggerAnalyzer() = default;
  ~TriggerAnalyzer();
  
  int decodeTriggers(PHCompositeNode *topNode);

  bool didTriggerFire(const std::string& triggername);
  bool didTriggerFire(int triggerbit);

  bool checkRawTrigger(const std::string& triggername);
  bool checkRawTrigger(int triggerbit);

  int getTriggerPrescale(const std::string& triggername);
  int getTriggerPrescale(int triggerbit);

  std::string getTriggerName(int triggerbit);

  uint64_t getTriggerRawScalers(const std::string& triggername);
  uint64_t getTriggerRawScalers(int triggerbit);

  uint64_t getTriggerLiveScalers(const std::string& triggername);
  uint64_t getTriggerLiveScalers(int triggerbit);

  uint64_t getTriggerScalers(const std::string& triggername);
  uint64_t getTriggerScalers(int triggerbit);

 private:

  Gl1Packet *gl1packet{nullptr};
  TriggerRunInfo *triggerruninfo{nullptr};

  uint64_t gl1_scaledvec{0};
  uint64_t gl1_livevec{0};
  uint64_t gl1_bco{0};

};

#endif /* CALOTRIGGER_TRIGGERANALYZER_H */

