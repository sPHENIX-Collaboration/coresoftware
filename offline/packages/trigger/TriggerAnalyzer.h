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

  int getTriggerPrescale(const std::string& triggername);
  int getTriggerPrescale(int triggerbit);

  bool checkRawTrigger(const std::string& triggername);
  bool checkRawTrigger(int triggerbit);

 private:

  Gl1Packet *gl1packet{nullptr};
  TriggerRunInfo *triggerruninfo{nullptr};

  uint64_t gl1_scaledvec{0};
  uint64_t gl1_livevec{0};
  uint64_t gl1_bco{0};
};

#endif /* CALOTRIGGER_TRIGGERANALYZER_H */

