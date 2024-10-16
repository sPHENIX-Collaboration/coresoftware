#include "TriggerAnalyzer.h"
#include <phool/PHNode.h>
#include <phool/getClass.h>

int TriggerAnalyzer::decodeTriggers(PHCompositeNode *topNode)
{
  gl1packet = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
  if (!gl1packet) 
    {
      std::cout << " no gl1 packet" << std::endl;
      return 1;
    }
  triggerruninfo = findNode::getClass<TriggerRunInfo>(topNode, "TriggerRunInfo");
  if (!triggerruninfo) 
    {
      std::cout << " no triggerruninfo" << std::endl;
      return 1;
    }

  gl1_scaledvec = gl1packet->lValue(0, "ScaledVector");
  gl1_livevec = gl1packet->lValue(0, "TriggerVector");
  gl1_bco = gl1packet->lValue(0, "BCO");

  return 0;
}

bool TriggerAnalyzer::didTriggerFire(const std::string& triggername)
{
  uint32_t bit = triggerruninfo->getTriggerBitByName(triggername);
  return (((gl1_scaledvec >> bit) & 0x1U) == 0x1U);
}

bool TriggerAnalyzer::didTriggerFire(int triggerbit)
{
  uint32_t bit = (uint32_t) triggerbit;
  return (((gl1_scaledvec >> bit) & 0x1U) == 0x1U);

}

int TriggerAnalyzer::getTriggerPrescale(const std::string& triggername)
{
  return triggerruninfo->getPrescaleByName(triggername);
}

int TriggerAnalyzer::getTriggerPrescale(int triggerbit)
{
  return triggerruninfo->getPrescaleByBit(triggerbit);
}

bool TriggerAnalyzer::checkRawTrigger(const std::string& triggername)
{
  uint32_t bit = triggerruninfo->getTriggerBitByName(triggername);
  return (((gl1_livevec >> bit) & 0x1U) == 0x1U);
  
}

bool TriggerAnalyzer::checkRawTrigger(int triggerbit)
{
  uint32_t bit = (uint32_t) triggerbit;
  return (((gl1_livevec >> bit) & 0x1U) == 0x1U);
}

std::string TriggerAnalyzer::getTriggerName(int triggerbit)
{
  return triggerruninfo->getTriggerName(triggerbit);
}

uint64_t TriggerAnalyzer::getTriggerScalers(const std::string& triggername)
{
  return triggerruninfo->getScalersByName(triggername);
}

uint64_t TriggerAnalyzer::getTriggerScalers(int triggerbit)
{
  return triggerruninfo->getScalersByBit(triggerbit);
}

uint64_t TriggerAnalyzer::getTriggerLiveScalers(const std::string& triggername)
{
  return triggerruninfo->getLiveScalersByName(triggername);
}

uint64_t TriggerAnalyzer::getTriggerLiveScalers(int triggerbit)
{
  return triggerruninfo->getLiveScalersByBit(triggerbit);
}

uint64_t TriggerAnalyzer::getTriggerRawScalers(const std::string& triggername)
{
  return triggerruninfo->getRawScalersByName(triggername);
}

uint64_t TriggerAnalyzer::getTriggerRawScalers(int triggerbit)
{
  return triggerruninfo->getRawScalersByBit(triggerbit);
}
