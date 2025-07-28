#include "TriggerAnalyzer.h"

#include "LL1Out.h"
#include "TriggerRunInfo.h"

#include <ffarawobjects/Gl1Packet.h>

#include <phool/getClass.h>

#include <iostream>

int TriggerAnalyzer::decodeTriggers(PHCompositeNode* topNode)
{
  triggerruninfo = findNode::getClass<TriggerRunInfo>(topNode, "TriggerRunInfo");
  if (!triggerruninfo)
  {
    std::cout << PHWHERE << " no triggerruninfo" << std::endl;
    return 1;
  }

  if (m_useEmulator)
  {
    ll1out_photon = findNode::getClass<LL1Out>(topNode, "LL1OUT_PHOTON");
    if (!ll1out_photon)
    {
      std::cout << PHWHERE << " no trigger emulator" << std::endl;
      return 1;
    }

    ll1out_jet = findNode::getClass<LL1Out>(topNode, "LL1OUT_JET");
    if (!ll1out_jet)
    {
      std::cout << PHWHERE << " no trigger emulator" << std::endl;
      return 1;
    }

    fillTriggerVector();

    return 0;
  }

  gl1packet = findNode::getClass<Gl1Packet>(topNode, 14001);
  if (!gl1packet)
  {
    gl1packet = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    if (!gl1packet)
    {
      gl1packet = findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");  // Different term used in track production
      if (!gl1packet)
      {
	std::cout << PHWHERE << "no gl1 packet" << std::endl;
	return 1;
      }
    }
  }

  gl1_scaledvec = gl1packet->lValue(0, "ScaledVector");
  gl1_livevec = gl1packet->lValue(0, "TriggerVector");
  gl1_bco = gl1packet->lValue(0, "BCO");

  return 0;
}

void TriggerAnalyzer::fillTriggerVector()
{
  gl1_scaledvec = 0x000000000000;
  gl1_livevec = 0x000000000000;
  gl1_bco = 0x000000000000;

  for (int i = 0; i < 4; i++)
  {
    if (ll1out_photon->passesThreshold(i + 1))
    {
      unsigned int bit = i + 28;
      gl1_scaledvec |= (0x1U << bit);
    }
  }
  for (int i = 0; i < 4; i++)
  {
    if (ll1out_jet->passesThreshold(i + 1))
    {
      unsigned int bit = i + 20;
      gl1_scaledvec |= (0x1U << bit);
    }
  }
  gl1_scaledvec &= 0x00000000ffffffff;
  return;
}

bool TriggerAnalyzer::didTriggerFire(const std::string& triggername)
{
  uint32_t bit = triggerruninfo->getTriggerBitByName(triggername);
  return (((gl1_scaledvec >> bit) & 0x1U) == 0x1U);
}

bool TriggerAnalyzer::didTriggerFire(int triggerbit) const
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

bool TriggerAnalyzer::checkRawTrigger(int triggerbit) const
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

void TriggerAnalyzer::Print() const
{
  for (int i = 0; i < 64; i++)
  {
    if (didTriggerFire(i))
    {
      std::cout << " Trigger " << i << " fired" << std::endl;
    }
  }
}
