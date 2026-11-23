#include "KFParticle_triggerInfo.h"

#include <calotrigger/TriggerAnalyzer.h>
#include <calotrigger/TriggerRunInfo.h>
#include <ffarawobjects/Gl1Packet.h>
#include <phool/getClass.h>

#include <TTree.h>

#include <map>
#include <string>

KFParticle_triggerInfo::KFParticle_triggerInfo() = default;  // Constructor

KFParticle_triggerInfo::~KFParticle_triggerInfo() = default;  // Destructor

bool KFParticle_triggerInfo::initializeBranches(PHCompositeNode *topNode, TTree *m_tree)
{
  //First, check whether we actually have the right information
  auto *gl1packet =findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");
  if (!gl1packet)
  {
    gl1packet = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    if (!gl1packet)
    {
      return false;
    }
  }

  auto *triggerruninfo =findNode::getClass<TriggerRunInfo>(topNode, "TriggerRunInfo");
  if (!triggerruninfo)
  {
    return false;
  }

  size_t pos;
  const std::map<std::string, std::string> forbiddenStrings = {
      {" ", "_"}, {",", ""}, {"/", "_"}, {"&", "and"}, {"=", "eq"},
      {"<", "l"}, {">", "g"}, {"+", "plus"}, {"-", "minus"}, {"*", "star"}};

  std::string trigger_header = "trigger_";

  if (!triggeranalyzer)
  {
    triggeranalyzer = std::make_unique<TriggerAnalyzer>();
  }
  triggeranalyzer->decodeTriggers(topNode);

  for (int i = 0; i < nTriggerBits; ++i)
  {
    std::string triggerName = triggeranalyzer->getTriggerName(i);

    if (triggerName.find("unknown") != std::string::npos)
    {
      continue;
    }

    std::string branchName = trigger_header + triggerName;

    for (const auto& [badString, goodString] : forbiddenStrings)
    {
      while ((pos = branchName.find(badString)) != std::string::npos)
      {
        branchName.replace(pos, 1, goodString);
      }
    }

    std::string branchNameType = branchName + "/O";
    m_tree->Branch(branchName.c_str(), &m_trigger_bit[i], branchNameType.c_str());
  }

  m_available = true;
  return true;
}

void KFParticle_triggerInfo::fill(PHCompositeNode *topNode)
{
  triggeranalyzer->decodeTriggers(topNode);
  for (int i = 0; i < nTriggerBits; ++i)
  {
    std::string triggerName = triggeranalyzer->getTriggerName(i);

    if (triggerName.find("unknown") != std::string::npos)
    {
      continue;
    }

    m_trigger_bit[i] = triggeranalyzer->didTriggerFire(i);
  }
}

void KFParticle_triggerInfo::clear()
{
  for (bool &bit : m_trigger_bit)
  {
    bit = std::numeric_limits<bool>::quiet_NaN();
  }
}
