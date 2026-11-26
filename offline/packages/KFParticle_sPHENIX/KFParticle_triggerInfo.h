#ifndef KFPARTICLESPHENIX_KFPARTICLETRIGGERINFO_H
#define KFPARTICLESPHENIX_KFPARTICLETRIGGERINFO_H

#include <limits>
#include <memory>

class PHCompositeNode;
class TriggerAnalyzer;
class TTree;

// Event-level GL1 trigger bits (one bool branch per named trigger) for the nTuple
class KFParticle_triggerInfo
{
 public:
  KFParticle_triggerInfo();  // Constructor

  ~KFParticle_triggerInfo();  // Destructor; out of line because TriggerAnalyzer is only forward-declared here

  KFParticle_triggerInfo(const KFParticle_triggerInfo &) = delete;
  KFParticle_triggerInfo &operator=(const KFParticle_triggerInfo &) = delete;
  KFParticle_triggerInfo(KFParticle_triggerInfo &&) = delete;
  KFParticle_triggerInfo &operator=(KFParticle_triggerInfo &&) = delete;

  // Books the trigger branches; returns false (and books nothing) if GL1 or TriggerRunInfo is missing
  bool initializeBranches(PHCompositeNode *topNode, TTree *m_tree);
  void fill(PHCompositeNode *topNode);
  void clear();

  // True once initializeBranches has found the trigger information
  bool isAvailable() const { return m_available; }

 private:
  std::unique_ptr<TriggerAnalyzer> triggeranalyzer;
  bool m_available{false};

  static const int nTriggerBits = 64;
  bool m_trigger_bit[nTriggerBits] = {std::numeric_limits<bool>::quiet_NaN()};
};

#endif // KFPARTICLESPHENIX_KFPARTICLETRIGGERINFO_H
