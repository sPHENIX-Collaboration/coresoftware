#ifndef KFPARTICLESPHENIX_KFPARTICLETRIGGERINFO_H
#define KFPARTICLESPHENIX_KFPARTICLETRIGGERINFO_H

#include <calotrigger/TriggerAnalyzer.h>

#include <limits>

class Gl1Packet;
class TTree;

class KFParticle_triggerInfo
{
 public:
  KFParticle_triggerInfo();  // Constructor

  virtual ~KFParticle_triggerInfo();  // Destructor

  protected:
    TriggerAnalyzer *triggeranalyzer{nullptr};

    bool buildTriggerBranches(PHCompositeNode *topNode, TTree *m_tree);
    void fillTriggerBranches(PHCompositeNode *topNode);
    void resetTriggerBranches();

  private:
    static const int nTriggerBits = 64;
    bool m_trigger_bit[nTriggerBits] = {std::numeric_limits<bool>::quiet_NaN()};

};

#endif // KFPARTICLESPHENIX_KFPARTICLETRIGGERINFO_H
