// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef JETDSTSKIMMER_H
#define JETDSTSKIMMER_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <map>

class PHCompositeNode;

class JetDSTSkimmer : public SubsysReco
{
 public:

  JetDSTSkimmer(const std::string &name = "JetDSTSkimmer");

// please declare default dtor/ctors in the header file
  ~JetDSTSkimmer() override = default;

  int process_event(PHCompositeNode *topNode) override;

  void SetJetNodeThresholds(const std::string &jetNodeName, float threshold) { m_JetNodePts[jetNodeName] = threshold; }
  void SetClusterNodeThresholds(const std::string &clusterNodeName, float threshold) { m_ClusterNodePts[clusterNodeName] = threshold; }

  void SetJetNodeThresholds(const std::map<std::string, float> &jetNodePts) { m_JetNodePts = jetNodePts; }
  void SetClusterNodeThresholds(const std::map<std::string, float> &clusterNodePts) { m_ClusterNodePts = clusterNodePts; }

  void ResetJetNodeThresholds() { m_JetNodePts.clear(); }
  void ResetClusterNodeThresholds() { m_ClusterNodePts.clear(); }

 private:
    bool isBackgroundEvent();

    std::map<std::string, float> m_JetNodePts;
    std::map<std::string, float> m_ClusterNodePts;
};

#endif // JETDSTSKIMMER_H
