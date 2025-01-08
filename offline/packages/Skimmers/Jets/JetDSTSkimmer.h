// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef JETDSTSKIMMER_H
#define JETDSTSKIMMER_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

class JetDSTSkimmer : public SubsysReco
{
 public:

  JetDSTSkimmer(const std::string &name = "JetDSTSkimmer");

// please declare default dtor/ctors in the header file
  ~JetDSTSkimmer() override = default;

  int process_event(PHCompositeNode *topNode) override;

  void SetMinJetPt(float minJetPt) { m_minJetPt = minJetPt; }
  void SetMinClusterPt(float minClusterPt) { m_minClusterPt = minClusterPt; }

  void SetJetNodeName(const std::string &jetNodeName) { m_JetNodeName = jetNodeName; }
  void SetClusterNodeName(const std::string &clusterNodeName) { m_ClusterNodeName = clusterNodeName; }

 private:
    bool isBackgroundEvent();

    float m_minJetPt{10};
    float m_minClusterPt{5};
    
    std::string m_JetNodeName{"AntiKt_Tower_r04_Sub1"};
    std::string m_ClusterNodeName{"CLUSTERINFO_CEMC"};
};

#endif // JETDSTSKIMMER_H
