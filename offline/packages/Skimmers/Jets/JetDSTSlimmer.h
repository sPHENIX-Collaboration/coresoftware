// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef JETDSTSLIMMER_H
#define JETDSTSLIMMER_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

class JetDSTSlimmer : public SubsysReco
{
 public:

  JetDSTSlimmer(const std::string &name = "JetDSTSlimmer");

  ~JetDSTSlimmer() override;

  int process_event(PHCompositeNode *topNode) override;

  void SetMinJetPt(float minJetPt) { m_minJetPt = minJetPt; }
  void SetMinClusterPt(float minClusterPt) { m_minClusterPt = minClusterPt; }

 private:
    float m_minJetPt{10};
    float m_minClusterPt{5};
    
    std::string m_JetNodeName{"AntiKt_Tower_r04_Sub1"};
    std::string m_ClusterNodeName{"CEMC_CALO_CLUSTER"};

    bool isBackgroundEvent();
};

#endif // JETDSTSLIMMER_H
