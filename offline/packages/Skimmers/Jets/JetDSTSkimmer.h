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

  ~JetDSTSkimmer() override;

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

#endif // JETDSTSKIMMER_H
