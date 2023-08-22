// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef __BBCRECO_H__
#define __BBCRECO_H__

#include <fun4all/SubsysReco.h>

#include <memory>
#include <string>

class PHCompositeNode;
class BbcPmtContainer;
class BbcVertexMap;
class BbcOut;
class TF1;
class TH1;

class BbcReco : public SubsysReco
{
 public:
  BbcReco(const std::string &name = "BbcReco");

  ~BbcReco() override;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 private:
  int createNodes(PHCompositeNode *topNode);
  int getNodes(PHCompositeNode *topNode);
  std::unique_ptr<TF1> m_gaussian = nullptr;

  float m_tres = 0.05;
  TH1 *h_evt_bbct[2];
  BbcVertexMap *m_bbcvertexmap = nullptr;
  BbcOut *m_bbcout {nullptr};
  BbcPmtContainer *m_bbcpmts = nullptr;
};

#endif  // __BBCRECO_H__
