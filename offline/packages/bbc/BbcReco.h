// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BBCRECONSTRUCTION_H
#define BBCRECONSTRUCTION_H

#include <fun4all/SubsysReco.h>

#include <memory>
#include <string>

class PHCompositeNode;
class BbcPmtContainer;
class BbcVertexMap;
class TF1;
class TH1;

class BbcReconstruction : public SubsysReco
{
 public:
  BbcReconstruction(const std::string &name = "BbcReconstruction");

  ~BbcReconstruction() override;

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
  BbcPmtContainer *m_bbcpmts = nullptr;
};

#endif  // BBCRECONSTRUCTION_H
